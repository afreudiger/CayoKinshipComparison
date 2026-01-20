# Calculation of STR-based relatedness

# read in data
rm(list = ls())
str <- read.table("STR_n98.txt", header = F)
# add column names
colnames(str) <- c("animal_id", paste0("L",rep(seq(1,43, by = 1), each = 2), rep(c("_A","_B"),43)))

# load libraries
library(dplyr) # %>%
library(HardyWeinberg) # HWPerm.mult

# Get unique locus base names (assumes colnames like "L1_A", "L1_B", ...)
locus_names <- unique(sub("_[A-Z]$", "", colnames(str[ , -1])))

# create list and data frame to store results
hwe <- list()
results <- data.frame(Locus = character(0), p_value = numeric(0), stringsAsFactors = FALSE)
# loop through each loci
for (loc in locus_names) {
  # Construct column names for the two alleles at this locus
  colA <- paste0(loc, "_A")
  colB <- paste0(loc, "_B")
  
  # Skip this locus if either allele column is missing from the data frame
  if (!(colA %in% names(str)) || !(colB %in% names(str))) next
  
  # Extract allele vectors for this locus
  # Each vector contains one allele per individual
  a <- str[[colA]]
  b <- str[[colB]]
  
  # Treat allele value 0 as missing data
  # Convert these to NA so they are ignored in downstream calculations
  a[a == 0] <- NA
  b[b == 0] <- NA
  
  # Order alleles within individuals
  # This ensures genotypes like "150/152" and "152/150"
  # are treated as the same genotype
  a1 <- pmin(a, b, na.rm = FALSE)  # smaller allele
  a2 <- pmax(a, b, na.rm = FALSE)  # larger allele
  
  # Create genotype strings per individual
  # Format: "allele1/allele2"
  # If either allele is missing, genotype is set to NA
  genotype <- ifelse(
    is.na(a1) | is.na(a2),
    NA,
    paste0(a1, "/", a2)
  )
  
  # Determine all unique alleles observed at this locus
  # Uses ordered alleles (a1, a2) and removes missing values
  alleles <- sort(unique(na.omit(c(a1, a2))))
  n_alleles <- length(alleles)
  
  # Handle loci with no usable data
  if (n_alleles < 1) {
    # No alleles observed → cannot construct a genotype matrix
    hwe[[loc]] <- list(geno_matrix = NULL, hwe = NA)
    next
  }
  
  # Handle monomorphic loci (only one allele observed)
  # HWE test is undefined, but we still store the genotype count
  if (n_alleles == 1) {
    
    # Create a 1x1 genotype matrix for the single homozygous genotype
    mat <- matrix(
      0,
      nrow = 1,
      ncol = 1,
      dimnames = list(alleles, alleles)
    )
    
    # Count number of non-missing homozygous genotypes
    mat[1, 1] <- sum(
      !is.na(genotype) &
        genotype == paste0(alleles[1], "/", alleles[1])
    )
    
    # Store matrix and NA for HWE result
    hwe[[loc]] <- list(geno_matrix = mat, hwe = NA)
    
    # Record NA p-value in results table
    results <- rbind(
      results,
      data.frame(Locus = loc, p_value = NA, stringsAsFactors = FALSE)
    )
    
    next
  }
  
  # Count observed genotypes (excluding missing)
  # Result is a named table like:
  # "150/150" "150/152" "152/152"
  geno_counts <- table(na.omit(genotype))
  
  # Initialize a lower-triangular genotype count matrix
  # Rows and columns correspond to allele values
  # Only the lower triangle (i >= j) is filled
  mat <- matrix(
    0,
    nrow = n_alleles,
    ncol = n_alleles,
    dimnames = list(alleles, alleles)
  )
  
  # Fill genotype matrix from genotype count table
  for (g in names(geno_counts)) {
    
    # Split genotype string into its two alleles
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    
    # Identify row/column indices for each allele
    idx_i <- match(parts[1], as.character(alleles))
    idx_j <- match(parts[2], as.character(alleles))
    
    # Skip if allele matching fails (should rarely occur)
    if (is.na(idx_i) || is.na(idx_j)) next
    
    # Store count in lower triangle to avoid double counting
    if (idx_i >= idx_j) {
      mat[idx_i, idx_j] <- geno_counts[g]
    } else {
      mat[idx_j, idx_i] <- geno_counts[g]
    }
  }
  
  # Store genotype matrix and perform Hardy–Weinberg test
  # Uses permutation-based test with 10,000 permutations
  hwe[[loc]] <- list(geno_matrix = mat)
  
  test <- try(
    HWPerm.mult(mat, B = 10000),
    silent = TRUE
  )
  
  # Handle possible errors from HWPerm.mult
  if (inherits(test, "try-error")) {
    
    # Store NA if the test failed
    hwe[[loc]][["hwe"]] <- NA
    
    results <- rbind(
      results,
      data.frame(Locus = loc, p_value = NA, stringsAsFactors = FALSE)
    )
    
  } else {
    
    # Store full test result and p-value
    hwe[[loc]][["hwe"]] <- test
    
    results <- rbind(
      results,
      data.frame(Locus = loc, p_value = test$pval, stringsAsFactors = FALSE)
    )
  }
}


# view results
head(results)
results
subset(results, results$p_value < 0.05)
# L17 and L27 are not in HWE --> remove

str_hwe <- str %>% select(-c("L17_A","L17_B","L27_A","L27_B"))
str(str_hwe)

# create vector with unique loci names
loci <- unique(sub("_[A-Z]$", "", names(str_hwe)[-1]))
# calculate the average number of markers per individual
str_hwe <- str_hwe %>%
  rowwise() %>%
  mutate(non_missing_loci = sum(
    sapply(loci, function(locus) {
      get(paste0(locus, "_A")) != 0 & get(paste0(locus, "_B")) != 0
    })
  )) %>%
  ungroup()
str(str_hwe)
mean(str_hwe$non_missing_loci)
# average marker per individual: 37.42
# work with 38, 20, 10 markers


set.seed(123)  # for reproducibility
# create random subsample of 20 loci
loci_20 <- sample(loci, 20)
str_hwe_20 <- str_hwe[,c("animal_id", paste0(rep(loci_20, each = 2), rep(c("_A","_B"),20)))]
# count average number of loci per individual
str_hwe_20 <- str_hwe_20 %>%
  rowwise() %>%
  mutate(non_missing_loci = sum(
    sapply(loci_20, function(locus) {
      get(paste0(locus, "_A")) != 0 & get(paste0(locus, "_B")) != 0
    })
  )) %>%
  ungroup()
mean(str_hwe_20$non_missing_loci) # 18.7

set.seed(123)  # for reproducibility
# create random subsample of 10 loci
loci_10 <- sample(loci, 10)
str_hwe_10 <- str_hwe[,c("animal_id", paste0(rep(loci_10, each = 2), rep(c("_A","_B"),10)))]
# count average number of loci per individual
str_hwe_10 <- str_hwe_10 %>%
  rowwise() %>%
  mutate(non_missing_loci = sum(
    sapply(loci_10, function(locus) {
      get(paste0(locus, "_A")) != 0 & get(paste0(locus, "_B")) != 0
    })
  )) %>%
  ungroup()
mean(str_hwe_10$non_missing_loci) # 8.78


# write.table(str_hwe,"STR_n98_HWE_ALL.txt", row.names = F, quote = F)
# write.table(str_hwe_20,"STR_n98_HWE_20.txt", row.names = F, quote = F)
# write.table(str_hwe_10,"STR_n98_HWE_10.txt", row.names = F, quote = F)







#####
# CALCULATE RELATEDNESS
rm(list = ls())
str_all <- read.table("STR_n98_HWE_ALL.txt", header = T)
str(str_all)
# remove last column with information on number of markers per individual
str_all <- str_all[,-84]


# compute pairwise relatedness using different estimators
# WARNING: HIGH MEMORY USAGE, CONSIDER OUTSOURCING ON CLUSTER
library(related) # coancestry
str_rel_all <- coancestry(
  str_all,
  lynchli = TRUE,     # Lynch & Li (1995)
  lynchrd = TRUE,     # Lynch & Ritland (1999)
  quellergt = TRUE,   # Queller & Goodnight (1989)
  ritland = TRUE,     # Ritland (1996)
  wang = TRUE,        # Wang (2002)
  trioml = TRUE,      # Maximum likelihood estimator
  dyadml = TRUE,
  ci95.num.bootstrap = 0 # skip bootstrap to reduce memory usage
)
# save output as data frame
rstr_all <- data.frame(str_rel_all[["relatedness"]])

# load function to create sorted dyads
cs.dyad <- function(ind1, ind2){
  dy = paste(ind1, ind2, sep = "_")   
  dy = unlist(lapply(as.character(dy), function(x){paste(sort(unlist(strsplit(x, split = "_"))), collapse = "_")}))
  return(dy)
}
str(rstr_all)
# create dyad ID
rstr_all$dyad_sorted <- cs.dyad(ind1 = rstr_all$ind1.id, ind2 = rstr_all$ind2.id)
# remove unnecessary columns
rstr_all <- rstr_all[,c(5:12)]
# name columns
colnames(rstr_all) <- c(paste0("str41_",colnames(rstr_all)[1:7]),"dyad_sorted")

# merge to data frame with other info
all <- read.table("n98_IBD_PED_PIHAT_RAB_2025-10-13.txt", 
                  sep = " ", header = T)
str(all)
all <- merge(all, rstr_all, by = "dyad_sorted")

# write.table(all, "n98_IBD_PED_PIHAT_RAB_STR_2025-11-04.txt", row.names = F, quote = F)







#####
# CALCULATE RELATEDNESS for n(STR) = 20
rm(list = ls())
str_20 <- read.table("STR_n98_HWE_20.txt", header = T)
str(str_20)
# remove column with info on number of available loci per individual
str_20 <- str_20[,-42]
# calculate relatedness with estimator that performed best with n(STR) = 41 
library(related)
str_rel_20 <- coancestry(
  str_20,      
  dyadml = TRUE,
  ci95.num.bootstrap = 0 # skip bootstrap for reduced memory usage
)
# save results as data frame
rstr_20 <- data.frame(str_rel_20[["relatedness"]])

# load function to create sorted dyads
cs.dyad <- function(ind1, ind2){
  dy = paste(ind1, ind2, sep = "_")   
  dy = unlist(lapply(as.character(dy), function(x){paste(sort(unlist(strsplit(x, split = "_"))), collapse = "_")}))
  return(dy)
}

str(rstr_20)
# create dyad ID
rstr_20$dyad_sorted <- cs.dyad(ind1 = rstr_20$ind1.id, ind2 = rstr_20$ind2.id)
# remove unnecessary information
rstr_20 <- rstr_20[,c(11:12)]
colnames(rstr_20) <- c("str20_dyadml","dyad_sorted")

# merge to data frame with other info
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)
all <- merge(all, rstr_20, by = "dyad_sorted")

# write.table(all, "n98_IBD_PED_PIHAT_RAB_STR_2025-11-04.txt", row.names = F, quote = F)





#####
# CALCULATE RELATEDNESS for n(STR) = 10
# rm(list = ls())
str_10 <- read.table("STR_n98_HWE_10.txt", header = T)
str(str_10)
# remove info on available markers per individual
str_10 <- str_10[,-22]
# calculate relatedness with estimator that performed best for n(STR) = 41
library(related)
str_rel_10 <- coancestry(
  str_10,
  dyadml = TRUE,
  ci95.num.bootstrap = 0 # skip bootstrap to reduce memory usage
)
# svae result as data frame
rstr_10 <- data.frame(str_rel_10[["relatedness"]])

# load function to create sorted dyads
cs.dyad <- function(ind1, ind2){
  dy = paste(ind1, ind2, sep = "_")   
  dy = unlist(lapply(as.character(dy), function(x){paste(sort(unlist(strsplit(x, split = "_"))), collapse = "_")}))
  return(dy)
}

str(rstr_10)
# create dyad ID
rstr_10$dyad_sorted <- cs.dyad(ind1 = rstr_10$ind1.id, ind2 = rstr_10$ind2.id)
# remove unnecessary columns
rstr_10 <- rstr_10[,c(11:12)]
colnames(rstr_10) <- c("str10_dyadml","dyad_sorted")

# merge to data frame with other info
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)
all <- merge(all, rstr_10, by = "dyad_sorted")

# write.table(all, "n98_IBD_PED_PIHAT_RAB_STR_2025-11-04.txt", row.names = F, quote = F)
