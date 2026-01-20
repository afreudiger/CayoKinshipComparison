# Compare the performance of different relatedness estimators

rm(list = ls())

# load function to create sorted dyads
cs.dyad <- function(ind1, ind2){
  dy = paste(ind1, ind2, sep = "_")   
  dy = unlist(lapply(as.character(dy), function(x){paste(sort(unlist(strsplit(x, split = "_"))), collapse = "_")}))
  return(dy)
}

# load ibd and ped relatedness
all <- read.csv("IBD_PED_wo-maxGP99-het-sampleswap_2025-09-17.csv", 
                header = T, sep = ",")
str(all)

# read in sequencing info of animal ID
# this is needed because SNP and PMR relatedness use different ID-codes for animals than STR, IBD and PED
seq_inf <- read.csv("IndividualsSequenced_2024-11-21.csv", 
                    header = T, sep = ";")
str(seq_inf)
# only keep the two ID-codes
seq_red <- seq_inf[,c("animal_id","id_ccg")]



#####
# load pi-hat data
library(data.table) # fread
library(stringr) # str_match
filenames <- list.files("./RelatednessEstimation", 
                        pattern = "*pi_hat*", full.names = TRUE)
# read all files in and save them in a list
pihat <- lapply(filenames, fread)
# rename the list according to the titles of the files
names(pihat) <- str_match(filenames, "RelatednessEstimation/*(.*?)_pi_hat.genome")[,2]


# extract ID from file name
library(stringr) # str_remove
# remove extensions from id columns --> only id should remain
for(i in names(pihat)){
  pihat[[i]] <- data.frame(pihat[[i]])
  for(j in 1:2){
    # remove file extensions to extract ID
    pihat[[i]][,paste0("IID",j)] <- str_remove(pihat[[i]][,paste0("IID",j)], pattern = ".RG.bam")
    # change col names of seq_red to match pihat
    colnames(seq_red) <- c(paste0("animal_id",j), paste0("IID",j))
    # merge animal ID to pihat
    pihat[[i]] <- merge(pihat[[i]], seq_red, by = paste0("IID",j), all.x = T, all.y = F)
    for(y in 1:nrow(pihat[[i]])){
      # if animal_id is NA
      if(is.na(pihat[[i]][y,paste0("animal_id",j)])){
        # replace with IID
        pihat[[i]][y,paste0("animal_id",j)] <- pihat[[i]][y,paste0("IID",j)]
      }
    }
  }
  # check if it worked well
  print(unique(c(pihat[[i]][,"IID1"], pihat[[i]][,"IID2"])))
  # create dyad ID
  pihat[[i]][,"dyad_sorted"] <- cs.dyad(ind1 = pihat[[i]][,"animal_id1"], ind2 = pihat[[i]][,"animal_id2"])
  # change column name of PI_HAT
  names(pihat[[i]])[names(pihat[[i]]) == "PI_HAT"] <- paste0("pihat_",i)
  all <- merge(all, pihat[[i]][,c("dyad_sorted",paste0("pihat_",i))], by = "dyad_sorted", all.x = F, all.y = T)
}
str(all)



#####
# load ngsrelate data
library(data.table) # fread
library(stringr) # str_match
filenames <- list.files("./RelatednessEstimation", 
                        pattern = "*ngsrelate.res", full.names = TRUE)
# read all files in and save them in a list
ngsrel <- lapply(filenames, fread)
# rename the list according to the titles of the files
names(ngsrel) <- str_match(filenames, "RelatednessEstimation/*(.*?)_ngsrelate.res")[,2]

# ngsrelate saves IDs as numbers --> translate into real IDs
ids <- read.table("ngsrelate_sample_ids.txt", 
                  header = F)
colnames(ids) <- "IID"
ids$a <- c(0:97)
str(ids)

# extract ID from file name
library(stringr) # str_remove
# remove extensions from id columns --> only id should remain
for(i in names(ngsrel)){
  ngsrel[[i]] <- data.frame(ngsrel[[i]])
  for(j in 1:2){
    # translate number into ID
    if(j == 1){
      colnames(ids) <- c(paste0("IID",j),"a")
      ngsrel[[i]] <- merge(ngsrel[[i]], ids, by = "a", all.x = T, all.y = F)
    } else if(j == 2){
      colnames(ids) <- c(paste0("IID",j),"b")
      ngsrel[[i]] <- merge(ngsrel[[i]], ids, by = "b", all.x = T, all.y = F)
    }
    # change col names of seq_red to match ngsrel
    colnames(seq_red) <- c(paste0("animal_id",j), paste0("IID",j))
    # merge animal ID to ngsrel
    ngsrel[[i]] <- merge(ngsrel[[i]], seq_red, by = paste0("IID",j), all.x = T, all.y = F)
    for(y in 1:nrow(ngsrel[[i]])){
      # if animal_id is NA
      if(is.na(ngsrel[[i]][y,paste0("animal_id",j)])){
        # replace with IID
        ngsrel[[i]][y,paste0("animal_id",j)] <- ngsrel[[i]][y,paste0("IID",j)]
      }
    }
  }
  # check if it worked well
  print(unique(c(ngsrel[[i]][,"IID1"], ngsrel[[i]][,"IID2"])))
  # create dyad ID
  ngsrel[[i]][,"dyad_sorted"] <- cs.dyad(ind1 = ngsrel[[i]][,"animal_id1"], ind2 = ngsrel[[i]][,"animal_id2"])
  # change column name of rab and 2of3_IBD
  names(ngsrel[[i]])[names(ngsrel[[i]]) == "X2of3_IDB"] <- paste0("twothree_",i)
  names(ngsrel[[i]])[names(ngsrel[[i]]) == "rab"] <- paste0("rab_",i)
  all <- merge(all, ngsrel[[i]][,c("dyad_sorted",paste0("rab_",i),paste0("twothree_",i))], by = "dyad_sorted", all.x = F, all.y = T)
}

str(all)



#####
# load PMR data
pmr <- read.table("pmr94.ds_r_peakPMRs.tsv", header = T, sep = "\t")
str(pmr)

# extract ID from file name
library(stringr) # str_remove
# remove extensions from id columns --> only id should remain
for(i in c("iid1", "iid2")){
  for(j in c(".RG.bam", ".bam", ".speedseq.mmul10.bam")){
    pmr[,i] <- str_remove(pmr[,i], pattern = j)
  }
}
# check if it worked
unique(c(pmr$iid1, pmr$iid2))

# add animal ID of individuals
colnames(seq_red) <- c("animal_id1", "iid1")
pmr <- merge(pmr, seq_red, by = "iid1", all.x = T, all.y = F)
colnames(seq_red) <- c("animal_id2", "iid2")
pmr <- merge(pmr, seq_red, by = "iid2", all.x = T, all.y = F)
str(pmr)


# some individuals have their real ID instead of ccg_id
# currently NA in data frame
# copy their iid to animal_id
for(i in 1:nrow(pmr)){
  if(is.na(pmr[i,"animal_id1"])){
    pmr[i,"animal_id1"] <- pmr[i,"iid1"]
  }
  if(is.na(pmr[i,"animal_id2"])){
    pmr[i,"animal_id2"] <- pmr[i,"iid2"]
  }
}
str(pmr)


# create dyad IDs
pmr$dyad_sorted <- cs.dyad(ind1 = pmr[,"animal_id1"], ind2 = pmr[,"animal_id2"])
str(pmr)

# create separate columns for PMR based on different coverages
for(i in unique(pmr$cov_ds)){
  temp <- subset(pmr[,c("dyad_sorted","r")], pmr$cov_ds == i)
  colnames(temp)[2] <- paste0("pmr_",i,"x")
  all <- merge(all, temp, by = "dyad_sorted", all.x = T, all.y = T)
}



# write.table(all, "n98_IBD_PED_PIHAT_RAB_STR_PMR_2025-11-04.txt", row.names = F, quote = F)




#####
# comparison STR estimators
rm(list = ls())
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_PMR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)

# plot accuracy and precision for STR-based relatedness
str <- with(all, data.frame(dyad_sorted = dyad_sorted,
                            ibd = prop_genome_shared_8cm,
                            kinlabel_plot = kinlabel_plot))
# extract relevant columns
str <- merge(str, all[,c(1,52:58)], by = "dyad_sorted")
str(str)
# calculate difference between IBD and each estimator
for(i in colnames(all[,c(52:58)])){
  str[,paste0(i,"_ibd")] <- str[,i] - str[,"ibd"]
}
str(str)
# create data frame for plotting
str_plot <- data.frame(dyad_sorted = str$dyad_sorted,
                       str = str$str41_trioml_ibd,
                       estimator = as.character("trioml"),
                       x_axis = 1)
library(stringr) # str_match
# loop through each estimator and add the data
for(i in c("str41_wang_ibd","str41_lynchli_ibd","str41_lynchrd_ibd","str41_ritland_ibd","str41_quellergt_ibd","str41_dyadml_ibd")){
  str_plot <- rbind(str_plot, data.frame(dyad_sorted = str$dyad_sorted,
                                         str = str[,i],
                                         estimator = as.character(str_match(i, "str41_*(.*?)_ibd")[,2]),
                                         x_axis = 1))
}

# order them in plot according to publication date
str_plot$x_axis <- ifelse(str_plot$estimator == "quellergt",1,
                          ifelse(str_plot$estimator == "lynchli",2,
                                 ifelse(str_plot$estimator == "ritland",3,
                                        ifelse(str_plot$estimator == "lynchrd",4,
                                               ifelse(str_plot$estimator == "wang",5,
                                                      ifelse(str_plot$estimator == "dyadml",6,7))))))

# create colours for plotting
vis_data <- data.frame(col = sapply(seq(0.8, 0.2, length.out = 7), function(g) rgb(red = 0, green = g, blue = 0)))
# plot it
pdf(paste0("Plot_IBD-STR_DiffEstimators_2025-11-04.pdf"), 
    width = 8, height = 5, useDingbats=FALSE)
par(mar = c(5,6,1,1))
# set up plot
plot(str_plot$str ~ str_plot$x_axis, xlim = c(0.7,7.3), 
     ylim = c(-0.6,0.4),
     type = "n", cex.axis = 1.4, xlab = "", 
     ylab = "", xaxt = "n", 
     yaxt = "n"
)
# add horizontal line for reference
abline(h = 0, lwd = 1, lty = 2, col = "black")

# loop through each estimator
for (i in c(1:7)) {
  dat <- subset(str_plot, x_axis == i)$str
  n <- length(dat)
  # Define whisker limits (like boxplot: Q1 - 1.5*IQR, Q3 + 1.5*IQR)
  q1 <- quantile(dat, 0.25, na.rm = TRUE)
  q3 <- quantile(dat, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- max(min(dat), q1 - 1.5 * iqr)
  upper <- min(max(dat), q3 + 1.5 * iqr)
  
  # Density only within whiskers
  d <- density(dat)
  keep <- d$x >= lower & d$x <= upper
  d$x <- d$x[keep]
  d$y <- d$y[keep]
  
  # Scale density width
  d$y <- d$y / max(d$y) * 0.3
  
  x_pos <- i
  
  # Draw violin (trimmed to whisker range)
  polygon(c(x_pos - d$y, rev(x_pos + d$y)),
          c(d$x, rev(d$x)),
          col = adjustcolor(vis_data[i, "col"], alpha.f = 0.5),
          border = vis_data[i, "col"])
  
  # Median line
  med <- median(dat, na.rm = TRUE)
  # find density value that's closest to position of median
  y_med <- d$y[which.min(abs(d$x - med))]
  # draw line as wide as violin is at this position
  # adjust by 0.007 to make sure if aligns well with the violin
  lines(c(x_pos - y_med + 0.007, x_pos + y_med - 0.007),
        c(med, med),
        col = vis_data[i, "col"], lwd = 2)
  
  # Plot outliers as points
  outliers <- dat[dat < lower | dat > upper]
  if (length(outliers) > 0) {
    x_jitter <- x_pos + runif(length(outliers), -0.02, 0.02)  # fixed jitter range
    points(x_jitter, outliers,
           pch = 16, col = adjustcolor(vis_data[i, "col"], alpha.f = 0.7), cex = 0.5)
  }
}

# Axes + labels
axis(side = 1, at = c(1:7), label = c("quellergt","lynchli","ritland","lynchrd","wang","dyadml","trioml"), 
     las = 1, cex.axis = 1.2)
axis(side = 2, c(-0.6,-0.4,-0.2,0,0.2,0.4), las = 2, cex.axis = 1.2)
mtext(side = 1, text = "estimator", line = 3.5, cex = 1.7)
mtext(side = 2, text = expression("difference to r"["IBD"]), line = 3.5, cex = 1.7)

dev.off()

# compare the precision and accuracy for the two best performing estimators
median(str$str41_dyadml_ibd); median(str$str41_trioml_ibd)
range(str$str41_dyadml_ibd); range(str$str41_trioml_ibd)

# compute the mean for each estimator
for(i in c("quellergt","lynchli","ritland","lynchrd","wang","dyadml","trioml")){
  print(paste(i, mean(subset(str_plot$str, str_plot$estimator == i)), sep = ": "))
}
# --> DyadML it is

# NOTE: To create Figure S3, this plot and the SNP plot were merged in INKSCAPE







#####
# comparison SNP estimators
rm(list = ls())
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_PMR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)


# plot accuracy and precision for SNP-based relatedness
snp <- with(all, data.frame(dyad_sorted = dyad_sorted,
                            ibd = prop_genome_shared_8cm,
                            kinlabel_plot = kinlabel_plot))
# extract relevant information
snp <- merge(snp, all[,c(1,34:51)], by = "dyad_sorted")

# calculate difference between IBD and each estimator
for(i in colnames(all[,c(34:51)])){
  snp[,paste0(i,"_ibd")] <- snp[,i] - snp[,"ibd"]
}
str(snp)

# set up data frame for plotting
snp_plot <- with(snp, data.frame(dyad_sorted = dyad_sorted,
                                 snp = pihat_filtered_vcf_ibd,
                                 estimator = as.character("pihat"),
                                 x_axis = 1))
snp_plot <- rbind(snp_plot, with(snp, data.frame(dyad_sorted = dyad_sorted,
                                                 snp = twothree_filtered_vcf_ibd,
                                                 estimator = as.character("twothree"),
                                                 x_axis = 2)))
snp_plot <- rbind(snp_plot, with(snp, data.frame(dyad_sorted = dyad_sorted,
                                                 snp = rab_filtered_vcf_ibd,
                                                 estimator = as.character("rab"),
                                                 x_axis = 3)))

# choose colours for plotting
vis_data <- data.frame(col = c("deepskyblue2","deepskyblue3","deepskyblue4"))
# plot it
pdf(paste0("Plot_IBD-SNP_DiffEstimators_2025-10-14.pdf"), 
    width = 4, height = 5) #width = 7, height = 7
par(mar = c(5,6,1,1))
# set up plot
plot(snp_plot$snp ~ snp_plot$x_axis, xlim = c(0.7,3.3), 
     ylim = c(-0.6,0.4),
     type = "n", cex.axis = 1.4, xlab = "", 
     ylab = "", xaxt = "n", 
     yaxt = "n"
)
# horizontal line for reference
abline(h = 0, lwd = 1, lty = 2, col = "black")
# loop through each estimator
for (i in c(1:3)) {
  dat <- subset(snp_plot, x_axis == i)$snp
  n <- length(dat)
  # Define whisker limits (like boxplot: Q1 - 1.5*IQR, Q3 + 1.5*IQR)
  q1 <- quantile(dat, 0.25, na.rm = TRUE)
  q3 <- quantile(dat, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- max(min(dat), q1 - 1.5 * iqr)
  upper <- min(max(dat), q3 + 1.5 * iqr)
  
  # Density only within whiskers
  d <- density(dat)
  keep <- d$x >= lower & d$x <= upper
  d$x <- d$x[keep]
  d$y <- d$y[keep]
  
  # Scale density width
  d$y <- d$y / max(d$y) * 0.3
  
  x_pos <- i
  
  # Draw violin (trimmed to whisker range)
  polygon(c(x_pos - d$y, rev(x_pos + d$y)),
          c(d$x, rev(d$x)),
          col = adjustcolor(vis_data[i, "col"], alpha.f = 0.5),
          border = vis_data[i, "col"])
  
  # Median line
  med <- median(dat, na.rm = TRUE)
  # find density value that's closest to position of median
  y_med <- d$y[which.min(abs(d$x - med))]
  # draw line as wide as violin is at this position
  # adjust by 0.007 to make sure if aligns well with the violin
  lines(c(x_pos - y_med + 0.007, x_pos + y_med - 0.007),
        c(med, med),
        col = vis_data[i, "col"], lwd = 2)
  
  # Plot outliers as points
  outliers <- dat[dat < lower | dat > upper]
  if (length(outliers) > 0) {
    x_jitter <- x_pos + runif(length(outliers), -0.02, 0.02)  # fixed jitter range
    points(x_jitter, outliers,
           pch = 16, col = adjustcolor(vis_data[i, "col"], alpha.f = 0.7), cex = 0.5)
  }
}

# Axes + labels
axis(side = 1, at = c(1,2,3), label = c("pi-hat",expression("r"["ab"]),"2-out-of-3"), 
     las = 1, cex.axis = 1.2)
axis(side = 2, c(-0.6,-0.4,-0.2,0,0.2,0.4), las = 2, cex.axis = 1.2)
mtext(side = 1, text = "estimator", line = 3.5, cex = 1.7)
mtext(side = 2, text = expression("difference to r"["IBD"]), line = 3.5, cex = 1.7)

dev.off()


# NOTE: To create Figure S3, this plot and the SNP plot were merged in INKSCAPE










#####
# Comparison datasets scatter plot
rm(list = ls())
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_PMR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)

# scatter plot STR/SNP/PMR ~ IBD
# plot them individually and merge them in INKSCAPE later
pdf(paste0("Plot_IBD-PIHAT_DYADML_PMR_Scatter_2025-11-04.pdf"), 
    width = 6, height = 6, useDingbats=FALSE)
# STR
par(mar = c(6,6,1,1))
# set up plot
plot(all$str41_dyadml ~ all$prop_genome_shared_8cm, xlim = c(0,0.6), 
     ylim = c(-0.1,0.7),
     cex.axis = 1.4, xlab = "", 
     ylab = "", xaxt = "n", 
     yaxt = "n", type = "n",
     col = adjustcolor(rgb(0,0.5,0), alpha = 0.4), pch = 19
)
# perfect correlation between STR and IBD for reference
abline(a = 0, b = 1, lwd = 1, lty = 2, col = "black")
# plot semi transparent points
points(all$str41_dyadml ~ all$prop_genome_shared_8cm,
       col = adjustcolor(rgb(0,0.5,0), alpha = 0.4), pch = 19)
# Axes + labels
axis(side = 1, c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), 
     las = 1, cex.axis = 1.2)
axis(side = 2, c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), las = 2, cex.axis = 1.2)
mtext(side = 1, text = expression("r"["IBD"]), line = 3.5, cex = 1.7)
mtext(side = 2, text = expression("r"["STR"]), line = 3.5, cex = 1.7)


# SNP
par(mar = c(6,6,1,1))
# set up plot
plot(all$pihat_filtered_vcf ~ all$prop_genome_shared_8cm, xlim = c(0,0.6), 
     ylim = c(-0.1,0.7),
     cex.axis = 1.4, xlab = "", 
     ylab = "", xaxt = "n", 
     yaxt = "n", type = "n",
     col = adjustcolor(rgb(0,0.5617647,0.75), alpha = 0.4), pch = 19
)
# perfect correlation between SNP and IBD for reference
abline(a = 0, b = 1, lwd = 1, lty = 2, col = "black")
# plot semi transparent points
points(all$pihat_filtered_vcf ~ all$prop_genome_shared_8cm,
       col = adjustcolor(rgb(0,0.5617647,0.75), alpha = 0.4), pch = 19)
# Axes + labels
axis(side = 1, c(0,0.1,0.2,0.3,0.4,0.5,0.6), 
     las = 1, cex.axis = 1.2)
axis(side = 2, c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), las = 2, cex.axis = 1.2)
mtext(side = 1, text = expression("r"["IBD"]), line = 3.5, cex = 1.7)
mtext(side = 2, text = expression("r"["SNP"]), line = 3.5, cex = 1.7)


# PMR
par(mar = c(6,6,1,1))
# set up plot
plot(all$pihat_filtered_vcf ~ all$prop_genome_shared_8cm, xlim = c(0,0.6), 
     ylim = c(-0.1,0.7),
     cex.axis = 1.4, xlab = "", 
     ylab = "", xaxt = "n", 
     yaxt = "n", type = "n",
     col = adjustcolor("#D59A0D", alpha = 0.4), pch = 19
)
# perfect correlation between STR and IBD for reference
abline(a = 0, b = 1, lwd = 1, lty = 2, col = "black")
# plot semi-transparent points
points(all$pmr_1x ~ all$prop_genome_shared_8cm,
       col = adjustcolor("#D59A0D", alpha = 0.4), pch = 19)
# Axes + labels
axis(side = 1, c(0,0.1,0.2,0.3,0.4,0.5,0.6), 
     las = 1, cex.axis = 1.2)
axis(side = 2, c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), las = 2, cex.axis = 1.2)
mtext(side = 1, text = expression("r"["IBD"]), line = 3.5, cex = 1.7)
mtext(side = 2, text = expression("r"["PMR"]), line = 3.5, cex = 1.7)

dev.off()







#####
# Comparison datasets violin plot
rm(list = ls())
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_PMR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)
# create data frame with info of interest for plotting
all_dataset <- with(all, data.frame(dyad_sorted = dyad_sorted,
                                    ibd = prop_genome_shared_8cm,
                                    kinlabel_plot = kinlabel_plot))
all_dataset <- merge(all_dataset, all[,c("dyad_sorted","str10_dyadml","str20_dyadml","str41_dyadml",
                                         "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
                                         "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")], by = "dyad_sorted")
str(all_dataset)
# calculate difference between IBD and each estimator
for(i in colnames(all_dataset[,c("str10_dyadml","str20_dyadml","str41_dyadml",
                                 "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
                                 "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")])){
  all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - all_dataset[,"ibd"]
}
str(all_dataset)

# Calculate R² for each estimator
r_squared <- numeric()
for(i in c("str10_dyadml","str20_dyadml","str41_dyadml",
           "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
           "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")){
  temp <- subset(all_dataset[,c("ibd",i)], !is.na(all_dataset[,i]))
  r_squared[i] <- cor(temp[, i], 
                      temp[, "ibd"])^2
}

# create data frame for plotting
all_plot <- data.frame(dyad_sorted = all_dataset$dyad_sorted,
                       r_dev = all_dataset$str10_dyadml_ibd,
                       dataset = as.character("str10_dyadml"),
                       x_axis = 1,
                       kinlabel_plot = all_dataset$kinlabel_plot)
for(i in c("str20_dyadml_ibd","str41_dyadml_ibd",
           "pihat_exons_1_ibd","pihat_exons_10_ibd","pihat_exons_100_ibd","pihat_filtered_vcf_ibd",
           "pmr_0.01x_ibd","pmr_0.05x_ibd","pmr_0.1x_ibd","pmr_1x_ibd")){
  temp <- data.frame(dyad_sorted = all_dataset$dyad_sorted,
                     r_dev = all_dataset[,i],
                     dataset = as.character(i),
                     x_axis = (all_plot[nrow(all_plot), x_axis] + 1),
                     kinlabel_plot = all_dataset$kinlabel_plot)
  all_plot <- rbind(all_plot, temp)
}

# define colours for STR data
vis_data <- data.frame(col = sapply(seq(0.8, 0.2, length.out = 3), function(g) rgb(red = 0, green = g, blue = 0)))
# define and add colours for SNP data
base <- col2rgb("deepskyblue") / 255
# Create 6 shades darker → lighter
vis_data <- rbind(vis_data, data.frame(
  col = sapply(seq(1, 0.5, length.out = 4), function(fac) {
    rgb(red   = base[1] * fac,
        green = base[2] * fac,
        blue  = base[3] * fac)
  })
))
# define and add colours for PMR data
base <- col2rgb("darkgoldenrod1") / 255
vis_data <- rbind(vis_data, data.frame(
  col = sapply(seq(1, 0.5, length.out = 4), function(fac) {
    rgb(red   = base[1] * fac,
        green = base[2] * fac,
        blue  = base[3] * fac)
  })
))

# plot it
pdf(paste0("Plot_IBD-STR_SNP_PMR_Datasets_Violin_2025-11-04.pdf"), 
    width = 13, height = 6, useDingbats=FALSE) #width = 7, height = 7
par(mar = c(4,6,3,1))
# set up plot
plot(all_plot$r_dev ~ all_plot$x_axis, xlim = c(0.9,11.1), 
     ylim = c(-0.5,1),
     type = "n", cex.axis = 1.4, xlab = "", 
     ylab = "", xaxt = "n", 
     yaxt = "n"
)
# add horizontal line for reference
abline(h = 0, lwd = 1, lty = 2, col = "black")
# add vertical lines for visual separation of estimation methods
abline(v = 3.5, lwd = 1, lty = 1, col = "black")
abline(v = 7.5, lwd = 1, lty = 1, col = "black")

# loop through each estimator and create a violin following boxplot conventions
for (i in c(1:11)) {
  dat <- subset(all_plot, x_axis == i & !is.na(r_dev))$r_dev
  n <- length(dat)
  # Define whisker limits (like boxplot: Q1 - 1.5*IQR, Q3 + 1.5*IQR)
  q1 <- quantile(dat, 0.25, na.rm = TRUE)
  q3 <- quantile(dat, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- max(min(dat), q1 - 1.5 * iqr)
  upper <- min(max(dat), q3 + 1.5 * iqr)
  
  # Density only within whiskers
  d <- density(dat)
  keep <- d$x >= lower & d$x <= upper
  d$x <- d$x[keep]
  d$y <- d$y[keep]
  
  # Scale density width
  d$y <- d$y / max(d$y) * 0.4
  
  x_pos <- i
  
  # Draw violin (trimmed to whisker range)
  polygon(c(x_pos - d$y, rev(x_pos + d$y)),
          c(d$x, rev(d$x)),
          col = adjustcolor(vis_data[i, "col"], alpha.f = 0.5),
          border = vis_data[i, "col"])
  
  # Median line
  med <- median(dat, na.rm = TRUE)
  # find density value that's closest to position of median
  y_med <- d$y[which.min(abs(d$x - med))]
  # draw line as wide as violin is at this position
  # adjust by 0.007 to make sure if aligns well with the violin
  lines(c(x_pos - y_med + 0.007, x_pos + y_med - 0.007),
        c(med, med),
        col = vis_data[i, "col"], lwd = 2)
  
  # Plot outliers as points
  outliers <- dat[dat < lower | dat > upper]
  if (length(outliers) > 0) {
    x_jitter <- x_pos + runif(length(outliers), -0.02, 0.02)  # fixed jitter range
    points(x_jitter, outliers,
           pch = 16, col = adjustcolor(vis_data[i, "col"], alpha.f = 0.7), cex = 0.5)
  }
  
  # Add R² values below x-axis labels
  est <- c("str10_dyadml","str20_dyadml","str41_dyadml",
           "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
           "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")[i]
  text(i, par("usr")[3],
       labels = sprintf("R²=%.3f", r_squared[est]),
       pos = 1, xpd = TRUE, cex = 0.9, offset = 2.5)
}

# Axes + labels
axis(side = 1, at = c(1:11), label = c("9 STR","19 STR","37 STR",
                                       "1% exom","10% exom","100% exom","WGS",
                                       "0.01x","0.05x","0.1x","1x"), 
     las = 1, cex.axis = 1.1)
axis(side = 2, c(-0.5,-0.25,0,0.25,0.5,0.75,1), las = 2, cex.axis = 1.2)
mtext(side = 2, text = expression("difference to r"["IBD"]), line = 3.5, cex = 1.7)

dev.off()








#####
# Comparison datasets violin plot for different degrees of freedom
rm(list = ls())
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_PMR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)
# extract relevant data
all_dataset <- with(all, data.frame(dyad_sorted = dyad_sorted,
                                    ibd = prop_genome_shared_8cm,
                                    kinlabel_plot = kinlabel_plot))
all_dataset <- merge(all_dataset, all[,c("dyad_sorted","str10_dyadml","str20_dyadml","str41_dyadml",
                                         "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
                                         "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")], by = "dyad_sorted")
str(all_dataset)


# calculate difference between IBD and each estimator
for(i in colnames(all_dataset[,c("str10_dyadml","str20_dyadml","str41_dyadml",
                                 "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
                                 "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")])){
  all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - all_dataset[,"ibd"]
}
str(all_dataset)
# assign each kin class to a degree of relatedness
all_dataset$kin_degree <- ifelse(all_dataset$kinlabel_plot == "parent-offspring","1st degree",
                                 ifelse(all_dataset$kinlabel_plot == "fullsib","1st degree",
                                        ifelse(all_dataset$kinlabel_plot == "matsib","2nd degree",
                                               ifelse(all_dataset$kinlabel_plot == "patsib","2nd degree",
                                                      ifelse(all_dataset$kinlabel_plot == "mat_grandparent-offspring","2nd degree",
                                                             ifelse(all_dataset$kinlabel_plot == "pat_grandparent-offspring","2nd degree",
                                                                    ifelse(all_dataset$kinlabel_plot == "aunt_uncle-niece_nephew","3rd degree",
                                                                           ifelse(all_dataset$kinlabel_plot == "cousins","4th degree",
                                                                                  ifelse(all_dataset$kinlabel_plot == "misc","misc","nonkin")))))))))

# create data frame for plotting
all_plot <- data.frame(dyad_sorted = all_dataset$dyad_sorted,
                       r_dev = all_dataset$str10_dyadml_ibd,
                       dataset = as.character("str10_dyadml"),
                       x_axis = 1,
                       kinlabel_plot = all_dataset$kinlabel_plot,
                       kin_degree = all_dataset$kin_degree)
for(i in c("str20_dyadml_ibd","str41_dyadml_ibd",
           "pihat_exons_1_ibd","pihat_exons_10_ibd","pihat_exons_100_ibd","pihat_filtered_vcf_ibd",
           "pmr_0.01x_ibd","pmr_0.05x_ibd","pmr_0.1x_ibd","pmr_1x_ibd")){
  temp <- data.frame(dyad_sorted = all_dataset$dyad_sorted,
                     r_dev = all_dataset[,i],
                     dataset = as.character(i),
                     x_axis = (all_plot[nrow(all_plot),"x_axis"] + 1),
                     kinlabel_plot = all_dataset$kinlabel_plot,
                     kin_degree = all_dataset$kin_degree)
  all_plot <- rbind(all_plot, temp)
}


# Calculate R² per estimator and per kinship degree
estimator_cols <- c("str10_dyadml","str20_dyadml","str41_dyadml",
                    "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
                    "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")
r_squared_by_degree <- data.frame()
for (deg in unique(all_plot$kin_degree)) {
  for (est in estimator_cols) {
    dat_sub <- subset(all_dataset[,c("ibd",est)], all_dataset$kin_degree == deg & !is.na(all_dataset[,est]))
    r2_val <- cor(dat_sub[, est], dat_sub[, "ibd"])^2
    r_squared_by_degree <- rbind(r_squared_by_degree,
                                 data.frame(kin_degree = deg,
                                            estimator = est,
                                            R2 = r2_val))
  }
}

# define colours for STR
vis_data <- data.frame(col = sapply(seq(0.8, 0.2, length.out = 3), function(g) rgb(red = 0, green = g, blue = 0)))
# define and add colours for SNP
base <- col2rgb("deepskyblue") / 255
# Create 6 shades darker → lighter
vis_data <- rbind(vis_data, data.frame(
  col = sapply(seq(1, 0.5, length.out = 4), function(fac) {
    rgb(red   = base[1] * fac,
        green = base[2] * fac,
        blue  = base[3] * fac)
  })
))
# define and add colours for PMR
base <- col2rgb("darkgoldenrod1") / 255
vis_data <- rbind(vis_data, data.frame(
  col = sapply(seq(1, 0.5, length.out = 4), function(fac) {
    rgb(red   = base[1] * fac,
        green = base[2] * fac,
        blue  = base[3] * fac)
  })
))

# plot it
pdf(paste0("Plot_IBD-STR_SNP_PMR_Datasets_KinClasses_2025-11-04.pdf"), 
    width = 13, height = 6, useDingbats=FALSE)
# loop through each degree of kinship and create one plot
for(y in c("1st degree","2nd degree", "3rd degree", "4th degree","nonkin","misc")){
  par(mar = c(6,6,3,1))
  # set up plot
  plot(all_plot$r_dev ~ all_plot$x_axis, xlim = c(0.9,11.1), 
       ylim = c(-0.5,1),
       type = "n", cex.axis = 1.4, xlab = "", 
       ylab = "", xaxt = "n", 
       yaxt = "n"
  )
  # add horizontal line for visual reference
  abline(h = 0, lwd = 1, lty = 2, col = "black")
  # add vertical lines to separate estimators
  abline(v = 3.5, lwd = 1, lty = 1, col = "black")
  abline(v = 7.5, lwd = 1, lty = 1, col = "black")
  
  # add violins for each estimator and dataset
  for (i in c(1:11)) {
    dat <- subset(all_plot, x_axis == i & kin_degree == y)$r_dev
    n <- length(dat)
    # Define whisker limits (like boxplot: Q1 - 1.5*IQR, Q3 + 1.5*IQR)
    q1 <- quantile(dat, 0.25, na.rm = TRUE)
    q3 <- quantile(dat, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower <- max(min(dat), q1 - 1.5 * iqr)
    upper <- min(max(dat), q3 + 1.5 * iqr)
    
    # Density only within whiskers
    d <- density(dat)
    keep <- d$x >= lower & d$x <= upper
    d$x <- d$x[keep]
    d$y <- d$y[keep]
    
    # Scale density width
    d$y <- d$y / max(d$y) * 0.4
    
    x_pos <- i
    
    # Draw violin (trimmed to whisker range)
    polygon(c(x_pos - d$y, rev(x_pos + d$y)),
            c(d$x, rev(d$x)),
            col = adjustcolor(vis_data[i, "col"], alpha.f = 0.5),
            border = vis_data[i, "col"])
    
    # Median line
    med <- median(dat, na.rm = TRUE)
    # find density value that's closest to position of median
    y_med <- d$y[which.min(abs(d$x - med))]
    # draw line as wide as violin is at this position
    # adjust by 0.007 to make sure if aligns well with the violin
    lines(c(x_pos - y_med + 0.007, x_pos + y_med - 0.007),
          c(med, med),
          col = vis_data[i, "col"], lwd = 2)
    
    # Plot outliers as points
    outliers <- dat[dat < lower | dat > upper]
    if (length(outliers) > 0) {
      x_jitter <- x_pos + runif(length(outliers), -0.02, 0.02)  # fixed jitter range
      points(x_jitter, outliers,
             pch = 16, col = adjustcolor(vis_data[i, "col"], alpha.f = 0.7), cex = 0.5)
    }
  }
  
  # Axes + labels
  axis(side = 1, at = c(1:11), label = c("9 STR","19 STR","37 STR",
                                         "1% exom","10% exom","100% exom","WGS",
                                         "0.01x","0.05x","0.1x","1x"), 
       las = 1, cex.axis = 1.1, line = 0)
  
  # Add R² values per kinship degree below x-axis labels
  for(i in 1:11){
    est_name <- estimator_cols[i]
    r2_val <- subset(r_squared_by_degree, estimator == est_name & kin_degree == y)$R2
    if (length(r2_val) == 0) r2_val <- NA
    text(i, par("usr")[3],
         labels = sprintf("R²=%.3f", r2_val),
         pos = 1, xpd = TRUE, cex = 0.9, offset = 2.5)
  }
  
  
  axis(side = 2, c(-0.5,-0.25,0,0.25,0.5,0.75,1), las = 2, cex.axis = 1.2)
  mtext(side = 2, text = expression("difference to r"["IBD"]), line = 3.5, cex = 1.7)
  # add title
  if(y %in% c("1st degree","2nd degree","3rd degree","4th degree")){
    suffix <- c("st","nd","rd","th")[match(y, c("1st degree","2nd degree","3rd degree","4th degree"))]
    number <- c("1","2","3","4")[match(y, c("1st degree","2nd degree","3rd degree","4th degree"))]
    mtext(side = 3, text = bquote(.(number)^.(suffix) ~ "degree" ~ "(n =" ~ .(length(dat)) * ")"), cex = 1.7)
  } else {
    mtext(side = 3, text = paste0(y," (n = ",length(dat),")"), cex = 1.7)
  }
}

dev.off()





#####
# calculate deviance of each estimator and dataset from IBD
rm(list = ls())
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_PMR_2025-11-04.txt", 
                  sep = " ", header = T)
str(all)
colnames(all)
# create suitable data frame
all_dataset <- with(all, data.frame(dyad_sorted = dyad_sorted,
                                    ibd = prop_genome_shared_8cm,
                                    kinlabel_plot = kinlabel_plot))
all_dataset <- merge(all_dataset, all[,c(1,34:67)], by = "dyad_sorted")
str(all_dataset)


# calculate difference between IBD and each estimator
for(i in colnames(all_dataset[,c(4:37)])){
  all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - all_dataset[,"ibd"]
}
str(all_dataset)

# print summary statistics for all
for(i in colnames(all_dataset[,c(38:71)])){
  print(i)
  print(summary(all_dataset[,i]))
  print("")
}
