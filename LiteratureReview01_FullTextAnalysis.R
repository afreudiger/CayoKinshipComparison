# Analysis of literature review data

# IMPORTANT: Some articles appear multiple times in data set because they used > 1 species or > 1 marker type.
#            Before you run any analysis, always create a unique subset of the data, only including variables
#            you are interested in. Always include publication ID.

# NOTE: This code only creates single-panel plots without labels. Merging single plots to multi-panel plots and adding labels was done using INKSCAPE.



#####
# phylogenetic class and order were added with the following code
# note that not all genera could be assigned to classes and orders by taxize
# the missing ones were searched in the internet and added manually 
library(taxize) # classification

# get unique genera to avoid redundant API calls
unique_genera <- unique(na.omit(data$genus))

# get classification for each genus (returns the full taxonomic hierarchy)
classifications <- classification(unique_genera, db = "ncbi")

# extract the class from each classification
get_class <- function(classification_result) {
  # check if result is NULL or not a data frame
  if (is.null(classification_result) || !is.data.frame(classification_result)) {
    return(NA)
  }
  # check if data frame is empty
  if (nrow(classification_result) == 0) {
    return(NA)
  }
  # look for the "class" rank in the classification
  class_row <- classification_result[classification_result$rank == "class", ]
  if (nrow(class_row) > 0) {
    return(class_row$name[1])
  } else {
    return(NA)
  }
}

# apply to all classifications
class_vector <- sapply(classifications, get_class)

# create a lookup data frame
class_lookup <- data.frame(
  genus = unique_genera,
  class = class_vector,
  stringsAsFactors = FALSE
)


# extract the order from classifications
get_order <- function(classification_result) {
  # check if input is NULL or not a data frame
  if (is.null(classification_result) || !is.data.frame(classification_result)) {
    return(NA)
  }
  # check if data frame is empty
  if (nrow(classification_result) == 0) {
    return(NA)
  }
  # filter for rows where rank equals "order"
  order_row <- classification_result[classification_result$rank == "order", ]
  # if order rank exists, return the first name; otherwise return NA
  if (nrow(order_row) > 0) {
    return(order_row$name[1])
  } else {
    return(NA)
  }
}
# apply to all
order_vector <- sapply(classifications, get_order)

# add order to class_lookup data frame
class_lookup$order <- order_vector




#####
# count publications per year
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)

# extract unique information
year_unique <- unique(data[,c("publication_id","year")])
# add frequency
year_unique$frequency <- as.numeric(1)
# count publications per year
pub_year <- data.frame(aggregate(year_unique$frequency, by = list(year_unique$year), FUN = sum))
colnames(pub_year) <- c("year","num_pub")

# plot total number publications
str(pub_year)
rownames(pub_year) <- pub_year$year

# predict number for 2025, as we only recorded 6 months here
pub_year["2025","num_pub"] <- pub_year["2025","num_pub"]*2

# plot it
pdf(paste0("Plot_NumbPublication-Year.pdf"), 
    width = 5, height = 5)
par(mar = c(4,4,1,1))
plot(num_pub ~ year, data = pub_year, type = "b", las = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     pch = 16, ylim = c(0,100))
points(num_pub ~ year, data = subset(pub_year, pub_year$year == "2025"),
       pch = 16, col = "grey70")
axis(1, at = c(2000,2005,2010,2015,2020,2025))
axis(2, at = c(20,40,60,80,100), las = 1)
mtext(side = 1, text = "year", line = 2.5, cex = 1.2)
mtext(side = 2, text = "articles published", line = 2.5, cex = 1.2)
legend("bottomleft", col = c("black","grey70"), legend = c("measured","predicted"), pch = 16, bty = "n")
box()
dev.off()



# count publications per year normalised by the total number of publications per year
# read in total counts per year
total <- read.csv("NumbPublications-Year.csv",
                  header = T, sep = ";")
# merge to data
pub_year <- merge(pub_year, total, by = "year")
# calculate percentage of relatedness studies for each year
pub_year$num_norm <- pub_year$num_pub / pub_year$pub_total
# multiply by 1e5 to improve plotting
pub_year$norm_plot <- pub_year$num_norm*1e5

# plot it
pdf(paste0("Plot_FractionTotalPublication-Year.pdf"), 
    width = 5, height = 5)
par(mar = c(4,4,1,1))
plot(norm_plot ~ year, data = pub_year, type = "b", las = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     pch = 16, ylim = c(0,10))
points(norm_plot ~ year, data = subset(pub_year, pub_year$year == "2025"),
       pch = 16, col = "grey70")
axis(1, at = c(2000,2005,2010,2015,2020,2025))
axis(2, at = c(0,2,4,6,8,10), las = 1)
mtext(side = 1, text = "year", line = 2.5, cex = 1.2)
mtext(side = 2, text = expression("fraction of total articles published [*"*10^{-5} * "]"),
  line = 2.2, cex = 1.2)
legend("bottomleft", col = c("black","grey70"), legend = c("measured","predicted"), pch = 16, bty = "n")
box()
dev.off()









#####
# marker type across years
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)
# create unique data
method <- data[,c("publication_id","year","seq_method")]

# pool RFLP, RAPD and minisatellites for plotting
method$seq_method <- replace(method$seq_method, 
                             method$seq_method %in% c("RFLP","minisatellite","RAPD"), 
                             "other")
unique(method$seq_method)
# create unique data with pooled seq method
method <- unique(method)
# check frequency of different methods
table(method$seq_method)

# define the order of seq methods for plotting (bottom to top)
method$seq_method <- factor(method$seq_method, 
                            levels = c("other","SNP","WGS","STR"))
# create count table for plotting
count_table <- table(method$year, method$seq_method)
percent_table <- prop.table(count_table, margin = 1) * 100
# transpose percent table
t_percent <- t(percent_table)
# get dimensions
n_methods <- nrow(t_percent)
n_years <- ncol(t_percent)
years <- as.numeric(colnames(t_percent))

# plot it!
pdf(paste0("Plot_Markers-Year.pdf"), 
    width = 5, height = 5)
par(mar = c(4,4,1,1))
plot(0, 0, type = "n", 
     xlim = c(min(years) - 0.1, max(years) + 0.1), 
     ylim = c(0.8, 99.2),
     xlab = "", 
     ylab = "",
     xaxt = "n", las = 1)
mtext(side = 1, text = "year", line = 2.5, cex = 1.2)
mtext(side = 2, text = "percentage [%]", line = 2.5, cex = 1.2)
axis(1, at = c(2000,2005,2010,2015,2020,2025))
colors <- c("#330000","#7F3136","#B8AC80","#D2D2CE")
# draw the stacked bars
bar_width <- 0.8
for (i in 1:n_years) {
  # initialise bottom position for stacking bars
  bottom <- 0
  # loop through each method to add segments to the stacked bar
  for (j in 1:n_methods) {
    # draw rectangle for this method's segment
    # x: centred on the year with specified bar width
    # y: from current bottom to bottom + percentage value
    rect(years[i] - bar_width/2, bottom, 
         years[i] + bar_width/2, bottom + t_percent[j, i],
         col = colors[j], border = NA)
    # update bottom position for next segment
    bottom <- bottom + t_percent[j, i]
  }
}

# add legend in new plot and manually merge using Inkscape
plot(0, 0, type = "n", 
     xlim = c(min(years) - 0.1, max(years) + 0.1), 
     ylim = c(0, 100),
     xlab = "", 
     ylab = "",
     xaxt = "n", las = 1)
legend("bottom", legend = rev(rownames(t_percent)), 
       col = rev(colors), pch = 15, pt.cex = 1.5, bty = "n")
dev.off()






#####
# number of STR markers across years
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)
# subset STR studies
str <- subset(data, data$seq_method == "STR")
# Convert no_markers to numeric
str$no_markers_num <- as.numeric(str$no_markers)

# create year index for statistical analysis
str$year_index <- str$year - 2000
# check dispersion parameter of model
summary(glm(no_markers_num ~ year_index, data = str, family = "poisson"))$dispersion
# run model and check output
summary(glm(no_markers_num ~ year_index, data = str, family = "poisson"))


# Get unique years
years <- sort(unique(str$year))
pdf(paste0("Plot_NoSTR-Year.pdf"), 
    width = 5, height = 5)
# Set up the plot
par(mar = c(4,4,1,1), mfrow = c(1,1))
plot(0, 0, type = "n",
     xlim = c(min(years) - 0.1, max(years) + 0.1),
     ylim = c(1, 100),        # must be > 0
     log = "y",               # <-- log scale
     las = 1,
     xlab = "",
     ylab = "",
     yaxt = "n")
axis(2, at = c(1, 2, 5, 10, 20, 50, 100), labels = c(1, 2, 5, 10, 20, 50, 100), las = 1)
mtext(side = 1, text = "year", line = 2.5, cex = 1.2)
mtext(side = 2, text = "number STR markers", line = 2.5, cex = 1.2)

# Create violin for each year
violin_width <- 0.3 # play around with violin width to find suitable value
# loop through years to create violins following boxplot conventions
for (i in 1:length(years)) {
  # get data for this year
  year_data <- str$no_markers_num[str$year == years[i] & !is.na(str$no_markers_num)]
  
  if (length(year_data) > 1) {
    # calculate boxplot statistics
    q1 <- quantile(year_data, 0.25)
    q3 <- quantile(year_data, 0.75)
    iqr <- q3 - q1
    lower_whisker <- max(min(year_data), q1 - 1.5 * iqr)
    upper_whisker <- min(max(year_data), q3 + 1.5 * iqr)
    
    # identify outliers
    outliers <- year_data[year_data < lower_whisker | year_data > upper_whisker]
    
    # calculate density
    dens <- density(year_data, adjust = 1)
    
    # trim density to whisker range
    keep <- dens$x >= lower_whisker & dens$x <= upper_whisker
    dens$x <- dens$x[keep]
    dens$y <- dens$y[keep]
    
    # scale density to fit within violin_width
    dens$y <- dens$y / max(dens$y) * violin_width
    
    # draw complete violin as single polygon (right side, then left side reversed)
    polygon(c(years[i] + dens$y, rev(years[i] - dens$y)),
            c(dens$x, rev(dens$x)),
            col = NA, border = "black")
    
    # add median line
    med <- median(year_data)
    segments(years[i] - violin_width, med, 
             years[i] + violin_width, med, 
             col = "#6E2023", lwd = 2)
    
    # add outliers as points
    if (length(outliers) > 0) {
      points(rep(years[i], length(outliers)), outliers, pch = 16, cex = 0.5, col = adjustcolor("black", alpha = 0.6))
    }
  }
}
dev.off()




#####
# number of SNP markers across years
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
# create subset of studies working with SNPs
snp <- subset(data, data$seq_method == "SNP")
# convert no_markers to numeric
snp$no_markers_num <- as.numeric(snp$no_markers)
# create year index for statistical testing
snp$year_index <- snp$year - 2009
# check dispersion factor
summary(glm(no_markers_num ~ year_index, data = snp, family = "poisson"))$dispersion
# run model and check results
summary(glm(no_markers_num ~ year_index, data = snp, family = "poisson"))



#####
# WGS coverage
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
# unique subset of studies using WGS
wgs <- subset(data, data$seq_method == "WGS")
# summary of coverages used
summary(wgs$WGS_coverage, na.rm = T)






#####
# ten most frequent journals
rm(list = ls())

data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",", na.strings = c("NA",""))
str(data)
# check how many journals are respresented in total
length(unique(data$journal))
# create unique data set
journal_unique <- unique(data[,c("publication_id","journal")])
# add count to each publication
journal_unique$frequency <- as.numeric(1)

# create new data frame with number of publications per journal
journal <- data.frame(aggregate(journal_unique$frequency, by = list(journal_unique$journal), FUN = sum))
colnames(journal) <- c("journal","num_pub")
# order data frame according to counts (highest to lowest)
journal <- journal[order(journal$num_pub, decreasing = T), ]
# calculate the percentage of publications per journal
journal$perc_pub <- (journal$num_pub/sum(journal$num_pub))*100
# create subset with 10 most frequent journals
journal <- journal[c(1:10),]
# create x-axis position index
journal$x_axt <- c(1:nrow(journal))
# check fraction of publications in top 10 journals
sum(journal$num_pub)/2034

# plot it!
pdf(paste0("Plot_Journals.pdf"), 
    width = 4.5, height = 5)
par(mar = c(4,1,1,1))
# create horizontal bar plot
bp <- barplot(
  journal$perc_pub,
  horiz = TRUE,
  col ="#AC8D61",
  border = NA,
  xlim = c(-0.5, 14.5),
  yaxt = "n"
)

text(
  x = 0.2,  # slightly shift text to the right
  y = bp,
  labels = c("Mol Ecol","Behav Ecol Sociobiol","Conserv Genet",
             "Ecol Evol","Plos One","Behav Ecol","J Hered",
             "Anim Behav","Proc R Soc B","Aquacult"),
  col = "black",
  cex = 0.8,
  adj = 0
)
box()
mtext(side = 1, text = "fraction publications [%]", line = 2.5, cex = 1.2)
dev.off()












#####
# phylogenetic classes
rm(list = ls())

data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",", na.strings = c("NA",""))
str(data)
# check number of phylogenetic classes
length(unique(data$class))
# create unique subset
taxon_unique <- unique(data[,c("publication_id","class")])
# add count
taxon_unique$frequency <- as.numeric(1)
# create new data frame with number of publications per phylo class
taxa <- data.frame(aggregate(taxon_unique$frequency, by = list(taxon_unique$class), FUN = sum))
colnames(taxa) <- c("class","num_pub")
# order data frame according to counts (highest to lowest)
taxa <- taxa[order(taxa$num_pub, decreasing = T), ]
# calculate the percentage of publications per class
taxa$perc_pub <- (taxa$num_pub/sum(taxa$num_pub))*100
# create subset with 10 most frequent journals
taxa <- taxa[c(1:10),]
# create x-axis position index
taxa$x_axt <- c(1:nrow(taxa))

748/2035 # percentage Mammalia

# plot it
pdf(paste0("Plot_PhyloClass.pdf"), 
    width = 5, height = 5)
par(mar = c(4,4,1,1))
# create horizontal barplot
bp <- barplot(
  taxa$num_pub,
  horiz = TRUE,
  col = "#AC8D61",
  border = NA,
  xlim = c(-22, 822),
  yaxt = "n"
)
text(
  x = 5,  # shift text slightly to the right
  y = bp,
  labels = taxa$class,
  col = "black",
  cex = 0.8,
  adj = 0
)
box()
mtext(side = 1, text = "number of publications", line = 2.5, cex = 1.2)
mtext(side = 2, text = "phylogenetic class", line = 2.5, cex = 1.2)
dev.off()

# check most common genus per class
rownames(taxa) <- taxa$class
for(i in rownames(taxa)){
  temp <- subset(data, data$class == i)
  print(i)
  print(which.max(table(temp$genus)))
}
# choose silhouette for each genus from https://www.phylopic.org/ and add in Inkscape




#####
# marker type per phylogenetic class
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)
# pool RFLP, RAPD and minisatellites
data$seq_method <- replace(data$seq_method, 
                           data$seq_method %in% c("RFLP","minisatellite","RAPD"), 
                             "other")
# create unique subset
method <- data[,c("publication_id","class","seq_method")]

# test association between class and seq method
tab <- table(method$class, method$seq_method)
# proportions within class
prop.table(tab, margin = 1)
# effect size
library(rstatix)
cramer_v(tab)


# define the order of seq method for plotting (bottom to top)
method$seq_method <- factor(method$seq_method, 
                            levels = c("other","SNP","WGS","STR"))

# select ten most common classes
top_class <- data.frame(class = c("Mammalia","Aves","Actinopterygii","Insecta","Reptilia","Lepidosauria",
                                  "Chondrichthyes","Amphibia","Malacostraca","Gastropoda"),
                        x_axis = 1:10)
# create subset of data with ten top classes
method <- merge(method, top_class, by = "class", all.x = F, all.y = T)

# create table with counts per method per class
count_table <- table(method$x_axis, method$seq_method)
# calculate proportions
percent_table <- prop.table(count_table, margin = 1) * 100
# transpose it
t_percent <- t(percent_table)
# get dimensions
n_methods <- nrow(t_percent)
n_taxa <- ncol(t_percent)
taxa <- as.numeric(colnames(t_percent))

# plot it
pdf(paste0("Plot_Markers-PhyloClass.pdf"), 
    width = 5, height = 5)
par(mar = c(4,4,1,1))
# set up the plot
plot(0, 0, type = "n", 
     xlim = c(min(taxa) - 0.25, max(taxa) + 0.25), 
     ylim = c(0.8, 99.2),
     xlab = "", 
     ylab = "",
     xaxt = "n", las = 1)
mtext(side = 1, text = "phylogenetic class", line = 2.5, cex = 1.2)
mtext(side = 2, text = "percentage [%]", line = 2.5, cex = 1.2)
# define colours
colors <- c("#330000","#7F3136","#B8AC80","grey91")

# draw stacked bars
bar_width <- 0.7 # play around with bar width to find a suitable value
for (i in 1:n_taxa) {
  # initialise bottom position for stacking bars
  bottom <- 0
  # draw rectangle for this method's segment
  # x: centred on the class with specified bar width
  # y: from current bottom to bottom + percentage value
  for (j in 1:n_methods) {
    rect(taxa[i] - bar_width/2, bottom, 
         taxa[i] + bar_width/2, bottom + t_percent[j, i],
         col = colors[j], border = NA)
    # update bottom position for next segment
    bottom <- bottom + t_percent[j, i]
  }
}

# add legend in new plot and manually merge using Inkscape
plot(0, 0, type = "n", 
     xlim = c(min(taxa) - 0.1, max(taxa) + 0.1), 
     ylim = c(0, 100),
     xlab = "", 
     ylab = "",
     xaxt = "n", las = 1)
legend("bottom", legend = rev(rownames(t_percent)), 
       col = rev(colors), pch = 15, pt.cex = 1.5, bty = "n")
dev.off()







#####
# sample type (invasive vs non-invasive)
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)

# classify if studies used invasive or non-invasive samples, or both types
# sample type is stored as comma separated string
data$sample_class <- sapply(data$sample_type, function(x) {
  # split and trim whitespace
  types <- trimws(unlist(strsplit(x, ",")))
  # indicate if study used invasive samples
  has_inv <- any(types %in% c("tissue", "blood", "saliva", "egg", "bones_teeth"))
  # indicate if study used non-invasive samples
  has_noninv <- any(types %in% c("fecal", "hair", "feather"))
  # if they used samples of both types --> both
  if (has_inv & has_noninv) {
    "both"
    # if they used only invasive samples --> inv
  } else if (has_inv) {
    "inv"
    # if they used non-invasive samples --> noninv
  } else if (has_noninv) {
    "noninv"
    # if sample type is unknown --> NA
  } else {
    NA
  }
})
# create unique dataset
sample_class <- unique(data[,c("publication_id","year","sample_class")])
# check how many studies used which sample type
table(sample_class$sample_class)
# claculate fraction of studies using at least partly non-invasive samples
(217+178)/2034

# check how many non-inv studies used SNP-sequencing
# create unique subset
noninv <- unique(subset(data[,c("publication_id","year","sample_class","seq_method")], sample_class %in% c("noninv","both")))
table(noninv$seq_method)
24/400 # fraction SNP
374/400 # fraction STR
# check WGS studies on non-inv samples in detail
View(subset(data, data$seq_method == "WGS" & data$sample_class %in% c("noninv","both")))





##### 
# type of relationship inferred
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)
# create unique subset
# publications that have several entries differing in the type of relationship inferred are merged 
# each type of relationship that's inferred at least once in the study is set == "yes"
data <- aggregate(
  data[, c("individual_identification",
           "parent_offspring",
           "full_sib",
           "half_sib",
           "other_degree_rel")],
  by = list(publication_id = data$publication_id),
  FUN = function(x) ifelse(any(x == "yes", na.rm = TRUE), "yes", "no")
)

# fraction parent-offspring
nrow(subset(data, parent_offspring == "yes"))/nrow(data)
# fraction full siblings
nrow(subset(data, full_sib == "yes"))/nrow(data)
# fraction half siblings
nrow(subset(data, half_sib == "yes"))/nrow(data)
# fraction other degree of relatedness
nrow(subset(data, other_degree_rel == "yes"))/nrow(data)
# fraction individual identification
nrow(subset(data, individual_identification == "yes"))/nrow(data)





#####
# sample size across years
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)
# create year index
data$year_index <- data$year - 2000
# check dispersion
summary(glm(sample_size ~ year_index, data = data, family = "poisson"))$dispersion
# run model and check results
summary(glm(sample_size ~ year_index, data = data, family = "poisson"))

# get unique years
years <- sort(unique(data$year))
pdf(paste0("Plot_SampleSize.pdf"), 
    width = 5, height = 5)
# set up the plot
par(mar = c(4,4,1,1), mfrow = c(1,1))
plot(0, 0, type = "n",
     xlim = c(min(years) - 0.1, max(years) + 0.1),
     ylim = c(1, 50000),        
     log = "y", # plot on log scale
     las = 1,
     xlab = "",
     ylab = "",
     yaxt = "n")
axis(2, at = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000), 
     labels = c("1", "5", "10", "50", "100", "500", "1,000", "5,000", "10,000", "50,000"), las = 1)
mtext(side = 1, text = "year", line = 2.5, cex = 1.2)
mtext(side = 2, text = "sample size", line = 3, cex = 1.2)

# create violin for each year
violin_width <- 0.3 # play around with width to find a suitable value
for (i in 1:length(years)) {
  # get data for this year
  year_data <- data$sample_size[data$year == years[i] & !is.na(data$sample_size)]
  
  if (length(year_data) > 1) {
    # calculate boxplot statistics
    q1 <- quantile(year_data, 0.25)
    q3 <- quantile(year_data, 0.75)
    iqr <- q3 - q1
    lower_fence <- q1 - 1.5 * iqr
    upper_fence <- q3 + 1.5 * iqr
    lower_whisker <- min(year_data[year_data >= lower_fence])
    upper_whisker <- max(year_data[year_data <= upper_fence])
    outliers <- year_data[year_data < lower_fence | year_data > upper_fence]
    # identify outliers
    outliers <- year_data[year_data < lower_whisker | year_data > upper_whisker]
    # calculate density
    dens <- density(year_data, adjust = 1)
    # trim density to whisker range (not full data range)
    keep <- dens$x >= lower_whisker & dens$x <= upper_whisker
    dens$x <- dens$x[keep]
    dens$y <- dens$y[keep]
    # Add whisker points 
    dens$x <- c(lower_whisker, dens$x, upper_whisker)
    dens$y <- c(0, dens$y, 0)  # violin width = 0 at whiskers
    dens$y <- dens$y / max(dens$y) * violin_width
    # draw violin as polygon
    polygon(c(years[i] + dens$y, rev(years[i] - dens$y)),
            c(dens$x, rev(dens$x)),
            col = NA, border = "black")
    # add median line
    med <- median(year_data)
    segments(years[i] - violin_width, med, 
             years[i] + violin_width, med, 
             col = "#6E2023", lwd = 2)
    
    # add outliers as points
    if (length(outliers) > 0) {
      points(rep(years[i], length(outliers)), outliers, pch = 16, cex = 0.5, col = adjustcolor("black", alpha = 0.6))
    }
  }
}
dev.off()







#####
# sample size per phylogenetic class
rm(list = ls())
data <- read.csv("Data_LiteratureReview.csv",
                 header = T, sep = ",")
str(data)

# check association between sample size and phylogenetic class
# effect size
tab <- table(na.omit(data$class, data$sample_size))
library(rstatix)
cramer_v(tab)

# plot ten most frequent phylogenetic classes only
top_class <- data.frame(class = c("Mammalia","Aves","Actinopterygii","Insecta","Reptilia","Lepidosauria",
                                  "Chondrichthyes","Amphibia","Malacostraca","Gastropoda"),
                        x_axis = 1:10)
# subset data of top ten classes
data <- merge(data, top_class, by = "class", all.x = F, all.y = T)

# get unique x_axis indices
x_axis <- sort(unique(data$x_axis))
pdf(paste0("Plot_SampleSize-Taxon.pdf"), 
    width = 5, height = 5)
# set up the plot
par(mar = c(4,4,1,1), mfrow = c(1,1))
plot(0, 0, type = "n",
     xlim = c(min(x_axis) - 0.05, max(x_axis) + 0.05),
     ylim = c(1, 50000),        # must be > 0
     log = "y",               # <-- log scale
     las = 1,
     xlab = "",
     ylab = "",
     yaxt = "n", xaxt = "n")
axis(2, at = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000), 
     labels = c("1", "5", "10", "50", "100", "500", "1,000", "5,000", "10,000", "50,000"), las = 1)
mtext(side = 1, text = "phylogenetic class", line = 2.5, cex = 1.2)
mtext(side = 2, text = "sample size", line = 3, cex = 1.2)

# create violin for each year
violin_width <- 0.3
for (i in 1:length(x_axis)) {
  # get data for this year
  taxon_data <- data$sample_size[data$x_axis == x_axis[i] & !is.na(data$sample_size)]
  if (length(taxon_data) > 1) {
    # calculate boxplot statistics
    q1 <- quantile(taxon_data, 0.25)
    q3 <- quantile(taxon_data, 0.75)
    iqr <- q3 - q1
    lower_fence <- q1 - 1.5 * iqr
    upper_fence <- q3 + 1.5 * iqr
    lower_whisker <- min(taxon_data[taxon_data >= lower_fence])
    upper_whisker <- max(taxon_data[taxon_data <= upper_fence])
    outliers <- taxon_data[taxon_data < lower_fence | taxon_data > upper_fence]
    # identify outliers
    outliers <- taxon_data[taxon_data < lower_whisker | taxon_data > upper_whisker]
    # calculate density
    dens <- density(taxon_data, adjust = 1)
    # trim density to whisker range (not full data range)
    keep <- dens$x >= lower_whisker & dens$x <= upper_whisker
    dens$x <- dens$x[keep]
    dens$y <- dens$y[keep]
    # add whisker points
    dens$x <- c(lower_whisker, dens$x, upper_whisker)
    dens$y <- c(0, dens$y, 0)  # violin width = 0 at whiskers
    dens$y <- dens$y / max(dens$y) * violin_width
    # draw violin as polygon
    polygon(c(x_axis[i] + dens$y, rev(x_axis[i] - dens$y)),
            c(dens$x, rev(dens$x)),
            col = NA, border = "black")
    # add median line
    med <- median(taxon_data)
    segments(x_axis[i] - violin_width, med, 
             x_axis[i] + violin_width, med, 
             col = "#6E2023", lwd = 2)
    # add outliers as points
    if (length(outliers) > 0) {
      points(rep(x_axis[i], length(outliers)), outliers, pch = 16, cex = 0.5, col = adjustcolor("black", alpha = 0.6))
    }
  }
}
dev.off()




