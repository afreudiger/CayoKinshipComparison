# Compare the performance of different relatedness estimators

# NOTE: This code only creates single-panel plots without labels. Merging single plots to multi-panel plots and adding labels was done using INKSCAPE.

#####
# comparison STR and SNP estimators
rm(list = ls())
all <- read.table("n98_IBD_PED_PIHAT_RAB_STR_PMR_2026-02-04.txt", 
                  sep = "\t", header = T)
str(all)

# plot accuracy and precision for STR-based relatedness
str_snp <- with(all, data.frame(dyad_sorted = dyad_sorted,
                                ibd = prop_genome_shared_8cm,
                                kinlabel_plot = kinlabel_plot))
# extract relevant columns
str_snp <- merge(str_snp, all[,c("dyad_sorted","str41_trioml","str41_wang","str41_lynchli","str41_lynchrd",
                                 "str41_ritland","str41_quellergt","str41_dyadml","pihat_filtered_vcf","rab_filtered_vcf",
                                 "twothree_filtered_vcf")], by = "dyad_sorted")
str(str_snp)
# calculate difference between IBD and each estimator
for(i in colnames(all[,c("str41_trioml","str41_wang","str41_lynchli","str41_lynchrd",
                         "str41_ritland","str41_quellergt","str41_dyadml","pihat_filtered_vcf","rab_filtered_vcf",
                         "twothree_filtered_vcf")])){
  str_snp[,paste0(i,"_ibd")] <- str_snp[,i] - str_snp[,"ibd"]
}
str(str_snp)
# create data frame for plotting
str_snp_plot <- data.frame(dyad_sorted = str_snp$dyad_sorted,
                           r = str_snp$str41_trioml_ibd,
                           estimator = as.character("str41_trioml_ibd"))
library(stringr) # str_match
# loop through each estimator and add the data
for(i in c("str41_wang_ibd","str41_lynchli_ibd","str41_lynchrd_ibd","str41_ritland_ibd","str41_quellergt_ibd",
           "str41_dyadml_ibd","pihat_filtered_vcf_ibd","rab_filtered_vcf_ibd","twothree_filtered_vcf_ibd")){
  str_snp_plot <- rbind(str_snp_plot, data.frame(dyad_sorted = str_snp$dyad_sorted,
                                                 r = str_snp[,i],
                                                 estimator = as.character(i)))
}

# order str in plot according to publication date
ax <- data.frame(estimator = c("str41_wang_ibd","str41_lynchli_ibd","str41_lynchrd_ibd","str41_ritland_ibd","str41_quellergt_ibd",
                               "str41_dyadml_ibd","str41_trioml_ibd","pihat_filtered_vcf_ibd","rab_filtered_vcf_ibd","twothree_filtered_vcf_ibd"),
                 x_axis = as.numeric(c(5,2,4,3,1,6,7,8,9,10)))

str_snp_plot <- merge(str_snp_plot, ax, by = "estimator")



# create colours for plotting
vis_data <- data.frame(col = c(sapply(seq(0.8, 0.2, length.out = 7), function(g) rgb(red = 0, green = g, blue = 0)),
                               "deepskyblue2","deepskyblue3","deepskyblue4"))
# plot it
pdf(paste0("Plot_IBD-STR-SNP_Estimators_2026-02-04.pdf"), 
    width = 9.5, height = 5, useDingbats=FALSE)
par(mar = c(5,6,1,1))
# set up plot
plot(str_snp_plot$r ~ str_snp_plot$x_axis, xlim = c(0.7,10.3), 
     ylim = c(-0.6,0.4),
     type = "n", cex.axis = 1.4, xlab = "", 
     ylab = "", xaxt = "n", 
     yaxt = "n"
)
# add horizontal line for reference
abline(h = 0, lwd = 1, lty = 2, col = "black")
abline(v = 7.5, lwd = 1, lty = 1, col = "black")

# loop through each estimator
for (i in c(1:10)) {
  dat <- subset(str_snp_plot, x_axis == i)$r
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
axis(1, at = 1:10, labels = FALSE)
text(
  x = 1:10,
  y = par("usr")[3] - 0.07,  # same height for all
  labels = c("quellergt","lynchli","ritland","lynchrd","wang",
             "dyadml","trioml","pi-hat",expression("r"["ab"]),"2-out-of-3"),
  xpd = TRUE
)
axis(side = 2, c(-0.6,-0.4,-0.2,0,0.2,0.4), las = 2, cex.axis = 1.2)
mtext(side = 1, text = "estimator", line = 3, cex = 1.7)
mtext(side = 2, text = expression("difference to r"["IBD"]), line = 3.5, cex = 1.7)

dev.off()

# compare the precision and accuracy for the two best performing estimators
median(str_snp$str41_dyadml_ibd); median(str_snp$str41_trioml_ibd)
range(str_snp$str41_dyadml_ibd); range(str_snp$str41_trioml_ibd)

# compute the mean for each estimator
for(i in c("str41_wang_ibd","str41_lynchli_ibd","str41_lynchrd_ibd","str41_ritland_ibd","str41_quellergt_ibd",
           "str41_dyadml_ibd","str41_trioml_ibd","pihat_filtered_vcf_ibd","rab_filtered_vcf_ibd","twothree_filtered_vcf_ibd")){
  print(paste(i, mean(subset(str_snp_plot$r, str_snp_plot$estimator == i)), sep = ": "))
}
# --> DyadML & PIHAT it is










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
# plot semi-transparent points (correct for background relatedness)
all$ibd_adjusted <- (all$prop_genome_shared_8cm - mean(subset(all$prop_genome_shared_8cm, all$kinlabel_plot == "nonkin")))
points(all$pmr_1x ~ all$ibd_adjusted,
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
  # correct PMR for background relatedness
  if(i %in% c("pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")){
    all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - (all_dataset[,"ibd"] - mean(subset(all$prop_genome_shared_8cm, all$kinlabel_plot == "nonkin")))
  } else {
    all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - all_dataset[,"ibd"]
  }
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
  # correct PMR for background relatedness
  if(i %in% c("pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")){
    all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - (all_dataset[,"ibd"] - mean(subset(all$prop_genome_shared_8cm, all$kinlabel_plot == "nonkin")))
  } else {
    all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - all_dataset[,"ibd"]
  }
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
all_dataset <- merge(all_dataset, all[,c("dyad_sorted","str10_dyadml","str20_dyadml","str41_dyadml",
                                         "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
                                         "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")], by = "dyad_sorted")
str(all_dataset)


# calculate difference between IBD and each estimator
for(i in colnames(all_dataset[,c("str10_dyadml","str20_dyadml","str41_dyadml",
                                 "pihat_exons_1","pihat_exons_10","pihat_exons_100","pihat_filtered_vcf",
                                 "pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")])){
  # correct PMR for background relatedness
  if(i %in% c("pmr_0.01x","pmr_0.05x","pmr_0.1x","pmr_1x")){
    all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - (all_dataset[,"ibd"] - mean(subset(all$prop_genome_shared_8cm, all$kinlabel_plot == "nonkin")))
  } else {
    all_dataset[,paste0(i,"_ibd")] <- all_dataset[,i] - all_dataset[,"ibd"]
  }
}
str(all_dataset)

# print summary statistics for all
for(i in colnames(all_dataset[,c(38:71)])){
  print(i)
  print(summary(all_dataset[,i]))
  print("")
}

