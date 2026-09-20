install.packages("BiocManager")
BiocManager::install("phyloseq")
##load(url("https://joey711.github.io/phyloseq-demo/HMPv35.RData"))
BiocManager::install("HMP16SData")
library(HMP16SData)
V35_data <- V35()
V35_data
# the abundance data (bacteria x samples)
abundance <- assay(V35_data)
dim(abundance)

# the sample metadata (includes body site)
metadata <- as.data.frame(colData(V35_data))
head(metadata)
table(metadata$HMP_BODY_SUBSITE)
skin_sites <- c("Left Antecubital Fossa", "Right Antecubital Fossa",
                "Left Retroauricular Crease", "Right Retroauricular Crease")

metadata_skin <- metadata[metadata$HMP_BODY_SUBSITE %in% skin_sites, ]
nrow(metadata_skin)
abundance_skin <- abundance[, rownames(metadata_skin)]
dim(abundance_skin)
metadata_skin$site_simple <- ifelse(grepl("Antecubital", metadata_skin$HMP_BODY_SUBSITE),
                                    "Antecubital Fossa", "Retroauricular Crease")
table(metadata_skin$site_simple)
abundance_t <- t(abundance_skin)
dim(abundance_t)
##attaching the label i am trying to predict
abundance_df <- as.data.frame(abundance_t)
abundance_df$site <- metadata_skin$site_simple
##installing and loading the random forset package
install.packages("randomForest")
library(randomForest)
##training and testing of model with random data
set.seed(42)  # makes the random split reproducible, so we get the same result each time
train_index <- sample(1:nrow(abundance_df), 0.7 * nrow(abundance_df))
train_data <- abundance_df[train_index, ]
test_data <- abundance_df[-train_index, ]
##looking at the raw abundance shapes
hist(abundance_t, main = "Distribution of Raw Abundance Values", xlab = "Abundance")
##sequencing depth
sample_totals <- rowSums(abundance_t)
hist(sample_totals, main = "Total Reads per Sample", xlab = "Total Abundance")
#how many samples each bact type shows up in (prevalence)
bacteria_prevalence <- colSums(abundance_t > 0)
hist(bacteria_prevalence, main = "Number of Samples Each Bacteria Appears In", xlab = "Sample Count")
##pca
pca_result <- prcomp(abundance_t, scale. = FALSE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$site <- metadata_skin$site_simple
head(pca_df)
library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = site)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA of Skin Microbiome Samples", x = "PC1", y = "PC2") +
  theme_minimal()
##retrying pca with log transform
abundance_log <- log1p(abundance_t)  # log1p = log(x + 1), handles zeros safely
pca_result_log <- prcomp(abundance_log, scale. = FALSE)
pca_df_log <- as.data.frame(pca_result_log$x[, 1:2])
pca_df_log$site <- metadata_skin$site_simple

ggplot(pca_df_log, aes(x = PC1, y = PC2, color = site)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA of Skin Microbiome Samples (log-transformed)", x = "PC1", y = "PC2") +
  theme_minimal()
##saving all the plots because i forgot to save
# 1. Raw PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = site)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA of Skin Microbiome Samples", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("pca_raw_plot.png", width = 8, height = 6)
# 2. Log-transformed PCA
ggplot(pca_df_log, aes(x = PC1, y = PC2, color = site)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA of Skin Microbiome Samples (log-transformed)", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("pca_log_plot.png", width = 8, height = 6)
# 3. Raw abundance histogram
png("histogram_raw_abundance.png")
hist(abundance_t, main = "Distribution of Raw Abundance Values", xlab = "Abundance")
dev.off()
# 4. Total reads per sample histogram
png("histogram_sample_totals.png")
hist(sample_totals, main = "Total Reads per Sample", xlab = "Total Abundance")
dev.off()
# 5. Bacteria prevalence histogram
png("histogram_bacteria_prevalence.png")
hist(bacteria_prevalence, main = "Number of Samples Each Bacteria Appears In", xlab = "Sample Count")
dev.off()
abundance_filtered <- abundance_t[, keep_cols]
abundance_log_filtered <- log1p(abundance_filtered)
dim(abundance_log_filtered)
##asking for 2 clusters
set.seed(42)
km_result <- kmeans(abundance_log_filtered, centers = 2)
table(km_result$cluster)
##cross tabulate
table(km_result$cluster, metadata_skin$site_simple)
(475 + 291) / 990
install.packages("mclust")
library(mclust)
adjustedRandIndex(km_result$cluster, metadata_skin$site_simple)
##which bacteria drives each cluster
# average abundance of each bacteria within each cluster
cluster_means <- aggregate(abundance_log_filtered, by = list(cluster = km_result$cluster), FUN = mean)
dim(cluster_means)
##top drivers
# difference between cluster 1 and cluster 2 for each bacteria
cluster_diff <- cluster_means[1, -1] - cluster_means[2, -1]
sorted_diff <- sort(cluster_diff, decreasing = TRUE)

# top 10 bacteria more abundant in cluster 1
head(sorted_diff, 10)

# top 10 bacteria more abundant in cluster 2
tail(sorted_diff, 10)
##converting to plain numeric vector first
cluster_diff <- as.numeric(cluster_means[1, -1] - cluster_means[2, -1])
names(cluster_diff) <- colnames(cluster_means)[-1]

sorted_diff <- sort(cluster_diff, decreasing = TRUE)

# top 10 bacteria more abundant in cluster 1
head(sorted_diff, 10)

# top 10 bacteria more abundant in cluster 2
tail(sorted_diff, 10)

library(SummarizedExperiment)
taxonomy <- as.data.frame(rowData(V35_data))
head(taxonomy)
top_cluster1_otus <- names(head(sorted_diff, 10))
top_cluster2_otus <- names(tail(sorted_diff, 10))

taxonomy[top_cluster1_otus, c("GENUS", "FAMILY")]
taxonomy[top_cluster2_otus, c("GENUS", "FAMILY")]

