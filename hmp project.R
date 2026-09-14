install.packages("BiocManager")
BiocManager::install("phyloseq")
load(url("https://joey711.github.io/phyloseq-demo/HMPv35.RData"))
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

