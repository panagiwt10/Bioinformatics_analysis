library(ggplot2)
library(vegan)
library(ggpubr)

install.packages("ggplot2")

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)  

if (!require("dplyr")) install.packages("dplyr")
if (!require("magrittr")) install.packages("magrittr")
if (!require("patchwork", quietly = TRUE)){
  install.packages("patchwork")  
}

library(dplyr)
library(magrittr)  # αυτό φέρνει το %>%

# 1) Links
url_meta <- "https://raw.githubusercontent.com/borenstein-lab/microbiome-metabolome-curated-data/main/data/processed_data/YACHIDA_CRC_2019/metadata.tsv"
url_counts <- "https://raw.githubusercontent.com/borenstein-lab/microbiome-metabolome-curated-data/main/data/processed_data/YACHIDA_CRC_2019/genera.counts.tsv"

meta <- read.delim("metadata.tsv")
counts <- read.delim("genera.counts.tsv")

# 1) Φτιάχνουμε σωστά το count matrix
rownames(counts) <- counts$Sample
count_mat <- as.matrix(counts[, -1])

# 2) Καθαρίζουμε πιθανά κενά
meta$Sample <- trimws(meta$Sample)
rownames(count_mat) <- trimws(rownames(count_mat))

# 3) Επιλέγουμε τις 2 ομάδες
meta2 <- meta %>% filter(Study.Group %in% c("Healthy", "Stage_III_IV"))

# 4) Κρατάμε ΜΟΝΟ τα κοινά Samples (ΑΥΤΟ είναι το κλειδί)
common_samples <- intersect(meta2$Sample, rownames(count_mat))

# (προαιρετικό) έλεγχος πόσα χάνονται
cat("meta2 samples:", nrow(meta2), "\n")
cat("count_mat samples:", nrow(count_mat), "\n")
cat("common samples:", length(common_samples), "\n")

# 5) Φιλτράρουμε και τα δύο ώστε να ταιριάζουν 
meta2 <- meta2 %>% filter(Sample %in% common_samples)
count_mat2 <- count_mat[common_samples, , drop = FALSE]

##############################################
# ALPHA DIVERSITY ANALYSIS

# Richness (Observed genera) per sample
observed_richness <- apply(count_mat2, 1, function(x) sum(x > 0))

# Shannon index per sample
shannon <- diversity(count_mat2, index = "shannon")

# Create dataframe with metadata + metrics
alpha_df <- data.frame(
  Sample = rownames(count_mat2),
  group = meta2$Study.Group[match(rownames(count_mat2), meta2$Sample)],
  Richness = observed_richness,
  Shannon = shannon
)

# --------------------------
# 2.2 Mean values per condition
# Means
alpha_df %>%
  group_by(group) %>%
  summarise(
    mean_Richness = mean(Richness),
    mean_Shannon = mean(Shannon),
    n = n()
  )

print("Mean alpha-diversity values per group:")
print(alpha_df)

# --------------------------
# 2.3 Wilcoxon rank-sum tests
wilcox_richness <- wilcox.test(Richness ~ group, data = alpha_df)
wilcox_shannon  <- wilcox.test(Shannon  ~ group, data = alpha_df)

print("Wilcoxon test - Richness:")
print(wilcox_richness)

print("Wilcoxon test - Shannon:")
print(wilcox_shannon)


###########################
# Boxplot for Richness
p1 <- ggplot(alpha_df, aes(x = group, y = Richness, fill = group)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(c("Breast milk", "Standard infant formula")), # Συγκρίσεις για τις 2 ομάδες
    method = "wilcox.test",  # Χρήση του Wilcoxon test
    label = "p.format",      # Προβολή του p-value
    label.x = 1.5,           # Προσαρμογή της θέσης του p-value
    hide.ns = FALSE          # Προβολή και όταν το p-value είναι μη σημαντικό
  ) +
  labs(title = "Alpha diversity: Richness (Observed genera)", x = NULL, y = "Richness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Boxplot for Shannon Index
p2 <- ggplot(alpha_df, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(c("Breast milk", "Standard infant formula")),
    method = "wilcox.test",
    label = "p.format", 
    label.x = 1.5,  # Adjust position of p-value label if needed
    hide.ns = FALSE
  ) +
  labs(title = "Shannon Index", x = NULL, y = "Shannon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show the plots
p1
p2

#Save the plots
ggsave("Richness_Boxplot.png", plot = p1, width = 6, height = 4, dpi = 300)
ggsave("Shannon_Boxplot.png", plot = p2, width = 6, height = 4, dpi = 300)

#Optional: save mean table
write.csv(alpha_df, "Q2_means_Healthy_vs_StageIIIIV.csv", row.names = FALSE)


#3 υποερωτημα

#allignement 
meta2 <- meta2 %>% arrange(match(Sample, rownames(count_mat2)))
stopifnot(all(meta2$Sample == rownames(count_mat2)))

# 1) jaccard Dissimilarity 
# Μετατρέπουμε counts -> 0/1
count_pa <- (count_mat2 > 0) * 1

dist_jaccard <- vegdist(count_pa, method= "jaccard")


# 2) PCoA
pcoa_res <- cmdscale(dist_jaccard, eig = TRUE, k = 2)

# % variance explained
var_exp <- round(100 * pcoa_res$eig[1:2] / sum(pcoa_res$eig), 1)

beta_df <- data.frame(
  PC1 = pcoa_res$points[, 1],
  PC2 = pcoa_res$points[, 2],
  Sample = rownames(pcoa_res$points)
) %>%
  left_join(meta2 %>% select(Sample, Study.Group), by = "Sample") %>%
  mutate(Group_short = recode(Study.Group,
                              "Healthy" = "Healthy",
                              "Stage_III_IV" = "Stage III-IV"))


# 3) PERMANOVA
permanova <- adonis2(dist_jaccard ~ Study.Group, data = meta2, permutations = 999)
r2 <- round(permanova$R2[1], 3)
pval <- signif(permanova$`Pr(>F)`[1], 3)

print(permanova)

# 4) Plot (PCoA) + ellipse + annotation
p_beta <- ggplot(beta_df, aes(PC1, PC2, color = Group_short)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(level = 0.68, linewidth = 1) +
  annotate(
    "text", x = -Inf, y = -Inf,
    hjust = -0.1, vjust = -0.2,
    label = paste0("PERMANOVA\nR² = ", r2, "\np = ", pval),
    size = 3.5
  ) +
  labs(
    title = "PCoA (Jaccard dissimilarity)",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "Group"
  ) +
  theme_bw()

p_beta
ggsave("Q3_Jaccard_PCoA_Healthy_vs_StageIIIIV.png", plot = p_beta, width = 7, height = 5, dpi = 300)