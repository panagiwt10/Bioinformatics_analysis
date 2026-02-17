## Load packages
if (!require("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
}

if (!require("vegan", quietly = TRUE)){
  install.packages('vegan', repos = c('https://vegandevs.r-universe.dev','https://cloud.r-project.org'))
}
if (!require("patchwork", quietly = TRUE)){
  install.packages("patchwork")  
}
if (!require("ggpubr", quietly = TRUE)){
  install.packages("ggpubr")  
}

library("tidyverse")
library("vegan")      # for diversity()
library("patchwork")  # to combine plots, optional
library("ggpubr")   #for the statistical analysis
library("tidyverse")

## Read data -------------------------------------------------------------
# Dataset: YACHIDA_CRC_2019

options(timeout = 300) 

url_metadata <- "https://raw.githubusercontent.com/borenstein-lab/microbiome-metabolome-curated-data/refs/heads/main/data/processed_data/YACHIDA_CRC_2019/metadata.tsv"
meta <- read.table(
  url_metadata,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

url_genera_counts <- "https://raw.githubusercontent.com/borenstein-lab/microbiome-metabolome-curated-data/main/data/processed_data/YACHIDA_CRC_2019/genera.counts.tsv"
genera_counts <- read.table("./genera.counts.tsv",
                            header = TRUE, sep = "\t",
                            check.names = FALSE,
                            stringsAsFactors = FALSE)

## -------------------------------
## ΕΡΩΤΗΜΑ 1.1 – Study design
## -------------------------------

# Πόσα δείγματα περιλαμβάνει το dataset
num_samples <- nrow(meta)
num_samples

# Σε ποιες ομάδες (συνθήκες) χωρίζονται τα δείγματα
group_distribution <- table(meta$Study.Group)
group_distribution


# Μετατροπή σε dataframe (για plots & πίνακες)
group_df <- as.data.frame(group_distribution)
colnames(group_df) <- c("Group", "Number_of_samples")


# Προβολή πίνακα ομάδων (για αναφορά)
group_df

## -------------------------------
## Plot 1: Barplot κατανομής δειγμάτων ανά ομάδα
## -------------------------------

ggplot(group_df, aes(x = Group, y = Number_of_samples, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Distribution of samples across study groups",
    x = "Study Group",
    y = "Number of samples"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = "Q1_Group_Distribution.png",
  width = 7,
  height = 4,
  dpi = 300
)

## Plot 2: Density Plot of Age by Group ---------------------------------
p_age <- ggplot(meta, aes(x = Age, fill = Study.Group)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density Plot of Age by Group", x = "Age", y = "Density") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p_age
ggsave("Age_Density_by_Group.png", plot = p_age, width = 7, height = 4, dpi = 300)