# Step 1: Load all required libraries at the top
library(tidyverse)
library(GEOquery)

# Step 2: Fetch the data
gset <- getGEO("GSE25724", GSEMatrix = TRUE, getGPL = FALSE)
data_input <- gset[[1]]
print(dim(data_input))

# Step 3: Clean the metadata
# Note: Check if "characteristics_ch1.1" actually contains the disease state for GSE25724!
metadata_modified <- pData(data_input) %>%
  select(title, characteristics_ch1.1) %>%
  rename(disease_state = characteristics_ch1.1) %>%
  mutate(disease_state = gsub("disease state: ", "", disease_state)) %>%
  rownames_to_column(var = "sample")
head(metadata_modified)

# Step 4: Extract and reshape the expression matrix (Wide to Long)
ex <- exprs(data_input)
dat.long <- as.data.frame(ex) %>%
  rownames_to_column(var = "gene") %>%
  gather(key = 'samples', value = 'expression', -gene) %>%
  left_join(metadata_modified, by = c("samples" = "sample"))
head(dat.long)

# Step 5: Visualizations

# A. Global Distribution (This will work immediately)
dat.long %>%
  ggplot(aes(x = expression, fill = disease_state)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Global Gene Expression Distribution")
head(unique(dat.long$gene), 10)

dat.long %>%
  filter(gene == "7892501") %>%  # <--- REAL ID GOES HERE
  ggplot(aes(x = disease_state, y = expression, fill = disease_state)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Target Gene Expression")


# 1. Automatically grab the very first gene in your data
first_gene <- unique(dat.long$gene)[1]

# 2. Plot the Boxplot
dat.long %>%
  filter(gene == first_gene) %>%
  ggplot(aes(x = disease_state, y = expression, fill = disease_state)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste("Target Gene Expression:", first_gene))

# 3. Automatically grab the first 5 genes
first_five_genes <- unique(dat.long$gene)[1:5]

# 4. Plot the Heatmap
dat.long %>%
  filter(gene %in% first_five_genes) %>%
  ggplot(aes(x = samples, y = gene, fill = expression)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
