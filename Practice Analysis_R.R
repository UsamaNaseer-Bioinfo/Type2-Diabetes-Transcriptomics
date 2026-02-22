# Gene Expression Analysis 
library(tidyverse)
library(GEOquery)
#Get the Data Fresh
gset <- getGEO("GSE25724", GSEMatrix = TRUE, getGPL = FALSE)
data_input<- gset[[1]]
print(dim(data_input))
meta_data <- pData(data_input)
# 2. Clean it up using the tidyverse method (%>%)
metadata_modefied <- meta_data %>%
  select(title, characteristics_ch1.1) %>%
  rename(disease_state = characteristics_ch1.1) %>%
  mutate(disease_state = gsub("disease state: ","",disease_state))
head(metadata_modefied)
# Extract the Expression Matrix
ex <- exprs(data_input)
#Convert to a dataframe and make the row names a 'gene' column
dat <- as.data.frame(ex) %>%
  rownames_to_column(var = "gene")
#Reshape from wide to long (just like your script)
data.modified <- dat %>%
  gather(key = 'samples' , value = 'expression',-gene )
head(data.modified)
#Turn the metadata row names into a 'samples' colum
metadata_ready<- metadata_modefied %>%
  rownames_to_column(var = "sample")
#Join the expression data with the metadata
dat.long <- data.modified %>%
  left_join(metadata_ready, by = c("samples" = "sample") )
head(dat.long)

library(ggplot2)
# Start with clean & join database
dat.long %>%
  filter(gene== "1007_s_at") %>%
  ggplot(aes(x = disease_state, y = expression, fill = disease_state)) +
  geom_boxplot() +
  theme_minimal()

#Select 5 specific genes from our dataset
gene_of_intrest <- c("1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at")
dat.long %>%
  filter(gene %in% gene_of_intrest) %>%
  ggplot(aes(x = samples, y = gene, fill = expression)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


dat.long %>%
  ggplot(aes(x = expression, fill = disease_state)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Global Gene Expression Distribution")