library(tidyverse)
library(ggprism)
library(umap)
library(Rtsne)
library(ggpubr)
library(neuralnet)
#library(factoextra)

brain_cancer_GSE <- read_csv("./data/Brain_GSE50161.csv")
type <- brain_cancer_GSE[ , "type"]
x <- brain_cancer_GSE %>%
  select(where(is.numeric)) %>%
  column_to_rownames("samples")

#cleaning----
x <- x %>% select(-contains("AFFX-")) #removing affy control probes
x <- x %>% select(-contains("BGP-")) #removes background control probes 

#PCA----
x_mean <- cbind(type, x)
x_mean <- x_mean %>%
  group_by(type) %>%
  summarise(across(everything(), mean))
pca_fit <- prcomp(x_mean %>% 
                    select(where(is.numeric)), 
                  center = TRUE, 
                  scale = FALSE)
var_explained_df <- data.frame(PC= paste0("PC",1:length(pca_fit$sdev)),
                               Variance=(pca_fit$sdev)^2/sum((pca_fit$sdev)^2))
var_explained_df$label <- round(var_explained_df$Variance, 3)

ggplot(data = var_explained_df, aes(x=PC,y=Variance, group = 1, label = label)) +
  geom_point() + 
  geom_line() + 
  geom_text(hjust = 0, nudge_x = 0.05) +
  theme_prism()

pc1 <- data.frame(pca_fit$rotation[, 1]) %>% 
  arrange(desc(pca_fit.rotation...1.))
top5 <- head(pc1, 5)
top5$gene <- c("GABRG2", "STMN2", "MYT1L", "NEFL", "MAP4")

#tSNE----
tSNE_fit <- Rtsne(x)
tSNE_df <- data.frame(tSNE1 = tSNE_fit$Y[,1], #1st tSNE dimension
                      tSNE2 = tSNE_fit$Y[,2], #2nd tSNE dimension
                      type = type)

ggplot(tSNE_df, aes(tSNE1, tSNE2, colour = type)) +
  geom_point() +
  theme_prism()

#UMAP----
umap_fit <- umap(x)
umap_df <- data.frame(UMAP1 = umap_fit$layout[,1], #1st UMAP dimension
                      UMAP2 = umap_fit$layout[,2], #2nd UMAP dimension
                      type = type)

ggplot(umap_df, aes(UMAP1, UMAP2, colour = type)) +
  geom_point() +
  theme_prism()

#k-means----
res.km <- kmeans(x, 
                 5, 
                 iter.max = 10, 
                 nstart = 25)

fviz_cluster(res.km, 
             data = brain_cancer_GSE %>%
               select(where(is.numeric)) %>%
               column_to_rownames("samples"),
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_prism()
)

#experimental----
training_data_rows <- floor(0.70 * nrow(brain_cancer_GSE))          
set.seed(123) 
training_indices <- sample(c(1:nrow(brain_cancer_GSE)), training_data_rows)

training_data <- brain_cancer_GSE[training_indices,] 
test_data <- brain_cancer_GSE[-training_indices,]

names(training_data) <- make.unique(names(training_data))
names(test_data) <- make.unique(names(test_data))
#names(training_data) <- gsub("-", "", names(training_data), fixed = TRUE)
#formula <- paste(names(training_data), collapse='+')
#formula <- substr(formula, 1, 13)

nn = neuralnet(type ~ .,  
             data=training_data, 
             #hidden=c(25000, 10000, 2000, 100, 10), 
             hidden = c(10), #for testing purposes
             linear.output = FALSE,
             lifesign = "full")
