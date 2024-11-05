# Load necessary libraries
library(limma)
library(caret)
library(glmnet)
library(vip)
library(pROC)

# Set working directory
setwd("directory") 

# Load expression and phenotype data (Extracted from Wallinhan et al. GSE)
Expresion <- readRDS("exprs_myco_filt.RDS") 
Pheno <- readRDS("pheno_myco_filt.RDS")  

# Assign row names and check data structure
rownames(Pheno) <- colnames(Expresion)
View(Expresion)
View(Pheno)

# Prepare data for analysis
GSE103119_Myco <- list("pheno" = Pheno, "genes" = Expresion)
condition <- as.factor(GSE103119_Myco$pheno$Bacterial)
age <- as.numeric(GSE103119_Myco$pheno$age)
data <- as.data.frame(GSE103119_Myco$genes) 
design <- model.matrix(~ 0 + condition + age) 

# Fit linear model
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(status = conditionMycoplasma - conditionViral, levels = colnames(design))
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)

# Get results and apply multiple testing correction
results <- topTable(fit3, coef = 1, n = Inf, adjust = "BH")  
results_0.05 <- results[results$adj.P.Val < 0.05, ] 
results_0.01 <- results[results$adj.P.Val < 0.01, ] 


# Filter results based on expression threshold
threshold <- 1
results_best <- results_0.01[results_0.01$AveExpr > 2 & abs(results_0.01$logFC) > threshold,]
View(results_best)

# Prepare for training and test splits
set.seed(123) # Set seed for reproducibility
unique_splits <- list()
list_train <- list()
list_test <- list()

# Create 1000 stratified random samples (Bacterial is the colum where Viral/Mycoplasma information is stored)
for(i in 1:1000) {
  trainIndex <- createDataPartition(y = Pheno$Bacterial, p = 0.7, list = FALSE)
  if (!identical(trainIndex, unique_splits)) {
    unique_splits[[i]] <- trainIndex  
    train_data <- Pheno[trainIndex, ]
    test_data <- Pheno[-trainIndex, ]
    list_train[[i]] <- train_data
    list_test[[i]] <- test_data
  }
} 

list_GSEexprs_70 <- list()
list_GSEexprs_30 <- list()

# Filter expression data for training and test sets
for (i in 1:1000) {
  rownames_train_data <- rownames(list_train[[i]])
  rownames_test_data <- rownames(list_test[[i]])
  
  GSEexprs_70 <- Expresion[, colnames(Expresion) %in% rownames_train_data]
  GSEexprs_30 <- Expresion[, colnames(Expresion) %in% rownames_test_data]
  
  list_GSEexprs_70[[i]] <- GSEexprs_70
  list_GSEexprs_30[[i]] <- GSEexprs_30 
} 

# Prepare datasets for models
list_GSE103119_70 <- list() 
list_GSE103119_30 <- list()

for (i in 1:1000) {
  GSE103119_70 <- list("pheno" = list_train[[i]], "genes" = list_GSEexprs_70[[i]])
  GSE103119_30 <- list("pheno" = list_test[[i]], "genes" = list_GSEexprs_30[[i]])
  
  list_GSE103119_70[[i]] <- GSE103119_70 
  list_GSE103119_30[[i]] <- GSE103119_30
}   

# Prepare matrices for the 1000 models with the same genes
list_matrices <- list()

for (i in 1:1000) {
  data <- as.data.frame(list_GSE103119_70[[i]]$genes)
  matriz <- t(data) 
  matriz <- as.data.frame(matriz) 
  matriz <- matriz[, colnames(matriz) %in% rownames(results_best)]  
  
  pheno_reduced <- list_GSE103119_70[[i]]$pheno[, c(48, 47)]
  matriz <- merge(matriz, pheno_reduced, by.x = 0, by.y = 0)
  
  rownames(matriz) <- matriz$Row.names
  matriz <- matriz[, -1] 
  matriz <- matriz[, !names(matriz) %in% "GENDER"] # Remove 'GENDER' column
  
  list_matrices[[i]] <- matriz
} 

# Fit LASSO coefficients for the 1000 datasets 
lasso_results <- list()
vip_results <- list()
set.seed(164) # Set seed for reproducibility

for (i in 1:1000) {
  current_matriz <- list_matrices[[i]]
  
  # Fit LASSO model
  cv.lassoModel <- cv.glmnet(
    x = data.matrix(current_matriz[, -which(names(current_matriz) == "Bacterial")]),
    y = current_matriz$Bacterial,
    standardize = TRUE,
    alpha = 1,
    nfolds = 10,
    family = "multinomial",
    parallel = TRUE
  )
  
  lasso_results[[i]] <- cv.lassoModel
  
  # Store variable importance results
  num_variables_totales <- ncol(current_matriz) - 1
  vip_results[[i]] <- vip(cv.lassoModel, num_features = num_variables_totales)
} 

# Extract coefficients from each LASSO model
lasso_coefs <- list() 
selected_mycoplasma_predictors <- list()
selected_mycoplasma_predictors_df <- list()

for (i in 1:1000) {
  cv.lassoModel <- lasso_results[[i]]
  idealLambda <- cv.lassoModel$lambda.min
  coefs <- coef(cv.lassoModel, s = idealLambda, exact = TRUE)
  
  lasso_coefs[[i]] <- coefs
  predictors_mycoplasma <- rownames(coefs$Mycoplasma)[which(coefs$Mycoplasma != 0)]
  selected_mycoplasma_predictors[[i]] <- predictors_mycoplasma  
  selected_mycoplasma_predictors_df[[i]] <- as.data.frame(selected_mycoplasma_predictors[[i]])
  
  # Save results to a unique file for each model
  file_name <- paste0("resultados_lasso_1000/modelo_", i, ".txt")
  write.table(selected_mycoplasma_predictors_df[[i]], file = file_name, sep = "\t", quote = FALSE)
}

# Filter out the intercept from the selected predictors
firma <- vector("list", length = length(selected_mycoplasma_predictors))

for (i in seq_along(selected_mycoplasma_predictors)) {
  current_predictors <- selected_mycoplasma_predictors[[i]]
  filtered_firma <- setdiff(current_predictors, "(Intercept)")
  firma[[i]] <- filtered_firma
}

# Prepare matrices for the LASSO models
list_matrices_3 <- list()

for (i in 1:1000) {
  data <- as.data.frame(list_GSE103119_70[[i]]$genes)
  matriz <- t(data)
  matriz <- as.data.frame(matriz) 
  matriz <- matriz[, colnames(matriz) %in% rownames(Expresion)]  
  
  pheno_reduced <- list_GSE103119_70[[i]]$pheno[, c(48, 47)]
  matriz <- merge(matriz, pheno_reduced, by.x = 0, by.y = 0)
  
  rownames(matriz) <- matriz$Row.names
  matriz <- matriz[, -1] 
  matriz <- matriz[, !names(matriz) %in% "GENDER"] # Remove 'GENDER' column
  
  list_matrices_3[[i]] <- matriz
} 

# Fit LASSO models using glm
lasso_results_3 <- vector("list", length = length(list_matrices_3))

for (i in seq_along(list_matrices_3)) {
  tryCatch({
    current_matriz <- list_matrices_3[[i]]
    current_firma <- firma[[i]]
    
    # Convert bacterial labels to numeric
    current_matriz$Bacterial <- gsub("Mycoplasma", "1", current_matriz$Bacterial)
    current_matriz$Bacterial <- gsub("Viral", "0", current_matriz$Bacterial)
    current_matriz$Bacterial <- as.numeric(current_matriz$Bacterial) 
    
    # Adjust column names if necessary
    if ("NKX3-1" %in% colnames(current_matriz)) {
      colnames(current_matriz)[colnames(current_matriz) == "NKX3-1"] <- "NKX3_1"
    }
    
    # Select only columns in the signature
    genes_en_firma <- intersect(current_firma, colnames(current_matriz))
    genes_en_firma_quoted <- lapply(genes_en_firma, function(x) if (grepl("[-]", x)) paste0("`", x, "`") else x)
    
    # Build the formula for the model
    formula <- as.formula(paste("Bacterial ~", paste(genes_en_firma_quoted, collapse = " + ")))
    
    # Fit the logistic regression model
    finalLasso <- glm(formula, data = current_matriz, family = binomial(link = "logit"))
    
    lasso_results_3[[i]] <- finalLasso
  }, error = function(e) {
    cat("Error in dataset", i, ":", conditionMessage(e), "\n")
    lasso_results_3[[i]] <- NULL
  })
}  

# Create test matrices with selected genes
list_matrices_test <- list()

for (i in 1:1000) {
  data_2 <- as.data.frame(list_GSE103119_30[[i]]$genes)
  matriz_test <- t(data_2)
  matriz_test <- as.data.frame(matriz_test) 
  
  pheno_reduced <- list_GSE103119_30[[i]]$pheno[, c(48, 47)]
  matriz_test <- merge(matriz_test, pheno_reduced, by.x = 0, by.y = 0)
  
  rownames(matriz_test) <- matriz_test$Row.names
  matriz_test <- matriz_test[, -1] 
  matriz_test$Bacterial <- gsub("Mycoplasma", "1", matriz_test$Bacterial)
  matriz_test$Bacterial <- gsub("Viral", "0", matriz_test$Bacterial)
  matriz_test$Bacterial <- as.numeric(matriz_test$Bacterial) 
  
  list_matrices_test[[i]] <- matriz_test 
}

# ROC analysis on training models
roc_list <- list()

for (i in seq_along(list_matrices_3)) {
  current_matriz <- list_matrices_3[[i]]
  finalLasso <- lasso_results_3[[i]]
  
  roc_result <- roc(current_matriz$Bacterial, fitted(finalLasso), smooth = FALSE)
  ci_results <- ci(roc_result)
  
  roc_list[[i]] <- list(roc_results = roc_result, ci = ci_results)
}

all_auc_values <- numeric(length(roc_list))
all_ci_values <- vector("list", length(roc_list))

for (i in seq_along(roc_list)) {
  all_auc_values[i] <- roc_list[[i]]$roc_results$auc
  all_ci_values[[i]] <- roc_list[[i]]$ci
}

# ROC analysis on test datasets
roc_list_test <- list()

for (i in seq_along(list_matrices_test)) {
  matriz_test <- list_matrices_test[[i]]
  finalLasso <- lasso_results_3[[i]]
  
  predictions <- predict(finalLasso, newdata = matriz_test, type = "response")
  
  roc_result_test <- roc(matriz_test$Bacterial, predictions, smooth = FALSE)
  ci_results_test <- ci(roc_result_test)
  
  roc_list_test[[i]] <- list(roc_results = roc_result_test, ci = ci_results_test)
} 

all_auc_values_test <- numeric(length(roc_list_test))
all_ci_values_test <- vector("list", length(roc_list_test))


for (i in seq_along(roc_list_test)) {
  all_auc_values_test[i] <- roc_list_test[[i]]$roc_results$auc
  # all_ci_values_test[[i]] <- roc_list_test[[i]]$ci  # Confidence intervals can be extracted here if needed
}

# Create a dataframe with AUC values for training and test
df_aucs <- data.frame(x = all_auc_values, y = all_auc_values_test)

# Create a new dataframe with counts of genes
df_counts <- data.frame(x = df_aucs$x, y = df_aucs$y, genes = numeric(length(df_aucs$x)))
View(df_counts)

# Add confidence intervals for training and test
df_counts$ci_train <- all_ci_values
df_counts$ci_test <- all_ci_values_test

# Iterate through the rows of df_counts to add a column with the number of genes
for (i in seq_along(df_counts$x)) {
  current_list <- firma[[i]]
  genes <- length(current_list)
  df_counts$genes[i] <- genes
}  

View(df_counts)

# Filter all rows that have a count of genes < 10
df_check <- df_counts[df_counts$genes < 11, ]
df_check$dataset <- rownames(df_check)

View(df_check) 

# Write the filtered results to a text file
write.table(df_check, "table_signatures_10_genes.txt", sep = "\t", dec = ",")

#### Load RNA-Seq expression data from EUCLIDS
pheno_merged <- readRDS("pheno_merged_euclids.rds")
counts_merged <- readRDS("counts_merged_euclids.rds")
cts <- readRDS("cts_euclids_b_v.rds")

# Prepare the matrix for EUCLIDS data
mat2 <- t(cts)
matriz_euclids <- mat2[, colnames(mat2) %in% rownames(counts_merged)] 
matriz_euclids <- as.data.frame(matriz_euclids)
View(matriz_euclids)  

# Merge phenotype data with expression matrix
pheno_merged_filtered <- pheno_merged[, c(2, 57)]
matriz_euclids <- merge(matriz_euclids, pheno_merged_filtered, by.x = 0, by.y = 0) 
rownames(matriz_euclids) <- matriz_euclids$Row.names
matriz_euclids <- matriz_euclids[, -c(1, 40647)] # Adjust position as needed
View(matriz_euclids)

# Select models with 10 genes or fewer to test in EUCLIDS data 
index_selected <- as.numeric(df_check$dataset)
lasso_results_selected <- vector("list", length = length(index_selected))

# Loop to select the appropriate LASSO results
for (i in seq_along(index_selected)) {
  indice_actual <- index_selected[i]
  if (indice_actual <= length(lasso_results_3)) {
    lasso_results_selected[[i]] <- lasso_results_3[[indice_actual]]
  } else {
    lasso_results_selected[[i]] <- NULL
  }
}

roc_list_EUCLIDS <- list()

# Iterate through the selected models to create ROC curves
for (i in seq_along(lasso_results_selected)) {
  predictions <- predict(lasso_results_selected[[i]], newdata = matriz_euclids, type = "response")
  roc_result <- roc(matriz_euclids$pneumonia.type, predictions)
  ci_results <- ci(roc_result)
  
  # Store the ROC results in the list
  roc_list_EUCLIDS[[i]] <- list(roc_result = roc_result, ci = ci_results)
}

# Initialize AUC and CI for EUCLIDS
all_auc_values_EUCLIDS <- numeric(length(roc_list_EUCLIDS))
all_ci_values_EUCLIDS <- vector("list", length(roc_list_EUCLIDS))

# Loop through roc_list_EUCLIDS to extract AUC and CI
for (i in seq_along(roc_list_EUCLIDS)) {
  all_auc_values_EUCLIDS[i] <- roc_list_EUCLIDS[[i]]$roc_result$auc
  all_ci_values_EUCLIDS[[i]] <- roc_list_EUCLIDS[[i]]$ci
}

# Add EUCLIDS AUC and CI to df_check
df_check$euclids <- all_auc_values_EUCLIDS
df_check$ci_euclids <- all_ci_values_EUCLIDS

# Calculate the combined weighting score
df_check$ponderation <- df_check$x * 0.4 + df_check$y * 0.4 - df_check$euclids * 0.2
View(df_check)  

# Select the best signatures based on specific indices
selected_signatures <- c(543, 294, 897, 82, 778, 780, 741, 744) 
df_seleccion_2 <- df_check[rownames(df_check) %in% selected_signatures, ]
View(df_seleccion_2) 

# Assign new indices for the selected signatures
df_seleccion_2$nuevos_indices <- 1:8
View(df_seleccion_2) 

# Select LASSO results for the chosen index
lasso_results_selected <- vector("list", length = length(selected_signatures))

for (i in seq_along(selected_signatures)) {
  indice_actual <- selected_signatures[i]
  if (indice_actual <= length(lasso_results_3)) {
    lasso_results_selected[[i]] <- lasso_results_3[[indice_actual]]
  } else {
    lasso_results_selected[[i]] <- NULL
  }
} 

#### Create a list of matrices excluding the matrix from which the signature was derived
nueva_lista_matrices_3 <- list_matrices_3[-543]  # Exclude matrix 543
nueva_lista_matrices_4 <- list_matrices_3[-294]  # Exclude matrix 294
nueva_lista_matrices_5 <- list_matrices_3[-897]  # Exclude matrix 897
nueva_lista_matrices_6 <- list_matrices_3[-82]   # Exclude matrix 82
nueva_lista_matrices_7 <- list_matrices_3[-778]  # Exclude matrix 778
nueva_lista_matrices_8 <- list_matrices_3[-780]  # Exclude matrix 780
nueva_lista_matrices_9 <- list_matrices_3[-741]  # Exclude matrix 741
nueva_lista_matrices_10 <- list_matrices_3[-744] # Exclude matrix 744

# Initialize lists to store ROC results for each excluded matrix
roc_list_retest3 <- list()
roc_list_retest4 <- list()
roc_list_retest5 <- list()
roc_list_retest6 <- list()
roc_list_retest7 <- list()
roc_list_retest8 <- list()
roc_list_retest9 <- list()
roc_list_retest10 <- list()

# Iterate over each test dataset and tºhe chosen LASSO model
for (i in seq_along(nueva_lista_matrices_10)) {
  
  # Get the current test dataset
  matriz_test <- nueva_lista_matrices_10[[i]]
  finalLasso <- lasso_results_3[[744]]  # Use the model corresponding to the excluded matrix in each case  (744,741,780,778,82,897,294,543)
  
  # Convert 'Bacterial' values to numeric (0 and 1)
  matriz_test$Bacterial <- gsub("Mycoplasma", "1", matriz_test$Bacterial)
  matriz_test$Bacterial <- gsub("Viral", "0", matriz_test$Bacterial)
  matriz_test$Bacterial <- as.numeric(matriz_test$Bacterial)
  
  # Calculate predictions
  predictions <- predict(finalLasso, newdata = matriz_test, type = "response")
  
  # Calculate the ROC curve
  roc_result_test <- roc(matriz_test$Bacterial, predictions, smooth = FALSE)
  
  # Store the ROC curve in the list
  roc_list_retest10[[i]] <- roc_result_test
} 


# Extract AUC values from the ROC lists
all_auc_values_retest3 <- unlist(lapply(roc_list_retest3, function(roc) roc$auc)) 
all_auc_values_retest4 <- unlist(lapply(roc_list_retest4, function(roc) roc$auc))
all_auc_values_retest5 <- unlist(lapply(roc_list_retest5, function(roc) roc$auc))
all_auc_values_retest6 <- unlist(lapply(roc_list_retest6, function(roc) roc$auc)) 
all_auc_values_retest7 <- unlist(lapply(roc_list_retest7, function(roc) roc$auc)) 
all_auc_values_retest8 <- unlist(lapply(roc_list_retest8, function(roc) roc$auc)) 
all_auc_values_retest9 <- unlist(lapply(roc_list_retest9, function(roc) roc$auc)) 
all_auc_values_retest10 <- unlist(lapply(roc_list_retest10, function(roc) roc$auc)) 

#####  save all AUC values
save_auc_values <- function(roc_list, file_path) {
  all_auc_values <- unlist(lapply(roc_list, function(roc) roc$auc))
  saveRDS(all_auc_values, file_path)
}  

##### TEST the 8  signatures other phenotypes included in Wallinhan et. all dataset following the same script 
pheno_others<-readRDS("pheno_myco_others.RDS")
expresion_others<-readRDS("exprs_myco_others.RDS")

mat <- t(Expresion_myco_others)
matriz_myco_others <- as.data.frame(mat)
pheno_reduced <- Pheno_myco_others[, c(45, 46)]

# Merge with actual columns
matriz_myco_others <- merge(matriz_myco_others, pheno_reduced, by.x = 0, by.y = 0)
rownames(matriz_myco_others) <- matriz_myco_others$Row.names
matriz_myco_others <- matriz_myco_others[, -1]
View(matriz_myco_others)

col_to_remove <- which(names(matriz_myco_others) == "viral.organism.ch")

# Eliminate column not used 
if (length(col_to_remove) > 0) {
  matriz_myco_others <- matriz_myco_others[, -col_to_remove]
}
View(matriz_myco_others)

roc_list_others <- list()

# Iterate over each test dataset and the selected LASSO model
for (i in seq_along(lasso_results_seleccionados)) {
  
  # Set the test matrix for "others" from the array script
  matriz_test <- matriz_myco_others
  finalLasso <- lasso_results_seleccionados[[i]]
  
  # Change variable values to numeric (0 and 1)
  matriz_myco_others$Bacterial <- gsub("Mycoplasma", "1", matriz_myco_others$Bacterial)
  matriz_myco_others$Bacterial <- gsub("Others", "0", matriz_myco_others$Bacterial)
  matriz_myco_others$Bacterial <- as.numeric(matriz_myco_others$Bacterial)
  
  # Calculate predictions
  predictions <- predict(finalLasso, newdata = matriz_myco_others, type = "response")
  
  # Calculate ROC curve
  roc_result_test <- roc(matriz_myco_others$Bacterial, predictions, smooth = FALSE)
  ci_results <- ci(roc_result_test)
  
  # Store the ROC curve in the list
  roc_list_others[[i]] <- list(roc_results = roc_result_test, ci = ci_results)
}
