####Maxent in R for Azadirachta indica A. Juss#########
options(java.parameters = "-Xmx16g")  # or higher depending on your RAM

# Load required libraries
library(terra)
library(caret)
library(dismo)
library(rJava)
library(plotly)
library(dplyr)
library(tidyr)
library(htmlwidgets)  # For saving interactive plots

# --- Define Paths ---
shape_path <- "E:/Maxent in r/Agroclimatic_regions/Agroclimatic_regions.shp"
env_dir <- "E:/Maxent in r/envdata"
soil_dir <- "E:/Maxent in r/soil rasters"
output_dir <- "E:/Maxent in r/outputs"

# --- List raster files ---
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
soil_files <- list.files(soil_dir, pattern = "\\.tif$", full.names = TRUE)

# --- Load shapefile and reproject ---
agro_shape <- vect(shape_path)
agro_shape_4326 <- project(agro_shape, "EPSG:4326")

# --- Function to reproject raster to EPSG:4326 ---
reproject_raster_to_4326 <- function(raster_path) {
  r <- rast(raster_path)
  if (!compareCRS(r, "EPSG:4326")) {
    r <- project(r, "EPSG:4326")
  }
  r
}

# --- Reproject environmental and soil rasters ---
env_rasters_4326 <- lapply(env_files, reproject_raster_to_4326)
soil_rasters_4326 <- lapply(soil_files, reproject_raster_to_4326)

# --- Function to check overlap of rasters with shapefile extent ---
check_overlap <- function(raster_list, shape_vect, file_paths) {
  report <- data.frame(file = character(),
                       extent_overlap = logical(),
                       stringsAsFactors = FALSE)
  for (i in seq_along(raster_list)) {
    r <- raster_list[[i]]
    ext_shape <- ext(shape_vect)
    ext_r <- ext(r)
    overlap <- !(ext_shape$xmax < ext_r$xmin ||
                   ext_shape$xmin > ext_r$xmax ||
                   ext_shape$ymax < ext_r$ymin ||
                   ext_shape$ymin > ext_r$ymax)
    report <- rbind(report,
                    data.frame(file = basename(file_paths[i]),
                               extent_overlap = overlap,
                               stringsAsFactors = FALSE))
  }
  report
}

# --- Check overlap for environmental and soil rasters ---
env_overlap_4326 <- check_overlap(env_rasters_4326, agro_shape_4326, env_files)
cat("Environmental Raster Overlap Report (EPSG:4326):\n")
print(env_overlap_4326)

soil_overlap_4326 <- check_overlap(soil_rasters_4326, agro_shape_4326, soil_files)
cat("Soil Raster Overlap Report (EPSG:4326):\n")
print(soil_overlap_4326)

# --- Filter valid raster files with overlap = TRUE ---
valid_env_files_4326 <- env_overlap_4326$file[env_overlap_4326$extent_overlap == TRUE]
valid_env_paths_4326 <- file.path(env_dir, valid_env_files_4326)
valid_soil_files_4326 <- soil_overlap_4326$file[soil_overlap_4326$extent_overlap == TRUE]
valid_soil_paths_4326 <- file.path(soil_dir, valid_soil_files_4326)

# --- Function to align, crop, resample, and mask rasters based on shapefile ---
align_crop_mask_4326 <- function(raster_paths, shape_vect) {
  aligned_list <- list()
  for (f in raster_paths) {
    r <- rast(f)
    if (!compareCRS(r, shape_vect)) {
      r <- project(r, crs(shape_vect))
    }
    
    r_crop <- tryCatch({
      crop(r, ext(shape_vect))
    }, error = function(e) {
      message(paste("Crop failed for", basename(f), ":", e$message))
      return(NULL)
    })
    if (is.null(r_crop)) next
    
    template <- rast(ext(shape_vect),
                     resolution = res(r_crop),
                     crs = crs(shape_vect))
    
    r_resample <- resample(r_crop, template, method = "bilinear")
    r_mask <- mask(r_resample, shape_vect)
    
    aligned_list <- c(aligned_list, list(r_mask))
  }
  if (length(aligned_list) == 0) stop("No rasters were successfully aligned.")
  rast(aligned_list)
}

# --- Align environmental and soil rasters ---
env_stack_aligned_4326 <- align_crop_mask_4326(valid_env_paths_4326, agro_shape_4326)
soil_stack_aligned_4326 <- align_crop_mask_4326(valid_soil_paths_4326, agro_shape_4326)

# --- Determine combined extent and target resolution ---
combined_extent <- union(ext(env_stack_aligned_4326), ext(soil_stack_aligned_4326))
res_env <- res(env_stack_aligned_4326)
res_soil <- res(soil_stack_aligned_4326)
target_res <- c(min(res_env[1], res_soil), min(res_env, res_soil))

# --- Define template raster for resampling ---
template <- rast(ext=combined_extent,
                 resolution=target_res,
                 crs=crs(env_stack_aligned_4326))

# --- Resample to common template ---
env_resampled <- resample(env_stack_aligned_4326, template, method="bilinear")
soil_resampled <- resample(soil_stack_aligned_4326, template, method="bilinear")

# --- Combine raster stacks ---
combined_stack_4326 <- c(env_resampled, soil_resampled)

# --- Sample points for correlation ---
set.seed(123)
sample_points <- spatSample(combined_stack_4326, size = 10000, method = "random", as.points = TRUE)

# --- Extract raster values at sampled points, remove NA ---
sample_values <- terra::extract(combined_stack_4326, sample_points)
sample_df <- as.data.frame(sample_values)
df <- na.omit(sample_df)

# --- Correlation matrix and filtering ---
if ("ID" %in% colnames(df)) {
  df_corr <- df[, !(colnames(df) == "ID")]
} else {
  df_corr <- df
}
cor_matrix <- cor(df_corr, method = "pearson")
high_corr_indices <- findCorrelation(cor_matrix, cutoff = 0.7, verbose = TRUE)
vars_to_remove <- colnames(cor_matrix)[high_corr_indices]
vars_to_keep <- setdiff(names(combined_stack_4326), vars_to_remove)

cat("Variables to remove due to correlation > 0.7:\n"); print(vars_to_remove)
cat("Variables kept after correlation filtering:\n"); print(vars_to_keep)

# --- Create filtered raster stack ---
selected_stack <- combined_stack_4326[[vars_to_keep]]

# --- Save filtered raster stack ---
writeRaster(selected_stack,
            filename = file.path(output_dir, "filtered_predictors_0.7cutoff.tif"),
            overwrite = TRUE)
cat("Filtered raster stack (cutoff 0.7) saved.\n")

# --- Plot correlation heatmap using Plotly ---
plotly_cor_heatmap <- plot_ly(
  z = cor_matrix,
  x = colnames(cor_matrix),
  y = rownames(cor_matrix),
  type = "heatmap",
  colorscale = "RdBu",
  zmin = -1,
  zmax = 1,
  reversescale = TRUE,
  colorbar = list(title = "Correlation")
) %>%
  layout(
    title = "Correlation Matrix Heatmap (Plotly Interactive)",
    xaxis = list(title = "", tickangle = -45),
    yaxis = list(title = "")
  )

heatmap_html_path <- file.path(output_dir, "correlation_heatmap_plotly.html")
htmlwidgets::saveWidget(plotly_cor_heatmap, heatmap_html_path)
cat("Interactive correlation heatmap saved to:", heatmap_html_path, "\n")

# --- Alternative heatmap using heatmaply ---
if (!requireNamespace("heatmaply", quietly = TRUE)) {
  install.packages("heatmaply")
}
library(heatmaply)

heatmaply_cor(cor_matrix,
              xlab = "Variables",
              ylab = "Variables",
              main = "Correlation Matrix Heatmap",
              file = file.path(output_dir, "correlation_heatmap_heatmaply.html"))
cat("Interactive heatmaply correlation heatmap saved as HTML.\n")

# --- Load and process occurrence data ---
occ_data <- read.csv("E:/Maxent in r/CSV-trial.csv", stringsAsFactors = FALSE)
standardize_names <- function(names_vec) {
  names_vec <- trimws(names_vec)
  names_vec <- gsub(" ", ".", names_vec, fixed = TRUE)
  names_vec <- gsub("__", "_", names_vec, fixed = TRUE)
  names_vec
}
colnames(occ_data) <- standardize_names(colnames(occ_data))

occ_points <- vect(occ_data, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# --- Extract environmental values for occurrence points ---
occ_env_values <- terra::extract(selected_stack, occ_points)
occ_env_df <- cbind(occ_data, occ_env_values)
if (!"presence" %in% colnames(occ_env_df)) occ_env_df$presence <- 1

# --- Generate background points ---
num_bg <- 10000
bg_points <- spatSample(agro_shape_4326, size = num_bg, method = "random")
bg_env_values <- terra::extract(selected_stack, bg_points)
bg_df <- data.frame(bg_env_values)
bg_df$presence <- 0

colnames(bg_df) <- standardize_names(colnames(bg_df))
names(selected_stack) <- standardize_names(names(selected_stack))

# --- Select common predictors and prepare modeling data ---
vars_common <- intersect(names(selected_stack), intersect(colnames(occ_env_df), colnames(bg_df)))
occ_subset <- occ_env_df[, c(vars_common, "presence")]
bg_subset <- bg_df[, c(vars_common, "presence")]
sdm_data <- rbind(occ_subset, bg_subset)

# --- Split into training and test sets ---
set.seed(123)
train_index <- createDataPartition(sdm_data$presence, p = 0.75, list = FALSE)
train_data <- sdm_data[train_index, ]
test_data <- sdm_data[-train_index, ]

# --- Train Maxent model ---
maxent_model <- maxent(x = train_data[, vars_common],
                       p = train_data$presence,
                       args = c("responsecurves=true"))

# --- Evaluate model performance ---
eval_results <- evaluate(p = test_data[test_data$presence == 1, vars_common],
                         a = test_data[test_data$presence == 0, vars_common],
                         model = maxent_model)
cat("Test AUC:", eval_results@auc, "\n")

# --- Prepare raster stack with selected predictors ---
selected_predictors_stack <- selected_stack[[vars_common]]
stopifnot(all(names(selected_predictors_stack) == vars_common))

# --- Predict suitability and save raster ---
prediction_path <- file.path(output_dir, "maxent_predicted_suitability.tif")
prediction_raster <- predict(selected_predictors_stack, maxent_model,
                             filename = prediction_path,
                             overwrite = TRUE,
                             na.rm = TRUE)

# --- Save model results and correlation matrix ---
write.csv(as.data.frame(maxent_model@results), 
          file = file.path(output_dir, "maxent_model_results.csv"), 
          row.names = TRUE)
write.csv(cor_matrix, file = file.path(output_dir, "correlation_matrix.csv"), row.names = TRUE)
cat("All outputs saved in:", output_dir, "\n")

# ---------- 1. ROC Curve ----------
roc_df <- data.frame(
  FalsePositiveRate = eval_results@FPR,
  TruePositiveRate = eval_results@TPR
)

roc_plot <- plot_ly(roc_df, x = ~FalsePositiveRate, y = ~TruePositiveRate, type = 'scatter', mode = 'lines',
                    line = list(color = 'darkorange')) %>%
  layout(title = paste0('ROC Curve (AUC = ', round(eval_results@auc, 3), ')'),
         xaxis = list(title = "False Positive Rate (1 - Specificity)"),
         yaxis = list(title = "True Positive Rate (Sensitivity)"),
         shapes = list(
           list(type = "line", x0 = 0, x1 = 1, y0 = 0, y1 = 1,
                line = list(dash = 'dash', color = "gray"))
         ))

roc_html <- file.path(output_dir, "maxent_ROC_curve.html")
htmlwidgets::saveWidget(roc_plot, roc_html)
cat("ROC curve plot saved to:", roc_html, "\n")

#-------2. Variable importance and permutation---------
contrib_rows <- grep("\\.contribution$", rownames(maxent_model@results), value = TRUE)
var_imp_contrib <- data.frame(
  Variable = gsub("\\.contribution$", "", contrib_rows),
  Importance = as.numeric(maxent_model@results[contrib_rows, 1])
)
var_imp_contrib <- var_imp_contrib[!is.na(var_imp_contrib$Importance) & var_imp_contrib$Importance > 0, ]
var_imp_contrib <- var_imp_contrib[order(var_imp_contrib$Importance, decreasing = TRUE), ]

perm_rows <- grep("\\.permutation.importance$", rownames(maxent_model@results), value = TRUE)
var_imp_perm <- data.frame(
  Variable = gsub("\\.permutation.importance$", "", perm_rows),
  Importance = as.numeric(maxent_model@results[perm_rows, 1])
)
var_imp_perm <- var_imp_perm[!is.na(var_imp_perm$Importance) & var_imp_perm$Importance > 0, ]
var_imp_perm <- var_imp_perm[order(var_imp_perm$Importance, decreasing = TRUE), ]

importance_combined <- merge(var_imp_contrib, var_imp_perm, by = "Variable", all = TRUE,
                             suffixes = c("_Contribution", "_Permutation"))

write.csv(importance_combined, file.path(output_dir, "maxent_variable_importance_combined.csv"),
          row.names = FALSE)

# --- Plot Contribution Importance ---
plot_contrib <- plot_ly(var_imp_contrib,
                        x = ~Importance,
                        y = ~reorder(Variable, Importance),
                        type = 'bar',
                        orientation = 'h',
                        marker = list(color = 'steelblue')) %>%
  layout(title = "Maxent Variable Importance (% Contribution)",
         xaxis = list(title = "Contribution (%)"),
         yaxis = list(title = "Variable"))

htmlwidgets::saveWidget(plot_contrib, file.path(output_dir, "maxent_variable_importance_contribution.html"))
cat("Contribution importance plot saved.\n")

# --- Plot Permutation Importance ---
plot_perm <- plot_ly(var_imp_perm,
                     x = ~Importance,
                     y = ~reorder(Variable, Importance),
                     type = 'bar',
                     orientation = 'h',
                     marker = list(color = 'darkorange')) %>%
  layout(title = "Maxent Variable Importance (Permutation Importance)",
         xaxis = list(title = "Permutation Importance (%)"),
         yaxis = list(title = "Variable"))

htmlwidgets::saveWidget(plot_perm, file.path(output_dir, "maxent_variable_importance_permutation.html"))
cat("Permutation importance plot saved.\n")

# --- Save static PNG for Variable Importance ---
png(filename = file.path(output_dir, "maxent_variable_importance_contribution.png"),
    width = 900, height = 600)
ggplot(var_imp_contrib, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Maxent Variable Importance (% Contribution)",
       x = "Variable", y = "Contribution (%)") +
  theme_minimal()
dev.off()

png(filename = file.path(output_dir, "maxent_variable_importance_permutation.png"),
    width = 900, height = 600)
ggplot(var_imp_perm, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  coord_flip() +
  labs(title = "Maxent Variable Importance (Permutation Importance)",
       x = "Variable", y = "Permutation Importance (%)") +
  theme_minimal()
dev.off()

cat("Static PNG variable importance charts saved.\n")

# ---------- 3. Response Curves ----------
library(dplyr)
library(plotly)
library(htmlwidgets)

generate_response_curve <- function(model, vars_common, data, raster_stack, steps = 100) {
  response_list <- list()
  medians <- apply(data[, vars_common], 2, median, na.rm = TRUE)
  
  for (var in vars_common) {
    var_seq <- seq(min(data[[var]], na.rm = TRUE), max(data[[var]], na.rm = TRUE), length.out = steps)
    newdata <- as.data.frame(matrix(rep(medians, each = steps), nrow = steps))
    colnames(newdata) <- vars_common
    newdata[[var]] <- var_seq
    preds <- predict(model, newdata)
    
    response_list[[var]] <- data.frame(
      Variable = var,
      VarValue = var_seq,
      Suitability = preds
    )
  }
  bind_rows(response_list)
}

res_df <- generate_response_curve(maxent_model, vars_common, train_data, selected_predictors_stack)

res_curve_plot <- plot_ly()
for (varname in unique(res_df$Variable)) {
  dat <- filter(res_df, Variable == varname)
  res_curve_plot <- add_trace(res_curve_plot,
                              x = dat$VarValue,
                              y = dat$Suitability,
                              name = varname,
                              mode = "lines",
                              type = "scatter")
}
res_curve_plot <- res_curve_plot %>%
  layout(title = "Manual Maxent Response Curves",
         xaxis = list(title = "Variable Value"),
         yaxis = list(title = "Predicted Suitability"),
         legend = list(orientation = "h", x = 0, y = -0.2))

response_html_path <- file.path(output_dir, "maxent_manual_response_curves.html")
htmlwidgets::saveWidget(res_curve_plot, response_html_path)

cat("Manual response curves plot saved to:", response_html_path, "\n")

# ---------- 4. Predicted Suitability Histogram ----------
library(terra)
library(plotly)
library(htmlwidgets)

pred_vals <- terra::values(prediction_raster)
pred_vals <- pred_vals[!is.na(pred_vals)]
pred_df <- data.frame(pred_vals)
colnames(pred_df) <- "Suitability"

set.seed(123)
sample_size <- 100000
if (nrow(pred_df) > sample_size) {
  pred_df_sub <- pred_df[sample(nrow(pred_df), sample_size), , drop = FALSE]
} else {
  pred_df_sub <- pred_df
}

pred_hist_plot <- plot_ly(pred_df_sub, x = ~Suitability, type = "histogram", nbinsx = 50,
                          marker = list(color = 'forestgreen')) %>%
  layout(title = paste0("Histogram of Predicted Habitat Suitability (Sampled ", nrow(pred_df_sub), " values)"),
         xaxis = list(title = "Suitability"),
         yaxis = list(title = "Frequency"))

hist_html_path <- file.path(output_dir, "maxent_predicted_suitability_histogram.html")
htmlwidgets::saveWidget(pred_hist_plot, hist_html_path)

cat("Sampled predicted suitability histogram saved to:", hist_html_path, "\n")
