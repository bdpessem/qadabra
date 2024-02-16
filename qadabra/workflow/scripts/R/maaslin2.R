library(biomformat)
library(dplyr)
library(Maaslin2)

# Set logging information
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load the input table
print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

# Load the metadata
print("Loading metadata")
metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@params[[1]][["factor_name"]]
target <- snakemake@params[[1]][["target_level"]]
reference <- snakemake@params[[1]][["reference_level"]]
confounders <- snakemake@params[[1]][["confounders"]]
#comments: 
#1. extract random-effects subjectID and/or site from metadata
random.effects <- snakemake@params[[1]][["random_effects"]]

# Harmonize table and metadata samples
print("Harmonizing table and metadata samples...")
samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)

# Create results table and modify result names
row.names(table) <- paste0("F_", row.names(table)) # Append F_ to features to avoid R renaming
table <- t(table[, sample_order])
fixed.effects <- c(covariate, confounders)
specify.reference <- paste(covariate, reference, sep = ",", collapse = ",")



# Run maAsLin
print("Running MaAsLin...")
#comments: 
# 1. random_effects can include the random effects for the model fi. random_effects="site,subject" --> request user to specify the subjectID in the metadata
# 2. normalization is default set to TSS
# 3. transformation is default set to log
# 4. correction method for computing q-value (or adjusted p-value) is set to BH, you can extract the q-value and don't need to perform the calculation later

fit.data <- Maaslin2::Maaslin2(
    input_data=table,
    input_metadata=metadata,
    output=snakemake@output[["out_dir"]],
    fixed_effects=fixed.effects,
    random_effects=random.effects,
    reference=specify.reference,
    min_prevalence=0,
    min_abundance=0,
    plot_scatter=F,
    plot_heatmap=F
)
results <- fit.data$results %>% dplyr::filter(metadata==covariate)
results <- results %>% distinct(feature, .keep_all = TRUE)

row.names(results) <- gsub("^F_", "", results$feature)
results <- results %>% select(-c("feature"))

# Extract coefficient, p-value, and adjusted p-value results
adjusted_p_values <- p.adjust(results$pval, method = "BH")
results$pval_BH_adj <- adjusted_p_values
# results$coef <- results$coef

#comment 4. alternative method 
results$qval <-adjusted_p_values

# Save results to output file
write.table(results, file=snakemake@output[["diff_file"]], sep="\t")
print("Saved differentials!")
