# Install necessary packages from CRAN or GitHub
# install.packages("remotes")
# remotes::install_github("privefl/bigsnpr")

library(bigsnpr)  # Loads bigsnpr for large-scale SNP data analysis
# install.packages('bigreadr')
# library(bigreadr)  
library(tidyverse)

# Load the LD reference map data
map_ldref <- readRDS('lassosum2_tutorial_data/map_hm3_plus.rds')  # Read in the LD reference data for SNP matching

library(data.table)
# Load the summary statistics for breast cancer
sumstats <- fread('lassosum2_tutorial_data/breast_cancer_summary_statistics.txt')

# Compute the effective sample size based on the provided sample sizes
sumstats$n_eff <- 4 / (1 / 137045 + 1 / 119078)

# Match the SNPs in the summary statistics with the LD reference map
info_snp <- snp_match(sumstats, map_ldref)
info_snp <- as_tibble(info_snp)  # Convert to tibble format

# Add columns for standard deviations based on LD reference allele frequencies and summary statistics
info_snp <- info_snp %>%
  mutate(sd_ldref = sqrt(2 * af_UKBB * (1 - af_UKBB)),  # Standard deviation from the reference allele frequency
         sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2))  # Standard deviation from summary statistics

# Add a column indicating whether a SNP should be removed based on certain conditions
info_snp <- info_snp %>%
  mutate(is_bad = case_when(
    sd_ss < 0.5 * sd_ldref ~ TRUE,
    sd_ss > sd_ldref + 0.1 ~ TRUE,
    sd_ss < 0.05 ~ TRUE,
    sd_ldref < 0.05 ~ TRUE,
    TRUE ~ FALSE
  ))

# Set the ggplot2 theme
library(ggplot2)
theme_set(theme_minimal())  # Use minimal theme
theme_update(axis.title = element_text(size = 13, face = 'bold'),
             axis.text = element_text(size = 12),
             legend.position = 'top',
             legend.justification = 'center',
             legend.title = element_text(size = 12, face = 'bold'),
             legend.text = element_text(size = 11),
             panel.grid = element_line(color = "gray93", size = 0.3))

# Create a scatter plot to visualize SNPs with problematic characteristics
g1 <- ggplot(info_snp, aes(x = sd_ldref, y = sd_ss, color = is_bad)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c('#DB5A42', '#E3A587')) +  # Color bad SNPs in a different color
  geom_abline(linetype = 2, color = "black") +  # Add a reference line
  geom_vline(xintercept = 0.05, linetype = 2, color = "black") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "black") +
  coord_equal() +  # Make the axes equal
  labs(x = "Standard deviations derived from \n allele frequencies of the LD reference",
       y = "Standard deviations derived from \n the summary statistics",
       color = "To remove")

# Save the plot as a TIFF file
tiff('plot_isbad_snps.tif', units = 'in', res = 500, compression = 'lzw', width = 5.5, height = 5.5)
g1
dev.off()

# Filter out the "bad" SNPs
df_betas <- info_snp %>% filter(is_bad == FALSE)

# Set chromosome number for further analysis
chr = 21

#### Intersect with SNPs in the test set #####
# Read from bed/bim/fam, it generates .bk and .rds files.
plink_files <- snp_readBed('lassosum2_tutorial_data/data_binary.bed')

# Load test data
test_data <- readRDS('data_binary.rds')

# Prepare SNPs from test data
snps_in_test <- test_data$map %>%
  as_tibble() %>%
  rename(rsid = marker.ID, a0 = allele1, a1 = allele2)
 
dim(snps_in_test) #811 6

# Join betas with SNPs in the test set by (a0, a1, rsid)
df_betas <- df_betas %>%
  inner_join(snps_in_test %>% select(rsid, a0, a1))

dim(df_betas) #408 18

#### Select LD of the SNPs ######

# Get the list of SNPs from the betas dataframe
ind_dfbetas <- df_betas %>% pull(`_NUM_ID_`)  # SNP indices

# Get the list of SNPs from the LD reference map
ind_ldref <- map_ldref %>% mutate(A = row_number()) %>% pull(A) 

# Find common SNPs between betas and LD reference map
ind_snps <- intersect(ind_dfbetas, ind_ldref)

# Read the LD correlation matrix
corr <- readRDS("LD_with_blocks_chr21.rds") #from https://github.com/comorment/ldpred2_ref
corr <- corr[ind_snps, ind_snps]  # Subset correlation matrix based on common SNPs
corr <- as_SFBM(corr, compact = TRUE)  # Convert to a compact format for memory efficiency

# Run the lassosum2 function to compute beta coefficients for polygenic risk scores
#delta parameter by default: L2-regularization c(0.001, 0.01, 0.1, 1)
#nlambda parameter by default: 30 (L1-regularization)
#df_betas dataframe needs $beta, $beta_se, $n_eff
beta_lassosum2 <- snp_lassosum2(corr, df_betas)#, ncores = 5) 

# Extract grid parameters from the lassosum2 result
params2 <- attr(beta_lassosum2, "grid_param")

######## Genotypes ########

# Load genotype data
test_data <- readRDS('lassosum2_tutorial_data/data_binary.rds')

# Obtain genotypes matrix

G <- test_data$genotypes
y <- test_data$fam$affection  # Affection status (case/control)

# Prepare SNP identifiers for selection
rsid <- snps_in_test %>%
  select(rsid) %>%
  mutate(order = row_number())

# Select the SNPs to use from betas
rsid_toselect <- rsid %>% inner_join(df_betas %>% select(rsid, `_NUM_ID_`))

# Copy selected SNPs from genotypes matrix
G <- big_copy(G, ind.col = rsid_toselect$order)
G[]

# Compute the polygenic risk score (PRS) for each individual based on SNPs
pred_grid2 <- big_prodMat(G, beta_lassosum2)
pred_grid2[] #rows= SNP coefficients per grid parameters (columns) (in params2)

# Function to compute the AUC (Area Under the Curve) for each predicted score based on params2
library('pROC')
obtain_AUC <- function(i) {
  config <- pred_grid2[, i]
  AUC <- auc(y, config) %>%
    as_tibble()  # Compute AUC for logistic regression
  result <- tibble('AUC' = AUC$value[1], 'config' = i)
  return(result)
}

# Compute AUC for all comb of lambda delta (params2)
AUC_est <- lapply(1:ncol(pred_grid2), obtain_AUC)
AUC_est <- do.call(rbind, AUC_est) %>%
  as_tibble()

# Add AUC to parameter grid
params2 <- cbind(params2, AUC_est)

# Create AUC plot
g2 <- ggplot(params2, aes(x = lambda, y = AUC, color = as.factor(delta))) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  labs(y = "Area Under the ROC Curve", color = "delta")

# Save the AUC plot as a TIFF file
tiff('plot_AUC.tif', units = 'in', res = 500, compression = 'lzw', width = 5, height = 5)
g2
dev.off()

# Select the best grid of parameters based on AUC
best_grid_lassosum2 <- params2 %>%
  mutate(id = row_number()) %>%
  slice_max(AUC, n = 1, with_ties = FALSE)

# Get the best polygenic risk score (PRS) prediction based on the best grid
best_prs <- pred_grid2[, best_grid_lassosum2$id]

# Prepare data for boxplot of PRS by case/control status
data <- tibble('PRS' = best_prs, 
               'case_control' = y %>% as_tibble() %>%
                 mutate(value = ifelse(value == 2, 'Case', 'Control')) %>% 
                 pull(value)
               )

# Create boxplot comparing PRS between cases and controls
library(ggpubr)
my_comparisons <- list(c("Case", 'Control'))
g3 <- ggplot(data, aes(x = case_control, y = PRS, colour = case_control, fill = case_control)) +
  geom_boxplot(alpha = 0.4) +
  scale_colour_manual(values = c('#F7B378', '#D78170')) +
  scale_fill_manual(values = c('#F7B378', '#D78170')) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  labs(x = 'Status', y = 'Polygenic Risk Score \n (69 SNPs)') +
  theme(legend.position = 'none')

# Save the boxplot as a TIFF file
tiff('boxplot.tif', units = 'in', res = 500, compression = 'lzw', width = 4, height = 4)
g3
dev.off()

# Get the SNPs with non-zero coefficients for PRS calculation
beta_prs <- beta_lassosum2[, best_grid_lassosum2$id]
selected_snps <- beta_prs %>% 
  as_tibble() %>% 
  filter(value != 0) %>% 
  row_number()

# Get the effect sizes (coefficients) of the selected SNPs
coefs <- beta_prs %>% 
  as_tibble() %>% 
  filter(value != 0) %>% 
  pull(value)

# Create a tibble with selected SNP details and their effect sizes
selected_snps_details <- df_betas[selected_snps, ] %>% 
  as_tibble() %>%
  mutate(effect_size = beta_prs[selected_snps]) %>%
  mutate(coef_PRS = coefs)

# Display the top 10 SNPs with the highest effect sizes
selected_snps_details %>% 
  slice_max(abs(coef_PRS), n = 10) %>% 
  select(chr, pos, a0, a1, rsid, coef_PRS)

