library(rrBLUP)

# Set seed for reproducible results
#set.seed(42)

# ==========================================
# Simulation parameters
# ==========================================
n_animals <- 100
p_markers <- 500
target_h2 <- 0.45  # Define your target heritability here

# Genotypes (-1, 0, 1)
M_raw <- matrix(sample(c(-1, 0, 1), size = n_animals * p_markers, replace = TRUE), nrow = n_animals)
rownames(M_raw) <- paste0("Animal_", 1:n_animals)
colnames(M_raw) <- paste0("SNP_", 1:p_markers)

# Simulate underlying biological components
true_marker_effects <- rnorm(p_markers, mean = 0, sd = 0.5)
true_breeding_values <- as.numeric(M_raw %*% true_marker_effects)

# Calculate empirical genetic variance, then enforce target Ve
true_Vg <- var(true_breeding_values)
target_Ve <- true_Vg * ((1 - target_h2) / target_h2)

# Generate noise. Note: scale() is used here to force the sample variance to be EXACTLY 
# target_Ve, eliminating any random sampling fluctuation in the true_h2 calculation.
raw_noise <- rnorm(n_animals)
environmental_noise <- as.numeric(scale(raw_noise) * sqrt(target_Ve))

# Verify the true simulated heritability (will now perfectly match target_h2)
true_Ve <- var(environmental_noise)
true_h2 <- true_Vg / (true_Vg + true_Ve)

# Final Phenotype = Mean(100) + Genetics + Environment
y_complete <- 100 + true_breeding_values + environmental_noise

# Mask the last 20 animals to create our "New/Unphenotyped" validation group
y <- y_complete
y[(n_animals - 19):n_animals] <- NA 

# 2. Calculate the scaling factor (S)
allele_freq <- (colMeans(M_raw) + 1) / 2
S <- sum(2 * allele_freq * (1 - allele_freq))

# ==========================================
# Method A: G-BLUP, using K
# ==========================================
G <- A.mat(M_raw)
fit_G <- mixed.solve(y = y, K = G)
u_GBLUP <- as.numeric(fit_G$u)

h2_GBLUP <- fit_G$Vu / (fit_G$Vu + fit_G$Ve)

# ==========================================
# Method B: RR-BLUP (aka SNP-BLUP), using Z
# ==========================================
Z_centered <- sweep(M_raw, 2, colMeans(M_raw))
fit_RR <- mixed.solve(y = y, Z = Z_centered)
u_RRBLUP <- as.numeric(Z_centered %*% fit_RR$u)

# Heritability requires scaling by S
Vg_RRBLUP <- S * fit_RR$Vu 
h2_RRBLUP <- Vg_RRBLUP / (Vg_RRBLUP + fit_RR$Ve)

# ==========================================
# Verification of two model equivalence
# ==========================================
cat("--- Heritability validation ---\n")
cat(sprintf("TRUE Simulated h2 : %f\n", true_h2))
cat(sprintf("Estimated G-BLUP  : %f\n", h2_GBLUP))
cat(sprintf("Estimated RR-BLUP : %f\n\n", h2_RRBLUP))

cat("--- Prediction equivalence on training animals (first 20) ---\n")
comparison_train <- data.frame(
  Animal = rownames(M_raw)[1:20],
  True_BV = round(true_breeding_values[1:20], 1),
  GBLUP_EBV = u_GBLUP[1:20],
  RRBLUP_EBV = u_RRBLUP[1:20]
)
print(comparison_train)

cat("\n")

cat("--- Prediction equivalence on validation animals ---\n")
comparison_val <- data.frame(
  Animal = rownames(M_raw)[(n_animals - 19):n_animals],
  True_BV = round(true_breeding_values[(n_animals - 19):n_animals], 1),
  GBLUP_EBV = u_GBLUP[(n_animals - 19):n_animals],
  RRBLUP_EBV = u_RRBLUP[(n_animals - 19):n_animals]
)
print(comparison_val)

cat("\n")

# ==========================================
# Accuracy (correlation b/w TBV and EBV)
# ==========================================
cat("--- Accuracy (correlation b/w TBV and EBV) ---\n")
cat(sprintf("Training data accuracy: %f\n", cor(true_breeding_values[1:(n_animals - 20)], u_GBLUP[1:(n_animals - 20)])))
cat(sprintf("Validation data accuracy: %f\n", cor(true_breeding_values[(n_animals - 19):n_animals], u_GBLUP[(n_animals - 19):n_animals])))
