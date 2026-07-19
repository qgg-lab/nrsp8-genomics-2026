# --- 0. Clear Environment ---
rm(list = ls())

# --- 1. Pedigree Matrix (A) - 4x4 ---
# Rows represent animals 1, 2, 3, 4.
A <- matrix(c(1.0,  0.5,  0.0,  0.25,
              0.5,  1.0,  0.0,  0.5,
              0.0,  0.0,  1.0,  0.5,
              0.25, 0.5,  0.5,  1.0),
            nrow = 4, byrow = TRUE)

# --- 2. Design Matrix (Z) - 2x4 ---
# Explicitly defines 2 rows (phenotypes) and 4 columns (animals).
# Animal 2 gets Record 1; Animal 3 gets Record 2.
Z <- matrix(c(0, 1, 0, 0,  # Row 1: Phenotype for Animal 2
              0, 0, 1, 0), # Row 2: Phenotype for Animal 3
            nrow = 2, byrow = TRUE)

# --- 3. Phenotypes (y) and k ---
y <- matrix(c(10, 20), nrow = 2) # 2x1 Column Matrix
k <- 2 # k = lambda in the diagram

# --- 4. The MME Setup ---
# Henderson's Mixed Model Equations: (Z'Z + A^-1 * k) * u = Z'y
A_inv <- solve(A)
ZZ    <- t(Z) %*% Z         # Result: 4x4
Zy    <- t(Z) %*% y         # Result: 4x1

# --- 5. Solving for Breeding Values (u_hat) ---
LHS   <- ZZ + (A_inv * k)
u_hat <- solve(LHS, Zy)

# --- 6. Results ---
cat("--- MME Solutions (u_hat) ---\n")
print(as.vector(u_hat))
