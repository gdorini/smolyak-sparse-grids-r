# Smolyak Sparse Grids in R via Clenshaw-Curtis Quadrature

A lightweight, robust, and highly optimized R implementation of Smolyak's algorithm for generating sparse grids based on Clenshaw-Curtis quadrature. 

This tool is designed to mitigate the **Curse of Dimensionality** in high-dimensional numerical integration, making it exceptionally useful for Statistical Computing, Bayesian Inference, and evaluating complex criteria in pseudo-Bayesian Optimal Designs.

## 🚀 Features
* **Memory Efficient:** Uses a custom recursive algorithm to generate valid multi-indices without relying on massive, memory-heavy permutation matrices.
* **Flexible Domains:** Automatically scales the standard $[-1, 1]^d$ domain to any arbitrary hyperrectangle $[a, b]^d$, adjusting weights via the Jacobian determinant.
* **Zero Dependencies:** Written in Base R. No need to install external packages.

## 🧠 Why Sparse Grids?
Standard tensor-product quadrature rules suffer from the **Curse of Dimensionality**. The number of integration points grows exponentially with the number of dimensions ($O(N^d)$). A full tensor grid guarantees exact integration for polynomials up to a maximum degree in each dimension, which forces the evaluation of highly complex cross-terms that usually contribute very little to the integral of smooth functions. Smolyak's algorithm mitigates this by constructing a sparse subset of these points. Instead of preserving the full tensor product space, it preserves exactness for polynomials of a given **total degree**. This effectively discards computationally expensive, high-order cross-terms, reducing the complexity to $O(N \log(N)^{d-1})$ while maintaining robust accuracy for smooth functions.

To illustrate this massive efficiency gain, consider a **5-dimensional** integration problem ($d=5$) evaluating up to 9 integration points per marginal dimension:
* A standard full tensor grid requires $9^5$ points (**59,049 evaluations**).
* The corresponding Smolyak Sparse Grid (Level $k=8$, $d=5$) achieves total-degree exactness using only **241 points**.

## 💻 Quick Start & Example
Load the `sparse_grids.R` script into your environment. You can call the `SPARSE.GRID(k, d, a, b)` function where:
* `k`: Accuracy level (usually $k > d$).
* `d`: Number of dimensions.
* `a`, `b`: Vectors defining the lower and upper bounds of integration.

### Test: Integrating a 4D Multivariate Normal Density
Let's approximate the volume under a standard multivariate normal distribution bounded by the hypercube $[-1, 1]^4$. The theoretical exact probability is approximately **0.21706**.

```R
source("sparse_grids.R")

# 1. Generate the Sparse Grid (Level 7, 4 Dimensions)
# A full tensor grid reaching the same marginal polynomial exactness 
# would require 9^4 = 6561 points. Our sparse grid uses only 137 points!
grid <- SPARSE.GRID(k = 7, d = 4)

# 2. Define the 4D Normal Density Function
normal_4d <- function(x1, x2, x3, x4) {
  (1 / ((2 * pi)^2)) * exp(-0.5 * (x1^2 + x2^2 + x3^2 + x4^2))
}

# 3. Evaluate the function at the sparse grid points
f_values <- mapply(normal_4d, grid$X1, grid$X2, grid$X3, grid$X4)

# 4. Calculate the Integral (Sum of f(x) * weights)
integral_approx <- sum(f_values * grid$W)

print(integral_approx)
# Output: 0.2177273
# Exact Theoretical Value: 0.21706...
# Result: Highly accurate approximation using only ~2% of the computational cost!



