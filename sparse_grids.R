# GENERATING SMOLYAK SPARSE GRIDS VIA CLENSHAW-CURTIS QUADRATURE


# FUNCTION FOR 1D CLENSHAW-CURTIS QUADRATURE POINTS AND WEIGHTS
QUAD_C.C <- function(l, a=-1, b=1){
  
  # NUMBER OF POINTS
  if(l == 1){
    m <- 1
  }else{
    m <- 2^(l-1)+1
  }
  
  N <- m-1
  j <- 0:N
  
  # STANDARD QUADRATURE POINTS AND WEIGHTS
  if(N == 0){
    X <- 0
    W <- 2
  } else {
    X <- round(cos(j*pi/N), 10)
    
    # VECTORIZED WEIGHT CALCULATION
    N.round <- floor(N/2)
    k <- 1:N.round
    
    # Coefficients for k: 2 for all, except 1 if k == N/2
    b_k <- rep(2, N.round)
    if(N.round == N/2){b_k[N.round] <- 1}
    
    # The k-dependent term: b_k / (4*k^2 - 1)
    v_k <- b_k / (4*k^2 - 1)
    
    # Matrix of cosine terms using outer product: cos(2 * pi * k * j / N)
    cos_matrix <- cos(2 * pi * outer(k, j, "*") / N)
    
    # Sum over k for each j (crossprod is highly optimized in base R)
    SUM_VAL <- as.numeric(crossprod(v_k, cos_matrix)) 
    
    # c vector: 1 for j=0 and j=N, 2 otherwise
    c_j <- rep(2, N+1)
    c_j[1] <- 1
    c_j[N+1] <- 1
    
    # Final weights
    W <- (c_j / N) * (1 - SUM_VAL)
  }
  
  # SCALING POINTS AND WEIGHTS TO CUSTOM DOMAIN [a, b]
  if(a != -1 || b != 1){
    SCALE <- (b-a)/2
    MEAN_VAL <- (b+a)/2
    X <- SCALE * X + MEAN_VAL
    W <- SCALE * W
  }
  
  return(rbind(X, W))
}

# RECURSIVE FUNCTION TO GENERATE VALID MULTI-INDICES
INDICES <- function(d, target_sum){
  
  # BASE CASE: THE SUM IS THE INDEX ITSELF FOR 1D
  if(d == 1){return(matrix(target_sum, nrow = 1))}
  
  # EMPTY MATRIX TO STORE VALID COMBINATIONS
  result <- matrix(ncol=d, nrow=0) 
  
  # MAXIMUM VALUE EACH DIMENSION CAN TAKE (ensuring remaining dims get at least 1)
  max_val <- target_sum - (d - 1) 
  
  # RETURN EMPTY MATRIX IF IT'S IMPOSSIBLE TO ASSIGN AT LEAST 1 TO EACH DIMENSION
  if(max_val < 1){return(result)} 
  
  for(i in 1:max_val){
    # GENERATE COMBINATIONS FOR REMAINING DIMENSIONS AND SUM
    sub_result <- INDICES(d-1, target_sum-i) 
    # BIND RESULTS TOGETHER IF VALID
    if(nrow(sub_result) > 0){result <- rbind(result, cbind(i, sub_result))} 
  }
  return(result)
}

# GENERATING SMOLYAK SPARSE GRIDS
SPARSE.GRID <- function(k, d, a = rep(-1, d), b = rep(1, d)){
  
  # STEP 1: GENERATE ALL VALID MULTI-INDICES AND SMOLYAK COEFFICIENTS
  L <- matrix(ncol = d, nrow = 0)
  Cl <- numeric()
  
  for(sum_val in (k-d+1):k){
    valid_indices <- INDICES(d, sum_val) 
    L <- rbind(L, valid_indices)
    
    # Calculating the Smolyak coefficient for inclusion/exclusion
    coef <- ((-1)^(k-sum_val)) * choose(d-1, k-sum_val)
    Cl <- c(Cl, rep(coef, nrow(valid_indices)))
  }
  
  # STEP 2: TENSOR PRODUCT
  all_sub_grids <- list()
  for(i in 1:nrow(L)){
    list_X <- list()
    list_W <- list()
    
    for(j in 1:d){
      quad <- QUAD_C.C(L[i, j])
      list_X[[j]] <- quad[1, ]
      list_W[[j]] <- quad[2, ]
    }
    
    # Calculating the full tensor product for spatial coordinates and weights
    PT.X <- as.matrix(expand.grid(list_X))
    PT.W <- as.matrix(expand.grid(list_W))
    
    # Multiplying corresponding weights and the Smolyak coefficient
    P.W <- apply(PT.W, 1, prod) * Cl[i]
    
    # Binding X coordinates with corresponding final weights
    all_sub_grids[[i]] <- cbind(PT.X, W = P.W)
  }
  
  # STEP 3: AGGREGATION AND CLEANUP (Standard domain [-1, 1]^d)
  # Uniting all points and weights into a single data frame
  df_complete <- as.data.frame(do.call(rbind, all_sub_grids))
  
  # Renaming columns (X1, X2, ..., Xd, W)
  col_names <- c(paste0("X", 1:d), "W")
  colnames(df_complete) <- col_names
  
  # ROUNDING (Crucial for accurately identifying overlapping points)
  df_complete[, 1:d] <- round(df_complete[, 1:d], 10)
  
  # AGGREGATING WEIGHTS FOR OVERLAPPING POINTS
  final_grid <- aggregate(W ~ ., data = df_complete, FUN = sum)
  
  # REMOVING ZERO WEIGHTS (Sparsity filter)
  final_grid <- final_grid[abs(final_grid$W) > 1e-12, ]
  
  # STEP 4: DOMAIN MAPPING (From [-1, 1]^d to generic [a, b]^d)
  if (!all(a == -1) || !all(b == 1)) {
    scale_vec <- (b - a) / 2
    mean_vec  <- (b + a) / 2
    
    for(j in 1:d){
      final_grid[, j] <- scale_vec[j] * final_grid[, j] + mean_vec[j]
    }
    
    # Adjusting weights by the determinant of the Jacobian matrix
    jacobian <- prod(scale_vec)
    final_grid$W <- final_grid$W * jacobian
  }
  
  return(final_grid)
}
