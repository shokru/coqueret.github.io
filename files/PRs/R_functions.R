# SCRIPT that gathers R functions

# Simple Toeplitz covariance matrix (auocorrelation of AR(1))
S <- function(s, r, TT){
    v <- r^(0:TT)                       # First line of scaled cov. matrix
    toeplitz(v) * s^2/(1-r^2)           # Cov matrix, T+1 elements
}

# Covariance matrix Sigma_xy
Sig_xy <- function(r_x, r_y, TT, k){
    
    idx <- 1:TT
    
    y_idx <- matrix(rep(idx, TT), ncol=TT)
    x_idx <- t(y_idx)
    
    Sig_xy_T <- (r_x^(x_idx-y_idx-k))*((x_idx-y_idx-k)>0) + (r_y^(y_idx+k-x_idx))*((y_idx+k-x_idx)>=0)
    
    Sig_xy_last_col <- r_x^(TT-idx)
    Sig_xy_last_row <- r_y^(2*k+TT-idx)
    
    sig_xy_last <- r_y^k
    
    rbind(cbind(Sig_xy_T, Sig_xy_last_col), c(Sig_xy_last_row, sig_xy_last))  # Full cov matrix
}

# Covariance matrix Sigma_x or Sigma_y
Sig <- function(r, TT, k){
    
    Top_T <- toeplitz(r^(0:(TT-1)))                         # First line of scaled cov. matrix
    eps <- r^((k+TT-1):k)

    rbind(cbind(Top_T, eps), c(eps, 1))
}

# BELOW: functions for integral computation - They are only valid if cor = 0 !!!
steps <- function(n, STOP, TT){  # Integral steps
    alpha <- STOP^(2/(TT-1))/n
    (alpha * 0:n)^((TT-1)/2)
}

disc_int_ind <- function(func, n, STOP, lambda_T, p, q, Q, s_x_bar){ # This uses trapezoidal approx for integrals
    TT <- length(lambda_T)
    points <- steps(n, STOP, TT)
    values <- func(points, lambda_T, p, q, Q, s_x_bar)               # Calls term_x
    trapez <- values[1:(n-1)] + values[2:n]
    delta <- points[2:n] - points[1:(n-1)]
    return(sum(delta*trapez)/2)
}

Delta_func1_ind <- function(tt, lambda_T, p, q, Q, s_x_bar){
    
    det_Delta <- 1
    tr_AC_x2  <- 0
    
    for(j in 1:length(lambda_T)){
        det_Delta <- det_Delta / sqrt(1+2*tt*lambda_T[j])
        tr_AC_x2  <- tr_AC_x2  + p[j]*q[j]/(1+2*tt*lambda_T[j])
    }
    return(det_Delta*tr_AC_x2)
}

Delta_func2_ind <- function(tt, lambda_T, p, q, Q, s_x_bar){

    det_Delta <- 1
    tr_AA_x2  <- 0
    tr_D      <- 0
    tr_AAD_x4 <- 0
  
    for(j in 1:length(lambda_T)){
        det_Delta <- det_Delta / sqrt(1+2*tt*lambda_T[j])
        tr_AA_x2  <- tr_AA_x2 + Q[j,j]*lambda_T[j]/(1+2*tt*lambda_T[j])
        tr_D      <- tr_D + p[j]^2/(1+2*tt*lambda_T[j])
        
        for(i in 1:length(lambda_T)){
        tr_AAD_x4 <- tr_AAD_x4 + Q[i,j]*p[i]*p[j]/(1+2*tt*lambda_T[i])/(1+2*tt*lambda_T[j])
        }
    }
    tr_D <- s_x_bar - 2*tt*tr_D
    return(det_Delta*tt*(tr_AA_x2*tr_D+2*tr_AAD_x4))
}

tail_int_Delta_func1 <- function(tt, lambda_T, p, q){
    TT <- length(lambda_T)
    sum(p[1:(TT-1)]*q[1:(TT-1)]/lambda_T[1:(TT-1)])/(2*tt)^((TT-1)/2)/sqrt(prod(lambda_T[1:(TT-1)]))/(TT-1)
}

tail_int_Delta_func2 <- function(tt, lambda_T, p, Q, s_x_bar){
    TT <- length(lambda_T)
    (s_x_bar-sum(p[1:(TT-1)]^2/lambda_T[1:(TT-1)]))*sum(diag(Q[1:(TT-1), 1:(TT-1)]))/(2*tt)^((TT-3)/2)/(2*sqrt(prod(lambda_T[1:(TT-1)])))/(TT-3)
}

disc_int <- function(func, n, STOP, lambda, SS_bar, A, B, C, D){ # This uses trapezoidal approx for integrals
    TT <- length(lambda)
    points <- steps(n, STOP, TT)
    values <- func(points, lambda, SS_bar, A, B, C, D)               # Calls term_x
    trapez <- values[1:(n-1)] + values[2:n]
    delta <- points[2:n] - points[1:(n-1)]
    return(sum(delta*trapez)/2)
}

Delta_func1_func2 <- function(tt, lambda, SS_bar, A, B, C, D){
    
    det_Delta = c()
    f1 = c()
    f2 = c()
    
    A_SS_bar <- A %*% SS_bar
    B_SS_bar <- B %*% SS_bar
    C_SS_bar <- C %*% SS_bar
    D_SS_bar <- D %*% SS_bar
    
    LL <- length(lambda)
    for(j in 1:length(tt)){
        # Compute the determinant of matrix Delta
        det_Delta[j] <- prod((1+2*tt[j]*lambda)^-0.5)

        # traces involved in functions f1 and f2
        W <- solve(diag(LL)+2*tt[j]*B_SS_bar)
        A_star <- A_SS_bar %*% W
        C_star <- C_SS_bar %*% W
        D_star <- D_SS_bar %*% W
        
        tr_A   <- sum(diag(A_star))
        tr_C   <- sum(diag(C_star))
        tr_AC  <- sum(diag(A_star %*% C_star))
        tr_AA  <- sum(diag(A_star %*% A_star))
        tr_D   <- sum(diag(D_star))
        tr_AAD <- sum(diag(A_star %*% A_star %*% D_star))
        tr_AD  <- sum(diag(A_star %*% D_star))
        
        # functions f1 and f2
        f1[j] <- tr_A*tr_C + 2*tr_AC
        f2[j] <- 2*tr_AA*tr_D + 8*tr_AAD + tr_A^2*tr_D + 4*tr_A*tr_AD
    }
    return(det_Delta*(tt*f2-2*f1))
}

MSE_integral <- function(r_x, r_y, s_y, cor, TT, k, n_points, STOP){

    # Build the general covariance matrix of z=(x,y) with unit variances
    S_xy <- cor/(1-r_x*r_y) * Sig_xy(r_x, r_y, TT, k)
    S_x  <- 1/(1-r_x^2) * Sig(r_x, TT, k)
    S_y  <- 1/(1-r_y^2) * Sig(r_y, TT, k)
    SS   <- rbind(cbind(S_x, t(S_xy)), cbind(S_xy, S_y))  # Full cov matrix
    
    # Build the matrices used in the paper for Th 1
    M_T <- diag(TT) - matrix(1, ncol = TT, nrow = TT) / TT
    N_T <- rbind(cbind(M_T,0),c(rep(-1/TT,TT),1))
    zero_T <- matrix(0, ncol = TT+1, nrow = TT+1)       # Matrix filled with zeros
    N <- rbind(cbind(N_T,zero_T), cbind(zero_T,N_T))
    
    J_T <- rbind(cbind(diag(TT),0),0)
    K_T <- diag(TT+1) - J_T
    
    A <- 1/2*rbind(cbind(zero_T,J_T), cbind(J_T,zero_T))
    B <- rbind(cbind(J_T,zero_T), cbind(zero_T,zero_T))
    C <- 1/2*rbind(cbind(zero_T,K_T),cbind(K_T,zero_T))
    D <- rbind(cbind(K_T,zero_T), cbind(zero_T,zero_T))
    E <- rbind(cbind(zero_T, zero_T), cbind(zero_T,K_T))
    
    # Build SS_bar and W needed to apply Th1 in Coqueret & Deguest
    SS_bar <- N %*% SS %*% t(N)
    
    # Diagonalization of B*SS_bar*B
    lambda <- eigen(B %*% SS_bar %*% B)$values
    
    # First term
    tr_E <- sum(diag(E %*% SS_bar))
    
    # Integral approximations
    t2_t3 <- disc_int(Delta_func1_func2, n_points, STOP, lambda, SS_bar, A, B, C, D)
    return(s_y^2 * (tr_E + t2_t3))
}

MSE_integral_ind <- function(r_x, r_y, s_y, TT, k, n_points, STOP){ # cor = 0 !!!
    
    # Build the covariance matrices for x and y with unit variances
    S_x  <- 1/(1-r_x^2) * Sig(r_x, TT, k)
    S_y  <- 1/(1-r_y^2) * Sig(r_y, TT, k)
    
    # Build the matrices used in the paper for Prop 2 & 3
    M_T <- diag(TT) - matrix(1, ncol = TT, nrow = TT) / TT
    N_T <- rbind(cbind(M_T,0),c(rep(-1/TT,TT),1))

    Sig_x_bar <- N_T %*% S_x %*% t(N_T)
    Sig_y_bar <- N_T %*% S_y %*% t(N_T)
    
    Sig_x_bar_T <- Sig_x_bar[1:TT, 1:TT]
    Sig_y_bar_T <- Sig_y_bar[1:TT, 1:TT]
    
    eps_x_bar <- Sig_x_bar[1:TT, TT+1]
    eps_y_bar <- Sig_y_bar[1:TT, TT+1]

    s_x_bar <- Sig_x_bar[TT+1, TT+1]
    s_y_bar <- Sig_y_bar[TT+1, TT+1]
    
    # Diadonalization of Sig_x_bar_T
    lambda_T <- eigen(Sig_x_bar_T)$value
    P_T      <- eigen(Sig_x_bar_T)$vectors
    
    # Define vector p, q and matrix Q
    p <- t(P_T) %*% eps_x_bar
    q <- t(P_T) %*% eps_y_bar
    Q <- t(P_T) %*% Sig_y_bar_T %*% P_T   
    
    # Integral approximations
    t2 <- disc_int_ind(Delta_func1_ind, n_points, STOP, lambda_T, p, q, Q, s_x_bar) + tail_int_Delta_func1(STOP, lambda_T, p, q)
    t3 <- disc_int_ind(Delta_func2_ind, n_points, STOP, lambda_T, p, q, Q, s_x_bar) + tail_int_Delta_func2(STOP, lambda_T, p, Q, s_x_bar) 
    return(s_y^2 * (s_y_bar - 2 * t2 + t3))
}

# BELOW: MSE via large covariance matrix
MSE_simulation <- function(a_x, a_y, r_x, r_y, s_x, s_y, cor, TT, k, n_sim){
    e <- vector("numeric", length = n_sim) # Empty error vector
    mu <- c(rep(a_x/(1-r_x),TT+1), rep(a_y/(1-r_y),TT+1))
    
    # Construction of the covariance matrix of z=(x,y) with real variances
    S_xy <- cor*s_x*s_y/(1-r_x*r_y) * Sig_xy(r_x, r_y, TT, k)
    S_x  <- s_x^2/(1-r_x^2) * Sig(r_x, TT, k)
    S_y  <- s_y^2/(1-r_y^2) * Sig(r_y, TT, k)
    SS   <- rbind(cbind(S_x, t(S_xy)), cbind(S_xy, S_y))  # Full cov matrix
    
    # Below are defined the matrices used in the paper
    M_T <- diag(TT) - matrix(1, ncol = TT, nrow = TT) / TT
    N_T <- rbind(cbind(M_T,0),c(rep(-1/TT,TT),1))
    zero_T <- matrix(0, ncol = TT+1, nrow = TT+1)       # Matrix filled with zeros
    N <- rbind(cbind(N_T,zero_T), cbind(zero_T,N_T))
    
    J_T <- rbind(cbind(diag(TT),0),0)
    K_T <- diag(TT+1) - J_T
    
    A <- 1/2*rbind(cbind(zero_T,J_T), cbind(J_T,zero_T))
    B <- rbind(cbind(J_T,zero_T), cbind(zero_T,zero_T))
    C <- 1/2*rbind(cbind(zero_T,K_T),cbind(K_T,zero_T))
    D <- rbind(cbind(K_T,zero_T), cbind(zero_T,zero_T))
    E <- rbind(cbind(zero_T, zero_T), cbind(zero_T,K_T))
    
    # Below, all random variates are simulated and the loop begins
    set.seed(42) # Freezing the random number generation
    z <- rmvnorm(n=n_sim, mean=mu, sigma=SS, method="chol")
    val1_tmp <- 0  # Initializing...
    val2_tmp <- 0
    val3_tmp <- 0
    for(j in 1:n_sim){
        z_bar <- N %*% z[j,]
        ratio <- (t(z_bar) %*% A %*% z_bar) / (t(z_bar) %*% B %*% z_bar)
        val1_tmp <- (t(z_bar) %*% E %*% z_bar)
        val2_tmp <- (t(z_bar) %*% C %*% z_bar) * ratio
        val3_tmp <- (t(z_bar) %*% D %*% z_bar) * ratio^2
        e[j] <- (val1_tmp - 2*val2_tmp + val3_tmp)
    }
    
    return(c(mean(e), var(e)))
}
