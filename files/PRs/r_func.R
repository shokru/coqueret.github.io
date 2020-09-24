




rolling <- function(y, x, time, y_mat, TT, t_vec, data){
    #print(c(y,x,time, y_mat, TT))
    time_ind <- which(time == t_vec)                         # Where we are in time
    train <- data %>% 
        dplyr::filter(date <= t_vec[time_ind - y_mat],              # No forward-looking bias!
                      date >= t_vec[time_ind - y_mat - TT + 1])  %>%    
        select(y,x) %>%
        as.matrix()
    test <- data %>%  filter(date == time)
    x_test <- test %>% pull(x) %>% as.numeric()
    y_test <- test %>% pull(y) %>% as.numeric()   
    res <- c()
    if(sum(is.na(train)) > 0 | sd(train[,2]) == 0){               # Incomplete data => OUT!
        res$date <- as.character(data %>%  filter(date == time) %>% pull(date))
        res$rho_y <- NA
        res$sig2_y <- NA
        res$rho_x <- NA
        res$sig2_x <- NA
        res$rho <- NA
        # Realized returns below
        res$r01m <- test$r01m
        res$r03m <- test$r03m
        res$r06m <- test$r06m
        res$r12m <- test$r12m
        res$r24m <- test$r24m
        # Prediction data
        res$y_test <- NA
        res$x_test <- NA
        res$a <- NA
        res$b <- NA
        res$pred <- NA
        res$s_b <- NA
        res$time <- time
        res$TT <- TT
        res$x <- x
        res$y <- y
        res$sd_y <- sd(train[,1])
        res$sd_x <- sd(train[,2])
        res$cor <- NA
    } else {
        res$date <- as.character(data %>%  filter(date == time) %>% pull(date))
        L <- length(train[,1])
        ar_y <- reg_C(train[2:L,1], train[1:(L-1),1])
        ar_x <- reg_C(train[2:L,2], train[1:(L-1),2])
        res$rho_y <- ar_y[2]
        res$sig2_y <- var(train[2:L,1] - ar_y[1] - ar_y[2] * train[1:(L-1),1])
        res$rho_x <- ar_x[2]
        res$sig2_x <- var(train[2:L,2] - ar_x[1] - ar_x[2] * train[1:(L-1),2])
        res$rho <- cor(train[2:L,2] - ar_x[1] - ar_x[2] * train[1:(L-1),2], 
                       train[2:L,1] - ar_y[1] - ar_y[2] * train[1:(L-1),1])
        
        # Realized returns below
        res$r01m <- test$r01m
        res$r03m <- test$r03m
        res$r06m <- test$r06m
        res$r12m <- test$r12m
        res$r24m <- test$r24m
        # Prediction data
        res$y_test <- y_test
        res$x_test <- x_test
        fit <- reg_C(train[,1], train[,2])
        res$a <- fit[1]
        res$b <- fit[2]
        res$pred <- fit[1] + fit[2] * x_test
        resid <- train[,1] - res$pred
        res$s_b <- sqrt(sum(resid^2, na.rm=T)/(TT-2)/sum((train[,2]-mean(train[,2], na.rm=T))^2, na.rm=T))
        res$time <- time
        res$TT <- TT
        res$x <- x
        res$y <- y
        res$sd_y <- sd(train[,1])
        res$sd_x <- sd(train[,2])    
        res$cor <- cor(train[,1], train[,2])
    }
    return(res)
}






errors <- function(y, x, time, y_mat, TT, t_vec, data){
    time_ind <- which(time == t_vec)                         # Where we are in time
    train <- data %>% 
        dplyr::filter(date <= t_vec[time_ind - y_mat],              # No forward-looking bias!
                      date >= t_vec[time_ind - y_mat - TT + 1])  %>%    
        select(y,x) %>%
        as.matrix()
    L <- length(train[,1])
    res <- c()
    if(sum(is.na(train)) > 0 | sd(train[,2]) == 0){               # Incomplete data => OUT!
        res$err_x <- rep(NA, TT)
        res$err_y <- rep(NA, TT)
        res$y <- y
        res$x <- x
        res$TT <- TT
    } else {
        ar_y <- reg_C(train[2:L,1], train[1:(L-1),1])
        ar_x <- reg_C(train[2:L,2], train[1:(L-1),2])
        res$err_y <- train[2:L,1] - ar_y[1] - ar_y[2] * train[1:(L-1),1]
        res$err_x <- train[2:L,2] - ar_x[1] - ar_x[2] * train[1:(L-1),2]
        res$y <- y
        res$x <- x
        res$TT <- TT
    }
    return(res)
}






rolling_ml <- function(y, x, time, y_mat, TT, t_vec, data){
    time_ind <- which(time == t_vec)                         # Where we are in time
    train <- data %>% 
        dplyr::filter(date <= t_vec[time_ind - y_mat],              # No forward-looking bias!
                      date >= t_vec[time_ind - y_mat - TT + 1])  %>%    
        select(y,x) %>%
        as.matrix()
    test <- data %>%  filter(date == time)
    x_test <- test %>% pull(x) %>% as.numeric()
    y_test <- test %>% pull(y) %>% as.numeric()   
    res <- c()
    if(sum(is.na(train)) > 0 | sd(train[,2]) == 0){               # Incomplete data => OUT!
        res$rho_y <- NA
        res$sig2_y <- NA
        res$rho_x <- NA
        res$sig2_x <- NA
        res$rho <- NA
        # Realized returns below
        res$r01m <- test$r01m
        res$r03m <- test$r03m
        res$r06m <- test$r06m
        res$r12m <- test$r12m
        res$r24m <- test$r24m
        # Prediction data
        res$y_test <- NA
        res$x_test <- NA
        res$a <- NA
        res$b <- NA
        res$pred <- NA
        res$time <- time
        res$TT <- TT
        res$x <- x
        res$y <- y
        res$sd_y <- sd(train[,1])
        res$sd_x <- sd(train[,2])
        res$cor <- NA
    } else {
        ar_y <- try(arima(train[,1], c(1,0,0), method = "ML"), silent = T)
        ar_x <- try(arima(train[,2], c(1,0,0), method = "ML"), silent = T)
        if(class(ar_y)!="try-error"){
            res$rho_y <- ar_y$coef[1]
            res$sig2_y <- ar_y$sigma2
        } else {
            res$rho_y <- NA
            res$sig2_y <- NA
        }
        if(class(ar_x)!="try-error"){
            res$rho_x <- ar_x$coef[1]
            res$sig2_x <- ar_x$sigma2
        } else {
            res$rho_x <- NA
            res$sig2_x <- NA
        }        
        if(class(ar_y)!="try-error" & class(ar_x)!="try-error"){
            res$rho <- cor(ar_y$residuals, ar_x$residuals)
        } else {
            res$rho <- NA
        }
        
        # Realized returns below
        res$r01m <- test$r01m
        res$r03m <- test$r03m
        res$r06m <- test$r06m
        res$r12m <- test$r12m
        res$r24m <- test$r24m
        # Prediction data
        res$y_test <- y_test
        res$x_test <- x_test
        fit <- reg_C(train[,1], train[,2])
        res$a <- fit[1]
        res$b <- fit[2]
        res$pred <- fit[1] + fit[2] * x_test
        res$time <- time
        res$TT <- TT
        res$x <- x
        res$y <- y
        res$sd_y <- sd(train[,1])
        res$sd_x <- sd(train[,2])    
        res$cor <- cor(train[,1], train[,2])
    }
    return(res)
}

