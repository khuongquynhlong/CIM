library(tidyverse)
library(magrittr)
library(knitr)
library(kableExtra)

#-------------------

#---------- Project 2
#===============================================================================
# Question 1
library(glm2)
data(crabs)
names(crabs)

plot(crabs$Width, crabs$Satellites)

head(crabs)



# 1.
crab_mod_1 <- glm(Satellites ~ Width, family = poisson(), data = crabs)
crab_mod_2 <- glm(Satellites ~ Width + GoodSpine, family = poisson(), data = crabs)

# likelihood ratio test
lrt_obs <- -2*(logLik(crab_mod_1)[1] - logLik(crab_mod_2)[1])
# P=value
1 - pchisq(lrt_obs, df = 1)


# 2. 
# Parametric

lambda <- mean(crabs$Satellites)
n <- nrow(crabs)

B <- 1000


# Parametric
lrt_para <- c()
set.seed(2023)
for (i in 1:B) {
    Sat_boot <- rpois(n, lambda)
    crab_boot_1 <- glm(Sat_boot ~ crabs$Width, family = poisson())
    crab_boot_2 <- glm(Sat_boot ~ crabs$Width + crabs$GoodSpine, family = poisson())
    lrt_para[i] <- -2*(logLik(crab_boot_1)[1] - logLik(crab_boot_2)[1])
}



# Non-parametric
lrt_nonpara <- c()
set.seed(2023)

for (i in 1:B) {
    Sat_boot <- sample(crabs$Satellites, n, replace = T)
    crab_boot_1 <- glm(Sat_boot ~ crabs$Width, family = poisson())
    crab_boot_2 <- glm(Sat_boot ~ crabs$Width + crabs$GoodSpine, family = poisson())
    lrt_nonpara[i] <- -2*(logLik(crab_boot_1)[1] - logLik(crab_boot_2)[1])
}


par(mfrow = c(1, 2))
hist(lrt_para, main = "Parametric")
hist(lrt_nonpara, main = "Non-parametric")
par(mfrow = c(1, 1))

# P-values
(1 + sum(lrt_para > lrt_obs))/(1 + B) 
(1 + sum(lrt_nonpara > lrt_obs))/(1 + B)

# 3
# permutaions
lrt_permu <- c()
set.seed(2023)

for (i in 1:B) {
    Sat_boot <- sample(crabs$Satellites, n, replace = F)
    crab_boot_1 <- glm(Sat_boot ~ crabs$Width, family = poisson())
    crab_boot_2 <- glm(Sat_boot ~ crabs$Width + crabs$GoodSpine, family = poisson())
    lrt_permu[i] <- -2*(logLik(crab_boot_1)[1] - logLik(crab_boot_2)[1])
}

hist(lrt_permu)

(1 + sum(lrt_permu > lrt_obs))/(1 + B)


# Question 2

extra <- c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0, 
           1.9, 0.8, 1.1, 0.1, -0.1, 4.4, 5.5, 1.6, 4.6, 4.3)
group <- c(rep(1, 10), rep(2, 10))
ID <- c(1:20)

sleep <- data.frame(extra, group, ID)


# 1
t_obs <- t.test(sleep$extra, sleep$group)$statistic
t_obs

B <- 1000
t_boot <- c()

set.seed(2023)
for (i in 1:B) {
    extra_boot <- sample(sleep$extra, length(sleep$extra), replace = T)
    t_boot[i] <- t.test(extra_boot, sleep$group)$statistic
}

hist(t_boot)
# 2-sided p-value
(1 + sum(abs(t_boot) > abs(t_obs)))/(1 + B)


# 3
psych::describeBy(sleep$extra, group = sleep$group)

B <- 1000
t_boot_m <- c()

mean_gr1 <- mean(sleep$extra[sleep$group == 1])
sd_gr1 <- sd(sleep$extra[sleep$group == 1])
mean_gr2 <- mean(sleep$extra[sleep$group == 2])
sd_gr2 <- sd(sleep$extra[sleep$group == 2])


set.seed(2023)
for (i in 1:B) {
    extra_boot1 <- rnorm(10, mean_gr1, sd_gr1^2)
    extra_boot2 <- rnorm(10, mean_gr2, sd_gr2^2)
    med1 <- median(extra_boot1)
    med2 <- median(extra_boot2)
    SM <- sum(abs(extra_boot1 - med1)) + sum(abs(extra_boot2 - med2))
    t_boot_m[i] <- (med1 - med2)/SM
}


med_gr1 <- median(sleep$extra[sleep$group == 1])
med_gr2 <- median(sleep$extra[sleep$group == 2])
SM_obs <- sum(abs(sleep$extra[sleep$group == 1] - med_gr1)) + 
    sum(abs(sleep$extra[sleep$group == 2] - med_gr2))
tm_obs <- (med_gr1 - med_gr2)/SM_obs


hist(t_boot_m)


# 2-sided p-value
(1 + sum(abs(t_boot_m) > abs(tm_obs)))/(1 + B)


par(mfrow = c(1, 2))
hist(t_boot, main = "Based on t-test statistic")
hist(t_boot_m, main = "Based on robust t-test statistic")
par(mfrow = c(1, 1))

# 5

theta <- mean_gr1 - mean_gr2

theta_hat <- c()
bias <- c()

B <- 1000
set.seed(2023)

for (i in 1:B) {
    extra_boot1 <- rnorm(10, mean_gr1, sd_gr1^2)
    extra_boot2 <- rnorm(10, mean_gr2, sd_gr2^2)
    med1 <- median(extra_boot1)
    med2 <- median(extra_boot2)
    theta_hat[i] <- med1 - med2
}

bias <- theta_hat - theta

hist(bias)

psych::describe(bias)


# Question 3

ID <- c(1:10)
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)
df <- data.frame(ID, x1, x2)

# 1
(ration_obs <- mean(df$x1)/mean(df$x2))

# Bootstrap
B <- 1000
ratio_hat_boot <- c()
n <- length(df$x1)
index <- c(1:n)

set.seed(2023) 

for (i in 1:B) {
    index.b <- sample(index, n, replace = TRUE)
    x1_boot <- df$x1[index.b]
    x2_boot <- df$x2[index.b]
    ratio_hat_boot[i] <- mean(x1_boot)/mean(x2_boot)
}



# Jackknife
ratio_hat_jn <- c()
set.seed(2023)

for (i in 1:n) {
    x1_jn <- df$x1[-i]
    x2_jn <- df$x2[-i]
    ratio_hat_jn[i] <- mean(x1_jn)/mean(x2_jn)
}

# SE
sd(ratio_hat_boot) # Bootstrap
(sum((ratio_hat_jn - mean(ratio_hat_jn))^2) * ((n - 1)/n))^(0.5) # Jackknife

# 3
quantile(ratio_hat_boot, p = c(0.025, 0.975))


# 4
t.test(ratio_hat_boot ~ 1, alternative = "greater")
t.test(ratio_hat_boot ~ 1, alternative = "less")



# Question 4
# 1:  95% C.I for mu_d using the classical method
n <- length(df$x2)
x <- df$x1 - df$x2
mu_d <- mean(x)
se <- sd(x)/sqrt(10)
c(mu_d - qt(0.975, df = 9)*se, mu_d + qt(0.975, df = 9)*se)


# 2: 
# Non-parametric bootstrap
B <- 1000
mu_hat_boot <- c()
n <- length(df$x1)
index <- c(1:n)

set.seed(2023) 
for (i in 1:B) {
    index.b <- sample(index, n, replace = TRUE)
    x1_boot <- df$x1[index.b]
    x2_boot <- df$x2[index.b]
    mu_hat_boot[i] <- mean(x1_boot) - mean(x2_boot)
}

# 95%CI
quantile(mu_hat_boot, prob = c(0.025, 0.975))


# BCa CI
bcboot <- function(x, nboot, theta, alpha = 
                       c(0.025, 0.975)) { 
    n <- length(x)
    thetahat <- theta(x)
    bootsam <- matrix(sample(x, size = n*nboot, replace = TRUE), nrow = nboot)
    thetastar <- apply(bootsam, 1, theta)
    z0 <- qnorm(sum(thetastar < thetahat)/nboot)
    u <- rep(0, n)
    for(i in 1:n){
        u[i] <- theta(x[-i])
    }
    uu <- mean(u)-u
    acc <- sum(uu^3)/(6*(sum(uu^2))^1.5)
    zalpha <- qnorm(alpha)
    tt <- pnorm(z0 + (z0 + zalpha)/(1-acc*(z0 + zalpha)))
    confpoints <- quantile(x = thetastar, probs = tt, type = 1)
    return(confpoints)
}

x <- df$x1 - df$x2

set.seed(2023)
bcboot(x, 1000, mean, alpha=c(0.025, 0.975))


# 3
# Non-parametric bootstrap
B <- 1000
t_boot <- c()
n <- length(df$x1)
x <- c(df$x1, df$x2)

set.seed(2023) 
for (i in 1:B) {
    x_boot <- sample(x, 2*n, replace = T)
    x1 <- x_boot[1:n]
    x2 <- x_boot[(n+1):(2*n)]
    t_boot[i] <- t.test(x1, x2)$statistic
}

hist(t_boot)
t_obs <- t.test(df$x1, df$x2)$statistic
(1 + sum(abs(t_boot) > t_obs))/(1 + B)



#---------- Project 3
#===============================================================================
# Question 1
#------ 1
head(chickwts)

table(chickwts$feed)

lm_mod <- lm(weight ~ feed, data = chickwts)
anova(lm_mod)

Ftest_obs <- anova(lm_mod)$`F value`[1]

#------ 2

Ftest_boot <- c()
B <- 1000
n <- nrow(chickwts)
set.seed(2023)
for (i in 1:B) {
    weight_boot <- sample(chickwts$weight, n, replace = T)
    lm_boot <- lm(weight_boot ~ chickwts$feed)
    Ftest_boot[i] <- anova(lm_boot)$`F value`[1]
}
hist(Ftest_boot)

 # P-value
(1 + sum(Ftest_boot > Ftest_obs))/(B + 1)


#------ 3 
Ftest_per <- c()
B <- 1000
n <- nrow(chickwts)
# Index
nindex <- table(chickwts$feed)
set.seed(2023)
for (i in 1:B) {
    weight_per <- sample(chickwts$weight, n, replace = F)
    gr_feed <- rep(c("casein","horsebean","linseed",
                     "meatmeal","soybean","sunflower"),
                   times = nindex)
    lm_per<- lm(weight_per ~ gr_feed)
    Ftest_per[i] <- anova(lm_per)$`F value`[1]
}

hist(Ftest_per)

# P-value
(1 + sum(Ftest_per > Ftest_obs))/(B + 1)

#------ 4: Semi-parametric

chickwts2 <- chickwts %>% filter(feed %in% c("sunflower", "soybean")) %>%
    mutate(feed01 = ifelse(feed == "sunflower", 1, 0))

lm_mod2 <- lm(weight ~ feed01, data = chickwts2)
lm_res <- lm_mod2$residuals

B <- 1000
n <- nrow(chickwts2)
theta_boot <- c()
set.seed(2023)
for (i in 1:B) {
    e_boot <- sample(lm_res, n, replace = T)
    wt_boot <- coef(lm_mod2)[1] + coef(lm_mod2)[2]*chickwts2$feed01 + e_boot
    lm_boot <- lm(wt_boot ~ chickwts2$feed01)
    theta_boot[i] <- coef(lm_boot)[2]
}

hist(theta_boot)
# 90%CI
quantile(theta_boot, probs = c(0.05, 0.95))

# Question 2
library(Ecdat)

data("Computers")
names(Computers)

#---------- 1
ols_mod <- lm(price ~ hd, data = Computers)
summary(ols_mod)

#---------- 2

price_pred <- predict(ols_mod)
# Error
pred_error <- Computers$price - price_pred
# Mean squared error
(mse1 <- mean(pred_error^2))
# Root Mean squared error
(rmse1 <- sqrt(mse1))

#---------- 3, 4
n <- nrow(Computers)
price_pred_cv <- c() # Q3
beta1_cv <- c() # Q4

for (i in 1:n) {
    price_cv <- Computers$price[-i]
    hd_cv <- Computers$hd[-i]
    ols_cv <- lm(price_cv ~ hd_cv)
    price_pred_cv[i] <- coef(ols_cv)[1] + coef(ols_cv)[2]*Computers$hd[i]
    beta1_cv[i] <- coef(ols_cv)[2]
}


pred_error_cv <- Computers$price - price_pred_cv

(mse_cv <- mean(pred_error_cv^2)); mse1
(rmse_cv <- sqrt(mse_cv)); rmse1


#---------- 4
hist(beta1_cv)

plot(1:n, beta1_cv, type = "l", xlab = "Iteration", ylab = "Beta 1")

#---------- 5
mean(price_pred)
quantile(price_pred, probs = c(0.025, 0.975))
quantile(price_pred_cv, probs = c(0.025, 0.975))


# Question 3: 
#---------- 1: Do bootstrap of CV
n <- nrow(Computers)
B <- 100

beta0_cv <- c(1:n)
index <- c(1:n)
sd_beta0 <- c(1:B)

set.seed(2023)
for (j in 1:B) {
    index.b <- sample(index, n, replace = TRUE)
    price_boot <- Computers$price[index.b]
    hd_boot <- Computers$hd[index.b]
    for (i in 1:n) {
        price_cv <- price_boot[-i]
        hd_cv <- hd_boot[-i]
        ols_cv <- lm(price_cv ~ hd_cv)
        beta0_cv[i] <- coef(ols_cv)[1]
    }
    sd_beta0[j] <- sd(beta0_cv)
}

#---------- 2
price_mod <- lm(price ~ hd, data = Computers)
beta1_obs <- coef(price_mod)[2]
price_res <- residuals(price_mod)


n <- nrow(Computers)
B <- 1000
beta1_boot <- c()
set.seed(2023)
for (i in 1:B) {
    price_res_boot <- sample(price_res, n, replace = T)
    price_boot <- coef(price_mod)[1] + price_res_boot
    price_mod_boot <- lm(price_boot ~ Computers$hd)
    beta1_boot[i] <- coef(price_mod_boot)[2]
}

hist(beta1_boot)

(1 + sum(abs(beta1_boot) > beta1_obs)) / (B + 1)


# Question 4: 

x <- c(0.68446806, -0.02596037,-0.90015774,0.72892605,-0.45612255, 
     0.19311847, -0.13297109, -0.99845382, 0.37278006, -0.20371894, 
     -0.15468803,0.19298230, -0.42755534, -0.04704525, 0.15273726, 
     0.03655799, 0.01315016, -0.59121428,4.50955771, 2.87272653)

#---------- 1
mu_x <- mean(x)
med_x <- median(x)

#---------- 2
B <- 1000
mu_boot <- c()
med_boot <- c()
n <- length(x)

for (i in 1:B) {
    x_boot <- sample(x, n, replace =  T)
    mu_boot[i] <- mean(x_boot)
    med_boot[i] <- median(x_boot)
}

hist(mu_boot)
hist(med_boot)


#---------- 3
se_mu <- sd(mu_boot)
se_med <- sd(med_boot)

# Semiparametric bootstrap
mod_mean <- lm(x ~ 1)
mod_res <- residuals(mod_mean)

n <- length(x)
B <- 1000
mu_boot_semi <- c()
med_boot_semi <- c()
set.seed(2023)

for (i in 1:B) {
    res_boot <- sample(mod_res, n, replace = T)
    x_boot <- coef(mod_mean) + res_boot
    mu_boot_semi[i] <- mean(x_boot)
    med_boot_semi[i] <- median(x_boot)
}

hist(mu_boot_semi)
hist(med_boot_semi)

quantile(mu_boot_semi, probs = c(0.025, 0.975))
quantile(med_boot_semi, probs = c(0.025, 0.975))

#---------- 4
# Bootstrap
B <- 1000
mu_boot <- c()
med_boot <- c()
n <- length(x)

set.seed(2023)
for (i in 1:B) {
    x_boot <- sample(x, n, replace =  T)
    mu_boot[i] <- mean(x_boot)
    med_boot[i] <- median(x_boot)
}

# MSE for mean
mu_error <- mu_x - mu_boot
(mu_sme <- mean(mu_error^2))

# MSE for median
med_error <- med_x - med_boot
(med_sme <- mean(med_error^2))

# Jackknife

B <- 1000
mu_jn <- c()
med_jn <- c()
n <- length(x)

for (i in 1:n) {
    x_jn <- x[-i]
    mu_jn[i] <- mean(x_jn)
    med_jn[i] <- median(x_jn)
}

# MSE for mean (adjust for inflation factor)
mu_er_jn <- mu_x - mu_jn
(mu_sme_jn <- mean(mu_er_jn^2))

# MSE for median (adjust for inflation factor)
med_er_jn <- med_x - med_jn
(med_sme_jn <- mean(med_er_jn^2))


#---------- 5
# Non-parametric Bootstrap
t_obs <- t.test(x, mu = 0.35, alternative="less")$statistic


B <- 1000
t_boot <- c()
xtidle <- x - mean(x) + 0.35

for (i in 1:B) {
    x_boot <- sample(xtidle, n, replace = T)
    t_boot[i] <- t.test(x_boot, mu = 0.35, alternative="less")$statistic
}

# p-values
(1 + sum(t_boot < t_obs))/(1 + B)



# Parametric bootstrap
B <- 1000
t_boot_par <- c()

for (i in 1:B) {
    x_boot <- rnorm(n, mean(xtidle), var(xtidle))
    t_boot_par[i] <- t.test(x_boot, mu = 0.35, alternative="less")$statistic
}

# p-values
(1 + sum(t_boot_par < t_obs))/(1 + B)
















