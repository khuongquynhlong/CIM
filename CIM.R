library(tidyverse)
library(magrittr)
library(knitr)
library(kableExtra)

#-------------------

#---------- Project 2
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


















