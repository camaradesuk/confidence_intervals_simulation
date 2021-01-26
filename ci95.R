# R code using simulation to estimate 95% confidence intervals for estimates of the # prevalence of reporting of study attributes which have been determined using NLP # approaches

# author; Malcolm Macleod, University of Edinburgh
# contact; malcolm.macleod@ed.ac.uk
# date: 20th January 2021
# version: 1.0.0

# input variables
# number of studies in corpus 
# number of simulations 
# performance of tool
# observed prevalence Po for which you want to calculate CI

# libraries
library(ggplot2)
set.seed(42)

# number of studies in corpus 
n <- 500 

# number of simulations 
nsim <- 10000

#performance of tool
sens <- 0.927 #sensitivity of test
spec <- 0.922 #specificity of test

#observed prevalence Po for which you want to calculate CI
obsprev <- 0.40

# create matrix for Pt
perf <- matrix(nrow=nsim,ncol=101)
colnames(perf) <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1)

# simulate studies from population with prevalence Pt [0.01 … 1.00], then simulate 
# performance of tool from sensitivity and specificity

for (i in 1:101) {
  prev <- (i-1)/100
  for (j in 1:nsim) {
    n3 <- 3*n
    mat <- matrix(runif(n3),n)
    simat <- matrix(nrow=n, ncol=4)
    for (k in 1:n)
    {
      if (mat[k,1] < prev) {
        simat[k,1] <- 1
        { if (mat[k,2] < sens) {
          simat[k,4] <- 1
        }}}
      {
        if (mat[k,1] > prev) {
          { if (mat[k,3] > spec) {
            simat[k,4] <- 1
          }}}}}
    
    # calculate unadjusted prevalence
    unadj <- sum(simat[,4], na.rm = TRUE)
    
    # calculate adjusted simulated prevalence Ps
    adj <- (unadj/n + spec -1)/(sens + spec -1)
    if (adj < 0) {
      adj <- 0
    }
    perf[j,i] <- adj
  }
}

# create bin to capture outputs where Ps = Po +/- 0.5%

opmin <- obsprev-0.005
opmax <- obsprev+0.005

# collect values of Pt where Ps = Po

for (p in 1:nsim) {
  { for (q in 1:101) { 
    if (perf[p,q] > opmin & perf[p,q] <  opmax) {
      perf[p,q] <- 1 
    } else {
      perf[p,q] <- 0
    }
  }}}

# create matrix to store results
perf <- rbind(perf, colSums(perf, na.rm = TRUE))
res <- rbind(c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
               0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 
               0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 
               0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 
               0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 
               0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 
               0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1))

nsim1 <- nsim + 1
res <- rbind(res, perf[nsim1,])

# count instances of Pt where Ps = Po 
denom <- rowSums(res)

# identify boundaries for median, l95 and u95% CIs
l95 <- (denom[2]*0.025)
u95 <- (denom[2]*0.975)
med <- (denom[2]*0.5)

res1 <- matrix(nrow = 2, ncol = 101)
for (i in 1:101) {
  res1[1,i] <- sum(res[2,1:i])
  res1[2,i] <- res[1,i] * res[2,i]
}
res <- rbind(res, res1)
# create transposed table
tres <- t(res)
tres <- as.data.frame(tres)

# create matrix for results
bsci <- matrix(nrow = 3, ncol = 5)
colnames(bsci) <- c("N", "lower95", "med", "mean", "upper95")

bsci[1,1] <- "N = "
bsci[1,2] <- "lower 95% = " 
bsci[1,3] <- "Median = " 
bsci[1,4] <- "Mean = " 
bsci[1,5] <- "Upper 95% = "

bsci[2,1] <- denom[2]
bsci[2,4] <- sum(res[4,])/denom[2]

# identify boundary crossing points
for (i in 1:101){
  j <- i-1
  srj <- sum(res[2,1:j])
  sri <- sum(res[2,1:i])
  if (srj < l95 & sri > l95) {
    bsci[2,2] <- res[1,i]
  }
  if (srj < u95 & sri > u95) {
    bsci[2,5] <- res[1,i]
  }
  if (srj < med & sri > med) {
    bsci[2,3] <- res[1,i]
  }}

for (i in 1:5) {
  bsci[3,i] <- paste(bsci[1,i],bsci[2,i])
}

#identify max value of y axis for graphical display

ylen <- max(tres[,2])

# plot distribution
ggplot(data=tres, aes(x=V1, y=V2)) +
  geom_line()+
  geom_point()+
  geom_vline(mapping=aes(xintercept=as.numeric(bsci[2,2])), color="blue") +
  geom_text(mapping=aes(x=as.numeric(bsci[2,2]),y=ylen/1.5,label = bsci[3,2]), size=4, angle=90, vjust=-0.4, hjust=0) +
  geom_vline(mapping=aes(xintercept=as.numeric(bsci[2,3])), color="blue") +
  geom_text(mapping=aes(x=as.numeric(bsci[2,3]),y=ylen/2, label=bsci[3,3]), size=4, angle=90, vjust=-0.4, hjust=0) +
  geom_vline(mapping=aes(xintercept=as.numeric(bsci[2,4])), color="blue") +
  geom_text(mapping=aes(x=as.numeric(format(bsci[2,4], digits = 4)),y=0, label=bsci[3,4]), size=4, angle=90, vjust=-0.4, hjust=0) +
  geom_vline(mapping=aes(xintercept=as.numeric(bsci[2,5])), color="blue") +
  geom_text(mapping=aes(x=as.numeric(bsci[2,5]),y=ylen/1.5, label=bsci[3,5]), size=4, angle=90, vjust=-0.4, hjust=0)
