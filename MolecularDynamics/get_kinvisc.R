#!/usr/bin/env Rscript
dir <- "/home/juan/MolecularDynamics/Proyect/VISCOSITY/"
setwd(dir)

lbox = 7.368063

GAF1 <- read.table(file = "GAF_k_n2.txt", header = TRUE, sep = "\t")
GAF2 <- read.table(file = "GAF_k_n4.txt", header = TRUE, sep = "\t")
GAF3 <- read.table(file = "GAF_k_n6.txt", header = TRUE, sep = "\t")
GAF4 <- read.table(file = "GAF_k_n8.txt", header = TRUE, sep = "\t")
GAF5 <- read.table(file = "GAF_k_n10.txt", header = TRUE, sep = "\t")


##---------------------------------------------------------------------
GAF1$GAF <- GAF1$GAF / (GAF1$GAF[1])
GAF2$GAF <- GAF2$GAF / (GAF2$GAF[1])
GAF3$GAF <- GAF3$GAF / (GAF3$GAF[1])
GAF4$GAF <- GAF4$GAF / (GAF4$GAF[1])
GAF5$GAF <- GAF5$GAF / (GAF5$GAF[1])

plot_col = c("blue4","blue1","cyan4","deepskyblue3","cyan3")

plot(GAF1$GAF, type = "l", ylim = c(0,1), lwd = 2, col = plot_col[1] , xlab = expression('Time'[LJ]), ylab = expression(paste("<g"[k]*"(t).g"[k]*"(0)>/<g"[k]*"(0)"^2*">")) , main = expression(bold("G"[k]*" Autocorrelation Function")))
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 
lines(GAF2$GAF, type = "l", lwd = 2, col = plot_col[2] , xlab = expression('Time'[LJ]), ylab = expression(paste("<g"[k]*"(t).g"[k]*"(0)>/<g"[k]*"(0)"^2*">")) , main = expression(bold("G"[k]*" Autocorrelation Function")))
lines(GAF3$GAF, type = "l", lwd = 2, col = plot_col[3] , xlab = expression('Time'[LJ]), ylab = expression(paste("<g"[k]*"(t).g"[k]*"(0)>/<g"[k]*"(0)"^2*">")) , main = expression(bold("G"[k]*" Autocorrelation Function")))
lines(GAF4$GAF, type = "l", lwd = 2, col = plot_col[4] , xlab = expression('Time'[LJ]), ylab = expression(paste("<g"[k]*"(t).g"[k]*"(0)>/<g"[k]*"(0)"^2*">")) , main = expression(bold("G"[k]*" Autocorrelation Function")))
lines(GAF5$GAF, type = "l", lwd = 2, col = plot_col[5] , xlab = expression('Time'[LJ]), ylab = expression(paste("<g"[k]*"(t).g"[k]*"(0)>/<g"[k]*"(0)"^2*">")) , main = expression(bold("G"[k]*" Autocorrelation Function")))
lines(legend(x = 350,y = 1, legend = c('n = 2','n = 4','n = 6','n = 8','n = 10'), fill = plot_col))

iters = c(1:60) 

lnGAF1 = -log(GAF1$GAF[1:60])
#lnGAF1 <- scale(lnGAF1)
lnGAF2 = -log(GAF2$GAF[1:60])
#lnGAF2 <- scale(lnGAF2)
lnGAF3 = -log(GAF3$GAF[1:60])
#lnGAF3 <- scale(lnGAF3)
lnGAF4 = -log(GAF4$GAF[1:60])
#lnGAF4 <- scale(lnGAF4)
lnGAF5 = -log(GAF5$GAF[1:60])
#lnGAF5 <- scale(lnGAF5)


plot(iters,lnGAF1, lwd = 2, type = "l", col = plot_col[1], ylim = c(-1.5,2), xlab = expression('Time'[LJ]), ylab = expression("t*viscosity*k"^2), main = "-ln(GACF)") 
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 
lines(iters,lnGAF2, lwd = 2, type = "l", col = plot_col[2], xlab = expression('Time'[LJ]), ylab = expression("viscosity*k"^2*"t"), main = "ln(GACF)") 
lines(iters,lnGAF3, lwd = 2, type = "l", col = plot_col[3], xlab = expression('Time'[LJ]), ylab = expression("viscosity*k"^2*"t"), main = "ln(GACF)") 
lines(iters,lnGAF4, lwd = 2, type = "l", col = plot_col[4], xlab = expression('Time'[LJ]), ylab = expression("viscosity*k"^2*"t"), main = "ln(GACF)") 
lines(iters,lnGAF5, lwd = 2, type = "l", col = plot_col[5], xlab = expression('Time'[LJ]), ylab = expression("viscosity*k"^2*"t"), main = "ln(GACF)") 
lines(legend(x = 5,y = 1.5, legend = c('n = 2','n = 4','n = 6','n = 8','n = 10'), fill = plot_col))


# We fit the ln(GACF) with time to obtain a slope
visc_fit1 <- lsfit(iters,lnGAF1) 
summary(visc_fit1)
visc_fit2 <- lsfit(iters,lnGAF2) 
summary(visc_fit2)
visc_fit3 <- lsfit(iters,lnGAF3) 
summary(visc_fit3)
visc_fit4 <- lsfit(iters,lnGAF4) 
summary(visc_fit4)
visc_fit5 <- lsfit(iters,lnGAF5) 
summary(visc_fit5)

# From here we extract the R^2 value
visc_fit1_2 <- lm( lnGAF1 ~ iters) 
summary(visc_fit1_2)
visc_fit2_2 <- lm( lnGAF2 ~ iters) 
summary(visc_fit2_2)
visc_fit3_2 <- lm( lnGAF3 ~ iters) 
summary(visc_fit3_2)
visc_fit4_2 <- lm( lnGAF4 ~ iters) 
summary(visc_fit4_2)
visc_fit5_2 <- lm( lnGAF5 ~ iters) 
summary(visc_fit5_2)

# We add the fitting line to the previous plot
slope1 = visc_fit1$coefficients[2]
cut1 = visc_fit1$coefficients[1]
plot(iters,slope1*iters+cut1, lwd = 2, type = "l", col = plot_col[1], ylim = c(-2.6,2), xlab = expression('Time'[LJ]), ylab = expression("predicted t*viscosity*k"^2), main = "Least Squares Fitting" )
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 
slope2 = visc_fit2$coefficients[2]
cut2 = visc_fit2$coefficients[1]
lines(iters,slope2*iters+cut2, lwd = 2, type = "l", col = plot_col[2])
slope3 = visc_fit3$coefficients[2]
cut3 = visc_fit3$coefficients[1]
lines(iters,slope3*iters+cut3, lwd = 2, type = "l", col = plot_col[3])
slope4 = visc_fit4$coefficients[2]
cut4 = visc_fit4$coefficients[1]
lines(iters,slope4*iters+cut4, lwd = 2, type = "l", col = plot_col[4])
slope5 = visc_fit5$coefficients[2]
cut5 = visc_fit5$coefficients[1]
lines(iters,slope5*iters+cut5, lwd = 2, type = "l", col = plot_col[5])
lines(legend(x = 10,y = 1.5, legend = c(expression("n = 2 [R"^2*" = 0.9596]"),expression("n = 4 [R"^2*" = 0.9608]"),expression("n = 6 [R"^2*" = 0.9625]"),expression("n = 8 [R"^2*" = 0.9636]"),expression("n = 10 [R"^2*" = 0.9746]")), fill = plot_col))



# Viscosity calculation
n <- c(2,4,6,8,10)
k <- 2*pi*n/lbox
slopes = c(slope1,slope2,slope3,slope4,slope5)
kin_visc = slopes/(k^2)

plot(k^2,kin_visc, type = "p", col = "blue", pch = 15, xlab = expression("k"[n]^2), ylab = expression(nu[k]), main = expression(bold("Relation between k"[n]^2*" and "*nu[k])))
grid(nx = NULL, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 

kinvisc_fit <- lsfit(k[1:4]^2,kin_visc[1:4])
lines(k^2,kinvisc_fit$coefficients[2]*k^2+kinvisc_fit$coefficients[1], type = "l", lwd = 2, col = 'red')
kinvisc_fit$coefficients[1]

