#!/usr/bin/env Rscript
dir <- "/home/juan/MolecularDynamics/Proyect"
setwd(dir)

#plot_col = rainbow(5)
file1 <- read.table(file = "file1.out", header = TRUE, sep = "\t")
file2 <- read.table(file = "file2.out", header = TRUE, sep = "\t")
file3 <- read.table(file = "file3.out", header = TRUE, sep = "\t")
file4 <- read.table(file = "file4.out", header = TRUE, sep = "\t")
file5 <- read.table(file = "file5.out", header = TRUE, sep = "\t")

plot(file1$TotEn,type = "l", col = "blue", ylim = c(-2,4), main = expression(bold(paste("Energy [10"^5*" steps, dt = 0.001, ",rho," = 0.5]"))), ylab = expression(paste("Energy [",epsilon,"]")), xlab = expression('Time'[LJ]))
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
#lines(file2$TotEn,type = "l", col = "red", ylim = c(0,5))
#lines(file3$TotEn,type = "l", col = "green", ylim = c(0,5))
#lines(file4$TotEn,type = "l", col = "orange", ylim = c(0,5))
#lines(file5$TotEn,type = "l", col = "cyan", ylim = c(0,5))
lines(file1$KinEn,type = "l", col = "green", ylim = c(-2,4))
lines(file1$PotEn,type = "l", col = "red", ylim = c(-2,4))
lines(legend(x = 30000,y = 1, legend = c('Total Energy','Kinetic Energy','Potential Energy'), fill = c('blue','green','red')))

#par(mfrow=c(1,1))

plot(file1$Temperature, type = "l", col = "blue",ylim = c(0,3) , ylab = expression(paste("Temperature [K"[b]*"T/",epsilon,"]")), xlab = expression('Time'[LJ]) , main = expression(bold(paste("Temperature [10"^5*" steps, dt = 0.001, ",rho," = 0.5]"))))
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0

plot(file1$MomentumX, ylim = c(-17,11) ,ylab = expression(paste("Momentum [m*",sigma,"/t"[LJ]*"]")), xlab = expression('Time'[LJ]) ,type = "p", col = "blue", main = expression(bold(paste("Momentum Conservation [10"^5*" steps, dt = 0.001, ",rho," = 0.5]"))))
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
lines(file1$MomentumY, type = 'p', col = 'green') 
lines(file1$MomentumZ, type = 'p', col = 'red') 
lines(legend(x = 30000,y = -5, legend = c('Mom. in X','Mom. in Y','Mom. in Z'), fill = c('blue','green','red')))


# Lennard-Jones potential plot

r = seq(from = 1, to = 2.5, by = 0.01)
LJ = 4*((1/r)^12 - (1/r)^6)
LJ_trunc = LJ + 0.01631689

plot(r,LJ, type = 'l', lwd = 2, col = "blue", xlab = expression(paste("Dist [",sigma,"]")), ylab =  expression(paste("Potential [",epsilon,"]")) , main = "Lennard-Jones potential" )
lines(r,LJ_trunc, type = 'l', lwd = 2, col = "red")
lines(legend(x = 1.6,y = -0.5, legend = c('Lennard-Jones potential','Truncated Lennard-Jones potential'), fill = c('blue','red')))
