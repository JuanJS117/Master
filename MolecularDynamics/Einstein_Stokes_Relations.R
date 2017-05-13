#!/usr/bin/env Rscript
dir <- "/home/juan/MolecularDynamics/Proyect/"
setwd(dir)

VAF1 <- read.table(file = "VAF1.txt", header = TRUE, sep = "\t")
VAF2 <- read.table(file = "VAF2.txt", header = TRUE, sep = "\t")
VAF3 <- read.table(file = "VAF3.txt", header = TRUE, sep = "\t")
VAF4 <- read.table(file = "VAF4.txt", header = TRUE, sep = "\t")
VAF5 <- read.table(file = "VAF5.txt", header = TRUE, sep = "\t")

FAF1 <- read.table(file = "FAF1.txt", header = TRUE, sep = "\t")
FAF2 <- read.table(file = "FAF2.txt", header = TRUE, sep = "\t")
FAF3 <- read.table(file = "FAF3.txt", header = TRUE, sep = "\t")
FAF4 <- read.table(file = "FAF4.txt", header = TRUE, sep = "\t")
FAF5 <- read.table(file = "FAF5.txt", header = TRUE, sep = "\t")

file1 <- read.table(file = "VISCOSITY/file_k_n2.out", header = TRUE, sep = "\t")
file2 <- read.table(file = "VISCOSITY/file_k_n4.out", header = TRUE, sep = "\t")
file3 <- read.table(file = "VISCOSITY/file_k_n6.out", header = TRUE, sep = "\t")
file4 <- read.table(file = "VISCOSITY/file_k_n8.out", header = TRUE, sep = "\t")
file5 <- read.table(file = "VISCOSITY/file_k_n10.out", header = TRUE, sep = "\t")

VAF1$Diff_coef <- scale(VAF1$Diff_coef)
VAF2$Diff_coef <- scale(VAF2$Diff_coef)
VAF3$Diff_coef <- scale(VAF3$Diff_coef)
VAF4$Diff_coef <- scale(VAF4$Diff_coef)
VAF5$Diff_coef <- scale(VAF5$Diff_coef)

FAF1$Fric_coef <- scale(FAF1$Fric_coef)
FAF2$Fric_coef <- scale(FAF2$Fric_coef)
FAF3$Fric_coef <- scale(FAF3$Fric_coef)
FAF4$Fric_coef <- scale(FAF4$Fric_coef)
FAF5$Fric_coef <- scale(FAF5$Fric_coef)



Diff_coefs <- c( max(VAF1$Diff_coef) , max(VAF2$Diff_coef) , max(VAF3$Diff_coef) , max(VAF4$Diff_coef) , max(VAF5$Diff_coef) )
Fric_coefs <- c( max(FAF1$Fric_coef) , max(FAF2$Fric_coef) , max(FAF3$Fric_coef) , max(FAF4$Fric_coef) , max(FAF5$Fric_coef) )
Temps <- c( mean(file1$Temperature) , mean(file2$Temperature) , mean(file3$Temperature) , mean(file4$Temperature) , mean(file5$Temperature) )

# Einstein Relation
plot(Diff_coefs, type = "p", ylim = c(0,2), pch = 15, col = "blue", main = "Einstein-Smoluchowski Relation", xlab = NA, ylab = "Diff. coefs")
grid(nx = NULL, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 
points(Temps/Fric_coefs, pch = 17, col = "red")
lines(legend(x = 2,y = 1.5, legend = c('Calculated Diffusion Coefficients','Predicted Diffusion Coefficients'), fill = c("blue","red")))

# Stokes Relation
kin_visc <- 0.01789644
pred_fric <- 6*pi*0.5*kin_visc
pred_frics <- c(pred_fric,pred_fric,pred_fric,pred_fric,pred_fric)
plot(Fric_coefs, type = "p", ylim = c(0,4) , pch = 15, col = "blue", main = "Stokes Relation", xlab = NA, ylab = "Fric. coefs")
grid(nx = NULL, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 
points(pred_frics, pch = 17, col = "red")
lines(legend(x = 2,y = 2, legend = c('Calculated Friction Coefficients','Predicted Friction Coefficient'), fill = c("blue","red")))
