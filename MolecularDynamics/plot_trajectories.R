#!/usr/bin/env Rscript
dir <- "/home/juan/MolecularDynamics/Proyect/"
setwd(dir)


lbox = 7.368063






file <- read.table(file = "file1.out", header = TRUE, sep = "\t")
plot(file$TotEn,type = "l", col = "blue", ylim = c(0,10))
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
plot(file$KinEn,type = "l", col = "blue", ylim = c(0,10))
plot(file$PotEn,type = "l", col = "red", ylim = c(-2,5))
plot(file$Temperature, type = "l", col = "blue",ylim = c(0,6))
plot(file$Pressure, type = "l", col = "blue")


VAF <- read.table(file = "VAF1.txt", header = TRUE, sep = "\t")
VAF$VAF <- VAF$VAF / (VAF$VAF[1])
plot(VAF$VAF,type = "l", lwd = 2  ,  col = "blue", xlab = expression('Time'[LJ]), ylab = "VACF", main = "Velocity Autocorrelation Function")
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
#VAF_norm <- scale(VAF$VAF) # Normalization
#plot(VAF_norm,type = "l")
VAF$Diff_coef <- scale(VAF$Diff_coef)
plot(VAF$Diff_coef, type = "l", lwd = 2, col = "red", xlab = expression('Time'[LJ]), ylab = "Diff. coef.", main = "Diffusion Coefficient")
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
Diffusion_Coefficient = max(VAF$Diff_coef)
Diffusion_Coefficient

FAF <- read.table(file = "FAF2.txt", header = TRUE, sep = "\t")
FAF$FAF <- FAF$FAF / (FAF$FAF[1])
plot(FAF$FAF, type = "l", lwd = 2, col = "blue", xlab = expression('Time'[LJ]), ylab = "FACF", main = "Force Autocorrelation Function")
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
#FAF_norm <- scale(FAF$FAF) # Normalization
#plot(FAF_norm,type = "l")
FAF$Fric_coef <- scale(FAF$Fric_coef)
plot(FAF$Fric_coef, type = "l", lwd = 2, col = "red", xlab = expression('Time'[LJ]), ylab = "Fric. coef.", main = "Friction Coefficient")
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
Friction_Coefficient = max(FAF$Fric_coef)
Friction_Coefficient

## Einstein-Smoluchowski relation
# If the diffusion coefficient is equal to the temperature 
# divided by the friction coefficient, then the relation is accomplished
mean(file$Temperature)
mean(file$Temperature)/Friction_Coefficient


GAF <- read.table(file = "GAF2.txt", header = TRUE, sep = "\t")
plot(GAF$GAF, type = "l", col = "blue")
GAF$GAF <- GAF$GAF / (GAF$GAF[1])
plot(GAF$GAF, type = "l", lwd = 2, col = "blue", xlab = expression('Time'[LJ]), ylab = expression(paste("<g"[k]*"(t).g"[k]*"(0)>/<g"[k]*"(0)"^2*">")) , main = expression(bold("G"[k]*" Autocorrelation Function")))
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 
iters = c(1:500) 
lnGAF = -log(GAF$GAF)
plot(iters,lnGAF, lwd = 2, type = "l", col = "red", xlab = expression('Time'[LJ]), ylab = expression("viscosity*k"^2*"t"), main = "ln(GACF)") 
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0 

# We fit the ln(GACF) with time to obtain a slope
visc_fit <- lsfit(iters,-log(GAF$GAF[1:500])) 
summary(visc_fit)

# From here we extract the R^2 value
visc_fit2 <- lm( lnGAF ~ iters) 
summary(visc_fit2)

# We add the fitting line to the previous plot
slope = visc_fit$coefficients[2]
lines(iters,slope*iters+visc_fit$coefficients[1], col = "blue")
lines(legend(x = 250,y = -2.5, legend = c(expression(paste("ln[<g"[k]*"(t).g"[k]*"(0)>/<g"[k]*"(0)"^2*">]")),expression('Fitting line [R'^2*' = 0.9968]')), fill = c('red','blue')))


#GAF_norm <- scale(GAF$GAF) # Normalization
#plot(GAF_norm,type = "l")



## Decay of VACF
iters <- c(1:1000)
iters <- iters/1000
VAF <- read.table(file = "VAF1.txt", header = TRUE, sep = "\t")
VAF$VAF <- VAF$VAF / (VAF$VAF[1])
plot(log(iters),log(VAF$VAF),type = "l", lwd = 2  ,  col = "blue", ylim = c(-10,2), xlab = expression('ln(Time'[LJ]*')'), ylab = expression("ln(VACF)"), main = "VACF Decay with time")
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0
lines(log(iters),log(iters^(-3/2)) - 5, type = "l", lwd = 2, col = "red")
lines(log(iters),log(exp(-iters)) , type = "l", lwd = 2, col = "green")
lines(legend(x = -6,y = -2.5, legend = c( "Real VACF", expression("t"^"-3/2"), expression("e"^-t) ) , fill = c('blue','red','green')))


library(rgl)
lbox = 3.8 #We expect LBOX to be higher than original value, so as to represent all particles after MC run
VAF <- read.table(file = "VAF2.txt", header = TRUE, sep = "\t")
VAF$VAF <- VAF$VAF / (VAF$VAF[1])
plot(VAF$VAF,type = "l", lwd = 2  ,  col = "blue", xlab = expression('Time'[LJ]), ylab = "VACF", main = "Velocity Autocorrelation Function")
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0VAF <- read.table(file = "VAF2.txt", header = TRUE, sep = "\t")
VAF$VAF <- VAF$VAF / (VAF$VAF[1])
plot(VAF$VAF,type = "l", lwd = 2  ,  col = "blue", xlab = expression('Time'[LJ]), ylab = "VACF", main = "Velocity Autocorrelation Function")
grid(nx = NA, ny = NULL, col = "black", lty = "dotted",lwd = par("lwd"), equilogs = TRUE) # Add a grid marking y = 0


frames = 99

for(i in 0:frames){
  # creating a name for each plot file with leading zeros
  if (i < 10) {name = paste('00000',i,'plot.png',sep='')}
  if (i < 100 && i >= 10) {name = paste('0000',i,'plot.png', sep='')}
  if (i < 1000 && i >= 100) {name = paste('000',i,'plot.png', sep='')}
  if (i < 10000 && i >= 1000) {name = paste('00',i,'plot.png', sep='')}
  if (i < 100000 && i >= 10000) {name = paste('0',i,'plot.png', sep='')}


  png(name,width = 800,height = 800)

  filename = paste((i)*10,'_traj.dat',sep='')
  traj <- read.table(file = filename, header = TRUE, sep = "\t")
  scatterplot3d::scatterplot3d(traj$PosX,traj$PosY,traj$PosZ,xlim = c(-lbox,lbox),ylim = c(-lbox,lbox),zlim = c(-lbox,lbox), xlab = 'x', ylab = 'y', zlab = 'z')

  dev.off()
}
