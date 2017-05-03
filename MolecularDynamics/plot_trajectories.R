#!/usr/bin/env Rscript
dir <- "/home/juan/MolecularDynamics/Proyect"
setwd(dir)


#file <- read.table(file = "file.out", header = TRUE, sep = "\t")
#plot(file$TotEn,type = "l", col = "blue")
#plot(file$KinEn,type = "l", col = "blue", ylim = c(-3.5,1.5))
#lines(file$PotEn,type = "l", col = "red", ylim = c(-3.5,1.5))
#plot(file$Temperature, type = "l")


library(rgl)
lbox = 3.8 #We expect LBOX to be higher than original value, so as to represent all particles after MC run

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
  scatterplot3d::scatterplot3d(traj$PosX,traj$PosY,traj$PosZ,xlim = c(-lbox,lbox),ylim = c(-lbox,lbox),zlim = c(-lbox,lbox))

  dev.off()
}
