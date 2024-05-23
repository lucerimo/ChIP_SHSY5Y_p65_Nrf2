########################################################################################################
######    Analisis de resultados del experimento de ciclo celular en celulas SH-SY5Y     ###############
######    Respuesta: PI (rojo)                                                           ###############
######    Variables: Experimento (Bloque)                                                ###############
######        Sustancia (14.2.2;Negativo;Positivo;Resveratrol;withaferina A)             ###############
######        Incubacion (4h, 24h)                                                       ###############
######        Concentracion de la sustancia (depende de la sustancia)                    ###############
######        Rotenona (Si, No)                                                          ###############
########################################################################################################
library(ggprism)
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(psych)
library(scales)
theme_set(theme_pubclean())
library(multcomp)
library(mvtnorm)
library(survival)
library(lmtest)
library(graphics)
library(car)
library(MASS)

#1. Cargar funciones


range <- function(x) {
  data.frame(ymin=min(x),
             ymax=max(x))
}

graf.barras<-function(v.respuesta,l.incub,compuesto,etiq.y,sd.respuesta,n.vector)
{
  
  barplot(v.respuesta,width=1,col=c("#f2f2f2ff","#4d4d4dff", "#aaccffff", "#5599ffff", "#0055d4ff"),main=paste(paste(l.incub, "h",collapse=""),collapse=""),
          xlab=paste("(µM)",collapse=""), ylab=etiq.y,xlim=c(0,6), cex.main =2,
          ylim=c(0,50))
  arrows(x0=c(0.7,1.9,3.1,4.3,5.5),y0=v.respuesta-(sd.respuesta/sqrt(n.vector)),x1=c(0.7,1.9,3.1,4.3,5.5),y1=v.respuesta+(sd.respuesta/sqrt(n.vector)),
         length= 0.1,angle= 90, code= 3)
}

graf.barras.cc<-function(v.respuesta,l.incub,compuesto,etiq.y,sd.respuesta,n.vector)
{
  barplot(v.respuesta,width=0.9,col=c("#f2f2f2ff","#4d4d4dff", "#aaccffff", "#5599ffff", "#0055d4ff","#f2f2f2ff","#4d4d4dff", "#aaccffff", "#5599ffff", "#0055d4ff","#f2f2f2ff","#4d4d4dff", "#aaccffff", "#5599ffff", "#0055d4ff"),
          main=paste(paste(l.incub, "h",collapse=""),collapse=""), cex.main = 1.75,
          xlab="", 
          ylab=etiq.y,xlim=c(0,16),
          ylim=c(0,100),
          space = c(0,0,0,0,0,1,0,0,0,0,1,0,0,0,0))
  #title(xlab=paste("", "\n", paste(compuesto,"(µM)",collapse=""), collapse = ""), line = 3.75)
  axis(1,pos=-15,at=c(0,1,2,3,4.5),col="black",line=2.5,tick=T,labels=rep("",5),lwd=2,lwd.ticks=0)
  axis(1,pos=-15,at=5.4+c(0,1,2,3,4.5),col="black",line=2.5,tick=T,labels=rep("",5),lwd=2,lwd.ticks=0)
  axis(1,pos=-15,at=10.85+c(0,1,2,3,4.5),col="black",line=2.5,tick=T,labels=rep("",5),lwd=2,lwd.ticks=0)
  text(c(2.25), -20, paste('G0/G1'), xpd=NA)
  text(c(7.65), -20, paste("S"), xpd=NA)
  text(c(13), -20, paste("G2/M"), xpd=NA)
  text(c(7.65), -27, paste("(µM)"), xpd=NA)
  arrows(x0=c(0.5,1.35,2.25,3.15,4.05,5.9,6.75,7.65,8.55,9.45,11.25,12.15,13.05,14.0,14.92),y0=v.respuesta-(sd.respuesta/sqrt(n.vector)),x1=c(0.5,1.35,2.25,3.15,4.05,5.9,6.75,7.65,8.55,9.45,11.25,12.15,13.05,14.0,14.92),y1=v.respuesta+(sd.respuesta/sqrt(n.vector)),
         length= 0.05,angle= 90, code= 3)
}
#par(mar=c(6,4,4,4))
#2. Cargar datos
{
datos<-read.csv("Resultados_2022_cellcycle.csv", sep = ";")
is.factor(datos$tratam)
is.factor(datos$Experimento)
levels(as.factor(datos$tratam))
levels(as.factor(datos$Experimento))
datos$tratam<-relevel(as.factor(datos$tratam),ref="CP")
datos$Experimento<-as.factor(datos$Experimento)
datos
datos$tratam<-relevel(as.factor(datos$tratam),ref="CP")
datos
attach(datos)
}

#3. Resumen de los datos
{
  #SubG1
subg1.mean<-tapply(SubG1,list(Sustancia,Cn,Incubacion,Rotenona),FUN=mean)
subg1.sd<-tapply(SubG1,list(Sustancia,Cn,Incubacion,Rotenona),FUN=sd)
#G0
g0.mean<-tapply(G0,list(Sustancia,Cn,Incubacion,Rotenona),FUN=mean)
g0.sd<-tapply(G0,list(Sustancia,Cn,Incubacion,Rotenona),FUN=sd)
#S
s.mean<-tapply(S,list(Sustancia,Cn,Incubacion,Rotenona),FUN=mean)
s.sd<-tapply(S,list(Sustancia,Cn,Incubacion,Rotenona),FUN=sd)
#G2
g2.mean<-tapply(G2,list(Sustancia,Cn,Incubacion,Rotenona),FUN=mean)
g2.sd<-tapply(G2,list(Sustancia,Cn,Incubacion,Rotenona),FUN=sd)
#conteo
n.cell<-tapply(SubG1,list(Sustancia,Cn,Incubacion,Rotenona),FUN=length)
}

#4. Construir tablas con datos por comupesto
{
#Subg1
#14.2.2
m.14_4h_rot.sg<-c(subg1.mean[2,1,1,1],subg1.mean[3,1,1,2],subg1.mean[1,2:4,1,2])
m.14_24h_rot.sg<-c(subg1.mean[2,1,2,1],subg1.mean[3,1,2,2],subg1.mean[1,2:4,2,2])
m.14_4h_srot.sg<-c(subg1.mean[2,1,1,1],subg1.mean[3,1,1,2],subg1.mean[1,2:4,1,1])
m.14_24h_srot.sg<-c(subg1.mean[2,1,2,1],subg1.mean[3,1,2,2],subg1.mean[1,2:4,2,1])
m.14.subg1<-as.matrix(rbind(m.14_4h_rot.sg,m.14_24h_rot.sg, m.14_4h_srot.sg,m.14_24h_srot.sg));colnames(m.14.subg1)<-c("CN","CP","0.5","1","2")
sd.14_4h_rot.sg<-c(subg1.sd[2,1,1,1],subg1.sd[3,1,1,2],subg1.sd[1,2:4,1,2])
sd.14_24h_rot.sg<-c(subg1.sd[2,1,2,1],subg1.sd[3,1,2,2],subg1.sd[1,2:4,2,2])
sd.14_4h_srot.sg<-c(subg1.sd[2,1,1,1],subg1.sd[3,1,1,2],subg1.sd[1,2:4,1,1])
sd.14_24h_srot.sg<-c(subg1.sd[2,1,2,1],subg1.sd[3,1,2,2],subg1.sd[1,2:4,2,1])
sd.14.sugb1<-as.matrix(rbind(sd.14_4h_rot.sg,sd.14_24h_rot.sg, sd.14_4h_srot.sg,sd.14_24h_srot.sg));colnames(sd.14.sugb1)<-c("CN","CP","0.5","1","2")
n.14.subg1<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[1,2:4,2,2]))
#Resveratrol
m.r_4h_rot.sg<-c(subg1.mean[2,1,1,1],subg1.mean[3,1,1,2],subg1.mean[5,5:7,1,2])
m.r_24h_rot.sg<-c(subg1.mean[2,1,2,1],subg1.mean[3,1,2,2],subg1.mean[5,5:7,2,2])
m.r_4h_srot.sg<-c(subg1.mean[2,1,1,1],subg1.mean[3,1,1,2],subg1.mean[5,5:7,1,1])
m.r_24h_srot.sg<-c(subg1.mean[2,1,2,1],subg1.mean[3,1,2,2],subg1.mean[5,5:7,2,1])
m.r.subg1<-as.matrix(rbind(m.r_4h_rot.sg,m.r_24h_rot.sg, m.r_4h_srot.sg,m.r_24h_srot.sg));colnames(m.r.subg1)<-c("CN","CP","0.5","1","2")
sd.r_4h_rot.sg<-c(subg1.sd[2,1,1,1],subg1.sd[3,1,1,2],subg1.sd[5,5:7,1,2])
sd.r_24h_rot.sg<-c(subg1.sd[2,1,2,1],subg1.sd[3,1,2,2],subg1.sd[5,5:7,2,2])
sd.r_4h_srot.sg<-c(subg1.sd[2,1,1,1],subg1.sd[3,1,1,2],subg1.sd[5,5:7,1,1])
sd.r_24h_srot.sg<-c(subg1.sd[2,1,2,1],subg1.sd[3,1,2,2],subg1.sd[5,5:7,2,1])
sd.r.subg1<-as.matrix(rbind(sd.r_4h_rot.sg,sd.r_24h_rot.sg, sd.r_4h_srot.sg,sd.r_24h_srot.sg));colnames(sd.r.subg1)<-c("CN","CP","0.5","1","2")
n.r.sugb1<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[5,5:7,2,2]))
#Withaferin A
m.WA_4h_rot.sg<-c(subg1.mean[2,1,1,1],subg1.mean[3,1,1,2],subg1.mean[6,2:4,1,2])
m.WA_24h_rot.sg<-c(subg1.mean[2,1,2,1],subg1.mean[3,1,2,2],subg1.mean[6,2:4,2,2])
m.WA_4h_srot.sg<-c(subg1.mean[2,1,1,1],subg1.mean[3,1,1,2],subg1.mean[6,2:4,1,1])
m.WA_24h_srot.sg<-c(subg1.mean[2,1,2,1],subg1.mean[3,1,2,2],subg1.mean[6,2:4,2,1])
m.WA.subG1<-as.matrix(rbind(m.WA_4h_rot.sg,m.WA_24h_rot.sg, m.WA_4h_srot.sg,m.WA_24h_srot.sg));colnames(m.WA.subG1)<-c("CN","CP","0.5","1","2")
sd.WA_4h_rot.sg<-c(subg1.sd[2,1,1,1],subg1.sd[3,1,1,2],subg1.sd[6,2:4,1,2])
sd.WA_24h_rot.sg<-c(subg1.sd[2,1,2,1],subg1.sd[3,1,2,2],subg1.sd[6,2:4,2,2])
sd.WA_4h_srot.sg<-c(subg1.sd[2,1,1,1],subg1.sd[3,1,1,2],subg1.sd[6,2:4,1,1])
sd.WA_24h_srot.sg<-c(subg1.sd[2,1,2,1],subg1.sd[3,1,2,2],subg1.sd[6,2:4,2,1])
sd.WA.subg1<-as.matrix(rbind(sd.WA_4h_rot.sg,sd.WA_24h_rot.sg, sd.WA_4h_srot.sg,sd.WA_24h_srot.sg));colnames(sd.WA.subg1)<-c("CN","CP","0.5","1","2")
n.WA.subg1<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[6,2:4,2,2]))

#G0
#14.2.2
m.14_4h_rot.g0<-c(g0.mean[2,1,1,1],g0.mean[3,1,1,2],g0.mean[1,2:4,1,2])
m.14_24h_rot.g0<-c(g0.mean[2,1,2,1],g0.mean[3,1,2,2],g0.mean[1,2:4,2,2])
m.14_4h_srot.g0<-c(g0.mean[2,1,1,1],g0.mean[3,1,1,2],g0.mean[1,2:4,1,1])
m.14_24h_srot.g0<-c(g0.mean[2,1,2,1],g0.mean[3,1,2,2],g0.mean[1,2:4,2,1])
m.14.g0<-as.matrix(rbind(m.14_4h_rot.g0,m.14_24h_rot.g0, m.14_4h_srot.g0,m.14_24h_srot.g0));colnames(m.14.g0)<-c("CN","CP","0.5","1","2")
sd.14_4h_rot.g0<-c(g0.sd[2,1,2,1],g0.sd[3,1,1,2],g0.sd[1,2:4,1,2])
sd.14_24h_rot.g0<-c(g0.sd[2,1,1,1],g0.sd[3,1,2,2],g0.sd[1,2:4,2,2])
sd.14_4h_srot.g0<-c(g0.sd[2,1,2,1],g0.sd[3,1,1,2],g0.sd[1,2:4,1,1])
sd.14_24h_srot.g0<-c(g0.sd[2,1,1,1],g0.sd[3,1,2,2],g0.sd[1,2:4,2,1])
sd.14.g0<-as.matrix(rbind(sd.14_4h_rot.g0,sd.14_24h_rot.g0, sd.14_4h_srot.g0,sd.14_24h_srot.g0));colnames(sd.14.sugb1)<-c("CN","CP","0.5","1","2")
n.14.g0<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[1,2:4,2,2]))
#Resveratrol
m.r_4h_rot.g0<-c(g0.mean[2,1,1,1],g0.mean[3,1,1,2],g0.mean[5,5:7,1,2])
m.r_24h_rot.g0<-c(g0.mean[2,1,2,1],g0.mean[3,1,2,2],g0.mean[5,5:7,2,2])
m.r_4h_srot.g0<-c(g0.mean[2,1,1,1],g0.mean[3,1,1,2],g0.mean[5,5:7,1,1])
m.r_24h_srot.g0<-c(g0.mean[2,1,2,1],g0.mean[3,1,2,2],g0.mean[5,5:7,2,1])
m.r.g0<-as.matrix(rbind(m.r_4h_rot.g0,m.r_24h_rot.g0, m.r_4h_srot.g0,m.r_24h_srot.g0));colnames(m.r.g0)<-c("CN","CP","0.5","1","2")
sd.r_4h_rot.g0<-c(g0.sd[2,1,1,1],g0.sd[3,1,1,2],g0.sd[5,5:7,1,2])
sd.r_24h_rot.g0<-c(g0.sd[2,1,2,1],g0.sd[3,1,2,2],g0.sd[5,5:7,2,2])
sd.r_4h_srot.g0<-c(g0.sd[2,1,1,1],g0.sd[3,1,1,2],g0.sd[5,5:7,1,1])
sd.r_24h_srot.g0<-c(g0.sd[2,1,2,1],g0.sd[3,1,2,2],g0.sd[5,5:7,2,1])
sd.r.g0<-as.matrix(rbind(sd.r_4h_rot.g0,sd.r_24h_rot.g0, sd.r_4h_srot.g0,sd.r_24h_srot.g0));colnames(sd.r.g0)<-c("CN","CP","0.5","1","2")
n.r.g0<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[5,5:7,2,2]))
#Withaferin A
m.WA_4h_rot.g0<-c(g0.mean[2,1,1,1],g0.mean[3,1,1,2],g0.mean[6,2:4,1,2])
m.WA_24h_rot.g0<-c(g0.mean[2,1,2,1],g0.mean[3,1,2,2],g0.mean[6,2:4,2,2])
m.WA_4h_srot.g0<-c(g0.mean[2,1,1,1],g0.mean[3,1,1,2],g0.mean[6,2:4,1,1])
m.WA_24h_srot.g0<-c(g0.mean[2,1,2,1],g0.mean[3,1,2,2],g0.mean[6,2:4,2,1])
m.WA.g0<-as.matrix(rbind(m.WA_4h_rot.g0,m.WA_24h_rot.g0, m.WA_4h_srot.g0,m.WA_24h_srot.g0));colnames(m.WA.g0)<-c("CN","CP","0.5","1","2")
sd.WA_4h_rot.g0<-c(g0.sd[2,1,1,1],g0.sd[3,1,1,2],g0.sd[6,2:4,1,2])
sd.WA_24h_rot.g0<-c(g0.sd[2,1,2,1],g0.sd[3,1,2,2],g0.sd[6,2:4,2,2])
sd.WA_4h_srot.g0<-c(g0.sd[2,1,1,1],g0.sd[3,1,1,2],g0.sd[6,2:4,1,1])
sd.WA_24h_srot.g0<-c(g0.sd[2,1,2,1],g0.sd[3,1,2,2],g0.sd[6,2:4,2,1])
sd.WA.g0<-as.matrix(rbind(sd.WA_4h_rot.g0,sd.WA_24h_rot.g0, sd.WA_4h_srot.g0,sd.WA_24h_srot.g0));colnames(sd.WA.g0)<-c("CN","CP","0.5","1","2")
n.WA.g0<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[6,2:4,2,2]))

#S
#14.2.2
m.14_4h_rot.s<-c(s.mean[2,1,1,1],s.mean[3,1,1,2],s.mean[1,2:4,1,2])
m.14_24h_rot.s<-c(s.mean[2,1,2,1],s.mean[3,1,2,2],s.mean[1,2:4,2,2])
m.14_4h_srot.s<-c(s.mean[2,1,1,1],s.mean[3,1,1,2],s.mean[1,2:4,1,1])
m.14_24h_srot.s<-c(s.mean[2,1,2,1],s.mean[3,1,2,2],s.mean[1,2:4,2,1])
m.14.s<-as.matrix(rbind(m.14_4h_rot.s,m.14_24h_rot.s, m.14_4h_srot.s,m.14_24h_srot.s));colnames(m.14.s)<-c("CN","CP","0.5","1","2")
sd.14_4h_rot.s<-c(s.sd[2,1,1,1],s.sd[3,1,1,2],s.sd[1,2:4,1,2])
sd.14_24h_rot.s<-c(s.sd[2,1,2,1],s.sd[3,1,2,2],s.sd[1,2:4,2,2])
sd.14_4h_srot.s<-c(s.sd[2,1,1,1],s.sd[3,1,1,2],s.sd[1,2:4,1,1])
sd.14_24h_srot.s<-c(s.sd[2,1,2,1],s.sd[3,1,2,2],s.sd[1,2:4,2,1])
sd.14.s<-as.matrix(rbind(sd.14_4h_rot.s,sd.14_24h_rot.s, sd.14_4h_srot.s,sd.14_24h_srot.s));colnames(sd.14.s)<-c("CN","CP","0.5","1","2")
n.14.s<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[1,2:4,2,2]))
#Resveratrol
m.r_4h_rot.s<-c(s.mean[2,1,1,1],s.mean[3,1,1,2],s.mean[5,5:7,1,2])
m.r_24h_rot.s<-c(s.mean[2,1,2,1],s.mean[3,1,2,2],s.mean[5,5:7,2,2])
m.r_4h_srot.s<-c(s.mean[2,1,1,1],s.mean[3,1,1,2],s.mean[5,5:7,1,1])
m.r_24h_srot.s<-c(s.mean[2,1,2,1],s.mean[3,1,2,2],s.mean[5,5:7,2,1])
m.r.s<-as.matrix(rbind(m.r_4h_rot.s,m.r_24h_rot.s, m.r_4h_srot.s,m.r_24h_srot.s));colnames(m.r.s)<-c("CN","CP","0.5","1","2")
sd.r_4h_rot.s<-c(s.sd[2,1,1,1],s.sd[3,1,1,2],s.sd[5,5:7,1,2])
sd.r_24h_rot.s<-c(s.sd[2,1,2,1],s.sd[3,1,2,2],s.sd[5,5:7,2,2])
sd.r_4h_srot.s<-c(s.sd[2,1,1,1],s.sd[3,1,1,2],s.sd[5,5:7,1,1])
sd.r_24h_srot.s<-c(s.sd[2,1,2,1],s.sd[3,1,2,2],s.sd[5,5:7,2,1])
sd.r.s<-as.matrix(rbind(sd.r_4h_rot.s,sd.r_24h_rot.s, sd.r_4h_srot.s,sd.r_24h_srot.s));colnames(sd.r.s)<-c("CN","CP","0.5","1","2")
n.r.s<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[5,5:7,2,2]))
#Withaferin A
m.WA_4h_rot.s<-c(s.mean[2,1,1,1],s.mean[3,1,1,2],s.mean[6,2:4,1,2])
m.WA_24h_rot.s<-c(s.mean[2,1,2,1],s.mean[3,1,2,2],s.mean[6,2:4,2,2])
m.WA_4h_srot.s<-c(s.mean[2,1,1,1],s.mean[3,1,1,2],s.mean[6,2:4,1,1])
m.WA_24h_srot.s<-c(s.mean[2,1,2,1],s.mean[3,1,2,2],s.mean[6,2:4,2,1])
m.WA.s<-as.matrix(rbind(m.WA_4h_rot.s,m.WA_24h_rot.s, m.WA_4h_srot.s,m.WA_24h_srot.s));colnames(m.WA.s)<-c("CN","CP","0.5","1","2")
sd.WA_4h_rot.s<-c(s.sd[2,1,1,1],s.sd[3,1,1,2],s.sd[6,2:4,1,2])
sd.WA_24h_rot.s<-c(s.sd[2,1,2,1],s.sd[3,1,2,2],s.sd[6,2:4,2,2])
sd.WA_4h_srot.s<-c(s.sd[2,1,1,1],s.sd[3,1,1,2],s.sd[6,2:4,1,1])
sd.WA_24h_srot.s<-c(s.sd[2,1,2,1],s.sd[3,1,2,2],s.sd[6,2:4,2,1])
sd.WA.s<-as.matrix(rbind(sd.WA_4h_rot.s,sd.WA_24h_rot.s, sd.WA_4h_srot.s,sd.WA_24h_srot.s));colnames(sd.WA.s)<-c("CN","CP","0.5","1","2")
n.WA.s<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[6,2:4,2,2]))

#G2
#14.2.2
m.14_4h_rot.g2<-c(g2.mean[2,1,1,1],g2.mean[3,1,1,2],g2.mean[1,2:4,1,2])
m.14_24h_rot.g2<-c(g2.mean[2,1,2,1],g2.mean[3,1,2,2],g2.mean[1,2:4,2,2])
m.14_4h_srot.g2<-c(g2.mean[2,1,1,1],g2.mean[3,1,1,2],g2.mean[1,2:4,1,1])
m.14_24h_srot.g2<-c(g2.mean[2,1,2,1],g2.mean[3,1,2,2],g2.mean[1,2:4,2,1])
m.14.g2<-as.matrix(rbind(m.14_4h_rot.g2,m.14_24h_rot.g2, m.14_4h_srot.g2,m.14_24h_srot.g2));colnames(m.14.g2)<-c("CN","CP","0.5","1","2")
sd.14_4h_rot.g2<-c(g2.sd[2,1,1,1],g2.sd[3,1,1,2],g2.sd[1,2:4,1,2])
sd.14_24h_rot.g2<-c(g2.sd[2,1,2,1],g2.sd[3,1,2,2],g2.sd[1,2:4,2,2])
sd.14_4h_srot.g2<-c(g2.sd[2,1,1,1],g2.sd[3,1,1,2],g2.sd[1,2:4,1,1])
sd.14_24h_srot.g2<-c(g2.sd[2,1,2,1],g2.sd[3,1,2,2],g2.sd[1,2:4,2,1])
sd.14.g2<-as.matrix(rbind(sd.14_4h_rot.g2,sd.14_24h_rot.g2, sd.14_4h_srot.g2,sd.14_24h_srot.g2));colnames(sd.14.g2)<-c("CN","CP","0.5","1","2")
n.14.g2<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[1,2:4,2,2]))
#Resveratrol
m.r_4h_rot.g2<-c(g2.mean[2,1,1,1],g2.mean[3,1,1,2],g2.mean[5,5:7,1,2])
m.r_24h_rot.g2<-c(g2.mean[2,1,2,1],g2.mean[3,1,2,2],g2.mean[5,5:7,2,2])
m.r_4h_srot.g2<-c(g2.mean[2,1,1,1],g2.mean[3,1,1,2],g2.mean[5,5:7,1,1])
m.r_24h_srot.g2<-c(g2.mean[2,1,2,1],g2.mean[3,1,2,2],g2.mean[5,5:7,2,1])
m.r.g2<-as.matrix(rbind(m.r_4h_rot.g2,m.r_24h_rot.g2, m.r_4h_srot.g2,m.r_24h_srot.g2));colnames(m.r.g2)<-c("CN","CP","0.5","1","2")
sd.r_4h_rot.g2<-c(g2.sd[2,1,1,1],g2.sd[3,1,1,2],g2.sd[5,5:7,1,2])
sd.r_24h_rot.g2<-c(g2.sd[2,1,2,1],g2.sd[3,1,2,2],g2.sd[5,5:7,2,2])
sd.r_4h_srot.g2<-c(g2.sd[2,1,1,1],g2.sd[3,1,1,2],g2.sd[5,5:7,1,1])
sd.r_24h_srot.g2<-c(g2.sd[2,1,2,1],g2.sd[3,1,2,2],g2.sd[5,5:7,2,1])
sd.r.g2<-as.matrix(rbind(sd.r_4h_rot.g2,sd.r_24h_rot.g2, sd.r_4h_srot.g2,sd.r_24h_srot.g2));colnames(sd.r.g2)<-c("CN","CP","0.5","1","2")
n.r.g2<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[5,5:7,2,2]))
#Withaferin A
m.WA_4h_rot.g2<-c(g2.mean[2,1,1,1],g2.mean[3,1,1,2],g2.mean[6,2:4,1,2])
m.WA_24h_rot.g2<-c(g2.mean[2,1,2,1],g2.mean[3,1,2,2],g2.mean[6,2:4,2,2])
m.WA_4h_srot.g2<-c(g2.mean[2,1,1,1],g2.mean[3,1,1,2],g2.mean[6,2:4,1,1])
m.WA_24h_srot.g2<-c(g2.mean[2,1,2,1],g2.mean[3,1,2,2],g2.mean[6,2:4,2,1])
m.WA.g2<-as.matrix(rbind(m.WA_4h_rot.g2,m.WA_24h_rot.g2, m.WA_4h_srot.g2,m.WA_24h_srot.g2));colnames(m.WA.g2)<-c("CN","CP","0.5","1","2")
sd.WA_4h_rot.g2<-c(g2.sd[2,1,1,1],g2.sd[3,1,1,2],g2.sd[6,2:4,1,2])
sd.WA_24h_rot.g2<-c(g2.sd[2,1,2,1],g2.sd[3,1,2,2],g2.sd[6,2:4,2,2])
sd.WA_4h_srot.g2<-c(g2.sd[2,1,1,1],g2.sd[3,1,1,2],g2.sd[6,2:4,1,1])
sd.WA_24h_srot.g2<-c(g2.sd[2,1,2,1],g2.sd[3,1,2,2],g2.sd[6,2:4,2,1])
sd.WA.g2<-as.matrix(rbind(sd.WA_4h_rot.g2,sd.WA_24h_rot.g2, sd.WA_4h_srot.g2,sd.WA_24h_srot.g2));colnames(sd.WA.g2)<-c("CN","CP","0.5","1","2")
n.WA.g2<-as.matrix(c(n.cell[2,1,2,1],n.cell[3,1,2,2],n.cell[6,2:4,2,2]))
}

#subg1
par(mar=c(4.2,4,1.5,1.5))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,3,1,3))
plot.new()
text(0.5,0.5, "SWitA + Rotenona (4h)",cex=2,font=2)
graf.barras(v.respuesta=m.14.subg1[1,],l.incub="4",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=sd.14.sugb1[1,],n.vector=n.14.subg1)
plot.new()
text(0.5,0.5,"SWitA",cex=2,font=2)
graf.barras(v.respuesta=m.14.subg1[3,],l.incub="4",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=sd.14.sugb1[3,],n.vector=n.14.subg1)
graf.barras(v.respuesta=m.14.subg1[2,],l.incub="24",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=sd.14.sugb1[2,],n.vector=n.14.subg1)
graf.barras(v.respuesta=m.14.subg1[4,],l.incub="24",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=sd.14.sugb1[4,],n.vector=n.14.subg1)

par(mar=c(4.2,4,1.75,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,3,1,3))
plot.new()
text(0.5,0.5, "Resveratrol + Rotenona (4h)",cex=2,font=2)
graf.barras(v.respuesta=m.r.subg1[1,],l.incub="4",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=sd.r.subg1[1,],n.vector=n.r.sugb1)
plot.new()
text(0.5,0.5,"Resveratrol",cex=2,font=2)
graf.barras(v.respuesta=m.r.subg1[3,],l.incub="4",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=sd.r.subg1[3,],n.vector=n.r.sugb1)
graf.barras(v.respuesta=m.r.subg1[2,],l.incub="24",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=sd.r.subg1[2,],n.vector=n.r.sugb1)
graf.barras(v.respuesta=m.r.subg1[4,],l.incub="24",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=sd.r.subg1[4,],n.vector=n.r.sugb1)

par(mar=c(4.2,4,1.75,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,3,1,3))
plot.new()
text(0.5,0.5, "WitA + Rotenona (4h)",cex=2,font=2)
graf.barras(v.respuesta=m.WA.subG1[1,],l.incub="4",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=sd.WA.subg1[1,],n.vector=n.WA.subg1)
plot.new()
text(0.5,0.5,"WitA",cex=2,font=2)
graf.barras(v.respuesta=m.WA.subG1[3,],l.incub="4",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=sd.WA.subg1[3,],n.vector=n.WA.subg1)
graf.barras(v.respuesta=m.WA.subG1[2,],l.incub="24",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=sd.WA.subg1[2,],n.vector=n.WA.subg1)
graf.barras(v.respuesta=m.WA.subG1[4,],l.incub="24",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=sd.WA.subg1[4,],n.vector=n.WA.subg1)

#cellcycle
par(mar=c(6,6,1.75,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1.3,3,1.3,3))
plot.new()
text(0.5,0.5, "SWitA + Rotenona (4h)",cex=2,font=2)
graf.barras.cc(v.respuesta=c(m.14.g0[1,],m.14.s[1,],m.14.g2[1,]),l.incub="4",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=c(sd.14.g0[1,],sd.14.s[1,],sd.14.g2[1,]),n.vector=c(n.14.g0, n.14.s, n.14.g2))
plot.new()
text(0.5,0.5,"SWitA",cex=2,font=2)
graf.barras.cc(v.respuesta=c(m.14.g0[3,],m.14.s[3,],m.14.g2[3,]),l.incub="4",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=c(sd.14.g0[3,],sd.14.s[3,],sd.14.g2[3,]),n.vector=c(n.14.g0, n.14.s, n.14.g2))
graf.barras.cc(v.respuesta=c(m.14.g0[2,],m.14.s[2,],m.14.g2[2,]),l.incub="24",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=c(sd.14.g0[2,],sd.14.s[2,],sd.14.g2[2,]),n.vector=c(n.14.g0, n.14.s, n.14.g2))
graf.barras.cc(v.respuesta=c(m.14.g0[4,],m.14.s[4,],m.14.g2[4,]),l.incub="24",compuesto="14.2.2",etiq.y="Percentage",sd.respuesta=c(sd.14.g0[4,],sd.14.s[4,],sd.14.g2[4,]),n.vector=c(n.14.g0, n.14.s, n.14.g2))

par(mar=c(6,6,1.75,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1.3,3,1.3,3))
plot.new()
text(0.5,0.5, "Resveratrol + Rotenona (4h)",cex=2,font=2)
graf.barras.cc(v.respuesta=c(m.r.g0[1,],m.r.s[1,],m.r.g2[1,]),l.incub="4",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=c(sd.r.g0[1,],sd.r.s[1,],sd.r.g2[1,]),n.vector=c(n.r.g0, n.r.s, n.r.g2))
plot.new()
text(0.5,0.5,"Resveratrol",cex=2,font=2)
graf.barras.cc(v.respuesta=c(m.r.g0[3,],m.r.s[3,],m.r.g2[3,]),l.incub="4",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=c(sd.r.g0[3,],sd.r.s[3,],sd.r.g2[3,]),n.vector=c(n.r.g0, n.r.s, n.r.g2))
graf.barras.cc(v.respuesta=c(m.r.g0[2,],m.r.s[2,],m.r.g2[2,]),l.incub="24",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=c(sd.r.g0[2,],sd.r.s[2,],sd.r.g2[2,]),n.vector=c(n.r.g0, n.r.s, n.r.g2))
graf.barras.cc(v.respuesta=c(m.r.g0[4,],m.r.s[4,],m.r.g2[4,]),l.incub="24",compuesto="Resveratrol",etiq.y="Percentage",sd.respuesta=c(sd.r.g0[4,],sd.r.s[4,],sd.r.g2[4,]),n.vector=c(n.r.g0, n.r.s, n.r.g2))

par(mar=c(6,6,1.75,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1.3,3,1.3,3))
plot.new()
text(0.5,0.5, "WitA + Rotenona (4h)",cex=2,font=2)
graf.barras.cc(v.respuesta=c(m.WA.g0[1,],m.WA.s[1,],m.WA.g2[1,]),l.incub="4",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=c(sd.WA.g0[1,],sd.WA.s[1,],sd.WA.g2[1,]),n.vector=c(n.WA.g0, n.WA.s, n.WA.g2))
plot.new()
text(0.5,0.5,"WitA",cex=2,font=2)
graf.barras.cc(v.respuesta=c(m.WA.g0[3,],m.WA.s[3,],m.WA.g2[3,]),l.incub="4",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=c(sd.WA.g0[3,],sd.WA.s[3,],sd.WA.g2[3,]),n.vector=c(n.WA.g0, n.WA.s, n.WA.g2))
graf.barras.cc(v.respuesta=c(m.WA.g0[2,],m.WA.s[2,],m.WA.g2[2,]),l.incub="24",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=c(sd.WA.g0[2,],sd.WA.s[2,],sd.WA.g2[2,]),n.vector=c(n.WA.g0, n.WA.s, n.WA.g2))
graf.barras.cc(v.respuesta=c(m.WA.g0[4,],m.WA.s[4,],m.WA.g2[4,]),l.incub="24",compuesto="Withaferin A",etiq.y="Percentage",sd.respuesta=c(sd.WA.g0[4,],sd.WA.s[4,],sd.WA.g2[4,]),n.vector=c(n.WA.g0, n.WA.s, n.WA.g2))



#6. Analisis estadistico

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#--------------------14.2.2 4 horas sin rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 4h SubG1",cex=2,font=2)
{
  datos.14.2.2.4h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="4"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.4h.sinrot.mfi.subg1<-lm(SubG1~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #verificar la interacción
  mod.14.2.2.4h.sinrot.mfi2.subg1<-lm(SubG1~tratam*Experimento,data=datos.14.2.2.4h.sinrot)
  anova(mod.14.2.2.4h.sinrot.mfi.subg1,mod.14.2.2.4h.sinrot.mfi2.subg1)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.sinrot.mfi.subg1)
  #se rechaza Ho cuando es heterocedastico
  bptest(SubG1~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.sinrot.mfi.subg1$residuals, "pnorm", 0, sd(mod.14.2.2.4h.sinrot.mfi.subg1$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.sinrot.mfi.subg1)
  #tky.14.2.2.4h.sinrot.mfi.subg1<-glht(mod.14.2.2.4h.sinrot.mfi2.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.14.2.2.4h.sinrot.mfi.subg1)
  interaction.plot(x.factor = datos.14.2.2.4h.sinrot$tratam,
                   trace.factor = datos.14.2.2.4h.sinrot$Experimento, 
                   response = datos.14.2.2.4h.sinrot$SubG1, fun = mean)
}

#--------------------14.2.2 24 horas sin rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 24h SubG1",cex=2,font=2)
{
  datos.14.2.2.24h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="24"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.24h.sinrot.mfi.subg1<-lm(SubG1~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #verificar la interacción
  mod.14.2.2.24h.sinrot.mfi2.subg1<-lm(SubG1~tratam*Experimento,data=datos.14.2.2.24h.sinrot)
  anova(mod.14.2.2.24h.sinrot.mfi.subg1,mod.14.2.2.24h.sinrot.mfi2.subg1)
  #se rechaza Ho cuando es heterocedastico
  bptest(SubG1~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.sinrot.mfi.subg1$residuals, "pnorm", 0, sd(mod.14.2.2.24h.sinrot.mfi.subg1$residuals))
  #verificar resultado
  anova(mod.14.2.2.24h.sinrot.mfi.subg1)
  tky.14.2.2.24h.sinrot.mfi.subg1<-glht(mod.14.2.2.24h.sinrot.mfi.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.24h.sinrot.mfi.subg1)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.sinrot.mfi.subg1)
  interaction.plot(x.factor = datos.14.2.2.24h.sinrot$tratam,
                   trace.factor = datos.14.2.2.24h.sinrot$Experimento, 
                   response = datos.14.2.2.24h.sinrot$SubG1, fun = mean)
}

#--------------------14.2.2 4 horas con rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 4h SubG1",cex=2,font=2)
{
  datos.14.2.2.4h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="4"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.4h.conrot.mfi.subg1<-lm(SubG1~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #verificar la interacción
  mod.14.2.2.4h.conrot.mfi2.subg1<-lm(SubG1~tratam*Experimento,data=datos.14.2.2.4h.conrot)
  anova(mod.14.2.2.4h.conrot.mfi.subg1,mod.14.2.2.4h.conrot.mfi2.subg1)
  #se rechaza Ho cuando es heterocedastico
  bptest(SubG1~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.conrot.mfi.subg1$residuals, "pnorm", 0, sd(mod.14.2.2.4h.conrot.mfi.subg1$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.conrot.mfi.subg1)
  #tky.14.2.2.4h.conrot.mfi.subg1<-glht(mod.14.2.2.4h.conrot.mfi.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.14.2.2.4h.conrot.mfi.subg1)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.conrot.mfi.subg1)
  interaction.plot(x.factor = datos.14.2.2.4h.conrot$tratam,
                   trace.factor = datos.14.2.2.4h.conrot$Experimento, 
                   response = datos.14.2.2.4h.conrot$SubG1, fun = mean)
}

#--------------------14.2.2 24 horas con rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 24h SubG1",cex=2,font=2)
{
  datos.14.2.2.24h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="24"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.24h.conrot.mfi.subg1<-lm(SubG1~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #verificar la interacción
  mod.14.2.2.24h.conrot.mfi2.subg1<-lm(SubG1~tratam*Experimento,data=datos.14.2.2.24h.conrot)
  anova(mod.14.2.2.24h.conrot.mfi.subg1,mod.14.2.2.24h.conrot.mfi2.subg1)
  #se rechaza Ho cuando es heterocedastico
  bptest(SubG1~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.conrot.mfi.subg1$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.mfi.subg1$residuals))
  #TRANSFORMACIÓN
  plot(boxcox(mod.14.2.2.24h.conrot.mfi.subg1))
  boxcox(mod.14.2.2.24h.conrot.mfi.subg1)$x[which.max(boxcox(mod.14.2.2.24h.conrot.mfi.subg1)$y)]
  mod.14.2.2.24h.conrot.mfit.subg1<-lm(SubG1^-0.46~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se rechaza Ho cuando es heterocedastico
  bptest(SubG1^-0.46~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.conrot.mfit.subg1$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.mfit.subg1$residuals))
  #verificar resultado
  anova(mod.14.2.2.24h.conrot.mfit.subg1)
  tky.14.2.2.24h.conrot.mfi.subg1<-glht(mod.14.2.2.24h.conrot.mfit.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.24h.conrot.mfi.subg1)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.conrot.mfi.subg1)
  interaction.plot(x.factor = datos.14.2.2.24h.conrot$tratam,
                   trace.factor = datos.14.2.2.24h.conrot$Experimento, 
                   response = datos.14.2.2.24h.conrot$SubG1, fun = mean)
}

#--------------------14.2.2 4 horas sin rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 4h G0",cex=2,font=2)
{
  datos.14.2.2.4h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="4"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.4h.sinrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #verificar la interacción
  mod.14.2.2.4h.sinrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.14.2.2.4h.sinrot)
  anova(mod.14.2.2.4h.sinrot.mfi.g0,mod.14.2.2.4h.sinrot.mfi2.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.sinrot.mfi.g0$residuals, "pnorm", 0, sd(mod.14.2.2.4h.sinrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.sinrot.mfi.g0)
  tky.14.2.2.4h.sinrot.mfi.g0<-glht(mod.14.2.2.4h.sinrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.4h.sinrot.mfi.g0)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.sinrot.mfi.g0)
  interaction.plot(x.factor = datos.14.2.2.4h.sinrot$tratam,
                   trace.factor = datos.14.2.2.4h.sinrot$Experimento, 
                   response = datos.14.2.2.4h.sinrot$G0, fun = mean)
}

#--------------------14.2.2 24 horas sin rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 24h G0",cex=2,font=2)
{
  datos.14.2.2.24h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="24"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.24h.sinrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #verificar la interacción
  mod.14.2.2.24h.sinrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.14.2.2.24h.sinrot)
  anova(mod.14.2.2.24h.sinrot.mfi.g0,mod.14.2.2.24h.sinrot.mfi2.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.sinrot.mfi.g0$residuals, "pnorm", 0, sd(mod.14.2.2.24h.sinrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.14.2.2.24h.sinrot.mfi.g0)
  #tky.14.2.2.24h.sinrot.mfi.g0<-glht(mod.14.2.2.24h.sinrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.14.2.2.24h.sinrot.mfi.g0)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.sinrot.mfi.g0)
  interaction.plot(x.factor = datos.14.2.2.24h.sinrot$tratam,
                   trace.factor = datos.14.2.2.24h.sinrot$Experimento, 
                   response = datos.14.2.2.24h.sinrot$G0, fun = mean)
}

#--------------------14.2.2 4 horas con rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 4h G0",cex=2,font=2)
{
  datos.14.2.2.4h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="4"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.4h.conrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #verificar la interacción
  mod.14.2.2.4h.conrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.14.2.2.4h.conrot)
  anova(mod.14.2.2.4h.conrot.mfi.g0,mod.14.2.2.4h.conrot.mfi2.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.conrot.mfi.g0$residuals, "pnorm", 0, sd(mod.14.2.2.4h.conrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.conrot.mfi.g0)
  tky.14.2.2.4h.conrot.mfi.g0<-glht(mod.14.2.2.4h.conrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.4h.conrot.mfi.g0)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.conrot.mfi.g0)
  interaction.plot(x.factor = datos.14.2.2.4h.conrot$tratam,
                   trace.factor = datos.14.2.2.4h.conrot$Experimento, 
                   response = datos.14.2.2.4h.conrot$G0, fun = mean)
}

#--------------------14.2.2 24 horas con rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 24h G0",cex=2,font=2)
{
  datos.14.2.2.24h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="24"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.24h.conrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #verificar la interacción
  mod.14.2.2.24h.conrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.14.2.2.24h.conrot)
  anova(mod.14.2.2.24h.conrot.mfi.g0,mod.14.2.2.24h.conrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.conrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.conrot.mfi.g0$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.14.2.2.24h.conrot.mfi.g0)
  tky.14.2.2.24h.conrot.mfi.g0<-glht(mod.14.2.2.24h.conrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.24h.conrot.mfi.g0)
  interaction.plot(x.factor = datos.14.2.2.24h.conrot$tratam,
                   trace.factor = datos.14.2.2.24h.conrot$Experimento, 
                   response = datos.14.2.2.24h.conrot$G0, fun = mean)
}

#--------------------14.2.2 4 horas sin rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 4h S",cex=2,font=2)
{
  datos.14.2.2.4h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="4"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.4h.sinrot.mfi.s<-lm(S~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #verificar la interacción
  mod.14.2.2.4h.sinrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.14.2.2.4h.sinrot)
  anova(mod.14.2.2.4h.sinrot.mfi.s,mod.14.2.2.4h.sinrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.sinrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.sinrot.mfi.s$residuals, "pnorm", 0, sd(mod.14.2.2.4h.sinrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.sinrot.mfi.s)
  #tky.14.2.2.4h.sinrot.mfi.s<-glht(mod.14.2.2.4h.sinrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.14.2.2.4h.sinrot.mfi.s)
  interaction.plot(x.factor = datos.14.2.2.4h.sinrot$tratam,
                   trace.factor = datos.14.2.2.4h.sinrot$Experimento, 
                   response = datos.14.2.2.4h.sinrot$S, fun = mean)
}

#--------------------14.2.2 24 horas sin rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 24h S",cex=2,font=2)
{
  datos.14.2.2.24h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="24"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.24h.sinrot.mfi.s<-lm(S~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #verificar la interacción
  mod.14.2.2.24h.sinrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.14.2.2.24h.sinrot)
  anova(mod.14.2.2.24h.sinrot.mfi.s,mod.14.2.2.24h.sinrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.sinrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.sinrot.mfi.s$residuals, "pnorm", 0, sd(mod.14.2.2.24h.sinrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.14.2.2.24h.sinrot.mfi.s)
  #tky.14.2.2.24h.sinrot.mfi.s<-glht(mod.14.2.2.24h.sinrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.14.2.2.24h.sinrot.mfi.s)
  interaction.plot(x.factor = datos.14.2.2.24h.sinrot$tratam,
                   trace.factor = datos.14.2.2.24h.sinrot$Experimento, 
                   response = datos.14.2.2.24h.sinrot$S, fun = mean)
}

#--------------------14.2.2 4 horas con rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 4h S",cex=2,font=2)
{
  datos.14.2.2.4h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="4"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.4h.conrot.mfi.s<-lm(S~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #verificar la interacción
  mod.14.2.2.4h.conrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.14.2.2.4h.conrot)
  anova(mod.14.2.2.4h.conrot.mfi.s,mod.14.2.2.4h.conrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.conrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.conrot.mfi.s$residuals, "pnorm", 0, sd(mod.14.2.2.4h.conrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.conrot.mfi.s)
  #tky.14.2.2.4h.conrot.mfi.s<-glht(mod.14.2.2.4h.conrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.14.2.2.4h.conrot.mfi.s)
  interaction.plot(x.factor = datos.14.2.2.4h.conrot$tratam,
                   trace.factor = datos.14.2.2.4h.conrot$Experimento, 
                   response = datos.14.2.2.4h.conrot$S, fun = mean)
}

#--------------------14.2.2 24 horas con rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 24h S",cex=2,font=2)
{
  datos.14.2.2.24h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="24"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.24h.conrot.mfi.s<-lm(S~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #verificar la interacción
  mod.14.2.2.24h.conrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.14.2.2.24h.conrot)
  anova(mod.14.2.2.24h.conrot.mfi.s,mod.14.2.2.24h.conrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.conrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.conrot.mfi.s$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.mfi.s$residuals))
  #TRANSFORMACIÓN
  plot(boxcox(mod.14.2.2.24h.conrot.mfi.s))
  boxcox(mod.14.2.2.24h.conrot.mfi.s)$x[which.max(boxcox(mod.14.2.2.24h.conrot.mfi.s)$y)]
  mod.14.2.2.24h.conrot.mfit.s<-lm(S^-0.06~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se rechaza Ho cuando es heterocedastico
  bptest(S^-0.06~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.conrot.mfit.s$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.mfit.s$residuals))
   #verificar resultado
  anova(mod.14.2.2.24h.conrot.mfit.s)
  #tky.14.2.2.24h.conrot.mfi2.s<-glht(mod.14.2.2.24h.conrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.14.2.2.24h.conrot.mfi2.s)
  interaction.plot(x.factor = datos.14.2.2.24h.conrot$tratam,
                   trace.factor = datos.14.2.2.24h.conrot$Experimento, 
                   response = datos.14.2.2.24h.conrot$S, fun = mean)
}

#--------------------14.2.2 4 horas sin rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 4h G2",cex=2,font=2)
{
  datos.14.2.2.4h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="4"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.4h.sinrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #verificar la interacción
  mod.14.2.2.4h.sinrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.14.2.2.4h.sinrot)
  anova(mod.14.2.2.4h.sinrot.mfi.g2,mod.14.2.2.4h.sinrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.sinrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.sinrot.mfi.g2$residuals, "pnorm", 0, sd(mod.14.2.2.4h.sinrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.sinrot.mfi.g2)
  tky.14.2.2.4h.sinrot.mfi.g2<-glht(mod.14.2.2.4h.sinrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.4h.sinrot.mfi.g2)
  interaction.plot(x.factor = datos.14.2.2.4h.sinrot$tratam,
                   trace.factor = datos.14.2.2.4h.sinrot$Experimento, 
                   response = datos.14.2.2.4h.sinrot$G2, fun = mean)
}

#--------------------14.2.2 24 horas sin rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 24h G2",cex=2,font=2)
{
  datos.14.2.2.24h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="24"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
  mod.14.2.2.24h.sinrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #verificar la interacción
  mod.14.2.2.24h.sinrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.14.2.2.24h.sinrot)
  anova(mod.14.2.2.24h.sinrot.mfi.g2,mod.14.2.2.24h.sinrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.sinrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.sinrot.mfi.g2$residuals, "pnorm", 0, sd(mod.14.2.2.24h.sinrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.14.2.2.24h.sinrot.mfi.g2)
  tky.14.2.2.24h.sinrot.mfi.g2<-glht(mod.14.2.2.24h.sinrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.24h.sinrot.mfi.g2)
  interaction.plot(x.factor = datos.14.2.2.24h.sinrot$tratam,
                   trace.factor = datos.14.2.2.24h.sinrot$Experimento, 
                   response = datos.14.2.2.24h.sinrot$G2, fun = mean)
}

#--------------------14.2.2 4 horas con rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 4h G2",cex=2,font=2)
{
  datos.14.2.2.4h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="4"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.4h.conrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #verificar la interacción
  mod.14.2.2.4h.conrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.14.2.2.4h.conrot)
  anova(mod.14.2.2.4h.conrot.mfi.g2,mod.14.2.2.4h.conrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.14.2.2.4h.conrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.14.2.2.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.4h.conrot.mfi.g2$residuals, "pnorm", 0, sd(mod.14.2.2.4h.conrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.14.2.2.4h.conrot.mfi.g2)
  tky.14.2.2.4h.conrot.mfi.g2<-glht(mod.14.2.2.4h.conrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.4h.conrot.mfi.g2)
  interaction.plot(x.factor = datos.14.2.2.4h.conrot$tratam,
                   trace.factor = datos.14.2.2.4h.conrot$Experimento, 
                   response = datos.14.2.2.4h.conrot$G2, fun = mean)
}

#--------------------14.2.2 24 horas con rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 24h G2",cex=2,font=2)
{
  datos.14.2.2.24h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="24"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
  mod.14.2.2.24h.conrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #verificar la interacción
  mod.14.2.2.24h.conrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.14.2.2.24h.conrot)
  anova(mod.14.2.2.24h.conrot.mfi.g2,mod.14.2.2.24h.conrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.14.2.2.24h.conrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.14.2.2.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.14.2.2.24h.conrot.mfi.g2$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.14.2.2.24h.conrot.mfi.g2)
  tky.14.2.2.24h.conrot.mfi.g2<-glht(mod.14.2.2.24h.conrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.14.2.2.24h.conrot.mfi.g2)
  interaction.plot(x.factor = datos.14.2.2.24h.conrot$tratam,
                   trace.factor = datos.14.2.2.24h.conrot$Experimento, 
                   response = datos.14.2.2.24h.conrot$G2, fun = mean)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------resv 4 horas sin rot--------------------------#
#--------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 4h subG1",cex=2,font=2)
{
datos.resv.4h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="R1"|tratam=="R2"|tratam=="R3")&Incubacion=="4"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
mod.resv.4h.sinrot.subg1<-lm(SubG1~tratam+Experimento,data=datos.resv.4h.sinrot)
#verificar la interacción
mod.resv.4h.sinrot.subg2<-lm(SubG1~tratam*Experimento,data=datos.resv.4h.sinrot)
anova(mod.resv.4h.sinrot.subg1,mod.resv.4h.sinrot.subg2)
#gráficos de supuestos
plot(mod.resv.4h.sinrot.subg1)
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.resv.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.sinrot.subg1$residuals, "pnorm", 0, sd(mod.resv.4h.sinrot.subg1$residuals))
#verificar resultado
anova(mod.resv.4h.sinrot.subg1)
res.tky.4h.sinrot.subg1<-glht(mod.resv.4h.sinrot.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.4h.sinrot.subg1)
interaction.plot(x.factor = datos.resv.4h.sinrot$tratam,
                 trace.factor = datos.resv.4h.sinrot$Experimento, 
                 response = datos.resv.4h.sinrot$SubG1, fun = mean)
}

#--------------------resv 24 horas sin rot--------------------------#
#--------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 24h subG1",cex=2,font=2)
{
  datos.resv.24h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="R1"|tratam=="R2"|tratam=="R3")&Incubacion=="24"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
mod.resv.24h.sinrot.subg1<-lm(SubG1~tratam+Experimento,data=datos.resv.24h.sinrot)
#verificar la interacción
mod.resv.24h.sinrot.subg2<-lm(SubG1~tratam*Experimento,data=datos.resv.24h.sinrot)
anova(mod.resv.24h.sinrot.subg1,mod.resv.24h.sinrot.subg2)
#gráficos de supuestos
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.resv.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.24h.sinrot.subg1$residuals, "pnorm", 0, sd(mod.resv.24h.sinrot.subg1$residuals))
#verificar resultado
anova(mod.resv.24h.sinrot.subg1)
res.tky.24h.sinrot.subg1<-glht(mod.resv.24h.sinrot.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.24h.sinrot.subg1)
plot(mod.resv.24h.sinrot.subg1)
interaction.plot(x.factor = datos.resv.24h.sinrot$tratam,
                 trace.factor = datos.resv.24h.sinrot$Experimento, 
                 response = datos.resv.24h.sinrot$SubG1, fun = mean)
}

#--------------------resv 4 horas con rot--------------------------#
#--------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 4h subG1",cex=2,font=2)
{
datos.resv.4h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="RR1"|tratam=="RR2"|tratam=="RR3")&Incubacion=="4"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
mod.resv.4h.conrot.subg1<-lm(SubG1~tratam+Experimento,data=datos.resv.4h.conrot)
#verificar la interacción
mod.resv.4h.conrot.subg2<-lm(SubG1~tratam*Experimento,data=datos.resv.4h.conrot)
anova(mod.resv.4h.conrot.subg1,mod.resv.4h.conrot.subg2)
#gráficos de supuestos
plot(mod.resv.4h.conrot.subg1)
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.resv.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.conrot.subg1$residuals, "pnorm", 0, sd(mod.resv.4h.conrot.subg1$residuals))
#verificar resultado
anova(mod.resv.4h.conrot.subg1)
#res.tky.4h.conrot.subg1<-glht(mod.resv.4h.conrot.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
#summary(res.tky.4h.conrot.subg1)
interaction.plot(x.factor = datos.resv.4h.conrot$tratam,
                 trace.factor = datos.resv.4h.conrot$Experimento, 
                 response = datos.resv.4h.conrot$SubG1, fun = mean)
}

#--------------------resv 24 horas con rot--------------------------#
#--------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 24h subG1",cex=2,font=2)
{
datos.resv.24h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="RR1"|tratam=="RR2"|tratam=="RR3")&Incubacion=="24"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
mod.resv.24h.conrot.subg1<-lm(SubG1~tratam+Experimento,data=datos.resv.24h.conrot)
#verificar la interacción
mod.resv.24h.conrot.subg2<-lm(SubG1~tratam*Experimento,data=datos.resv.24h.conrot)
anova(mod.resv.24h.conrot.subg1,mod.resv.24h.conrot.subg2)
#gráficos de supuestos
plot(mod.resv.24h.conrot.subg1)
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.resv.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.24h.conrot.subg1$residuals, "pnorm", 0, sd(mod.resv.24h.conrot.subg1$residuals))
#verificar resultado
anova(mod.resv.24h.conrot.subg1)
res.tky.24h.conrot.subg1<-glht(mod.resv.24h.conrot.subg1, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.24h.conrot.subg1)
interaction.plot(x.factor = datos.resv.24h.conrot$tratam,
                 trace.factor = datos.resv.24h.conrot$Experimento, 
                 response = datos.resv.24h.conrot$SubG1, fun = mean)
}


#--------------------resv 4 horas sin rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 4h G0",cex=2,font=2)
{
  mod.resv.4h.sinrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.resv.4h.sinrot)
  #verificar la interacción
  mod.resv.4h.sinrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.resv.4h.sinrot)
  anova(mod.resv.4h.sinrot.mfi.g0,mod.resv.4h.sinrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.resv.4h.sinrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.resv.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.4h.sinrot.mfi.g0$residuals, "pnorm", 0, sd(mod.resv.4h.sinrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.resv.4h.sinrot.mfi.g0)
  tky.resv.4h.sinrot.mfi.g0<-glht(mod.resv.4h.sinrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.4h.sinrot.mfi.g0)
  interaction.plot(x.factor = datos.resv.4h.sinrot$tratam,
                   trace.factor = datos.resv.4h.sinrot$Experimento, 
                   response = datos.resv.4h.sinrot$G0, fun = mean)
}

#--------------------resv 24 horas sin rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 24h G0",cex=2,font=2)
{
  mod.resv.24h.sinrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.resv.24h.sinrot)
  #verificar la interacción
  mod.resv.24h.sinrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.resv.24h.sinrot)
  anova(mod.resv.24h.sinrot.mfi.g0,mod.resv.24h.sinrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.resv.24h.sinrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.resv.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.24h.sinrot.mfi.g0$residuals, "pnorm", 0, sd(mod.resv.24h.sinrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.resv.24h.sinrot.mfi.g0)
  tky.resv.24h.sinrot.mfi.g0<-glht(mod.resv.24h.sinrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.24h.sinrot.mfi.g0)
  interaction.plot(x.factor = datos.resv.24h.sinrot$tratam,
                   trace.factor = datos.resv.24h.sinrot$Experimento, 
                   response = datos.resv.24h.sinrot$G0, fun = mean)
}

#--------------------resv 4 horas con rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 4h G0",cex=2,font=2)
{
  mod.resv.4h.conrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.resv.4h.conrot)
  #verificar la interacción
  mod.resv.4h.conrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.resv.4h.conrot)
  anova(mod.resv.4h.conrot.mfi.g0,mod.resv.4h.conrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.resv.4h.conrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.resv.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.4h.conrot.mfi.g0$residuals, "pnorm", 0, sd(mod.resv.4h.conrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.resv.4h.conrot.mfi.g0)
  #tky.resv.4h.conrot.mfi.g0<-glht(mod.resv.4h.conrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.resv.4h.conrot.mfi.g0)
  interaction.plot(x.factor = datos.resv.4h.conrot$tratam,
                   trace.factor = datos.resv.4h.conrot$Experimento, 
                   response = datos.resv.4h.conrot$G0, fun = mean)
}

#--------------------resv 24 horas con rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 24h G0",cex=2,font=2)
{
  mod.resv.24h.conrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.resv.24h.conrot)
  #verificar la interacción
  mod.resv.24h.conrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.resv.24h.conrot)
  anova(mod.resv.24h.conrot.mfi.g0,mod.resv.24h.conrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.resv.24h.conrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.resv.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.24h.conrot.mfi.g0$residuals, "pnorm", 0, sd(mod.resv.24h.conrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.resv.24h.conrot.mfi.g0)
  tky.resv.24h.conrot.mfi.g0<-glht(mod.resv.24h.conrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.24h.conrot.mfi.g0)
  interaction.plot(x.factor = datos.resv.24h.conrot$tratam,
                   trace.factor = datos.resv.24h.conrot$Experimento, 
                   response = datos.resv.24h.conrot$G0, fun = mean)
}

#--------------------resv 4 horas sin rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 4h S",cex=2,font=2)
{
  mod.resv.4h.sinrot.mfi.s<-lm(S~tratam+Experimento,data=datos.resv.4h.sinrot)
  #verificar la interacción
  mod.resv.4h.sinrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.resv.4h.sinrot)
  anova(mod.resv.4h.sinrot.mfi.s,mod.resv.4h.sinrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.resv.4h.sinrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.resv.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.4h.sinrot.mfi.s$residuals, "pnorm", 0, sd(mod.resv.4h.sinrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.resv.4h.sinrot.mfi.s)
  tky.resv.4h.sinrot.mfi.s<-glht(mod.resv.4h.sinrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.4h.sinrot.mfi.s)
  interaction.plot(x.factor = datos.resv.4h.sinrot$tratam,
                   trace.factor = datos.resv.4h.sinrot$Experimento, 
                   response = datos.resv.4h.sinrot$S, fun = mean)
}

#--------------------resv 24 horas sin rot--------------------------#
#---------------------------S-------------------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 24h S",cex=2,font=2)
{
  mod.resv.24h.sinrot.mfi.s<-lm(S~tratam+Experimento,data=datos.resv.24h.sinrot)
  #verificar la interacción
  mod.resv.24h.sinrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.resv.24h.sinrot)
  anova(mod.resv.24h.sinrot.mfi.s,mod.resv.24h.sinrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.resv.24h.sinrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.resv.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.24h.sinrot.mfi.s$residuals, "pnorm", 0, sd(mod.resv.24h.sinrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.resv.24h.sinrot.mfi.s)
  tky.resv.24h.sinrot.mfi.s<-glht(mod.resv.24h.sinrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.24h.sinrot.mfi.s)
  interaction.plot(x.factor = datos.resv.24h.sinrot$tratam,
                   trace.factor = datos.resv.24h.sinrot$Experimento, 
                   response = datos.resv.24h.sinrot$S, fun = mean)
}

#--------------------resv 4 horas con rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 4h S",cex=2,font=2)
{
  mod.resv.4h.conrot.mfi.s<-lm(S~tratam+Experimento,data=datos.resv.4h.conrot)
  #verificar la interacción
  mod.resv.4h.conrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.resv.4h.conrot)
  anova(mod.resv.4h.conrot.mfi.s,mod.resv.4h.conrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.resv.4h.conrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.resv.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.4h.conrot.mfi.s$residuals, "pnorm", 0, sd(mod.resv.4h.conrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.resv.4h.conrot.mfi.s)
  tky.resv.4h.conrot.mfi.s<-glht(mod.resv.4h.conrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.4h.conrot.mfi.s)
  interaction.plot(x.factor = datos.resv.4h.conrot$tratam,
                   trace.factor = datos.resv.4h.conrot$Experimento, 
                   response = datos.resv.4h.conrot$S, fun = mean)
}

#--------------------resv 24 horas con rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 24h S",cex=2,font=2)
{
  mod.resv.24h.conrot.mfi.s<-lm(S~tratam+Experimento,data=datos.resv.24h.conrot)
  #verificar la interacción
  mod.resv.24h.conrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.resv.24h.conrot)
  anova(mod.resv.24h.conrot.mfi.s,mod.resv.24h.conrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.resv.24h.conrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.resv.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.24h.conrot.mfi.s$residuals, "pnorm", 0, sd(mod.resv.24h.conrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.resv.24h.conrot.mfi.s)
  tky.resv.24h.conrot.mfi.s<-glht(mod.resv.24h.conrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.24h.conrot.mfi.s)
  interaction.plot(x.factor = datos.resv.24h.conrot$tratam,
                   trace.factor = datos.resv.24h.conrot$Experimento, 
                   response = datos.resv.24h.conrot$S, fun = mean)
}

#--------------------resv 4 horas sin rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 4h G2",cex=2,font=2)
{
  mod.resv.4h.sinrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.resv.4h.sinrot)
  #verificar la interacción
  mod.resv.4h.sinrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.resv.4h.sinrot)
  anova(mod.resv.4h.sinrot.mfi.g2,mod.resv.4h.sinrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.resv.4h.sinrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.resv.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.4h.sinrot.mfi.g2$residuals, "pnorm", 0, sd(mod.resv.4h.sinrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.resv.4h.sinrot.mfi.g2)
  tky.resv.4h.sinrot.mfi.g2<-glht(mod.resv.4h.sinrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.4h.sinrot.mfi.g2)
  interaction.plot(x.factor = datos.resv.4h.sinrot$tratam,
                   trace.factor = datos.resv.4h.sinrot$Experimento, 
                   response = datos.resv.4h.sinrot$G2, fun = mean)
}

#--------------------resv 24 horas sin rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 24h G2",cex=2,font=2)
{
  mod.resv.24h.sinrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.resv.24h.sinrot)
  #verificar la interacción
  mod.resv.24h.sinrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.resv.24h.sinrot)
  anova(mod.resv.24h.sinrot.mfi.g2,mod.resv.24h.sinrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.resv.24h.sinrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.resv.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.24h.sinrot.mfi.g2$residuals, "pnorm", 0, sd(mod.resv.24h.sinrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.resv.24h.sinrot.mfi.g2)
  tky.resv.24h.sinrot.mfi.g2<-glht(mod.resv.24h.sinrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.24h.sinrot.mfi.g2)
  interaction.plot(x.factor = datos.resv.24h.sinrot$tratam,
                   trace.factor = datos.resv.24h.sinrot$Experimento, 
                   response = datos.resv.24h.sinrot$G2, fun = mean)
}

#--------------------resv 4 horas con rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 4h G2",cex=2,font=2)
{
  mod.resv.4h.conrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.resv.4h.conrot)
  #verificar la interacción
  mod.resv.4h.conrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.resv.4h.conrot)
  anova(mod.resv.4h.conrot.mfi.g2,mod.resv.4h.conrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.resv.4h.conrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.resv.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.4h.conrot.mfi.g2$residuals, "pnorm", 0, sd(mod.resv.4h.conrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.resv.4h.conrot.mfi.g2)
  tky.resv.4h.conrot.mfi.g2<-glht(mod.resv.4h.conrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.4h.conrot.mfi.g2)
  interaction.plot(x.factor = datos.resv.4h.conrot$tratam,
                   trace.factor = datos.resv.4h.conrot$Experimento, 
                   response = datos.resv.4h.conrot$G2, fun = mean)
}

#--------------------resv 24 horas con rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 24h G2",cex=2,font=2)
{
  mod.resv.24h.conrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.resv.24h.conrot)
  #verificar la interacción
  mod.resv.24h.conrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.resv.24h.conrot)
  anova(mod.resv.24h.conrot.mfi.g2,mod.resv.24h.conrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.resv.24h.conrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.resv.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.resv.24h.conrot.mfi.g2$residuals, "pnorm", 0, sd(mod.resv.24h.conrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.resv.24h.conrot.mfi.g2)
  tky.resv.24h.conrot.mfi2.g2<-glht(mod.resv.24h.conrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.resv.24h.conrot.mfi2.g2)
  interaction.plot(x.factor = datos.resv.24h.conrot$tratam,
                   trace.factor = datos.resv.24h.conrot$Experimento, 
                   response = datos.resv.24h.conrot$G2, fun = mean)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------WitA 4 horas sin rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 4h SubG1",cex=2,font=2)
{
datos.wita.4h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="WA1"|tratam=="WA2"|tratam=="WA3")&Incubacion=="4"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
mod.wita.4h.sinrot.mfi<-lm(SubG1~tratam+Experimento,data=datos.wita.4h.sinrot)
#verificar la interacción
mod.wita.4h.sinrot.mfi2<-lm(SubG1~tratam*Experimento,data=datos.wita.4h.sinrot)
anova(mod.wita.4h.sinrot.mfi,mod.wita.4h.sinrot.mfi2)
#gráficos de supuestos
plot(mod.wita.4h.sinrot.mfi)
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.wita.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.wita.4h.sinrot.mfi$residuals))
#verificar resultado
anova(mod.wita.4h.sinrot.mfi)
#tky.wita.4h.sinrot.mfi<-glht(mod.wita.4h.sinrot.mfi, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
#summary(tky.wita.4h.sinrot.mfi)
interaction.plot(x.factor = datos.wita.4h.sinrot$tratam,
                 trace.factor = datos.wita.4h.sinrot$Experimento, 
                 response = datos.wita.4h.sinrot$SubG1, fun = mean)
}

#--------------------WitA 24 horas sin rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 24h SubG1",cex=2,font=2)
{
datos.wita.24h.sinrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="WA1"|tratam=="WA2"|tratam=="WA3")&Incubacion=="24"&(Experimento=="1"|Experimento=="2"|Experimento=="3"|Experimento=="4"),]
mod.wita.24h.sinrot.mfi<-lm(SubG1~tratam+Experimento,data=datos.wita.24h.sinrot)
#verificar la interacción
mod.wita.24h.sinrot.mfi2<-lm(SubG1~tratam*Experimento,data=datos.wita.24h.sinrot)
anova(mod.wita.24h.sinrot.mfi,mod.wita.24h.sinrot.mfi2)
#gráficos de supuestos
plot(mod.wita.24h.sinrot.mfi)
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.wita.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.24h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.wita.24h.sinrot.mfi$residuals))
#TRANSFORMACIÓN
#plot(boxcox(mod.wita.24h.sinrot.mfi))
#boxcox(mod.wita.24h.sinrot.mfi)$x[which.max(boxcox(mod.wita.24h.sinrot.mfi)$y)]
#mod.wita.24h.sinrot.mfit<-lm(Media.geometrica^0.02~tratam+Experimento,data=datos.wita.24h.sinrot)
#se rechaza Ho cuando es heterocedastico
#bptest(Media.geometrica^0.02~tratam+Experimento,data=datos.wita.24h.sinrot)
#se acepta Ho cuando es normal
#ks.test(mod.wita.24h.sinrot.mfit$residuals, "pnorm", 0, sd(mod.resv.24h.conrot.mfit$residuals))
#verificar resultado
anova(mod.wita.24h.sinrot.mfi)
tky.wita.24h.sinrot.mfi<-glht(mod.wita.24h.sinrot.mfi, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.24h.sinrot.mfi)
interaction.plot(x.factor = datos.wita.24h.sinrot$tratam,
                 trace.factor = datos.wita.24h.sinrot$Experimento, 
                 response = datos.wita.24h.sinrot$SubG1, fun = mean)

}

#--------------------WitA 4 horas con rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 4h SubG1",cex=2,font=2)
{
  datos.wita.4h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="WAR1"|tratam=="WAR2"|tratam=="WAR3")&Incubacion=="4"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
mod.wita.4h.conrot.mfi<-lm(SubG1~tratam+Experimento,data=datos.wita.4h.conrot)
#verificar la interacción
mod.wita.4h.conrot.mfi2<-lm(SubG1~tratam*Experimento,data=datos.wita.4h.conrot)
anova(mod.wita.4h.conrot.mfi,mod.wita.4h.conrot.mfi2)
#gráficos de supuestos
plot(mod.wita.4h.conrot.mfi)
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.wita.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.conrot.mfi$residuals, "pnorm", 0, sd(mod.wita.4h.conrot.mfi$residuals))
#verificar resultado
anova(mod.wita.4h.conrot.mfi)
#tky.wita.4h.conrot.mfi<-glht(mod.wita.4h.conrot.mfi, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
#summary(tky.wita.4h.conrot.mfi)
interaction.plot(x.factor = datos.wita.4h.conrot$tratam,
                 trace.factor = datos.wita.4h.conrot$Experimento, 
                 response = datos.wita.4h.conrot$SubG1, fun = mean)

}

#--------------------WitA 24 horas con rot--------------------------#
#---------------------------SubG1-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 24h SubG1",cex=2,font=2)
{
datos.wita.24h.conrot<-datos[(tratam=="CN"|tratam=="CP"|tratam=="WAR1"|tratam=="WAR2"|tratam=="WAR3")&Incubacion=="24"&(Experimento=="5"|Experimento=="6"|Experimento=="7"|Experimento=="8"),]
mod.wita.24h.conrot.mfi<-lm(SubG1~tratam+Experimento,data=datos.wita.24h.conrot)
#verificar la interacción
mod.wita.24h.conrot.mfi2<-lm(SubG1~tratam*Experimento,data=datos.wita.24h.conrot)
anova(mod.wita.24h.conrot.mfi,mod.wita.24h.conrot.mfi2)
#gráficos de supuestos
plot(mod.wita.24h.conrot.mfi)
#se rechaza Ho cuando es heterocedastico
bptest(SubG1~tratam+Experimento,data=datos.wita.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.24h.conrot.mfi$residuals, "pnorm", 0, sd(mod.wita.24h.conrot.mfi$residuals))
#verificar resultado
anova(mod.wita.24h.conrot.mfi)
tky.wita.24h.conrot.mfi<-glht(mod.wita.24h.conrot.mfi, linfct = mcp(tratam = "Dunnet"), alternative = "two.sided")
summary(tky.wita.24h.conrot.mfi)
interaction.plot(x.factor = datos.wita.24h.conrot$tratam,
                 trace.factor = datos.wita.24h.conrot$Experimento, 
                 response = datos.wita.24h.conrot$SubG1, fun = mean)
}

#--------------------wita 4 horas sin rot--------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 4h G0",cex=2,font=2)
{
  mod.wita.4h.sinrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.wita.4h.sinrot)
  #verificar la interacción
  mod.wita.4h.sinrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.wita.4h.sinrot)
  anova(mod.wita.4h.sinrot.mfi.g0,mod.wita.4h.sinrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.wita.4h.sinrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.wita.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.4h.sinrot.mfi.g0$residuals, "pnorm", 0, sd(mod.wita.4h.sinrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.wita.4h.sinrot.mfi.g0)
  tky.wita.4h.sinrot.mfi.g0<-glht(mod.wita.4h.sinrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.4h.sinrot.mfi.g0)
  interaction.plot(x.factor = datos.wita.4h.sinrot$tratam,
                   trace.factor = datos.wita.4h.sinrot$Experimento, 
                   response = datos.wita.4h.sinrot$G0, fun = mean)
}

#--------------------wita 24 horas sin rot--------------------------#
#---------------------------G0-------------------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 24h G0",cex=2,font=2)
{
  mod.wita.24h.sinrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.wita.24h.sinrot)
  #verificar la interacción
  mod.wita.24h.sinrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.wita.24h.sinrot)
  anova(mod.wita.24h.sinrot.mfi.g0,mod.wita.24h.sinrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.wita.24h.sinrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.wita.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.24h.sinrot.mfi.g0$residuals, "pnorm", 0, sd(mod.wita.24h.sinrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.wita.24h.sinrot.mfi.g0)
  tky.wita.24h.sinrot.mfi.g0<-glht(mod.wita.24h.sinrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.24h.sinrot.mfi.g0)
  interaction.plot(x.factor = datos.wita.24h.sinrot$tratam,
                   trace.factor = datos.wita.24h.sinrot$Experimento, 
                   response = datos.wita.24h.sinrot$G0, fun = mean)
}

#--------------------wita 4 horas con rot--------------------------#
#---------------------------G0-------------------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 4h G0",cex=2,font=2)
{
  mod.wita.4h.conrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.wita.4h.conrot)
  #verificar la interacción
  mod.wita.4h.conrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.wita.4h.conrot)
  anova(mod.wita.4h.conrot.mfi.g0,mod.wita.4h.conrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.wita.4h.conrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.wita.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.4h.conrot.mfi.g0$residuals, "pnorm", 0, sd(mod.wita.4h.conrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.wita.4h.conrot.mfi.g0)
  tky.wita.4h.conrot.mfi.g0<-glht(mod.wita.4h.conrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.4h.conrot.mfi.g0)
  interaction.plot(x.factor = datos.wita.4h.conrot$tratam,
                   trace.factor = datos.wita.4h.conrot$Experimento, 
                   response = datos.wita.4h.conrot$G0, fun = mean)
}

#--------------------wita 24 horas con rot--------------------------#
#---------------------------G0-------------------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 24h G0",cex=2,font=2)
{
  mod.wita.24h.conrot.mfi.g0<-lm(G0~tratam+Experimento,data=datos.wita.24h.conrot)
  #verificar la interacción
  mod.wita.24h.conrot.mfi2.g0<-lm(G0~tratam*Experimento,data=datos.wita.24h.conrot)
  anova(mod.wita.24h.conrot.mfi.g0,mod.wita.24h.conrot.mfi2.g0)
  #gráficos de supuestos
  plot(mod.wita.24h.conrot.mfi.g0)
  #se rechaza Ho cuando es heterocedastico
  bptest(G0~tratam+Experimento,data=datos.wita.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.24h.conrot.mfi.g0$residuals, "pnorm", 0, sd(mod.wita.24h.conrot.mfi.g0$residuals))
  #verificar resultado
  anova(mod.wita.24h.conrot.mfi.g0)
  tky.wita.24h.conrot.mfi2.g0<-glht(mod.wita.24h.conrot.mfi.g0, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.24h.conrot.mfi2.g0)
  interaction.plot(x.factor = datos.wita.24h.conrot$tratam,
                   trace.factor = datos.wita.24h.conrot$Experimento, 
                   response = datos.wita.24h.conrot$G0, fun = mean)
}

#--------------------wita 4 horas sin rot--------------------------#
#---------------------------S-------------------------------------#
#---------------------------G0-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 4h S",cex=2,font=2)
{
  mod.wita.4h.sinrot.mfi.s<-lm(S~tratam+Experimento,data=datos.wita.4h.sinrot)
  #verificar la interacción
  mod.wita.4h.sinrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.wita.4h.sinrot)
  anova(mod.wita.4h.sinrot.mfi.s,mod.wita.4h.sinrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.wita.4h.sinrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.wita.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.4h.sinrot.mfi.s$residuals, "pnorm", 0, sd(mod.wita.4h.sinrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.wita.4h.sinrot.mfi.s)
  tky.wita.4h.sinrot.mfi.s<-glht(mod.wita.4h.sinrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.4h.sinrot.mfi.s)
  interaction.plot(x.factor = datos.wita.4h.sinrot$tratam,
                   trace.factor = datos.wita.4h.sinrot$Experimento, 
                   response = datos.wita.4h.sinrot$S, fun = mean)
}

#--------------------wita 24 horas sin rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 24h S",cex=2,font=2)
{
  mod.wita.24h.sinrot.mfi.s<-lm(S~tratam+Experimento,data=datos.wita.24h.sinrot)
  #verificar la interacción
  mod.wita.24h.sinrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.wita.24h.sinrot)
  anova(mod.wita.24h.sinrot.mfi.s,mod.wita.24h.sinrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.wita.24h.sinrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.wita.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.24h.sinrot.mfi.s$residuals, "pnorm", 0, sd(mod.wita.24h.sinrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.wita.24h.sinrot.mfi.s)
  tky.wita.24h.sinrot.mfi.s<-glht(mod.wita.24h.sinrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.24h.sinrot.mfi.s)
  interaction.plot(x.factor = datos.wita.24h.sinrot$tratam,
                   trace.factor = datos.wita.24h.sinrot$Experimento, 
                   response = datos.wita.24h.sinrot$S, fun = mean)
}

#--------------------wita 4 horas con rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 4h S",cex=2,font=2)
{
  mod.wita.4h.conrot.mfi.s<-lm(S~tratam+Experimento,data=datos.wita.4h.conrot)
  #verificar la interacción
  mod.wita.4h.conrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.wita.4h.conrot)
  anova(mod.wita.4h.conrot.mfi.s,mod.wita.4h.conrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.wita.4h.conrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.wita.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.4h.conrot.mfi.s$residuals, "pnorm", 0, sd(mod.wita.4h.conrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.wita.4h.conrot.mfi.s)
  #tky.wita.4h.conrot.mfi.s<-glht(mod.wita.4h.conrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  #summary(tky.wita.4h.conrot.mfi.s)
  interaction.plot(x.factor = datos.wita.4h.conrot$tratam,
                   trace.factor = datos.wita.4h.conrot$Experimento, 
                   response = datos.wita.4h.conrot$S, fun = mean)
}

#--------------------wita 24 horas con rot--------------------------#
#---------------------------S-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 24h S",cex=2,font=2)
{
  mod.wita.24h.conrot.mfi.s<-lm(S~tratam+Experimento,data=datos.wita.24h.conrot)
  #verificar la interacción
  mod.wita.24h.conrot.mfi2.s<-lm(S~tratam*Experimento,data=datos.wita.24h.conrot)
  anova(mod.wita.24h.conrot.mfi.s,mod.wita.24h.conrot.mfi2.s)
  #gráficos de supuestos
  plot(mod.wita.24h.conrot.mfi.s)
  #se rechaza Ho cuando es heterocedastico
  bptest(S~tratam+Experimento,data=datos.wita.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.24h.conrot.mfi.s$residuals, "pnorm", 0, sd(mod.wita.24h.conrot.mfi.s$residuals))
  #verificar resultado
  anova(mod.wita.24h.conrot.mfi.s)
  tky.wita.24h.conrot.mfi2.s<-glht(mod.wita.24h.conrot.mfi.s, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.24h.conrot.mfi2.s)
  interaction.plot(x.factor = datos.wita.24h.conrot$tratam,
                   trace.factor = datos.wita.24h.conrot$Experimento, 
                   response = datos.wita.24h.conrot$S, fun = mean)
}

#--------------------wita 4 horas sin rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 4h G2",cex=2,font=2)
{
  mod.wita.4h.sinrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.wita.4h.sinrot)
  #verificar la interacción
  mod.wita.4h.sinrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.wita.4h.sinrot)
  anova(mod.wita.4h.sinrot.mfi.g2,mod.wita.4h.sinrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.wita.4h.sinrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.wita.4h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.4h.sinrot.mfi.g2$residuals, "pnorm", 0, sd(mod.wita.4h.sinrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.wita.4h.sinrot.mfi.g2)
  tky.wita.4h.sinrot.mfi.g2<-glht(mod.wita.4h.sinrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.4h.sinrot.mfi.g2)
  interaction.plot(x.factor = datos.wita.4h.sinrot$tratam,
                   trace.factor = datos.wita.4h.sinrot$Experimento, 
                   response = datos.wita.4h.sinrot$G2, fun = mean)
}

#--------------------wita 24 horas sin rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 24h G2",cex=2,font=2)
{
  mod.wita.24h.sinrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.wita.24h.sinrot)
  #verificar la interacción
  mod.wita.24h.sinrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.wita.24h.sinrot)
  anova(mod.wita.24h.sinrot.mfi.g2,mod.wita.24h.sinrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.wita.24h.sinrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.wita.24h.sinrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.24h.sinrot.mfi.g2$residuals, "pnorm", 0, sd(mod.wita.24h.sinrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.wita.24h.sinrot.mfi.g2)
  tky.wita.24h.sinrot.mfi.g2<-glht(mod.wita.24h.sinrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.24h.sinrot.mfi.g2)
  interaction.plot(x.factor = datos.wita.24h.sinrot$tratam,
                   trace.factor = datos.wita.24h.sinrot$Experimento, 
                   response = datos.wita.24h.sinrot$G2, fun = mean)
}

#--------------------wita 4 horas con rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAr 4h G2",cex=2,font=2)
{
  mod.wita.4h.conrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.wita.4h.conrot)
  #verificar la interacción
  mod.wita.4h.conrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.wita.4h.conrot)
  anova(mod.wita.4h.conrot.mfi.g2,mod.wita.4h.conrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.wita.4h.conrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.wita.4h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.4h.conrot.mfi.g2$residuals, "pnorm", 0, sd(mod.wita.4h.conrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.wita.4h.conrot.mfi.g2)
  tky.wita.4h.conrot.mfi.g2<-glht(mod.wita.4h.conrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.4h.conrot.mfi.g2)
  interaction.plot(x.factor = datos.wita.4h.conrot$tratam,
                   trace.factor = datos.wita.4h.conrot$Experimento, 
                   response = datos.wita.4h.conrot$G2, fun = mean)
}

#--------------------wita 24 horas con rot--------------------------#
#---------------------------G2-------------------------------------#
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 24h G2",cex=2,font=2)
{
  mod.wita.24h.conrot.mfi.g2<-lm(G2~tratam+Experimento,data=datos.wita.24h.conrot)
  #verificar la interacción
  mod.wita.24h.conrot.mfi2.g2<-lm(G2~tratam*Experimento,data=datos.wita.24h.conrot)
  anova(mod.wita.24h.conrot.mfi.g2,mod.wita.24h.conrot.mfi2.g2)
  #gráficos de supuestos
  plot(mod.wita.24h.conrot.mfi.g2)
  #se rechaza Ho cuando es heterocedastico
  bptest(G2~tratam+Experimento,data=datos.wita.24h.conrot)
  #se acepta Ho cuando es normal
  ks.test(mod.wita.24h.conrot.mfi.g2$residuals, "pnorm", 0, sd(mod.wita.24h.conrot.mfi.g2$residuals))
  #verificar resultado
  anova(mod.wita.24h.conrot.mfi.g2)
  tky.wita.24h.conrot.mfi2.g2<-glht(mod.wita.24h.conrot.mfi.g2, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
  summary(tky.wita.24h.conrot.mfi2.g2)
  interaction.plot(x.factor = datos.wita.24h.conrot$tratam,
                   trace.factor = datos.wita.24h.conrot$Experimento, 
                   response = datos.wita.24h.conrot$G2, fun = mean)
}



