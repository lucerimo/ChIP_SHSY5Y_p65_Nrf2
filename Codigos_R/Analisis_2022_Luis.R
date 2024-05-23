########################################################################################################
######    Analisis de resultados del experimento de produccion de ROS en celulas SH-SY5Y ###############
######    Respuesta: oxidacion de dihidrorhodamina 123 (verde)                           ###############
######    Variables: Experimento (Bloque)                                                ###############
######        Sustancia (14.2.2;Negativo;Positivo;Resveratrol;withaferina A)             ###############
######        Incubacion (4h, 24h)                                                       ###############
######        Concentracion de la sustancia (depende de la sustancia)                    ###############
######        Rotenona (Si, No)                                                          ###############
########################################################################################################
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
library(MASS)
library(gridExtra)

#1. Cargar funciones

#median <- function(x) {
#  data.frame(y = mean(x))
#}  
  
range <- function(x) {
    data.frame(ymin=min(x),
               ymax=max(x))
  }

graf.barras<-function(v.respuesta,l.incub,compuesto,etiq.y,sd.respuesta,n.vector)
{
barplot(v.respuesta,width=1,main=paste(l.incub,"h", paste(collapse="")),
         xlab=paste("(µM)",collapse=""), ylab=etiq.y,xlim=c(0,7),
        ylim=c(0,max(v.respuesta+(sd.respuesta/sqrt(n.vector)))),
        col = c("#f2f2f2ff","#4d4d4dff", "#999999ff", "#aaccffff", "#5599ffff", "#0055d4ff"))
arrows(x0=c(0.7,1.9,3.1,4.3,5.5,6.7),y0=v.respuesta-(sd.respuesta/sqrt(n.vector)),x1=c(0.7,1.9,3.1,4.3,5.5,6.7),y1=v.respuesta+(sd.respuesta/sqrt(n.vector)),
       length= 0.1,angle= 90, code= 3)

}
#2. Cargar datos
datos_ros<-read.csv("Resultados_2022_ROS.csv", sep = ";")
is.factor(datos_ros$tratam)
is.factor(datos_ros$Experimento)
levels(as.factor(datos_ros$tratam))
levels(as.factor(datos_ros$Experimento))
datos_ros$tratam<-relevel(as.factor(datos_ros$tratam),ref="CP")
datos_ros$Experimento<-as.factor(datos_ros$Experimento)
datos_ros
attach(datos_ros)

#3. Resumen de los datos
  #MFI
medias.mfi<-tapply(Media.geometrica,list(Sustancia,Cn,Incubacion,Rotenona),FUN=mean)
sd.mfi<-tapply(Media.geometrica,list(Sustancia,Cn,Incubacion,Rotenona),FUN=sd)
  #viabilidad
viab.mfi<-tapply(Viabilidad,list(Sustancia,Cn,Incubacion,Rotenona),FUN=mean)
Viabilidadsd.mfi<-tapply(Viabilidad,list(Sustancia,Cn,Incubacion,Rotenona),FUN=sd)
  #conteo
n.mfi<-tapply(Media.geometrica,list(Sustancia,Cn,Incubacion,Rotenona),FUN=length)
#4. Construir tablas con datos por comupesto
  #MFI
  #14.2.2
{
m.14_4h_rot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[1,2:4,1,2])
m.14_24h_rot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[1,2:4,2,2])
m.14_4h_srot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[1,2:4,1,1])
m.14_24h_srot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[1,2:4,2,1])
m.14<-as.matrix(rbind(m.14_4h_rot,m.14_24h_rot, m.14_4h_srot,m.14_24h_srot));colnames(m.14)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
sd.14_4h_rot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[1,2:4,1,2])
sd.14_24h_rot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[1,2:4,2,2])
sd.14_4h_srot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[1,2:4,1,1])
sd.14_24h_srot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[1,2:4,2,1])
sd.14<-as.matrix(rbind(sd.14_4h_rot,sd.14_24h_rot, sd.14_4h_srot,sd.14_24h_srot));colnames(sd.14)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
n.14<-as.matrix(c(n.mfi[2,1,2,1],n.mfi[3,1,2,2],n.mfi[4,1,2,1],n.mfi[1,2:4,2,2]))
}

  #Resveratrol
{
m.r_4h_rot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[5,5:7,1,2])
m.r_24h_rot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[5,5:7,2,2])
m.r_4h_srot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[5,5:7,1,1])
m.r_24h_srot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[5,5:7,2,1])
m.r<-as.matrix(rbind(m.r_4h_rot,m.r_24h_rot, m.r_4h_srot,m.r_24h_srot));colnames(m.r)<-c("Untreated","Rotenone","DMSO","10","20","40")
sd.r_4h_rot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[5,5:7,1,2])
sd.r_24h_rot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[5,5:7,2,2])
sd.r_4h_srot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[5,5:7,1,1])
sd.r_24h_srot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[5,5:7,2,1])
sd.r<-as.matrix(rbind(sd.r_4h_rot,sd.r_24h_rot, sd.r_4h_srot,sd.r_24h_srot));colnames(sd.r)<-c("Untreated","Rotenone","DMSO","10","20","40")
n.r<-as.matrix(c(n.mfi[2,1,2,1],n.mfi[3,1,2,2],n.mfi[4,1,2,1],n.mfi[5,5:7,2,2]))
}

  #Withaferin A
{
m.WA_4h_rot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[6,2:4,1,2])
m.WA_24h_rot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[6,2:4,2,2])
m.WA_4h_srot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[6,2:4,1,1])
m.WA_24h_srot<-c(medias.mfi[2,1,1,1],medias.mfi[3,1,1,2],medias.mfi[4,1,1,1],medias.mfi[6,2:4,2,1])
m.WA<-as.matrix(rbind(m.WA_4h_rot,m.WA_24h_rot, m.WA_4h_srot,m.WA_24h_srot));colnames(m.WA)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
sd.WA_4h_rot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[6,2:4,1,2])
sd.WA_24h_rot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[6,2:4,2,2])
sd.WA_4h_srot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[6,2:4,1,1])
sd.WA_24h_srot<-c(sd.mfi[2,1,1,1],sd.mfi[3,1,1,2],sd.mfi[4,1,1,1],sd.mfi[6,2:4,2,1])
sd.WA<-as.matrix(rbind(sd.WA_4h_rot,sd.WA_24h_rot, sd.WA_4h_srot,sd.WA_24h_srot));colnames(sd.WA)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
n.WA<-as.matrix(c(n.mfi[2,1,2,1],n.mfi[3,1,2,2],n.mfi[4,1,2,1],n.mfi[6,2:4,2,2]))
}
  #VIABILIADAD
#14.2.2
{
m.14_4h_rot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[1,2:4,1,2])
m.14_24h_rot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[1,2:4,2,2])
m.14_4h_srot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[1,2:4,1,1])
m.14_24h_srot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[1,2:4,2,1])
m.14.vib<-as.matrix(rbind(m.14_4h_rot.vib,m.14_24h_rot.vib, m.14_4h_srot.vib,m.14_24h_srot.vib));colnames(m.14.vib)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
sd.14_4h_rot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[1,2:4,1,2])
sd.14_24h_rot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[1,2:4,2,2])
sd.14_4h_srot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[1,2:4,1,1])
sd.14_24h_srot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[1,2:4,2,1])
sd.14.vib<-as.matrix(rbind(sd.14_4h_rot.vib,sd.14_24h_rot.vib, sd.14_4h_srot.vib,sd.14_24h_srot.vib));colnames(sd.14.vib)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
n.14.vib<-as.matrix(c(viab.mfi[2,1,2,1],viab.mfi[3,1,2,2],viab.mfi[4,1,2,1],viab.mfi[1,2:4,2,2]))
}

  #resveratrol
{
m.r_4h_rot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[5,5:7,1,2])
m.r_24h_rot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[5,5:7,2,2])
m.r_4h_srot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[5,5:7,1,1])
m.r_24h_srot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[5,5:7,2,1])
m.r.vib<-as.matrix(rbind(m.r_4h_rot.vib,m.r_24h_rot.vib, m.r_4h_srot.vib,m.r_24h_srot.vib));colnames(m.r.vib)<-c("Untreated","Rotenone","DMSO","10","20","40")
sd.r_4h_rot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[5,5:7,1,2])
sd.r_24h_rot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[5,5:7,2,2])
sd.r_4h_srot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[5,5:7,1,1])
sd.r_24h_srot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[5,5:7,2,1])
sd.r.vib<-as.matrix(rbind(sd.r_4h_rot.vib,sd.r_24h_rot.vib, sd.r_4h_srot.vib,sd.r_24h_srot.vib));colnames(sd.r.vib)<-c("Untreated","Rotenone","DMSO","10","20","40")
n.r.vib<-as.matrix(c(viab.mfi[2,1,2,1],viab.mfi[3,1,2,2],viab.mfi[4,1,2,1],viab.mfi[5,5:7,2,2]))
}

#Withaferin A
{
m.WA_4h_rot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[6,2:4,1,2])
m.WA_24h_rot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[6,2:4,2,2])
m.WA_4h_srot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[6,2:4,1,1])
m.WA_24h_srot.vib<-c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[6,2:4,2,1])
m.WA.vib<-as.matrix(rbind(m.WA_4h_rot.vib,m.WA_24h_rot.vib, m.WA_4h_srot.vib,m.WA_24h_srot.vib));colnames(m.WA.vib)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
sd.WA_4h_rot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[6,2:4,1,2])
sd.WA_24h_rot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[6,2:4,2,2])
sd.WA_4h_srot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[6,2:4,1,1])
sd.WA_24h_srot.vib<-c(Viabilidadsd.mfi[2,1,1,1],Viabilidadsd.mfi[3,1,1,2],Viabilidadsd.mfi[4,1,1,1],Viabilidadsd.mfi[6,2:4,2,1])
sd.WA.vib<-as.matrix(rbind(sd.WA_4h_rot.vib,sd.WA_24h_rot.vib, sd.WA_4h_srot.vib,sd.WA_24h_srot.vib));colnames(sd.WA.vib)<-c("Untreated","Rotenone","DMSO","0.5","1","2")
n.WA.vib<-as.matrix(c(viab.mfi[2,1,1,1],viab.mfi[3,1,1,2],viab.mfi[4,1,1,1],viab.mfi[6,2:4,2,2]))
}

#5. Graficos de barras
  #MFI
#par(mfrow = c(2, 2))
par(mar=c(4.5,4,1,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1.2,3,1.2,3))
plot.new()
text(0.5,0.5, "SWitA + Rotenona (4h)",cex=1.5,font=2)
graf.barras(v.respuesta=m.14[1,],l.incub="4",compuesto="14.2.2",etiq.y="MFI",sd.respuesta=sd.14[1,],n.vector=n.14)
plot.new()
text(0.5,0.5,"SWitA",cex=1.5,font=2)
graf.barras(v.respuesta=m.14[3,],l.incub="4",compuesto="14.2.2",etiq.y="MFI",sd.respuesta=sd.14[3,],n.vector=n.14)
graf.barras(v.respuesta=m.14[2,],l.incub="24",compuesto="14.2.2",etiq.y="MFI",sd.respuesta=sd.14[2,],n.vector=n.14)
graf.barras(v.respuesta=m.14[4,],l.incub="24",compuesto="14.2.2",etiq.y="MFI",sd.respuesta=sd.14[4,],n.vector=n.14)

par(mar=c(4.5,4,1,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1.2,3,1.2,3))
plot.new()
text(0.5,0.5, "Resveratrol + Rotenona (4h)",cex=1.5,font=2)
graf.barras(v.respuesta=m.r[1,],l.incub="4",compuesto="Resveratrol",etiq.y="MFI",sd.respuesta=sd.r[1,],n.vector=n.r)
plot.new()
text(0.5,0.5,"Resveratrol",cex=1.5,font=2)
graf.barras(v.respuesta=m.r[3,],l.incub="4",compuesto="Resveratrol",etiq.y="MFI",sd.respuesta=sd.r[3,],n.vector=n.r)
graf.barras(v.respuesta=m.r[2,],l.incub="24",compuesto="Resveratrol",etiq.y="MFI",sd.respuesta=sd.r[2,],n.vector=n.r)
graf.barras(v.respuesta=m.r[4,],l.incub="24",compuesto="Resveratrol",etiq.y="MFI",sd.respuesta=sd.r[4,],n.vector=n.r)


par(mar=c(4.5,4,1,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1.2,3,1.2,3))
plot.new()
text(0.5,0.5, "WitA + Rotenona (4h)",cex=1.5,font=2)
graf.barras(v.respuesta=m.WA[1,],l.incub="4",compuesto="Withaferin A",etiq.y="MFI",sd.respuesta=sd.WA[1,],n.vector=n.WA)
plot.new()
text(0.5,0.5,"WitA",cex=1.5,font=2)
graf.barras(v.respuesta=m.WA[3,],l.incub="4",compuesto="Withaferin A",etiq.y="MFI",sd.respuesta=sd.WA[3,],n.vector=n.WA)
graf.barras(v.respuesta=m.WA[2,],l.incub="24",compuesto="Withaferin A",etiq.y="MFI",sd.respuesta=sd.WA[2,],n.vector=n.WA)
graf.barras(v.respuesta=m.WA[4,],l.incub="24",compuesto="Withaferin A",etiq.y="MFI",sd.respuesta=sd.WA[4,],n.vector=n.WA)
 
#%Viabilidad
par(mar=c(4.5,4.5,1,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1.2,3,1.2,3))
plot.new()
text(0.5,0.5, "SWitA + Rotenona (4h)",cex=1.5,font=2)
graf.barras(v.respuesta=m.14.vib[1,],l.incub="4",compuesto="14.2.2",etiq.y="Viabilidad (%)",sd.respuesta=sd.14.vib[1,],n.vector=n.14.vib)
plot.new()
text(0.5,0.5,"SWitA",cex=1.5,font=2)
graf.barras(v.respuesta=m.14.vib[3,],l.incub="4",compuesto="14.2.2",etiq.y="Viabilidad (%)",sd.respuesta=sd.14.vib[3,],n.vector=n.14.vib)
graf.barras(v.respuesta=m.14.vib[2,],l.incub="24",compuesto="14.2.2",etiq.y="Viabilidad (%)",sd.respuesta=sd.14.vib[2,],n.vector=n.14.vib)
graf.barras(v.respuesta=m.14.vib[4,],l.incub="24",compuesto="14.2.2",etiq.y="Viabilidad (%)",sd.respuesta=sd.14.vib[4,],n.vector=n.14.vib)


par(mar=c(4.5,4.5,1.75,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,3,1,3))
plot.new()
text(0.5,0.5, "Resveratrol + Rotenona (4h)",cex=1.5,font=2)
graf.barras(v.respuesta=m.r.vib[1,],l.incub="4",compuesto="Resveratrol",etiq.y="Viabilidad (%)",sd.respuesta=sd.r.vib[1,],n.vector=n.r.vib)
plot.new()
text(0.5,0.5,"Resveratrol",cex=1.5,font=2)
graf.barras(v.respuesta=m.r.vib[3,],l.incub="4",compuesto="Resveratrol",etiq.y="Viabilidad (%)",sd.respuesta=sd.r.vib[3,],n.vector=n.r.vib)
graf.barras(v.respuesta=m.r.vib[2,],l.incub="24",compuesto="Resveratrol",etiq.y="Viabilidad (%)",sd.respuesta=sd.r.vib[2,],n.vector=n.r.vib)
graf.barras(v.respuesta=m.r.vib[4,],l.incub="24",compuesto="Resveratrol",etiq.y="Viabilidad (%)",sd.respuesta=sd.r.vib[4,],n.vector=n.r)

par(mar=c(4.5,4.5,1.75,1.75))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,3,1,3))
plot.new()
text(0.5,0.5, "WitA + Rotenona (4h)",cex=1.5,font=2)
graf.barras(v.respuesta=m.WA.vib[1,],l.incub="4",compuesto="Withaferin A",etiq.y="Viabilidad (%)",sd.respuesta=sd.WA.vib[1,],n.vector=n.WA.vib)
plot.new()
text(0.5,0.5,"WitA",cex=1.5,font=2)
graf.barras(v.respuesta=m.WA.vib[3,],l.incub="4",compuesto="Withaferin A",etiq.y="Viabilidad (%)",sd.respuesta=sd.WA.vib[3,],n.vector=n.WA.vib)
graf.barras(v.respuesta=m.WA.vib[2,],l.incub="24",compuesto="Withaferin A",etiq.y="Viabilidad (%)",sd.respuesta=sd.WA.vib[2,],n.vector=n.WA.vib)
graf.barras(v.respuesta=m.WA.vib[4,],l.incub="24",compuesto="Withaferin A",etiq.y="Viabilidad (%)",sd.respuesta=sd.WA.vib[4,],n.vector=n.WA.vib)

#6. Analisis estadistico
#--------------------resv 4 horas sin rot--------------------------#
#--------------------------MFI-------------------------------------#
datos.resv.4h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="R1"|tratam=="R2"|tratam=="R3")&Incubacion=="4",]
mod.resv.4h.sinrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.resv.4h.sinrot)
#verificar la interacción
mod.resv.4h.sinrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.resv.4h.sinrot)
anova(mod.resv.4h.sinrot.mfi2,mod.resv.4h.sinrot.mfi) #Se escoje el modelo completo
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.resv.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.resv.4h.sinrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.resv.4h.sinrot.mfi))
boxcox(mod.resv.4h.sinrot.mfi)$x[which.max(boxcox(mod.resv.4h.sinrot.mfi)$y)]
mod.resv.4h.sinrot.mfit<-lm(Media.geometrica^-0.46~tratam+Experimento,data=datos.resv.4h.sinrot)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica^-0.46~tratam+Experimento,data=datos.resv.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.sinrot.mfit$residuals, "pnorm", 0, sd(mod.resv.4h.sinrot.mfit$residuals))
#verificar resultado
anova(mod.resv.4h.sinrot.mfit)
res.tky.4h.sinrot.mfi<-glht(mod.resv.4h.sinrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.4h.sinrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 4h",cex=2,font=2)
plot(mod.resv.4h.sinrot.mfit)
interaction.plot(x.factor = datos.resv.4h.sinrot$tratam,
                 trace.factor = datos.resv.4h.sinrot$Experimento, 
                 response = datos.resv.4h.sinrot$Media.geometrica, fun = mean)


#--------------------resv 24 horas sin rot--------------------------#
#--------------------------MFI-------------------------------------#
datos.resv.24h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="R1"|tratam=="R2"|tratam=="R3")&Incubacion=="24",]
mod.resv.24h.sinrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.resv.24h.sinrot)
#verificar la interacción
mod.resv.24h.sinrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.resv.24h.sinrot)
anova(mod.resv.24h.sinrot.mfi2,mod.resv.24h.sinrot.mfi)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.resv.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.24h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.resv.24h.sinrot.mfi$residuals))
#verificar resultado
anova(mod.resv.24h.sinrot.mfi)
res.tky.24h.sinrot.mfi<-glht(mod.resv.24h.sinrot.mfi, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.24h.sinrot.mfi)
#gráficos de supuestos
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv 24h",cex=2,font=2)
plot(mod.resv.24h.sinrot.mfi)
interaction.plot(x.factor = datos.resv.24h.sinrot$tratam,
                 trace.factor = datos.resv.24h.sinrot$Experimento, 
                 response = datos.resv.24h.sinrot$Media.geometrica, fun = mean)


#--------------------resv 4 horas con rot--------------------------#
#--------------------------MFI-------------------------------------#
datos.resv.4h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="RR1"|tratam=="RR2"|tratam=="RR3")&Incubacion=="4",]
mod.resv.4h.conrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.resv.4h.conrot)
#verificar la interacción
mod.resv.4h.conrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.resv.4h.conrot)
anova(mod.resv.4h.conrot.mfi,mod.resv.4h.conrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.resv.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.conrot.mfi$residuals, "pnorm", 0, sd(mod.resv.4h.conrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.resv.4h.conrot.mfi))
boxcox(mod.resv.4h.conrot.mfi)$x[which.max(boxcox(mod.resv.4h.conrot.mfi)$y)]
mod.resv.4h.conrot.mfit<-lm(Media.geometrica^-0.18~tratam+Experimento,data=datos.resv.4h.conrot)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica^-0.18~tratam+Experimento,data=datos.resv.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.conrot.mfit$residuals, "pnorm", 0, sd(mod.resv.4h.conrot.mfit$residuals))
#verificar resultado
anova(mod.resv.4h.conrot.mfit)
res.tky.4h.conrot.mfi<-glht(mod.resv.4h.conrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.4h.conrot.mfi)
#gráficos de supuestos
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 4h",cex=2,font=2)
plot(mod.resv.4h.conrot.mfit)
interaction.plot(x.factor = datos.resv.4h.conrot$tratam,
                 trace.factor = datos.resv.4h.conrot$Experimento, 
                 response = datos.resv.4h.conrot$Media.geometrica, fun = mean)


#--------------------resv 24 horas con rot--------------------------#
#--------------------------MFI-------------------------------------#
datos.resv.24h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="RR1"|tratam=="RR2"|tratam=="RR3")&Incubacion=="24",]
mod.resv.24h.conrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.resv.24h.conrot)
#verificar la interacción
mod.resv.24h.conrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.resv.24h.conrot)
anova(mod.resv.24h.conrot.mfi,mod.resv.24h.conrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.resv.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.24h.conrot.mfi$residuals, "pnorm", 0, sd(mod.resv.24h.conrot.mfi$residuals))
#verificar resultado
anova(mod.resv.24h.conrot.mfi)
res.tky.24h.conrot.mfi<-glht(mod.resv.24h.conrot.mfi, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.24h.conrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 24h",cex=2,font=2)
plot(mod.resv.24h.conrot.mfi)
interaction.plot(x.factor = datos.resv.24h.conrot$tratam,
                 trace.factor = datos.resv.24h.conrot$Experimento, 
                 response = datos.resv.24h.conrot$Media.geometrica, fun = mean)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#--------------------14.2.2 4 horas sin rot--------------------------#
#---------------------------MFI-------------------------------------#
datos.14.2.2.4h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="4",]
mod.14.2.2.4h.sinrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
#verificar la interacción
mod.14.2.2.4h.sinrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.14.2.2.4h.sinrot)
anova(mod.14.2.2.4h.sinrot.mfi,mod.14.2.2.4h.sinrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.4h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.14.2.2.4h.sinrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.14.2.2.4h.sinrot.mfi))
boxcox(mod.14.2.2.4h.sinrot.mfi)$x[which.max(boxcox(mod.14.2.2.4h.sinrot.mfi)$y)]
mod.14.2.2.4h.sinrot.mfit<-lm(Media.geometrica^-0.74~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica^-0.74~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.4h.sinrot.mfit$residuals, "pnorm", 0, sd(mod.14.2.2.4h.sinrot.mfit$residuals))
#verificar resultado
anova(mod.14.2.2.4h.sinrot.mfit)
tky.14.2.2.4h.sinrot.mfi<-glht(mod.14.2.2.4h.sinrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.4h.sinrot.mfi)
#gráficos de supuestos
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 4h",cex=2,font=2)
plot(mod.14.2.2.4h.sinrot.mfit)
interaction.plot(x.factor = datos.14.2.2.4h.sinrot$tratam,
                 trace.factor = datos.14.2.2.4h.sinrot$Experimento, 
                 response = datos.14.2.2.4h.sinrot$Media.geometrica, fun = mean)


#--------------------14.2.2 24 horas sin rot--------------------------#
#---------------------------MFI-------------------------------------#
#14.2.2 4 horas sin rot
datos.14.2.2.24h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="24",]
mod.14.2.2.24h.sinrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
#verificar la interacción
mod.14.2.2.24h.sinrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.14.2.2.24h.sinrot)
anova(mod.14.2.2.24h.sinrot.mfi,mod.14.2.2.24h.sinrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.24h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.14.2.2.24h.sinrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.14.2.2.24h.sinrot.mfi))
boxcox(mod.14.2.2.24h.sinrot.mfi)$x[which.max(boxcox(mod.14.2.2.24h.sinrot.mfi)$y)]
mod.14.2.2.24h.sinrot.mfit<-lm(Media.geometrica^0.42~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica^0.42~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.24h.sinrot.mfit$residuals, "pnorm", 0, sd(mod.14.2.2.24h.sinrot.mfit$residuals))
#verificar resultado
anova(mod.14.2.2.24h.sinrot.mfit)
tky.14.2.2.24h.sinrot.mfi<-glht(mod.14.2.2.24h.sinrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.24h.sinrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA 24h",cex=2,font=2)
plot(mod.14.2.2.24h.sinrot.mfit)
interaction.plot(x.factor = datos.14.2.2.24h.sinrot$tratam,
                 trace.factor = datos.14.2.2.24h.sinrot$Experimento, 
                 response = datos.14.2.2.24h.sinrot$Media.geometrica, fun = mean)


#--------------------14.2.2 4 horas con rot--------------------------#
#---------------------------MFI-------------------------------------#
#14.2.2 4 horas sin rot
datos.14.2.2.4h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="4",]
mod.14.2.2.4h.conrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.14.2.2.4h.conrot)
#verificar la interacción
mod.14.2.2.4h.conrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.14.2.2.4h.conrot)
anova(mod.14.2.2.4h.conrot.mfi,mod.14.2.2.4h.conrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.14.2.2.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.4h.conrot.mfi$residuals, "pnorm", 0, sd(mod.14.2.2.4h.conrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.14.2.2.4h.conrot.mfi))
boxcox(mod.14.2.2.4h.conrot.mfi)$x[which.max(boxcox(mod.14.2.2.4h.conrot.mfi)$y)]
mod.14.2.2.4h.conrot.mfit<-lm(Media.geometrica^-0.06~tratam+Experimento,data=datos.14.2.2.4h.conrot)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica^-0.06~tratam+Experimento,data=datos.14.2.2.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.4h.conrot.mfit$residuals, "pnorm", 0, sd(mod.14.2.2.4h.conrot.mfit$residuals))
#verificar resultado
anova(mod.14.2.2.4h.conrot.mfit)
tky.14.2.2.4h.conrot.mfi<-glht(mod.14.2.2.4h.conrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.4h.conrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 4h",cex=2,font=2)
plot(mod.14.2.2.24h.sinrot.mfit)
interaction.plot(x.factor = datos.14.2.2.4h.conrot$tratam,
                 trace.factor = datos.14.2.2.4h.conrot$Experimento, 
                 response = datos.14.2.2.4h.conrot$Media.geometrica, fun = mean)


#--------------------14.2.2 24 horas con rot--------------------------#
#---------------------------MFI-------------------------------------#
datos.14.2.2.24h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="24",]
mod.14.2.2.24h.conrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.14.2.2.24h.conrot)
#verificar la interacción
mod.14.2.2.24h.conrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.14.2.2.24h.conrot)
anova(mod.14.2.2.24h.conrot.mfi,mod.14.2.2.24h.conrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.14.2.2.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.24h.conrot.mfi$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.mfi$residuals))
#verificar resultado
anova(mod.14.2.2.24h.conrot.mfi)
tky.14.2.2.24h.conrot.mfi<-glht(mod.14.2.2.24h.conrot.mfi, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.24h.conrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR 24h",cex=2,font=2)
plot(mod.14.2.2.24h.conrot.mfi)
interaction.plot(x.factor = datos.14.2.2.24h.conrot$tratam,
                 trace.factor = datos.14.2.2.24h.conrot$Experimento, 
                 response = datos.14.2.2.24h.conrot$Media.geometrica, fun = mean)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------WitA 4 horas sin rot--------------------------#
#---------------------------MFI-------------------------------------#
datos.wita.4h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WA1"|tratam=="WA2"|tratam=="WA3")&Incubacion=="4",]
mod.wita.4h.sinrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.wita.4h.sinrot)
#verificar la interacción
mod.wita.4h.sinrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.wita.4h.sinrot)
anova(mod.wita.4h.sinrot.mfi,mod.wita.4h.sinrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.wita.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.wita.4h.sinrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.wita.4h.sinrot.mfi))
boxcox(mod.wita.4h.sinrot.mfi)$x[which.max(boxcox(mod.wita.4h.sinrot.mfi)$y)]
mod.wita.4h.sinrot.mfit<-lm(Media.geometrica^-0.58~tratam+Experimento,data=datos.wita.4h.sinrot)
#se rechaza Ho cuando es heterocedastico
bptest(log(Media.geometrica) ~ tratam + Experimento, data = datos.wita.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.sinrot.mfit$residuals, "pnorm", 0, sd(mod.wita.4h.sinrot.mfit$residuals))
#verificar resultado
anova(mod.wita.4h.sinrot.mfit)
tky.wita.4h.sinrot.mfi<-glht(mod.wita.4h.sinrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.4h.sinrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 4h",cex=2,font=2)
plot(mod.wita.4h.sinrot.mfit)
interaction.plot(x.factor = datos.wita.4h.sinrot$tratam,
                 trace.factor = datos.wita.4h.sinrot$Experimento, 
                 response = datos.wita.4h.sinrot$Media.geometrica, fun = mean)


#--------------------WitA 24 horas sin rot--------------------------#
#---------------------------MFI-------------------------------------#
datos.wita.24h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WA1"|tratam=="WA2"|tratam=="WA3")&Incubacion=="24",]
mod.wita.24h.sinrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.wita.24h.sinrot)
#verificar la interacción
mod.wita.24h.sinrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.wita.24h.sinrot)
anova(mod.wita.24h.sinrot.mfi,mod.wita.24h.sinrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.wita.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.24h.sinrot.mfi$residuals, "pnorm", 0, sd(mod.wita.24h.sinrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.wita.24h.sinrot.mfi))
boxcox(mod.wita.24h.sinrot.mfi)$x[which.max(boxcox(mod.wita.24h.sinrot.mfi)$y)]
mod.wita.24h.sinrot.mfit<-lm(Media.geometrica^0.46~tratam+Experimento,data=datos.wita.24h.sinrot)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica^0.46~tratam+Experimento,data=datos.wita.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.24h.sinrot.mfit$residuals, "pnorm", 0, sd(mod.wita.24h.sinrot.mfit$residuals))
#verificar resultado
anova(mod.wita.24h.sinrot.mfit)
tky.wita.24h.sinrot.mfi<-glht(mod.wita.24h.sinrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.24h.sinrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA 24h",cex=2,font=2)
plot(mod.wita.24h.sinrot.mfit)
interaction.plot(x.factor = datos.wita.24h.sinrot$tratam,
                 trace.factor = datos.wita.24h.sinrot$Experimento, 
                 response = datos.wita.24h.sinrot$Media.geometrica, fun = mean)



#--------------------WitA 4 horas con rot--------------------------#
#---------------------------MFI-------------------------------------#
datos.wita.4h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WAR1"|tratam=="WAR2"|tratam=="WAR3")&Incubacion=="4",]
mod.wita.4h.conrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.wita.4h.conrot)
#verificar la interacción
mod.wita.4h.conrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.wita.4h.conrot)
anova(mod.wita.4h.conrot.mfi,mod.wita.4h.conrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.wita.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.conrot.mfi$residuals, "pnorm", 0, sd(mod.wita.4h.conrot.mfi$residuals))
#TRANSFORMACIÓN
plot(boxcox(mod.wita.4h.conrot.mfi))
boxcox(mod.wita.4h.conrot.mfi)$x[which.max(boxcox(mod.wita.4h.conrot.mfi)$y)]
mod.wita.4h.conrot.mfit<-lm(Media.geometrica^-0.14~tratam+Experimento,data=datos.wita.4h.conrot)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica^-0.14~tratam+Experimento,data=datos.wita.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.conrot.mfit$residuals, "pnorm", 0, sd(mod.wita.4h.conrot.mfit$residuals))
#verificar resultado
anova(mod.wita.4h.conrot.mfit)
tky.wita.4h.conrot.mfi<-glht(mod.wita.4h.conrot.mfit, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.4h.conrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR 4h",cex=2,font=2)
plot(mod.wita.4h.conrot.mfit)
interaction.plot(x.factor = datos.wita.4h.conrot$tratam,
                 trace.factor = datos.wita.4h.conrot$Experimento, 
                 response = datos.wita.4h.conrot$Media.geometrica, fun = mean)



#--------------------WitA 24 horas con rot--------------------------#
#---------------------------MFI-------------------------------------#
datos.wita.24h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WAR1"|tratam=="WAR2"|tratam=="WAR3")&Incubacion=="24",]
mod.wita.24h.conrot.mfi<-lm(Media.geometrica~tratam+Experimento,data=datos.wita.24h.conrot)
#verificar la interacción
mod.wita.24h.conrot.mfi2<-lm(Media.geometrica~tratam*Experimento,data=datos.wita.24h.conrot)
anova(mod.wita.24h.conrot.mfi,mod.wita.24h.conrot.mfi2)
#se rechaza Ho cuando es heterocedastico
bptest(Media.geometrica~tratam+Experimento,data=datos.wita.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.24h.conrot.mfi$residuals, "pnorm", 0, sd(mod.wita.24h.conrot.mfi$residuals))
#verificar resultado
anova(mod.wita.24h.conrot.mfi)
tky.wita.24h.conrot.mfi<-glht(mod.wita.24h.conrot.mfi, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.24h.conrot.mfi)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAr 24h",cex=2,font=2)
plot(mod.wita.24h.conrot.mfi)
interaction.plot(x.factor = datos.wita.24h.conrot$tratam,
                 trace.factor = datos.wita.24h.conrot$Experimento, 
                 response = datos.wita.24h.conrot$Media.geometrica, fun = mean)


################################################################################

#7. Analisis estadistico Incertidumbre
#--------------------resv 4 horas sin rot--------------------------#
#--------------------------Viabilidad-------------------------------------#
mod.resv.4h.sinrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.resv.4h.sinrot)
#verificar la interacción
mod.resv.4h.sinrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.resv.4h.sinrot)
anova(mod.resv.4h.sinrot.viab2,mod.resv.4h.sinrot.viab) #Se escoje el modelo completo
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.resv.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.sinrot.viab$residuals, "pnorm", 0, sd(mod.resv.4h.sinrot.viab$residuals))
#verificar resultado
anova(mod.resv.4h.sinrot.viab)
res.tky.4h.sinrot.viab<-glht(mod.resv.4h.sinrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.4h.sinrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv Viab 4h",cex=2,font=2)
plot(mod.resv.4h.sinrot.viab)
interaction.plot(x.factor = datos.resv.4h.sinrot$tratam,
                 trace.factor = datos.resv.4h.sinrot$Experimento, 
                 response = datos.resv.4h.sinrot$Viabilidad, fun = mean)


#--------------------resv 24 horas sin rot--------------------------#
#--------------------------viab-------------------------------------#
datos.resv.24h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="R1"|tratam=="R2"|tratam=="R3")&Incubacion=="24",]
mod.resv.24h.sinrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.resv.24h.sinrot)
#verificar la interacción
mod.resv.24h.sinrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.resv.24h.sinrot)
anova(mod.resv.24h.sinrot.viab2,mod.resv.24h.sinrot.viab)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.resv.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.24h.sinrot.viab$residuals, "pnorm", 0, sd(mod.resv.24h.sinrot.viab$residuals))
#verificar resultado
anova(mod.resv.24h.sinrot.viab)
res.tky.24h.sinrot.viab<-glht(mod.resv.24h.sinrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.24h.sinrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "Resv Viab 24h",cex=2,font=2)
plot(mod.resv.24h.sinrot.viab)
interaction.plot(x.factor = datos.resv.24h.sinrot$tratam,
                 trace.factor = datos.resv.24h.sinrot$Experimento, 
                 response = datos.resv.24h.sinrot$Viabilidad, fun = mean)


#--------------------resv 4 horas con rot--------------------------#
#--------------------------viab-------------------------------------#
datos.resv.4h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="RR1"|tratam=="RR2"|tratam=="RR3")&Incubacion=="4",]
mod.resv.4h.conrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.resv.4h.conrot)
#verificar la interacción
mod.resv.4h.conrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.resv.4h.conrot)
anova(mod.resv.4h.conrot.viab,mod.resv.4h.conrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.resv.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.4h.conrot.viab$residuals, "pnorm", 0, sd(mod.resv.4h.conrot.viab$residuals))
#verificar resultado
anova(mod.resv.4h.conrot.viab)
res.tky.4h.conrot.viab<-glht(mod.resv.4h.conrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.4h.conrot.viab)
#gráficos de supuestos
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR Viab 4h",cex=2,font=2)
plot(mod.resv.4h.conrot.viab)
interaction.plot(x.factor = datos.resv.4h.conrot$tratam,
                 trace.factor = datos.resv.4h.conrot$Experimento, 
                 response = datos.resv.4h.conrot$Viabilidad, fun = mean)


#--------------------resv 24 horas con rot--------------------------#
#--------------------------viab-------------------------------------#
datos.resv.24h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="RR1"|tratam=="RR2"|tratam=="RR3")&Incubacion=="24",]
mod.resv.24h.conrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.resv.24h.conrot)
#verificar la interacción
mod.resv.24h.conrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.resv.24h.conrot)
anova(mod.resv.24h.conrot.viab,mod.resv.24h.conrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.resv.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.resv.24h.conrot.viab$residuals, "pnorm", 0, sd(mod.resv.24h.conrot.viab$residuals))
#verificar resultado
anova(mod.resv.24h.conrot.viab)
res.tky.24h.conrot.viab<-glht(mod.resv.24h.conrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(res.tky.24h.conrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "ResvR 24h",cex=2,font=2)
plot(mod.resv.24h.conrot.viab)
interaction.plot(x.factor = datos.resv.24h.conrot$tratam,
                 trace.factor = datos.resv.24h.conrot$Experimento, 
                 response = datos.resv.24h.conrot$Viabilidad, fun = mean)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#--------------------14.2.2 4 horas sin rot--------------------------#
#---------------------------viab-------------------------------------#
datos.14.2.2.4h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="4",]
mod.14.2.2.4h.sinrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
#verificar la interacción
mod.14.2.2.4h.sinrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.14.2.2.4h.sinrot)
anova(mod.14.2.2.4h.sinrot.viab,mod.14.2.2.4h.sinrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.14.2.2.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.4h.sinrot.viab$residuals, "pnorm", 0, sd(mod.14.2.2.4h.sinrot.viab$residuals))
#verificar resultado
anova(mod.14.2.2.4h.sinrot.viab)
tky.14.2.2.4h.sinrot.viab<-glht(mod.14.2.2.4h.sinrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.4h.sinrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA Viab 4h",cex=2,font=2)
plot(mod.14.2.2.4h.sinrot.viab)
interaction.plot(x.factor = datos.14.2.2.4h.sinrot$tratam,
                 trace.factor = datos.14.2.2.4h.sinrot$Experimento, 
                 response = datos.14.2.2.4h.sinrot$Viabilidad, fun = mean)


#--------------------14.2.2 24 horas sin rot--------------------------#
#---------------------------viab-------------------------------------#
#14.2.2 4 horas sin rot
datos.14.2.2.24h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SW1"|tratam=="SW2"|tratam=="SW3")&Incubacion=="24",]
mod.14.2.2.24h.sinrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
#verificar la interacción
mod.14.2.2.24h.sinrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.14.2.2.24h.sinrot)
anova(mod.14.2.2.24h.sinrot.viab,mod.14.2.2.24h.sinrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.14.2.2.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.24h.sinrot.viab$residuals, "pnorm", 0, sd(mod.14.2.2.24h.sinrot.viab$residuals))
#verificar resultado
anova(mod.14.2.2.24h.sinrot.viab)
tky.14.2.2.24h.sinrot.viab<-glht(mod.14.2.2.24h.sinrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.24h.sinrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitA Viab 24h",cex=2,font=2)
plot(mod.14.2.2.24h.sinrot.viab)
interaction.plot(x.factor = datos.14.2.2.24h.sinrot$tratam,
                 trace.factor = datos.14.2.2.24h.sinrot$Experimento, 
                 response = datos.14.2.2.24h.sinrot$Viabilidad, fun = mean)


#--------------------14.2.2 4 horas con rot--------------------------#
#---------------------------viab-------------------------------------#
#14.2.2 4 horas sin rot
datos.14.2.2.4h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="4",]
mod.14.2.2.4h.conrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.14.2.2.4h.conrot)
#verificar la interacción
mod.14.2.2.4h.conrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.14.2.2.4h.conrot)
anova(mod.14.2.2.4h.conrot.viab,mod.14.2.2.4h.conrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.14.2.2.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.4h.conrot.viab$residuals, "pnorm", 0, sd(mod.14.2.2.4h.conrot.viab$residuals))
#verificar resultado
anova(mod.14.2.2.4h.conrot.viab)
tky.14.2.2.4h.conrot.viab<-glht(mod.14.2.2.4h.conrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.4h.conrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR Viab 4h",cex=2,font=2)
plot(mod.14.2.2.24h.sinrot.viab)
interaction.plot(x.factor = datos.14.2.2.4h.conrot$tratam,
                 trace.factor = datos.14.2.2.4h.conrot$Experimento, 
                 response = datos.14.2.2.4h.conrot$Viabilidad, fun = mean)


#--------------------14.2.2 24 horas con rot--------------------------#
#---------------------------viab-------------------------------------#
datos.14.2.2.24h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="SWR1"|tratam=="SWR2"|tratam=="SWR3")&Incubacion=="24",]
mod.14.2.2.24h.conrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.14.2.2.24h.conrot)
#verificar la interacción
mod.14.2.2.24h.conrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.14.2.2.24h.conrot)
anova(mod.14.2.2.24h.conrot.viab,mod.14.2.2.24h.conrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.14.2.2.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.14.2.2.24h.conrot.viab$residuals, "pnorm", 0, sd(mod.14.2.2.24h.conrot.viab$residuals))
#verificar resultado
anova(mod.14.2.2.24h.conrot.viab)
tky.14.2.2.24h.conrot.viab<-glht(mod.14.2.2.24h.conrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.14.2.2.24h.conrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "SWitAR Viab 24h",cex=2,font=2)
plot(mod.14.2.2.24h.conrot.viab)
interaction.plot(x.factor = datos.14.2.2.24h.conrot$tratam,
                 trace.factor = datos.14.2.2.24h.conrot$Experimento, 
                 response = datos.14.2.2.24h.conrot$Viabilidad, fun = mean)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#--------------------WitA 4 horas sin rot--------------------------#
#---------------------------viab-------------------------------------#
datos.wita.4h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WA1"|tratam=="WA2"|tratam=="WA3")&Incubacion=="4",]
mod.wita.4h.sinrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.wita.4h.sinrot)
#verificar la interacción
mod.wita.4h.sinrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.wita.4h.sinrot)
anova(mod.wita.4h.sinrot.viab,mod.wita.4h.sinrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.wita.4h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.sinrot.viab$residuals, "pnorm", 0, sd(mod.wita.4h.sinrot.viab$residuals))
#verificar resultado
anova(mod.wita.4h.sinrot.viab)
tky.wita.4h.sinrot.viab<-glht(mod.wita.4h.sinrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.4h.sinrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA Viab 4h",cex=2,font=2)
plot(mod.wita.4h.sinrot.viab)
interaction.plot(x.factor = datos.wita.4h.sinrot$tratam,
                 trace.factor = datos.wita.4h.sinrot$Experimento, 
                 response = datos.wita.4h.sinrot$Viabilidad, fun = mean)


#--------------------WitA 24 horas sin rot--------------------------#
#---------------------------viab-------------------------------------#
datos.wita.24h.sinrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WA1"|tratam=="WA2"|tratam=="WA3")&Incubacion=="24",]
mod.wita.24h.sinrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.wita.24h.sinrot)
#verificar la interacción
mod.wita.24h.sinrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.wita.24h.sinrot)
anova(mod.wita.24h.sinrot.viab,mod.wita.24h.sinrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.wita.24h.sinrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.24h.sinrot.viab$residuals, "pnorm", 0, sd(mod.wita.24h.sinrot.viab$residuals))
#verificar resultado
anova(mod.wita.24h.sinrot.viab)
tky.wita.24h.sinrot.viab<-glht(mod.wita.24h.sinrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.24h.sinrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitA Viab 24h",cex=2,font=2)
plot(mod.wita.24h.sinrot.viab)
interaction.plot(x.factor = datos.wita.24h.sinrot$tratam,
                 trace.factor = datos.wita.24h.sinrot$Experimento, 
                 response = datos.wita.24h.sinrot$Viabilidad, fun = mean)



#--------------------WitA 4 horas con rot--------------------------#
#---------------------------viab-------------------------------------#
datos.wita.4h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WAR1"|tratam=="WAR2"|tratam=="WAR3")&Incubacion=="4",]
mod.wita.4h.conrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.wita.4h.conrot)
#verificar la interacción
mod.wita.4h.conrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.wita.4h.conrot)
anova(mod.wita.4h.conrot.viab,mod.wita.4h.conrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.wita.4h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.4h.conrot.viab$residuals, "pnorm", 0, sd(mod.wita.4h.conrot.viab$residuals))
#verificar resultado
anova(mod.wita.4h.conrot.viab)
tky.wita.4h.conrot.viab<-glht(mod.wita.4h.conrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.4h.conrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR Viab 4h",cex=2,font=2)
plot(mod.wita.4h.conrot.viab)
interaction.plot(x.factor = datos.wita.4h.conrot$tratam,
                 trace.factor = datos.wita.4h.conrot$Experimento, 
                 response = datos.wita.4h.conrot$Viabilidad, fun = mean)



#--------------------WitA 24 horas con rot--------------------------#
#---------------------------viab-------------------------------------#
datos.wita.24h.conrot<-datos_ros[(tratam=="CN"|tratam=="CP"|tratam=="DMSO"|tratam=="WAR1"|tratam=="WAR2"|tratam=="WAR3")&Incubacion=="24",]
mod.wita.24h.conrot.viab<-lm(Viabilidad~tratam+Experimento,data=datos.wita.24h.conrot)
#verificar la interacción
mod.wita.24h.conrot.viab2<-lm(Viabilidad~tratam*Experimento,data=datos.wita.24h.conrot)
anova(mod.wita.24h.conrot.viab,mod.wita.24h.conrot.viab2)
#se rechaza Ho cuando es heterocedastico
bptest(Viabilidad~tratam+Experimento,data=datos.wita.24h.conrot)
#se acepta Ho cuando es normal
ks.test(mod.wita.24h.conrot.viab$residuals, "pnorm", 0, sd(mod.wita.24h.conrot.viab$residuals))
#verificar resultado
anova(mod.wita.24h.conrot.viab)
tky.wita.24h.conrot.viab<-glht(mod.wita.24h.conrot.viab, linfct = mcp(tratam = "Tukey"), alternative = "two.sided")
summary(tky.wita.24h.conrot.viab)
#gráficos de supuestos
par(mar=c(2.5,2.5,1.75,1.75))
layout(matrix(c(1,2,3,6,1,4,5,6),ncol=2),heights=c(1.5,3,3,3))
plot.new()
text(0.5,0.5, "WitAR Viab 24h",cex=2,font=2)
plot(mod.wita.24h.conrot.viab)
interaction.plot(x.factor = datos.wita.24h.conrot$tratam,
                 trace.factor = datos.wita.24h.conrot$Experimento, 
                 response = datos.wita.24h.conrot$Viabilidad, fun = mean)





