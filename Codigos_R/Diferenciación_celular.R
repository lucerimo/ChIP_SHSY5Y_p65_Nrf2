#############################################################################
###   Análisis de diferenciación #########
### Se miden 30 células de cada condición y se realiza una prueba z########
#############################################################################
library(car)

#mediciones celulares n=33
data <- list(
  size_und <- c(183.2,275.8,231.8,271.9,180.5,182.9,323.5,182.7,262.2,365.4,320.2,261.0,229.0,190.4,263.8,164.4,264.7,230.0,314.9,286.5,201.7,329.1,266.1,166.7,237.8,295.4,288.8,231.5,206.1,360.9,195.9,185.6,221.3),
  size_diff <-c(350.5,418.3,363.6,403.1,351.0,576.4,402.8,501.9,484.4,543.5,444.8,578.4,454.6,507.3,418.4,594.1,442.7,442.1,439.5,480.9,377.1,414.4,526.4,429.4,434.1,477.2,478.6,517.6,464.5,545.7,366.6,531.9,557.2)
)

boxplot(data, names = c("No diferenciada", "Diferenciada"), main="Longitud de las células por condición", 
        ylab="Longitud (px)", col = c("#4d4d4dff",  "#5599ffff"),  xlab="Condición celular")


shapiro.test(size_diff)
shapiro.test(size_und)
leveneTest(c(size_diff, size_und) ~ factor(rep(1:2, each= 33)))
t.test(size_diff, size_und)

