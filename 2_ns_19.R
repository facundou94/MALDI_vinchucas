################ MALDI-TOF ANALISIS VINCHUCA ###################################
################ 1) NO SUPERVISADO - 19 INDIVIDUOS  ############################
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################

library("readBrukerFlexData")
library("binda")
library("fs")
library("readxl")
library("MALDIquant")
library("MALDIquantForeign")
library("MALDIrppa")
library("tidyverse")
library("dplyr")
library("clValid")
library(cluster)
library(factoextra)
library(ggplot2)
library(gridExtra)

### CARGA DE ARCHIVOS ##########################################################
################################################################################

# Creación de la ruta relativa de los archivos
ruta_proyecto <- here("OneDrive/Documents/Proyectos/Vinchucas")
ruta_datos <- file.path(ruta_proyecto)

# Load the Rdata files using the relative path
load(file.path(ruta_datos, "matint_19_ind_dico.Rdata"))
load(file.path(ruta_datos, "matint_19_ind.Rdata"))


### SELECCIÓN DE PICOS #########################################################
################################################################################

# Selección de picos para binary discriminant analysis (BDA)
factor_tipo <- factor(df_metadata_prom_mue$estado)
is.binaryMatrix(matint_19_ind_dico) # TRUE
br <- binda.ranking(matint_19_ind_dico, factor_tipo, verbose = FALSE)

# Gráfico de picos vs score 
nueva_columna <- c()
matriz <- matrix(br, nrow = 244, ncol = 4) #244 es la cantidad de picos
for (i in 1:244) {
  nuevo_valor <- colnames(matint_19_ind_dico)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)
plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros")
# Crear un gradiente de colores (por ejemplo, de azul a rojo)
colores <- colorRampPalette(c("green4", "red2"))(244)
# Agregar puntos con colores en forma de gradiente
for (i in 1:244) {
  points(df_br$nueva_columna[i], df_br$V2[i], col = colores[i]) 
}
# Agregar puntos con relleno de colores en forma de gradiente
for (i in 1:244) {
  points(df_br$nueva_columna[i], df_br$V2[i], pch = 19, col = colores[i]) 
}

# Selección de picos mas preponderantes
top.b5 <- br[1:5]  ## primeros 5 picos
top.b10 <- br[1:10]  ## primeros 10 picos
top.b15 <- br[1:15]  ## primeros 15 picos
top.b20 <- br[1:20]  ## primeros 20 picos 
top_actual <- top.b20

# Elección de mejores algoritmos de clustering
comparacion <- clValid(
  obj        = matint_19_ind_dico[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)
summary(comparacion)
optimalScores(comparacion) #Se puede ir probando con distintos top picos


### ALGORITMO DE CLUSTERING ####################################################
################################################################################

# HKMEANS clustering con top20 y 2 clusters

top_actual <- top.b20
K.num <- 2 # clusters
var2 = 0.95

hkm.res20 <- hkmeans(matint_19_ind_dico[, top_actual], 
                   K.num)

cluster_hkmean20 <- fviz_cluster(hkm.res20, ellipse.type = "convex", 
                                data = matint_19_ind_dico[, top_actual],
                                ellipse.level = var2,
                                show.clust.cent = F, 
                                geom = "point", main = "INF VS NO INF - hkmeans - Top 20 - 2 cluster")

# Personalización del ploteo
cluster_hkmean20 <- cluster_hkmean20 + 
  geom_point(data = cluster_hkmean20$data, 
             aes(x = x, y = y, color = factor_tipo)) +
  scale_color_manual(values = c("maroon","steelblue4","steelblue4", "maroon" )) +
  scale_size_continuous(range = c(2, 4)) + 
  labs(color = "Cluster", size = "Sexo") +
  theme(legend.position = "right")

# Muestra el gráfico
print(cluster_hkmean20)
