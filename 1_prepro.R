################ MALDI-TOF ANALISIS VINCHUCA ###################################
################ 1) PREPROCESAMIENTO  ##########################################
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################


library(binda)
library(here)
library(dplyr)
library(readBrukerFlexData)
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(stringr)


### CARGA DE ESPECTROS #########################################################
################################################################################


# Creación de la ruta relativa de los archivos
ruta_proyecto <- "C:/Users/urtea/OneDrive/Documents/Proyectos/MALDI_Vinchucas/Datos"
ruta_datos <- file.path(ruta_proyecto)

# Importar espectros
Spectra_list <- importBrukerFlex(file.path(ruta_datos), verbose=FALSE)


### OBTENCIÓN DE METADATA DE ESPECTROS #########################################
################################################################################

# Creación de columnas vacías
col_num <- c()
col_estado <- c()
col_numero <- c()
col_sexo <- c()
col_rep_m <- c()
col_rep_t <- c()

# Patrones auxiliares para buscar el día de la muestra
patron_i <- "Vinchucas infectadas"
patron_ni <- "Vinchucas no infectadas"

# Ciclo que extrae dia, tipo, numero, well y réplica de cada muestra
for(i in 1:length(Spectra_list)) {
  nombre <- Spectra_list[[i]]@metaData$file
  
  # Encuentra la posición del patrón infectadas
  posicion_i <- str_locate(nombre, patron_i)[1, 2]
  # Encuentra la posición del patrón infectadas
  posicion_ni <- str_locate(nombre, patron_ni)[1, 2]

  # Verifica si posicion_ni es NA
  if (is.na(posicion_ni)) {
    posicion <- posicion_i
    estado <- "infectado"
  }
  else {
    posicion <- posicion_ni
    estado <- "no infectado"
  }
  
  # Agrega estado a la columna correspondiente
  col_estado <- c(col_estado, estado)
  
  # Extrae la parte restante del nombre después del estado
  resultado <- substr(nombre, (posicion+2), nchar(nombre))
  
  # Usa expresiones regulares para extraer los datos:
  # Ejemplo: "Vinchuca 1 hembra 23-05-24\\0_A1\\1\\1SLin\\fid"
  numero_vinchuca <- str_extract(resultado, "Vinchuca [0-9]+") %>% str_extract("[0-9]+")
  sexo <- str_extract(resultado, "hembra|macho")
  replica_muestra <- str_extract(resultado, "[A-Z][0-9]+")
  replica_tecnica <- str_extract(resultado, "\\\\[0-9]+\\\\1SLin") %>% str_extract("[0-9]+")
  
  # Almacena los valores extraídos en sus respectivas columnas
  col_numero <- c(col_numero, numero_vinchuca)
  col_sexo <- c(col_sexo, sexo)
  col_rep_m <- c(col_rep_m, replica_muestra)
  col_rep_t <- c(col_rep_t, replica_tecnica)
  
  # print(paste("Estado:", estado, "| Número:", numero_vinchuca, "| Sexo:", sexo, 
  #             "| Réplica Muestra:", replica_muestra, "| Réplica Técnica:", replica_tecnica))
  
  Spectra_list[[i]]@metaData$numero_vinchuca <- numero_vinchuca
  Spectra_list[[i]]@metaData$estado <- estado
  Spectra_list[[i]]@metaData$sexo <- sexo
  Spectra_list[[i]]@metaData$rep_m <- replica_muestra
  Spectra_list[[i]]@metaData$rep_t <- replica_tecnica
  
  
  # Data Frame con los datos limpios
  df_metadata <- data.frame(estado = col_estado, numero = col_numero, sexo = col_sexo, 
                            rep_m = col_rep_m, rep_t = col_rep_t)
  
  # Creación de factores de agrupamiento para su uso posterior
  df_metadata$factor_num <- paste0(df_metadata$estado, "_", df_metadata$numero)
  df_metadata$factor_mue <- paste0(df_metadata$estado, "_", df_metadata$rep_m)
  df_metadata$factor_sex <- paste0(df_metadata$estado, "_", df_metadata$sexo)
  
}

### CONTROL DE CALIDAD Y LIMPIEZA DE ESPECTROS #################################
################################################################################


# Screening inicial: Detección de espectros de baja calidad
sc.results <- screenSpectra(Spectra_list, meta = df_metadata)
summary(sc.results)
plot(sc.results, labels = TRUE)

# plot(Spectra_list[[253]]) # Ploteo de espectros ruidosos

# Descartamos espectros defectuosos 
Spectra_list_f1 <- sc.results$fspectra # Filtramos espectros
df_metadata_f1 <- sc.results$fmeta # Filtramos metadatos

### FILTRADO Y TRANSFORMACIÓN DE ESPECTROS #####################################
################################################################################

# Parámetros de procesamiento de espectros
thScale <- 10 # Smoothing
ite <- 105 # Baseline correction
SigNoi <- 2.5 # Peak extraction
hws <- 20 # Peak extraction
tol <- 0.03 # Peak binning

# Transformación/filtrado/corrección de espectros con parámetros definidos
# 1) Transformación de intensidad por medio de función sqrt
Spectra_list_f1 <- transfIntensity(Spectra_list_f1, fun = sqrt)
plot(Spectra_list_f1[[30]])
# 2) Suavizado del espectro mediante el método Wavelet
Spectra_list_f1 <- wavSmoothing(Spectra_list_f1, method = "Wavelet", n.levels = 4)
plot(Spectra_list_f1[[30]])
# Detección de la linea de base
baseline <- estimateBaseline(Spectra_list_f1[[30]], method = "SNIP",
                             iterations = ite)
plot(Spectra_list_f1[[30]])
lines(baseline, col="red", lwd=2)
# 3) Remoción de linea de base mediante método SNIP
Spectra_list_f2 <- removeBaseline(Spectra_list_f1, method = "SNIP", 
                                  iterations = ite)
plot(Spectra_list_f2[[30]])
# 4) Calibración de intensidad mediante método PQN
Spectra_list_f2 <- calibrateIntensity(Spectra_list_f2, method = "PQN")
plot(Spectra_list_f2[[30]])
# 5) Alineación de espectros
Spectra_list_f3 <- alignSpectra(Spectra_list_f2,
                                halfWindowSize=20,
                                SNR=2,
                                tolerance=0.02, # Parámetro sensible
                                warpingMethod="lowess")
plot(Spectra_list_f3[[30]])


### PROMEDIO DE LECTURAS DE UNA MISMA RÉPLICA TÉCNICA #################
################################################################################


# Promedio de lecturas de una misma well
Spectra_list_prom_rep <- averageMassSpectra(Spectra_list_f3,
                                            labels = factor(df_metadata_f1$factor_mue), 
                                            method = "mean")

# Creo la nueva metadata de los espectros promediados
df_metadata_prom_rep <- df_metadata_f1 %>% 
  distinct(factor_mue, .keep_all = TRUE)

# Promedio de wells de una misma muestra
Spectra_list_prom_muestra <- averageMassSpectra(Spectra_list_prom_rep,
                                                labels = factor(df_metadata_prom_rep$factor_num), 
                                                method = "mean")

# Creo la nueva metadata de los espectros promediados
df_metadata_prom_mue <- df_metadata_prom_rep %>% 
  distinct(factor_num, .keep_all = TRUE)


### EXTRACCIÓN DE PICOS Y ALINEACIÓN ###########################################
################################################################################

# A partir de acá probamos trabajar con Spectra_list_prom_rep

# Análisis de la SNR en espectros para chequear que utilizamos el valor correcto
noise <- estimateNoise(Spectra_list_prom_rep[[20]])
plot(Spectra_list_prom_rep[[20]], xlim=c(4000, 20000), ylim=c(0, 0.002))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue") # Se ve que es correcto el 2

# Detección de picos a partir del umbral definido de SNR
peaks <- detectPeaks(Spectra_list_prom_rep, 
                     SNR = SigNoi, 
                     halfWindowSize = 40)

# Ploteo de picos detectados en un espectro de ejemplo
plot(Spectra_list_prom_rep[[20]], xlim=c(4000, 20000), ylim=c(0, 0.002))
points(peaks[[20]], col="red", pch=4)

# Alineado de picos
peaks <- alignPeaks(peaks, 
                    minFreq = 0.8, 
                    tolerance = tol)

#summaryPeaks(peaks[1:10])  # resumen estadistico de picos (primeros 10)

# Conteo de picos por perfil
cP <- countPeaks(peaks)

# Gráfico de picos
plot(cP, type = "n")
text(cP, label = 1:length(cP))


# Patrones de picos
peakPatterns(peaks)

# Filtrado de picos de baja frecuencia de aparición
picos_filtrados <- filterPeaks(peaks, 
                               minFreq = 0.25, 
                               labels = df_metadata_prom_rep$factor_num ) #labels

# Patrones de picos
peakPatterns(picos_filtrados)

# Conteo de picos por perfil
cP2 <- countPeaks(picos_filtrados)

# Gráfico
plot(cP2, type = "n")
text(cP2, label = 1:length(cP2))

# Fusión de picos de la misma muestra
picos_fusion_muestra <- mergeMassPeaks(picos_filtrados, 
                                       labels = df_metadata_prom_rep$factor_num, 
                                       method = "median")

# Patrones de picos
peakPatterns(picos_fusion_muestra)



### CREACIÓN DE MATRIZ DE INTENSIDADES Y DICOTÓMICA ############################
################################################################################


# Matriz de intensidades 19 individuos
matint_19_ind <- intensityMatrix(picos_fusion_muestra, 
                              Spectra_list_prom_muestra) # sin valores NA

# Matriz de intensidades de 80 muestras
# matint_na_51 <- intensityMatrix(picos_fusion_muestra) # con valores NA
matint_80_mue <- intensityMatrix(picos_filtrados, 
                             Spectra_list_prom_rep) # sin valores NA

# Definición de umbrales
thr1 <- optimizeThreshold(matint_19_ind,
                          df_metadata_prom_mue$estado, 
                          verbose = T)
thr2 <- optimizeThreshold(matint_80_mue,
                          df_metadata_prom_rep$estado, 
                          verbose = T)

# Dicotomización
matint_19_ind_dico <- dichotomize(matint_19_ind, thr1)
matint_80_mue_dico <- dichotomize(matint_80_mue, thr2)

# Agrego nombres a las filas de cada df
rownames(matint_19_ind_dico) <- df_metadata_prom_mue$factor_num
rownames(matint_80_mue_dico) <- df_metadata_prom_rep$factor_mue
rownames(matint_19_ind) <- df_metadata_prom_mue$factor_num
rownames(matint_80_mue) <- df_metadata_prom_rep$factor_mue


### GUARDAR DATOS ##############################################################
################################################################################

# Establecer el directorio de trabajo en la ubicación del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Se puede mejorar, subcarpeta

# Guardar matrices y metadata asociada como archivo .Rdata
save(matint_19_ind_dico,df_metadata_prom_mue, file = "matint_19_ind_dico.Rdata")
save(matint_80_mue_dico, df_metadata_prom_rep, file = "matint_80_mue_dico.Rdata")
save(matint_19_ind,df_metadata_prom_mue, file = "matint_19_ind.Rdata")
save(matint_80_mue, df_metadata_prom_rep, file = "matint_80_mue.Rdata")

# Guardar matrices y metadata asociada como archivo .csv
write.csv(matint_19_ind_dico, "matint_19_ind_dico.csv", row.names = TRUE)
write.csv(matint_19_ind, "matint_19_ind.csv", row.names = TRUE)
write.csv(matint_80_mue_dico, "matint_80_mue_dico.csv", row.names = TRUE)
write.csv(matint_80_mue, "matint_80_mue.csv", row.names = TRUE)
write.csv(df_metadata_prom_mue, "df_19.csv", row.names = TRUE)
write.csv(df_metadata_prom_rep, "df_80.csv", row.names = TRUE)

#
#
#
### FIN ########################################################################
################################################################################