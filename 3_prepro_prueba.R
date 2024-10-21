################ MALDI-TOF ANALISIS VINCHUCA ###################################
################ 3) PREPROCESAMIENTO de datos prueba ###########################
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
#ruta_proyecto <- "C:/Users/urtea/OneDrive/Documents/Proyectos/MALDI_Vinchucas/Datos"
ruta_proyecto <- "C:/Users/Facundo/Documents/Proyectos/MALDI_Vinchucas/Datos_prueba"
ruta_datos <- file.path(ruta_proyecto)

# Importar espectros
Spectra_list <- importBrukerFlex(file.path(ruta_datos), verbose=FALSE)


# for(i in 1:length(Spectra_list)) {
#   print(Spectra_list[[i]]@metaData$file)
# }

### OBTENCIÓN DE METADATA DE ESPECTROS #########################################
################################################################################

# Creación de columnas vacías
col_num <- c()
col_estado <- c()
col_numero <- c()
col_sexo <- c()
col_rep_m <- c()
col_rep_t <- c()

# # Patrones auxiliares para buscar el día de la muestra
patron_h <- "embra"
patron_m <- "Macho"

# Ciclo que extrae dia, tipo, numero, well y réplica de cada muestra
for(i in 1:length(Spectra_list)) {
  nombre <- Spectra_list[[i]]@metaData$file

  # Encuentra la posición del patrón infectadas
  posicion_h <- str_locate(nombre, patron_h)[1, 2]
  # Encuentra la posición del patrón infectadas
  posicion_m <- str_locate(nombre, patron_m)[1, 2]
  
  # Verifica si posicion_ni es NA
  if (is.na(posicion_m)) {
    sexo <- "Hembra"
    
    # Extrae la parte restante del nombre después del estado
    resultado <- substr(nombre, (posicion_h+2), nchar(nombre))
    
    # Usa expresiones regulares para extraer los datos:
    # Ejemplo: "Vinchuca 1 hembra 23-05-24\\0_A1\\1\\1SLin\\fid"
    numero_vinchuca <- str_extract(resultado, "Hembra [0-9]+") %>% str_extract("[0-9]+")
    estado <- "infectado"
    replica_muestra <- str_extract(resultado, "[A-Z][0-9]+")
    replica_tecnica <- str_extract(resultado, "\\\\[0-9]+\\\\1SLin") %>% str_extract("[0-9]+")
    
  }
  else {
    posicion <- posicion_m
    sexo <- "Macho"
    
    # Extrae la parte restante del nombre después del estado
    resultado <- substr(nombre, (posicion_m+2), nchar(nombre))
    
    # Usa expresiones regulares para extraer los datos:
    # Ejemplo: "Vinchuca 1 hembra 23-05-24\\0_A1\\1\\1SLin\\fid"
    numero_vinchuca <- str_extract(resultado, "Macho [0-9]+") %>% str_extract("[0-9]+")
    estado <- "infectado"
    replica_muestra <- str_extract(resultado, "[A-Z][0-9]+")
    replica_tecnica <- str_extract(resultado, "\\\\[0-9]+\\\\1SLin") %>% str_extract("[0-9]+")
  }
  
  # Almacena los valores extraídos en sus respectivas columnas
  col_estado <- c(col_estado, estado)
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