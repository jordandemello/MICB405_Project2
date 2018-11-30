cat("\014")
library(KEGGREST)
library(tidyverse)
library(cowplot)
library(pathview)
library(magrittr)

raw_data <- read.csv("Saanich_TimeSeries_Chemical_DATA.csv", header=TRUE)


filtered_data <- raw_data %>%
  dplyr::filter(Depth == 100)

PO4_data <- filtered_data %>%
  dplyr::select(Date, PO4) %>%
  dplyr::filter(PO4 != "NAN" & PO4 != "ND" & PO4 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(PO4))) %>%
  dplyr::filter(Y_Bin < 50)

SI_data  <- filtered_data %>%
  dplyr::select(Date, SI) %>%
  dplyr::filter(SI != "NAN" & SI != "ND" & SI != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(SI)))

NO3_data <- filtered_data %>%
  dplyr::select(Date, NO3) %>%
  dplyr::filter(NO3 != "NAN" & NO3 != "ND" & NO3 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(NO3)))

NH4_Mean_data <- filtered_data %>%
  dplyr::select(Date, Mean_NH4) %>%
  dplyr::filter(Mean_NH4 != "NAN" & Mean_NH4 != "ND" & Mean_NH4 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(Mean_NH4))) %>%
  dplyr::filter(Y_Bin < 20)

NO2_Mean_data <- filtered_data %>%
  dplyr::select(Date, Mean_NO2) %>%
  dplyr::filter(Mean_NO2 != "NAN" & Mean_NO2 != "ND" & Mean_NO2 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(Mean_NO2)))

H2S_Mean_data <- filtered_data %>%
  dplyr::select(Date, Mean_H2S) %>%
  dplyr::filter(Mean_H2S != "NAN" & Mean_H2S != "ND" & Mean_H2S != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(Mean_H2S)))

N2_Mean_data  <- filtered_data %>%
  dplyr::select(Date, Mean_N2) %>%
  dplyr::filter(Mean_N2 != "NAN" & Mean_N2 != "ND" & Mean_N2 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(Mean_N2)))

O2_Mean_data  <- filtered_data %>%
  dplyr::select(Date, Mean_O2) %>%
  dplyr::filter(Mean_O2 != "NAN" & Mean_O2 != "ND" & Mean_O2 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(Mean_O2)))

CO2_Mean_data <- filtered_data %>%
  dplyr::select(Date, Mean_co2) %>%
  dplyr::filter(Mean_co2 != "NAN" & Mean_co2 != "ND" & Mean_co2 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(Mean_co2)))

N2O_Mean_data <- filtered_data %>%
  dplyr::select(Date, Mean_N2O) %>%
  dplyr::filter(Mean_N2O != "NAN" & Mean_N2O != "ND" & Mean_N2O != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = (as.numeric(as.character(Mean_N2O)))/1000)

CH4_Mean_data <- filtered_data %>%
  dplyr::select(Date, Mean_CH4) %>%
  dplyr::filter(Mean_CH4 != "NAN" & Mean_CH4 != "ND" & Mean_CH4 != "NaN") %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = (as.numeric(as.character(Mean_CH4)))/1000) %>%
  dplyr::filter(Y_Bin < 0.4)

Cell_Density_data <- filtered_data %>%
  dplyr::select(Date, c(16)) %>%
  dplyr::rename(Cell_Density = 2) %>%
  dplyr::filter(Cell_Density != "NAN" & Cell_Density != "ND" & Cell_Density != "NaN") %>%
  dplyr::arrange(Date)%>%
  dplyr::mutate(Date_Bin = as.Date(Date)) %>%
  dplyr::mutate(Y_Bin = as.numeric(as.character(Cell_Density))) %>%
  dplyr::filter(Y_Bin < 3000000)
  
  
CH4_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in CH4 Levels at 100m Over Time", x = "Date of Sample", y = "[CH4] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

CO2_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in CO2 Levels at 100m Over Time", x = "Date of Sample", y = "[CO2] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

H2S_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in H2S Levels at 100m Over Time", x = "Date of Sample", y = "[H2S] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

N2_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in N2 Levels at 100m Over Time", x = "Date of Sample", y = "[N2] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

N2O_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in N2O Levels at 100m Over Time", x = "Date of Sample", y = "[N2O] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

NH4_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in NH4 Levels at 100m Over Time", x = "Date of Sample", y = "[NH4] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

NO2_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in NO2 Levels at 100m Over Time", x = "Date of Sample", y = "[NO2] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

NO3_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in NO3 Levels at 100m Over Time", x = "Date of Sample", y = "[NO3] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

O2_Mean_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in O2 Levels at 100m Over Time", x = "Date of Sample", y = "[O2] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

PO4_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in PO4 Levels at 100m Over Time", x = "Date of Sample", y = "[PO4] in μM") +
  geom_smooth(method = loess) +
  theme_classic()
  
SI_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in SI Levels at 100m Over Time", x = "Date of Sample", y = "[SI] in μM") +
  geom_smooth(method = loess) +
  theme_classic()

Cell_Density_data %>%
  ggplot(aes(x=Date_Bin, y=Y_Bin)) + 
  geom_point() + 
  labs(title = "Fluctuations in Cell Density Levels at 100m Over Time", x = "Date of Sample", y = "Cell Density in Cells/mL") +
  geom_smooth(method = loess) +
  theme_classic()
