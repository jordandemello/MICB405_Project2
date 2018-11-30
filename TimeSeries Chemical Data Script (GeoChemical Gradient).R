cat("\014")
library(KEGGREST)
library(tidyverse)
library(cowplot)
library(pathview)
library(magrittr)

raw_data <- read.csv("Saanich_TimeSeries_Chemical_DATA.csv", header=TRUE)

gathered_data <- raw_data %>%
  dplyr::rename(Cell_Density=16) %>%
  dplyr::mutate(H2S=(as.numeric(as.character(Mean_H2S))), N2=(as.numeric(as.character(Mean_N2))), O2=(as.numeric(as.character(Mean_O2))), CO2=(as.numeric(as.character(Mean_co2))), NO2=(as.numeric(as.character(Mean_NO2))), NH4=(as.numeric(as.character(Mean_NH4))), CH4=(as.numeric(as.character(Mean_CH4)))/1000, N2O=(as.numeric(as.character(Mean_N2O)))/1000) %>%
  dplyr::select(Depth, PO4, SI, NH4, NO3, NO2, H2S, N2, CO2, N2O, O2, CH4) %>% 
  gather(key="Key", value="Value", na.rm=TRUE, -Depth) %>%
  dplyr::filter(Value >= 0 & Value != "NAN" & Value != "ND" & Value != "NaN") %>%
  dplyr::mutate(Depth_Bucket=(floor(as.numeric(as.character(Depth))/10)*10), Concentration=as.numeric(as.character(Value)), Chemical=as.character(Key)) %>%
  dplyr::select(Depth, Chemical, Concentration) %>%
  dplyr::arrange(Depth)
  
summarized_data <- gathered_data %>%
  group_by(Depth, Chemical) %>%
  summarize(Average_Concentration = mean(Concentration), Count=n()) %>%
  dplyr::filter(Count>1)
  
max_summarized_data <- summarized_data %>%
  group_by(Chemical) %>%
  dplyr::select(Chemical, Average_Concentration) %>%
  summarize(Mean_Concentration=(mean(Average_Concentration)))
  
relative_data <- summarized_data %>%
  left_join(max_summarized_data, by="Chemical") %>%
  dplyr::mutate(Relative_Level = ((Average_Concentration)/Mean_Concentration))

cell_data <- raw_data %>%
  dplyr::rename(Cell_Density=16) %>%
  dplyr::select(Depth, Cell_Density) %>%
  dplyr::filter(Cell_Density != "NAN" & Cell_Density != "ND" & Cell_Density != "NaN") %>%
  dplyr::mutate(Depth_Bucket=(floor(as.numeric(as.character(Depth))/10)*10), Concentration=as.numeric(as.character(Cell_Density))) %>%
  dplyr::filter(Concentration >= 0) %>%
  dplyr::select(Depth, Concentration) %>%
  dplyr::arrange(Depth)

summarized_data%>%
  ggplot(aes(x=Depth, y=Average_Concentration, color=Chemical), group=Chemical) +
  geom_smooth(method = loess, se = FALSE) +
  labs(title = "Absolute Geochemical Gradient over Depth", x = "Depth (m)", y = "Chemical Concentration (μM)") +
  coord_flip() + 
  geom_vline(xintercept=100, linetype="dashed", color = "red") +
  scale_x_reverse() +
  theme_classic()

relative_data%>%
  ggplot(aes(x=Depth, y=Relative_Level, color=Chemical), group=Chemical) +
  geom_smooth(method = loess, se = FALSE) +
  labs(title = "Relative Geochemical Gradient over Depth", x = "Depth (m)", y = "Fraction of Average Concentration") +
  coord_flip() + 
  geom_vline(xintercept=100, linetype="dashed", color = "red") +
  scale_x_reverse() +
  theme_classic()

cell_data%>%
  ggplot(aes(x=Depth, y=Concentration)) +
  geom_smooth(method = loess, se = FALSE) +
  labs(title = "Cell Concentration Gradient over Depth", x = "Depth (m)", y = "Cell Concentration (Cells/mL)") +
  coord_flip() + 
  geom_vline(xintercept=100, linetype="dashed", color = "red") +
  scale_x_reverse() +
  theme_classic()
