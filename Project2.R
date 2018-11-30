# Group 5 Project 2 R scripts
##Load packages used today
## R v3.4+
library(tidyverse)
library(dplyr)
library(pathview)

##Set the directory path
setwd("~/Desktop/MICB405_Proj2")

##Read the files
read_delim(file="SaanichInlet_100m_ORFs_ko.cleaned.txt", col_names=TRUE, delim=",")
df <- read.table(file="SaanichInlet_100m_ORFs_ko.cleaned.txt", header = FALSE)

###Prokka rRNA output
SI_100m_rRNA <- read.delim("~/Desktop/MICB405_Proj2/SI_100m_rRNA.tsv")
prokka_rrna <- subset(SI_100m_rRNA, select = -c(gene, EC_number, COG)) %>% 
  filter(ftype == "rRNA")
  
prokka_trna <- subset(SI_100m_rRNA, select = -c(gene, EC_number, COG)) %>% 
  filter(ftype == "tRNA")

z <- subset(korch74, select = c(mag, prokka_id, orf_id, ko, Completeness, Contamination)) %>% #korch74 data frame is generated later on in the script, by bubble plots
  filter(Completeness > 90 & Contamination < 5) %>% 
  filter(mag %in% c("SaanichInlet_100m.14"))

unique(z$ko) #looks at number of unique ko values in MAG 14

###Chemical Data
chemdat <- read.csv(file="Saanich_TimeSeries_Chemical_DATA.csv")

###Look at the gtdbk data
gtdb_dat %>% 
  group_by(Phylum) %>% 
  summarise(count = n_distinct(mag)) %>% 
  kable()

#Contamination vs. Completion
library(cowplot)

checkm <- read.csv("MetaBAT2_SaanichInlet_100m_min1500_checkM_stdout.tsv", sep = "\t")

newvars3 <- c("Bin.Id", "Completeness", "Contamination", "Strain.heterogeneity")
checkm3  <- checkm[newvars3]

###CheckM value distribution boxplot
p1 <- checkm3 %>% 
  ggplot(aes(x="", y=Contamination)) +
  geom_boxplot() +
  labs(x="Contamination", y="Percentage (%)")

p2 <- checkm3 %>% 
  ggplot(aes(x="", y=Completeness)) +
  geom_boxplot() +
  labs(x="Completeness", y="Percentage (%)")

p <- plot_grid(p1, p2, labels=c("A", "B"), align="h", axis="tb", rel_widths=c(1/2, 1/2))

p

#Other tables for CheckM data
contam_complete = checkm2 %>%
  group_by(Contamination, Completeness)

##Merge CheckM stdout with gtdbk data
mergedat <- merge(checkm, gtdb_dat, by.x="Bin.Id", by.y="mag")

###By Phylum
mergedat %>%
  mutate(Phylum=substring(Phylum,4)) %>%
  ggplot(aes(x = Completeness, y = Contamination)) +
  geom_point(aes(color = Phylum)) 

###With bin numbers
mergedat %>%
  mutate(Phylum=substring(Phylum,4)) %>%
  mutate(BinID=substring(`Bin.Id`,18)) %>%
  ggplot(aes(x = Completeness, y = Contamination, color = Phylum, size = Contamination)) +
  geom_point() +
  geom_text(aes(label=BinID),
            colour="black",
            size=3) +
  labs(x = "Completeness (%)", y = "Contamination (%)")

mergedat %>%
  mutate(Family=substring(Family,4)) %>%
  mutate(BinID=substring(`Bin.Id`,18)) %>%
  ggplot(aes(x = Completeness, y = Contamination, color = Family, size = Contamination)) +
  geom_point() +
  geom_text(aes(label=BinID),
            colour="black",
            size=3) +
  labs(x = "Completeness (%)", y = "Contamination (%)")

###By Kingdom
mergedat %>%
  ggplot(aes(x = Completeness, y = Contamination)) +
  geom_point(aes(color = Kingdom, size = Contamination)) 

#RPKM bubble-plot of each Nitrogen/Sulphur-cycling gene versus taxonomy

###Load all the rpkm data for each Cruise and join them with the ko table (scroll to Pathview to view the code to generate the ko table)
ko_rpkm42 <- left_join(ko, rpkm42, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% 
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag")

ko_rpkm48 <- left_join(ko, rpkm48, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% 
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag")

ko_rpkm72 <- left_join(ko, rpkm72, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% 
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag")

ko_rpkm74 <- left_join(ko, rpkm74, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% 
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag")

kr <- data.frame(ko_rpkm42, ko_rpkm48, ko_rpkm72, ko_rpkm74)

###Example Bubble Plot
ko_rpkm %>%
  mutate(Family=substring(Family,4)) %>%
  filter(ko %in% c("K00370", "K02567", "K00362", "K03385")) %>%
  ggplot(aes(x = Family, y = ko)) +
  geom_point(aes(color = Family, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###Actual Bubble Plot

Nitro_genes <- c("nirB", "nirA", "narB", "nirK", "narG", "nasA", "nosZ", "anfG", "napA", "nifD", "nrfA", "norB", "NR", "hao", "pmoA", "nirS", "NIT-6", "K20932", "hdh", "vnfD")


###For Cruise SI074 (Nitrogen Pathway)
newvars <- c("Bin.Id", "Completeness", "Contamination")
checkm2  <- checkm[newvars]
korch74 <- merge(ko_rpkm74, checkm2, by.x="mag", by.y="Bin.Id")
korch74 %>%
  mutate(Family=substring(Family,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00362", "K00366", "K00367", "K00368", "K00370", "K00372", "K00376", "K00531", "K02567", "K02586", "K03385", "K04561", "K10534", "K10535", "K10944", "K15864", "K17877", "K20932", "K20935", "K22896")) %>%
  filter(rpkm > 0) %>% 
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Family, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels= Ng74fam) +
  labs(x = "mags", y = "Nitrogen pathway gene") +
  scale_size_continuous(range = c(7, 12))

###The labels for the y axis were added after the plot was run. The gene names  match the KO values in order. 
Ng74fam <- c("nirB", "nirA", "narB", "nirK", "narG", "nasA", "napA", "nrfA", "norB", "NR")

suppnit <- c("nirA", "narB", "nirK", "narG", "napA", "nrfA", "norB", "NR")

korch74 %>%
  mutate(Phylum=substring(Phylum,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00362", "K00366", "K00367", "K00368", "K00370", "K00372", "K00376", "K00531", "K02567", "K02586", "K03385", "K04561", "K10534", "K10535", "K10944", "K15864", "K17877", "K20932", "K20935", "K22896")) %>%
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Phylum, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels= Ng74phy) +
  labs(x = "mags", y = "Nitrogen pathway gene") +
  scale_size_continuous(range = c(3, 7))

Ng74phy <- c("nirB", "nirA", "narB", "nirK", "narG", "nasA", "napA", "nrfA", "norB", "NR")

###For Cruise SI048 (For Nitrogen Pathway)
korch48 <- merge(ko_rpkm48, checkm2, by.x="mag", by.y="Bin.Id")
korch48 %>%
  mutate(Family=substring(Family,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00362", "K00366", "K00367", "K00368", "K00370", "K00372", "K00376", "K00531", "K02567", "K02586", "K03385", "K04561", "K10534", "K10535", "K10944", "K15864", "K17877", "K20932", "K20935", "K22896")) %>%
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Family, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels= Ng48fam) +
  labs(x = "mags", y = "Nitrogen pathway gene") +
  scale_size_continuous(range = c(1, 3))

Ng48fam <- c("nirB", "nirA", "narB", "nirK", "narG", "nasA", "napA", "nrfA", "norB", "NR")

korch48 %>%
  mutate(Phylum=substring(Phylum,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00362", "K00366", "K00367", "K00368", "K00370", "K00372", "K00376", "K00531", "K02567", "K02586", "K03385", "K04561", "K10534", "K10535", "K10944", "K15864", "K17877", "K20932", "K20935", "K22896")) %>%
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Phylum, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels= Ng48phy) +
  labs(x = "mags", y = "Nitrogen pathway gene") +
  scale_size_continuous(range = c(1, 3))

Ng48phy <- c("nirB", "nirA", "narB", "nirK", "narG", "nasA", "napA", "nrfA", "norB", "NR")

###For Cruise SI_074 (For Sulfur Pathway)

korch74 %>%
  mutate(Family=substring(Family,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00380", "K00390", "K00392", "K00394", "K00955", "K00958",  "K11180", "K13811", "K17222", "K17223", "K17224", "K17225", "K17226", "K17227", "K22622")) %>%
  filter(rpkm > 0) %>% 
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Family, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "mags", y = "Sulfur pathway gene") +
  scale_y_discrete(labels= suppsulf) +
  scale_size_continuous(range = c(7, 12)) 

suppsulf <- c("cysJ", "cysH", "sir", "aprA", "cysNC", "sat", "soxA", "soxB", "soxY", "soxY")  
  


korch74 <- merge(ko_rpkm48, checkm2, by.x="mag", by.y="Bin.Id")
korch74 %>%
  mutate(Phylum=substring(Phylum,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00380", "K00390", "K00392", "K00394", "K00955", "K00958",  "K11180", "K13811", "K17222", "K17223", "K17224", "K17225", "K17226", "K17227", "K22622")) %>%
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Phylum, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels= Sg74phy) +
  labs(x = "mags", y = "Sulfur pathway gene") +
  scale_size_continuous(range = c(7, 12))

Sg74phy <- c("cysJ", "cysH", "sir", "aprA", "cysNC", "sat", "soxA", "soxX", "soxB", "soxY", "soxZ")


###For Cruise SI048 (For Sulfur Pathway)
korch48 <- merge(ko_rpkm48, checkm2, by.x="mag", by.y="Bin.Id")
korch48 %>%
  mutate(Family=substring(Family,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00380", "K00390", "K00392", "K00394", "K00955", "K00958",  "K11180", "K13811", "K17222", "K17223", "K17224", "K17225", "K17226", "K17227", "K22622")) %>%
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Family, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels= Sg48fam) +
  labs(x = "mags", y = "Sulfur pathway gene") +
  scale_size_continuous(range = c(1, 3))

Sg48fam <- c("cysJ", "cysH", "sir", "aprA", "cysNC", "sat", "soxA", "soxX", "soxB", "soxY", "soxZ")

korch48 <- merge(ko_rpkm48, checkm2, by.x="mag", by.y="Bin.Id")
korch48 %>%
  mutate(Phylum=substring(Phylum,4)) %>%
  mutate(mag=substring(`mag`,18)) %>%
  filter(Completeness > 50 & Contamination < 10) %>%
  filter(ko %in% c("K00380", "K00390", "K00392", "K00394", "K00955", "K00958",  "K11180", "K13811", "K17222", "K17223", "K17224", "K17225", "K17226", "K17227", "K22622")) %>%
  ggplot(aes(x= mag, y = ko)) +
  geom_point(aes(color = Phylum, size = rpkm)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels= Sg48phy) +
  labs(x = "mags", y = "Sulfur pathway gene") +
  scale_size_continuous(range = c(1, 3))

Sg48phy <- c("cysJ", "cysH", "sir", "aprA", "cysNC", "sat", "soxA", "soxX", "soxB", "soxY", "soxZ")

#Connor's Method for bubble plot (Metagenomic Data)

metag_rpkm <- read.table("~/Desktop/MICB405_Proj2/SaanichInlet_100m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence = gsub('m_', 'm.', Sequence)) %>% 
  mutate(Sequence = gsub('Inlet_', 'Inlet.', Sequence)) %>% 
  separate(col=Sequence, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))

checkm_dat <- read.table("~/Desktop/MICB405_Proj2/MetaBAT2_SaanichInlet_100m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)

rpkm_dat <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  filter(Completeness > 50 & Contamination < 10) %>% # good quality MAGs
  filter(Phylum == "p__Proteobacteria") %>% # Proteobacteria only - we kept to this because it's the most dominant bacteria phylum at 100m
  mutate(mag = reorder(mag, Order, sort)) %>% # sort by their taxonomic Order so everything shows up together
  mutate(mag=substring(`mag`,18)) %>% 
  mutate(family=substring(Family,4))
  
  

ggplot(rpkm_dat, aes(x=Sample, y=mag, col=Family)) +
  geom_point(aes(size=g_rpkm)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 90),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(x = "Cruise", y = "mag (Saanich Inlet 100m)")


#Pathview

ko <- read.table("~/Desktop/MICB405_Proj2/SaanichInlet_100m_ORFs_ko.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)


rpkm <- read.table("~/Desktop/MICB405_Proj2/SI042_SaanichInlet_ORFs_RPKM.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)

prokka_mag_map <- read.table("Prokka_MAG_map.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)

arc_class <- read.table("gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))


ko_rpkm <- left_join(ko, rpkm, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% 
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag")

#by Family
t_rpkm <- ko_rpkm %>% 
  group_by(Family, ko) %>% 
  summarise(total = sum(rpkm)) %>% 
  spread(key = Family, value = total)
###We repeated this exact script for all five cruises at 100m depth. For the sake of simplicity, I'm not going to 

pv_mat <- dplyr::select(t_rpkm, -ko)
rownames(pv_mat) <- t_rpkm$ko

pv_mat %>% 
  filter(ko == "K00370", ko == "K02567", ko == "K00362", ko == "K03385") %>%
  ggplot(aes(y = ko)) +
  geom_point(aes(color = Family, size = rpkm)) 



# Nitrogen metabolism
pv.out <- pathview(gene.data = pv_mat,
                   species = "ko",
                   pathway.id="00910",
                   kegg.dir = "~/Desktop/MICB405_Proj2/pathview_output.dir/")


# Sulfer metabolism
pv.out <- pathview(gene.data = pv_mat,
                   species = "ko",
                   pathway.id="00920",
                   kegg.dir = "~/Desktop/MICB405_Proj2/pathview_output.dir/")

# Carbon fixation pathways in prokaryotes 00720
pv.out <- pathview(gene.data = pv_mat,
                   species = "ko",
                   pathway.id="00720",
                   kegg.dir = "~/Desktop/MICB405_Proj2/pathview_output.dir/")


