library(ARTofR)
# READ IN THE DATA------------------------------------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/Halls-Cave-Lizards/main/Fossil_lizard_15bin.csv")
Fossil_lizard_15bin <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = FALSE, na.strings=c("","NA")) # this is a matrix of measured specimens 
head(Fossil_lizard_15bin)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Visualize data                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(ggplot2)

age_depth <- function(Depth) { # Function to assign ages to depths based on linear relationship defined by Tomé et al. 2021
  68.61*(Depth) + 758.82
}
bins <- seq(from = 0, to = 295, by = 15) #VECTOR OF 15 CM BINS
age <- sapply(bins, age_depth) # ASSIGN AGES TO ALL BINS


LIZVIZ <- Fossil_lizard_15bin %>% group_by(Bin) %>% summarise(NISP = sum(NISP))
LIZVIZ$age <- age
LIZVIZ_plot<-ggplot(LIZVIZ, aes(age, NISP))+ geom_point()+
  geom_line(colour="blue")+
  scale_x_reverse(breaks =seq(0,20000,2000))+
  scale_y_continuous(breaks =seq(0,2000,300))+
  xlab("Age (cal. BP)")+ylab("NISP")+
  coord_flip() +
  theme_classic()
LIZVIZ_plot


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     PREPARING THE DATA   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)
levels(as.factor(Fossil_lizard_15bin$Family)) #find all levels in data
Fossil_lizard_15bin_fam <- filter(Fossil_lizard_15bin,!grepl('cf|Pleurodonta',Family)) #filter out uncertain identifications

#MAKE VARIABLES FACTORS
Fossil_lizard_15bin_fam$Bin <- as.factor(Fossil_lizard_15bin_fam$Bin)
Fossil_lizard_15bin_fam$Family <- as.factor(Fossil_lizard_15bin_fam$Family)
Fossil_lizard_15bin_fam$Details <- as.factor(Fossil_lizard_15bin_fam$Details)
Fossil_lizard_15bin_fam$Element <- as.factor(Fossil_lizard_15bin_fam$Element)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     PREPARING THE DATA (NISP BASED ON ALL ELEMENTS)   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SUMMARIZE INTO NISP DATA
LIZNISP <- Fossil_lizard_15bin_fam %>% 
  filter(!is.na(Details))  %>%
  filter(!is.na(Family)) %>% 
  group_by(Family,Bin, .drop=FALSE) %>% tally
LIZNISP <- LIZNISP %>% group_by(Family, Bin) %>% summarise(NISP = max(n))



LIZNISP2 <- LIZNISP %>% group_by(Bin) %>% summarise(NISP = sum(NISP))
LIZNISP2$age <- age
LIZVIZ2_plot<-ggplot(LIZNISP2, aes(age, NISP))+ geom_point()+
  geom_line(colour="blue")+
  scale_x_reverse(breaks =seq(0,20000,2000))+
  scale_y_continuous(breaks =seq(0,300,30))+
  xlab("Age (cal. BP)")+ylab("NISP")+
  coord_flip() +
  theme_classic()
LIZVIZ2_plot





# MAKE INTO LONG FORMAT WITH FAMILIES AS COLUMNS AND BINS AS ROWS
library(tidyverse)
library(tidyr)
LIZNISP_wide <- spread(LIZNISP, Family, NISP)
LIZNISP_wide <- LIZNISP_wide %>% remove_rownames %>% column_to_rownames(var="Bin")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  PREPARING THE DATA (NISP BASED ON ALL ELEMENTS EXCEPT DENTARY AND MAXILLA)   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fossil_lizard_15bin_other_ele <- filter(Fossil_lizard_15bin_fam, !grepl('dentary|\\bmaxilla\\b',Element)) #don't include counts of dentary and the maxilla

# SUMMARIZE INTO NISP DATA
LIZNISP_other <- Fossil_lizard_15bin_other_ele %>% 
  filter(!is.na(Details))  %>%
  filter(!is.na(Family)) %>% 
  group_by(Family,Bin, .drop=FALSE) %>% tally
LIZNISP_other <- LIZNISP_other %>% group_by(Family, Bin) %>% summarise(NISP = max(n))

# MAKE INTO LONG FORMAT WITH FAMILIES AS COLUMNS AND BINS AS ROWS
library(tidyverse)
library(tidyr)
LIZNISP_other_wide <- spread(LIZNISP_other, Family, NISP)
LIZNISP_other_wide <- LIZNISP_other_wide %>% remove_rownames %>% column_to_rownames(var="Bin")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     PREPARING THE DATA (MNI BASED ON THE DENTARY AND MAXILLA)   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fossil_lizard_15bin_MNI <- filter(Fossil_lizard_15bin_fam, grepl('dentary|\\bmaxilla\\b',Element)) #only include counts of dentary and the maxilla
Fossil_lizard_15bin_MNI<-droplevels(Fossil_lizard_15bin_MNI)

# SUMMARIZE INTO MNI DATA
LIZMNI <- Fossil_lizard_15bin_MNI %>% 
  filter(!is.na(Details))  %>%
  filter(!is.na(Family)) %>% 
  group_by(Family,Bin,Element,Details, .drop=FALSE) %>% tally
LIZMNI <- LIZMNI %>% group_by(Family, Bin) %>% summarise(MNI = max(n))

# MAKE INTO LONG FORMAT WITH FAMILIES AS COLUMNS AND BINS AS ROWS
library(tidyverse)
library(tidyr)
LIZMNI_wide <- spread(LIZMNI, Family, MNI)
LIZMNI_wide <- LIZMNI_wide %>% remove_rownames %>% column_to_rownames(var="Bin")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Pollen plots                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(grid)
library(analogue)
LIZMNI.pct <- data.frame(tran(LIZMNI_wide, method = 'percent')) # CONVERT MNI INTO PERCENTS

Stratiplot(age ~ ., LIZMNI.pct, sort = 'wa', type = 'poly',
           ylab ='Years Before Present')

MNIdf <- data.frame(yr=rep(age,ncol(LIZMNI.pct)),
              per=as.vector(as.matrix(LIZMNI.pct)),
              taxa=as.factor(rep(colnames(LIZMNI.pct),each=nrow(LIZMNI.pct))))

theme_new <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grids
                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                   strip.text.x = element_text(size=10, angle=90, vjust=0), # Taxa names
                   strip.background = element_blank(),
                   strip.text.y = element_text(angle = 0),
                   legend.position="none",panel.border = element_blank(),
                   axis.text.x=element_text(angle=45,hjust=1)) # Axis tick label angle

MNIplot<- ggplot(MNIdf)+
  geom_line(aes(yr,per))+
  geom_area(aes(yr,per))+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,100,10))+
  xlab("Age (cal. BP)")+ylab("%")+
  coord_flip()+
  theme_new+
  facet_grid(~MNIdf$taxa,scales = "free", space = "free")


LIZNISP.pct <- data.frame(tran(LIZNISP_wide, method = 'percent')) # CONVERT MNI INTO PERCENTS

Stratiplot(age ~ ., LIZNISP.pct, sort = 'wa', type = 'poly',
           ylab ='Years Before Present')

NISPdf <- data.frame(yr=rep(age,ncol(LIZNISP.pct)),
                    per=as.vector(as.matrix(LIZNISP.pct)),
                    taxa=as.factor(rep(colnames(LIZNISP.pct),each=nrow(LIZNISP.pct))))
NISPplot<-ggplot(NISPdf)+
  geom_line(aes(yr,per))+
  geom_area(aes(yr,per))+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,100,10))+
  xlab("Age (cal. BP)")+ylab("%")+
  coord_flip()+
  theme_new+
  facet_grid(~NISPdf$taxa,scales = "free", space = "free")

LIZNISP_OTHER.pct <- data.frame(tran(LIZNISP_other_wide, method = 'percent')) # CONVERT MNI INTO PERCENTS

NISPdf_other <- data.frame(yr=rep(age,ncol(LIZNISP_OTHER.pct)),
                     per=as.vector(as.matrix(LIZNISP_OTHER.pct)),
                     taxa=as.factor(rep(colnames(LIZNISP_OTHER.pct),each=nrow(LIZNISP_OTHER.pct))))
NISPother_plot<-ggplot(NISPdf_other)+
  geom_line(aes(yr,per))+
  geom_area(aes(yr,per))+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,100,10))+
  xlab("Age (cal. BP)")+ylab("%")+
  coord_flip()+
  theme_new+
  facet_grid(~NISPdf_other$taxa,scales = "free", space = "free")


library(ggpubr)
Pollen.plots <- ggarrange(MNIplot, NISPplot,
                    labels = c("MNI", "NISP", "NISP other elements"),
                    ncol = 1, nrow = 2)
Pollen.plots

famcolors <- c("#BDD9BF", "#FFC857", "#A997DF", "#929084", "#2E4052")
ggplot(MNIdf, aes(fill=taxa, y=per, x=yr)) + 
  geom_bar(position="fill", stat="identity") +theme_classic(base_size = 18) +
  scale_fill_manual(name = "Family", values=c(famcolors)) +
  ylab("Relative Abundace") + xlab("Years BP")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                                RAREFACTION                               ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chao et al. (2014) 

library(iNEXT)
library(ggplot2)

#   MNI
LIZMNI_wideT<-t(LIZMNI_wide) #transpose data

out<-iNEXT(LIZMNI_wideT, q=c(1), datatype="abundance") 

ggiNEXT(out, type=1, facet.var="site")

ggiNEXT(out, type=2)

Shannon_MNI_estimates<- ChaoShannon(LIZMNI_wideT) #Estimate Shannon div
Shannon_MNI_estimates <- tibble::rownames_to_column(Shannon_MNI_estimates, "Bin")
Shannon_MNI_estimates$Bin <- factor(Shannon_MNI_estimates$Bin, levels=unique(Shannon_MNI_estimates$Bin))


#   NISP
LIZNISP_wideT<-t(LIZNISP_wide) #transpose data

out2<-iNEXT(LIZNISP_wideT, q=c(1), datatype="abundance")

ggiNEXT(out2, type=1, facet.var="site")

Shannon_NISP_estimates <- ChaoShannon(LIZNISP_wideT) #Estimate Shannon div

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                           Shannon diversity plots                        ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HMNI <- data.frame(yr=rep(age,ncol(Shannon_MNI_estimates)), #make dataframe
                   Obs=as.vector(as.matrix(Shannon_MNI_estimates$Observed)),
                   Est=as.vector(as.matrix(Shannon_MNI_estimates$Estimator)),
                   L95=as.vector(as.matrix(Shannon_MNI_estimates$`95% Lower`)),
                   U95=as.vector(as.matrix(Shannon_MNI_estimates$`95% Upper`)))

HMNI_plot<-ggplot(HMNI, aes(yr, Est))+
  geom_line(colour="blue")+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,5,.5))+
  xlab("Age (cal. BP)")+ylab("H")+
  coord_flip() +
  theme_classic()+ geom_ribbon(aes(ymin = L95, ymax = U95), alpha = 0.2)+
  geom_line(aes(y = Obs), color = "black")


HNISP <- data.frame(yr=rep(age,ncol(Shannon_NISP_estimates)),#make dataframe
                    Obs=as.vector(as.matrix(Shannon_NISP_estimates$Observed)),
                    Est=as.vector(as.matrix(Shannon_NISP_estimates$Estimator)),
                    L95=as.vector(as.matrix(Shannon_NISP_estimates$`95% Lower`)),
                    U95=as.vector(as.matrix(Shannon_NISP_estimates$`95% Upper`)))

HNISP_plot<-ggplot(HNISP, aes(yr, Est))+
  geom_line(colour="blue")+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,5,.5))+
  xlab("Age (cal. BP)")+ylab("H")+
  coord_flip() +
  theme_classic() + geom_ribbon(aes(ymin = L95, ymax = U95), alpha = 0.2) +
  geom_line(aes(y = Obs), color = "black")

H.plots <- ggarrange(HMNI_plot, HNISP_plot,
                          labels = c("MNI", "NISP"),
                          ncol = 2, nrow = 1)
H.plots


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                           ELEMENT REPRESENTATION                         ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

levels(as.factor(Fossil_lizard_15bin$Element))
Fossil_lizard_15bin_element <- filter(Fossil_lizard_15bin,!grepl('tooth bearing element|jugal?|longbones|postorbitofrontal?|roofing bone',Element)) #filter out uncertain identifications
Fossil_lizard_15bin_element$Element <- ifelse(Fossil_lizard_15bin_element$Element == "ilium", # combine common elements
                        "pelvis", Fossil_lizard_15bin_element$Element)
Fossil_lizard_15bin_element$Element <- ifelse(Fossil_lizard_15bin_element$Element == "surangular", # combine common elements
                                              "surangular_articular", Fossil_lizard_15bin_element$Element)
Fossil_lizard_15bin_element$Element <- ifelse(Fossil_lizard_15bin_element$Element == "articular", # combine common elements
                                              "surangular_articular", Fossil_lizard_15bin_element$Element)
levels(as.factor(Fossil_lizard_15bin_element$Element))

Element_rep <- aggregate(data = Fossil_lizard_15bin_element,   # Applying aggregate to count # of unique elements per bin
                          Element ~ Bin,
                          function(x) length(unique(x)))

Element_rep <- data.frame(yr=rep(age,ncol(Element_rep)),#make dataframe
                    Bin=as.vector(as.matrix(Element_rep$Bin)),
                    Num=as.vector(as.matrix(Element_rep$Element)))

Element_rep_plot<-ggplot(Element_rep, aes(yr, Num))+ geom_point()+ #plot
  geom_line(colour="black")+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,100,5))+
  xlab("Age (cal. BP)")+ylab("# Unique Elements")+
  coord_flip() +
  theme_classic()
Element_rep_plot


p.plots <- ggarrange(LIZVIZ_plot, Element_rep_plot,
                     labels = c("NISP", "Elements"),
                     ncol = 2, nrow = 1)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                              ELEMENT BREAKAGE                            ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fossil_lizard_15bin_break <- filter(Fossil_lizard_15bin,!grepl('tooth bearing element|jugal?|longbones|postorbitofrontal?|roofing bone|atlas|axis|vertebra|metacarpals_metatarsals|osteoderm|phalanges|ribs',Element)) #filter out elements not scored for completeness
Fossil_lizard_15bin_break$Bin<- as.factor(Fossil_lizard_15bin_break$Bin)
Fossil_lizard_15bin_break$Element<- as.factor(Fossil_lizard_15bin_break$Element)

LizBreak <- Fossil_lizard_15bin_break %>% 
  filter(!is.na(Quadrant))  %>%
  group_by(Bin,Element, .drop=FALSE) %>% summarise(PctPresent = mean(Quadrant))

levels(Fossil_lizard_15bin_break$Element)

LizBreak <- LizBreak %>% 
  mutate(Position = recode(Element, "articular"="C"     ,       "braincase"="C" ,            "coronoid"  ="C" ,          
                          "dentary"="C" ,              "ectopterygoid" ="C" ,       "femur" ="P",               
                          "frontal" ="C" ,             "humerus" ="P",             "ilium"  ="P",             
                          "interclavicle" ="P",       "maxilla" ="C" ,             "nasal" ="C" ,              
                          "parietal"   ="C" ,          "pelvis" ="P",              "postorbital"="C" ,         
                          "prefrontal" ="C" ,          "premaxilla"="C" ,           "pterygoid" ="C" ,          
                          "quadrate" ="C" ,            "radius"   ="P",            "scapulocorocoid" ="P",    
                          "splenial" ="C" ,            "squamosal"="C" ,            "surangular" ="C" ,         
                          "surangular_articular"="C" , "tibia" ="P",               "ulna" ="P",               
                          "vomer"="C"  ))

BreakSum <- LizBreak %>% 
  filter(!is.na(PctPresent))  %>%
  group_by(Position) %>% summarise(PctCompl = mean(PctPresent))

BreakBinSum <- LizBreak %>% 
  filter(!is.na(PctPresent))  %>%
  group_by(Bin,Position) %>% summarise(PctCompl = mean(PctPresent))
BreakBinSum$Bin<-as.integer(BreakBinSum$Bin)

BreakBinSum_plot<-ggplot(BreakBinSum, aes(Bin, PctCompl))+
  geom_line(aes(color = Position))+
  scale_x_reverse(breaks =seq(0,20,1))+
  xlab("Bin")+ylab("Percent complete")+
  coord_flip() +
  theme_classic()
BreakBinSum_plot

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                                  CORROSION                               ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fossil_lizard_15bin$Bin<- as.factor(Fossil_lizard_15bin$Bin)
Fossil_lizard_15bin$TeethCorrosion<- as.factor(Fossil_lizard_15bin$TeethCorrosion)

Corrosion <- Fossil_lizard_15bin %>% 
  filter(!is.na(TeethCorrosion))  %>%
  group_by(Bin, TeethCorrosion, .drop=FALSE) %>% tally()
Corrosion$Bin<-as.integer(Corrosion$Bin)

Corrosion_wide <- spread(Corrosion, TeethCorrosion, n)
Corrosion_wide <- Corrosion_wide %>% remove_rownames %>% column_to_rownames(var="Bin")

Corrosion.pct <- data.frame(tran(Corrosion_wide, method = 'percent')) # CONVERT INTO PERCENTS

Corrosiondf <- data.frame(Bin=rep(1:20,ncol(Corrosion.pct)),
                     per=as.vector(as.matrix(Corrosion.pct)),
                     Type=as.factor(rep(colnames(Corrosion.pct),each=nrow(Corrosion.pct))))

Corrosion_plot<-ggplot(Corrosiondf, aes(Bin,per))+
  geom_line(aes(color = Type))+
  scale_x_reverse(breaks =seq(0,20,1))+
  scale_y_continuous(breaks =seq(0,100,10))+
  xlab("Bin")+ylab("%")+
  coord_flip()+
  theme_classic()
Corrosion_plot


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                        GENERALIZED ADDITIVE MODELS                       ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library("mgcv")
library("scam")
library("ggplot2")
library("cowplot")
library("tidyr")
#devtools::install_github("gavinsimpson/gratia")
library("gratia") # need to change the name of the package


## Default ggplot theme
theme_set(theme_bw())

small <- readRDS("small-water-isotope-data.rds") #load in dataset
head(small)














