library(ARTofR)
# READ IN THE DATA------------------------------------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/Halls-Cave-Lizards/main/Fossil_lizard_15bin.csv")
Fossil_lizard_15bin <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = FALSE, na.strings=c("","NA")) # this is a matrix of measured specimens 
head(Fossil_lizard_15bin)

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

# SUMMARIZE INTO MNI DATA
LIZNISP <- Fossil_lizard_15bin_fam %>% 
  filter(!is.na(Details))  %>%
  filter(!is.na(Family)) %>% 
  group_by(Family,Bin, .drop=FALSE) %>% tally
LIZNISP <- LIZNISP %>% group_by(Family, Bin) %>% summarise(NISP = max(n))

# MAKE INTO LONG FORMAT WITH FAMILIES AS COLUMNS AND BINS AS ROWS
library(tidyverse)
library(tidyr)
LIZNISP_wide <- spread(LIZNISP, Family, NISP)
LIZNISP_wide <- LIZNISP_wide %>% remove_rownames %>% column_to_rownames(var="Bin")

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

age_depth <- function(Depth) { # Function to assign ages to depths based on linear relationship defined by Tomé et al. 2021
  68.61*(Depth) + 758.82
}
bins <- seq(from = 15, to = 300, by = 15) #VECTOR OF 15 CM BINS
age <- sapply(bins, age_depth) # ASSIGN AGES TO ALL BINS

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

library(ggpubr)
Pollen.plots <- ggarrange(MNIplot, NISPplot,
                    labels = c("MNI", "NISP"),
                    ncol = 1, nrow = 2)
Pollen.plots


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

ggiNEXT(LIZMNI_wideT, type=1, facet.var="site")

ggiNEXT(LIZMNI_wideT, type=2)

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
Fossil_lizard_15bin_element$Element <- ifelse(Fossil_lizard_15bin_element$Element == "surangular",
                                              "surangular_articular", Fossil_lizard_15bin_element$Element)
Fossil_lizard_15bin_element$Element <- ifelse(Fossil_lizard_15bin_element$Element == "articular",
                                              "surangular_articular", Fossil_lizard_15bin_element$Element)
levels(as.factor(Fossil_lizard_15bin_element$Element))

Element_rep <- aggregate(data = Fossil_lizard_15bin_element,   # Applying aggregate to count # of unique elements per bin
                          Element ~ Bin,
                          function(x) length(unique(x)))

Element_rep <- data.frame(yr=rep(age,ncol(Element_rep)),#make dataframe
                    Bin=as.vector(as.matrix(Element_rep$Bin)),
                    Num=as.vector(as.matrix(Element_rep$Element)))

Element_rep_plot<-ggplot(Element_rep, aes(yr, Num))+
  geom_line(colour="black")+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,100,5))+
  xlab("Age (cal. BP)")+ylab("# Unique Elements")+
  coord_flip() +
  theme_classic()
Element_rep_plot












