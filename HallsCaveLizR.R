library(ARTofR)
# READ IN THE DATA------------------------------------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/Halls-Cave-Lizards/main/Fossil_lizard_15bin.csv")
Fossil_lizard_15bin <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = FALSE, na.strings=c("","NA")) # this is a matrix of measured specimens 
head(Fossil_lizard_15bin)

f3 <- curl("https://raw.githubusercontent.com/TIMAVID/Halls-Cave-Lizards/main/Ages.csv")
Ages <- read.csv(f3, header = TRUE, sep = ",", stringsAsFactors = FALSE, na.strings=c("","NA")) # this is a matrix of measured specimens 
head(Ages)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Visualize data                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(ggplot2)

# age_depth <- function(Depth) { # Function to assign ages to depths based on linear relationship defined by Tomé et al. 2021
#   68.61*(Depth) + 758.82
# }
# bins <- seq(from = 0, to = 295, by = 15) #VECTOR OF 15 CM BINS
# age <- sapply(bins, age_depth) # ASSIGN AGES TO ALL BINS


LIZVIZ <- Fossil_lizard_15bin %>% group_by(Bin) %>% summarise(NISP = sum(NISP))
LIZVIZ <- LIZVIZ[-c(21), ]
LIZVIZ$age <- Ages$MEDAge
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
LIZNISP2$age <- Ages$MEDAge
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

Stratiplot(Ages$MEDAge ~ ., LIZMNI.pct, sort = 'wa', type = 'poly',
           ylab ='Years Before Present')

MNIdf <- data.frame(yr=rep(Ages$MEDAge,ncol(LIZMNI.pct)),
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

Stratiplot(Ages$MEDAge ~ ., LIZNISP.pct, sort = 'wa', type = 'poly',
           ylab ='Years Before Present')

NISPdf <- data.frame(yr=rep(Ages$MEDAge,ncol(LIZNISP.pct)),
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

NISPdf_other <- data.frame(yr=rep(Ages$MEDAge,ncol(LIZNISP_OTHER.pct)),
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
  geom_bar(position="fill", stat="identity", width=500) +theme_classic(base_size = 18) +
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

HMNI <- data.frame(yr=(Ages$MEDAge), #make dataframe
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


HNISP <- data.frame(yr=(Ages$MEDAge),#make dataframe
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

 #load in dataset
head(HNISP)

### fit model using GCV -------
HNISP.gcv <- gam(Obs ~ s(yr, k = 7), data = HNISP)

N <- 500 # number of points at which to evaluate the smooth
## data to predict at
newHNISP <- with(HNISP, data.frame(yr = seq(min(yr), max(yr),
                                              length.out = N)))

## add GAM GCV results
fit_gcv <- predict(HNISP.gcv, newdata = newHNISP, se.fit = TRUE)
crit.t <- qt(0.975, df.residual(HNISP.gcv))
newGCV <- data.frame(yr = newHNISP[["yr"]],
                     fit = fit_gcv$fit,
                     se.fit = fit_gcv$se.fit)
newGCV <- transform(newGCV,
                    upper = fit + (crit.t * se.fit),
                    lower = fit - (crit.t * se.fit))


## plot GCV fits
braya_fitted <- ggplot(HNISP, aes(y = Obs, x = yr)) +
  geom_point() +
  geom_ribbon(data = newGCV,
              mapping = aes(x = yr, ymax = upper, ymin = lower),
              alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = newGCV,
            mapping = aes(y = fit, x = yr)) +
  labs(y = "H (Shannon div.)", x = "Cal Year BP") +
  scale_color_manual(values = c("#5e3c99", "#e66101")) +
  scale_fill_manual(values = c("#5e3c99", "#e66101")) +
  theme(legend.position = "right")
braya_fitted

# Checking if the size of the basis expansion is sufficient
HNISP_low_k <- gam(Obs ~ s(yr, k = 7), data = HNISP, method = "REML")

gam.check(HNISP_low_k)

### Accounting for heteroscedasticity due to time averaging -------

TimeAvgHNISP <- cbind(HNISP, Ages)


HNISP_reml <- gam(Obs ~ s(yr, k = 10), data = TimeAvgHNISP,
                  method = "REML",
                  weights = TimeInterval / mean(TimeInterval))
summary(HNISP_reml)

gam.check(HNISP_reml)


### Posterior simulation ------

set.seed(1) # set the random seed to make this reproducible
nsim <- 20 # how many simulations to draw

## data points to simulate at
newHNISP <- with(HNISP,
                 data.frame(yr = seq(min(yr), max(yr),
                                       length.out = N)))
HNISP_pred <- cbind(newHNISP,
                    data.frame(predict(HNISP_reml, newHNISP,
                                       se.fit = TRUE)))
## simulate
set.seed(1)
sims2 <- simulate(HNISP_reml, nsim = nsim, data = newHNISP,
                  unconditional = TRUE)
## rearrange the output into a long/tidy format
colnames(sims2) <- paste0("sim", seq_len(nsim))

sims2 <- setNames(stack(as.data.frame(sims2)),
                  c("simulated", "run"))
sims2 <- transform(sims2, yr = rep(newHNISP$yr, nsim),
                   simulated = simulated)
HNISPSim.plt <- ggplot(HNISP_pred, aes(x = yr, y = fit)) +
  geom_line(data = sims2,
            mapping = aes(y = simulated, x = yr, group = run),
            colour = "grey80") +
  geom_line(lwd = 1) +
  labs(y = "H (Shannon div.)", x = "Year CE")
HNISPSim.plt



HNISP.cint <- confint(HNISP_reml, parm = "yr", newdata = newHNISP,
                   type = "confidence", partial_match = TRUE)
# HNISP.sint <- confint(HNISP_reml, parm = "yr", newdata = newHNISP,
#                    type = "simultaneous", partial_match = TRUE) #not working for some reason


HNISPInt.plt <- ggplot(HNISP.cint, aes(x = yr, y = est)) +
  # geom_ribbon(data = bs.sint,
  #             mapping = aes(ymin = lower, ymax = upper, x = Year),
  #             fill = "grey80", inherit.aes = FALSE) +
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, x = yr),
              fill = "grey60", inherit.aes = FALSE) +
  geom_line(lwd = 1) +
  labs(y = "H (Shannon div.)", x = "Year CE")
HNISPInt.plt


### First derivative ----------

HNISP.d <- derivatives(HNISP_reml)

# HNISP.d <- fderiv(HNISP_reml, data = newHNISP, n = N)
# HNISP.sint <- with(newHNISP,
#                    cbind(confint(HNISP.d, nsim = nsim,
#                                  type = "simultaneous"),
#                          Year = yr))

HNISP_deriv_plt <- ggplot(HNISP.d, aes(x = data, y = derivative)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, fill = "black") +
  geom_line() +
  labs(x = "Cal. Year BP", y = "First derivative")
HNISP_deriv_plt






### Adaptive spline, weights as sampleInterval-----
mod_ad <- gam(Obs ~ s(yr, k = 10, bs = "ad"), data = TimeAvgHNISP,
              method = "REML",
              weights = TimeInterval / mean(TimeInterval))


## wrap this in a function that will return all the plots & derived objects
processGAM <- function(mod) {
  ## Predict from model
  N <- 500
  newYear <- with(TimeAvgHNISP,
                  data.frame(yr = seq(min(yr), max(yr),
                                        length.out = N)))
  newYear <- cbind(newYear,
                   data.frame(predict(mod, newYear, se.fit = TRUE)))
  out <- list(objects = newYear)
  out
}
plts_ad <- processGAM(mod = mod_ad) # Adaptive smooth with weights


pltData <- do.call("rbind", lapply(list(plts_ad),
                                   `[[`, "objects"))

pltData <- transform(pltData,
                     Model = rep(c("Adaptive"),
                                 each = nrow(plts_ad$objects)))

allFits <- ggplot(pltData, aes(x = yr, y = fit)) +
  geom_point(aes(x = yr, y = Obs), data = TimeAvgHNISP) +
  geom_line(aes(colour = Model)) + labs(y = "braya_ylabel", x = "Year") +
  theme(legend.position = "right") +
  scale_colour_manual(name = "",
                      values = c("#e66101", "#fdb863", "#5e3c99"))
allFits


### Accounting for age-model uncertainty------------

knots <- with(TimeAvgHNISP, list(Year = seq(min(yr), max(yr), length = 14))) #fix the knots at the extremes of the observed values of Year and spread the remaining 12 knots evenly inbetween
HNISP_reml <- gam(Obs ~ s(yr, k = 10), data = TimeAvgHNISP,
                  method = "REML",
                  weights = TimeInterval / mean(TimeInterval), knots = knots)

swAge <- read.csv("Radiocarbon_Dates_Cal.csv")
## monotonic spline age-depth model
#swAge$Error[1] <- 1.1
swAgeMod <- scam(mu ~ s(Depth, k = 7, bs = "mpi"), data = swAge,
                 weights = 1 / (swAge$sigma*2), gamma = 1.4)

## predict from the age model for a smooth set of points in `Depth`
newAge <- with(swAge, data.frame(Depth = seq(min(Depth), max(Depth),
                                             length.out = 200)))
newAge <- transform(newAge,
                    fitted = predict(swAgeMod, newdata = newAge,
                                     type = "response"))

newSims <- as.data.frame(simulate(swAgeMod, nsim = 25, data = newAge))
newSims <- cbind(Depth = newAge$Depth, newSims)
newSims <- gather(newSims, Simulation, Age, -Depth)


## simulate from age model; each column is a simulation
ageSims <- simulate(swAgeMod, nsim = 100, data = TimeAvgHNISP, seed = 42)

ageSims <- as.data.frame(ageSims)

fitSWModels <- function(x, y, knots, w) {
  dat <- data.frame(Obs = y, yr = x)
  m <- gam(Obs ~ s(yr, k = 10), data = dat,
           method = "REML",
           weights = w / mean(w), knots = knots)
}

## generate new trends using draws from age-model posterior
simTrendMods <- lapply(ageSims, fitSWModels, y = TimeAvgHNISP$Obs, knots = knots, w= TimeAvgHNISP$TimeInterval)


simTrendModsDeriv <- lapply(simTrendMods, fderiv, n = N)

## function wrapper to predict new trends at locations over the
## range of `Year`
predSWModels <- function(mod, newdata) {
  predict(mod, newdata, type = "response")
}

## predict from fitted model to produce a smooth trend for each posterior
## sample
simTrends <- lapply(simTrendMods, predSWModels, newdata = newGCV)
## arrange in a tidy format form plottings
simTrends <- data.frame(Year = with(newGCV, rep(yr, length(simTrends))),
                        Trend = unlist(simTrends),
                        Group = rep(seq_along(simTrends),
                                    times = lengths(simTrends)))

# For each of the models we just fitted, we
# simulate 50 draws from the model posterior distribution. We start with a wrapper function
# around the simulate() code we want to run on each model, then do the actual posterior draws
# for each model using lapply(). The final step just arranges data for plotting.


## wrapper to simulate from a fitted GAM with the required arguments
simulateSWModels <- function(mod, newdata, nsim, seed = 42) {
  sims <- simulate(mod, nsim = nsim, data = newdata, seed = seed)
  as.vector(sims)
}
## now do the posterior simulation
NSIM <- 50 # number of posterior samples *per* model
simSimulate <- lapply(simTrendMods, simulateSWModels, newdata = newGCV,
                      nsim = NSIM, seed = 42)

## arrange in a tidy format
simSimulate <-
  data.frame(yr = with(newGCV,
                         rep(yr, times = NSIM * length(simSimulate))),
             Trend = unlist(simSimulate),
             Group = rep(seq_len(NSIM * length(simSimulate)),
                         each = nrow(newGCV)))


## plot the estimated age model ad posterior simulations from it
plt1 <- ggplot(swAge, aes(y = mu, x = Depth)) +
  geom_line(data = newSims,
            mapping = aes(y = Age, x = Depth, group = Simulation),
            alpha = 1, colour = "grey80") +
  geom_line(data = newAge, mapping = aes(y = fitted, x = Depth)) +
  geom_point(size = 1.5, colour = "red") +
  geom_errorbar(aes(ymin =  mu - (sigma*2), ymax = mu + (sigma*2), width = 0),
                colour = "red") +
  labs(y = "Cal. Year BP", x = "Depth")

## plot the simulated trends showing the effect of age-model uncertainty
plt2 <- ggplot(simTrends, aes(x = Year, y = Trend, group = Group)) +
  geom_line(alpha = 0.1, colour = "grey80") +
  geom_line(data = newGCV,
            mapping = aes(x = yr, y = fit), inherit.aes = FALSE) +
  geom_point(data = TimeAvgHNISP,
             mapping = aes(x = yr, y = Obs),
             inherit.aes = FALSE, size = 0.7) +
  labs(x = "Year", y = "H (Shannon div.)")

## plot simulated trends showing the effect of age-model uncertainty and
## the effect of uncertainty in the estimated trend itself
plt3 <- ggplot(simSimulate, aes(x = yr, y = Trend, group = Group)) +
  geom_line(alpha = 0.2, colour = "grey80") +
  geom_point(data = TimeAvgHNISP,
             mapping = aes(x = yr, y = Obs),
             inherit.aes = FALSE,
             size = 0.7) +
  geom_line(data = newGCV,
            mapping = aes(x = yr, y = fit),
            inherit.aes = FALSE) +
  labs(x = "Year", y = "H (Shannon div.)")
## align all plots vertically
plots <- align_plots(plt1, plt2, plt3, align = 'v', axis = 'l')

## create the two rows of figures, from `plots`
top_row <- plot_grid(plots[[1]], NULL, ncol = 2, labels = "a")
bot_row <- plot_grid(plots[[2]], plots[[3]], ncol = 1, labels = c("b", "c"))
## comDepthe the two rows, top row has 1 plot row, bottom row has 2, hence
## the rel_heights to even this out
plot_grid(top_row, bot_row, ncol = 1, rel_heights = c(0.5, 1))




#derivative with age uncertainties------
plot(simSimulate$yr, simSimulate$Trend)

system.time(test.gcv <- gam(Trend ~ s(Year, k = 40), data = simSimulate))
gam.check(test.gcv)

test.d <- derivatives(test.gcv)

test_deriv_plt <- ggplot(test.d, aes(x = data, y = derivative)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, fill = "black") +
  geom_line() +
  labs(x = "Cal. Year BP", y = "First derivative")
test_deriv_plt





library(dplyr)
df <-simSimulate %>% group_by(Group) %>% filter(n() > 1) %>%
  mutate(first_d = Trend - lag(Trend))
df<-as.data.frame(df)  #Back to data frame

plt3 <- ggplot(df, aes(x = yr, y = first_d, group = Group)) +
  geom_line(alpha = 0.2, colour = "grey80") +
  # geom_point(data = TimeAvgHNISP,
  #            mapping = aes(x = yr, y = Obs),
  #            inherit.aes = FALSE,
  #            size = 0.7) +
  # geom_line(data = newGCV,
  #           mapping = aes(x = yr, y = fit),
  #           inherit.aes = FALSE) +
  labs(x = "Year", y = "H (Shannon div.)")



