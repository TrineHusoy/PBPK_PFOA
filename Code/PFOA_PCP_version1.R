####################################################################################################################################################################################
# In this project we want to estimate the PFOA exposure from PCPs  , using simmilar approach as Karrer et al 2020 for PCP
# Date 16.10.20
#######################################################################################

HOME <- "Your work directory"

setwd(HOME)

newday <- file.path('Your work directory/Results', Sys.Date())
dir.create(newday)


# Load required packages and libraries

library(ggplot2)
library(writexl)
library(triangle) # selected distribution
library(ggpubr)
library(Hmisc)
library(openxlsx)
library(reshape2)
library(trapezoid)
library(truncnorm)
library(tidyverse)


####  Read data #################################################

ProductUsePCP <- read_csv2("./Data/PCP_frequency_dummy.csv")

AmountsFemale <- read_csv("./Data/AmountsPerApplication_female.csv")
AmountsMale <- read_csv("./Data/AmountsPerApplication_male.csv")

EEFsFemale <-read_csv("./Data/EEF_distributions_female.csv")
EEFsMale <- read_csv("./Data/EEF_distributions_male.csv")

SexWeight <- read_csv2("./Data/EuroMix_dummy_sex_weight.csv")

PFOAconcPCP_LB<- read_csv2("./Data/2_SumPCPPFAS_LB_050122.csv")
PFOAconcPCP_MB<- read_csv2("./Data/2_SumPCPPFAS_MB_050122.csv")
PFOAconcPCP_UB<- read_csv2("./Data/2_SumPCPPFAS_UB_050122.csv")

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Clean data
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# merge PCP use with sex and weight ##

ProductUsePCP <- merge(SexWeight, ProductUsePCP, by = "IDkode", all = TRUE)


## Extract PCP use for males and females separately

ProductUsePCP_male <- ProductUsePCP[ProductUsePCP$cat_male == 2,]
ProductUsePCP_female <- ProductUsePCP[ProductUsePCP$cat_male == 1,]

# replace NAs by zero
ProductUsePCP_male[is.na(ProductUsePCP_male)] <- 0
ProductUsePCP_female[is.na(ProductUsePCP_female)] <- 0
PFOAconcPCP_LB[is.na(PFOAconcPCP_LB)]<- 0
PFOAconcPCP_MB[is.na(PFOAconcPCP_MB)]<- 0
PFOAconcPCP_UB[is.na(PFOAconcPCP_UB)]<- 0

# chengge from chaaracers to numeric
PFOAconcPCP_LB <- PFOAconcPCP_LB[,2:8] %>%  lapply(function(x) as.numeric(x))
PFOAconcPCP_LB <- as_tibble(PFOAconcPCP_LB)
PFOAconcPCP_MB <- PFOAconcPCP_MB[,2:8] %>%  lapply(function(x) as.numeric(x))
PFOAconcPCP_MB <- as_tibble(PFOAconcPCP_MB)
PFOAconcPCP_UB <- PFOAconcPCP_UB[,2:8] %>%  lapply(function(x) as.numeric(x))
PFOAconcPCP_UB <- as_tibble(PFOAconcPCP_UB)

## delete the categories not needed any longer ##

EEFsMale <- EEFsMale[-c(5),]
EEFsFemale <- EEFsFemale[-c(5),]

# also those products men don't use

ProductUsePCP_male <- ProductUsePCP_male %>% dplyr::select(-cat_male, -AntiWrincleCream, -Sunscreen, -LipGlossBalm, -Foundation, -IntimateSoap,
                                                           -HairTreatProd, -EyeMakeup, -RougePowder, -MakeupRemov, -Oils, -FootCream)

ProductUsePCP_female <- ProductUsePCP_female %>% dplyr::select(-cat_male)


### remove PCP conc for cat that is not used by males ###
# change from ng/ml to ug/ml



PFOAconcPCP_LB_male <- PFOAconcPCP_LB/1000 # from ng til ug
PFOAconcPCP_MB_male <- PFOAconcPCP_MB/1000
PFOAconcPCP_UB_male <- PFOAconcPCP_UB/1000
PFOAconcPCP_LB_male <- PFOAconcPCP_LB_male[-c(7,8,12,13,14,16,18,19,20,21,24),]
PFOAconcPCP_MB_male <- PFOAconcPCP_MB_male[-c(7,8,12,13,14,16,18,19,20,21,24),]
PFOAconcPCP_UB_male <- PFOAconcPCP_UB_male[-c(7,8,12,13,14,16,18,19,20,21,24),]

PFOAconcPCP_LB_female <- PFOAconcPCP_LB
PFOAconcPCP_MB_female <- PFOAconcPCP_MB
PFOAconcPCP_UB_female <- PFOAconcPCP_UB

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Functions
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Sum <-function(x){
  x %>%
    group_by(IDkode) %>%
    summarise(
      N = n(),
      mean = mean(value, na.rm=TRUE),
      sd=sd(value, na.rm=TRUE),
      min=min(value, na.rm=TRUE),
      P05=quantile(value, .05, na.rm=TRUE),
      P50=quantile(value, .50, na.rm=TRUE),
      P95=quantile(value, .95, na.rm=TRUE),
      max=max(value, na.rm=TRUE)
    )
}  


# function to calculate the amount (g) of each PCP for each ID from the frequency
pcp_amount <- function(x,y,z){
  set.seed(123)
  MC <- 1000
  for (u in 1:MC) {
    for (i in 1:nrow(x)) {
      y[,(i+2),u] <- z[,(i+2)] * rtruncnorm (nrow(z), 
                                             a=x$LB[i],b=x$UB[i], 
                                             mean = x$mean[i], 
                                            sd = x$SD[i])
    }
  }
  return(y)
}


# function to calculate the amount of chemical (ng) from amount PCP (g) used
compound_amount <- function(x,y,z){
  set.seed(123)
  MC <- 1000
  for (u in 1:MC) {
    for (i in 1:nrow(x)) {
      y[,(i+2),u] <- z[,(i+2),u] * rtriangle(nrow(z),
                                             a=x$PFOA_P05[i], 
                                             b=x$PFOA_P95[i],
                                             c=x$PFOA_P50[i])
    }
  } 
  return(y)
}


# function to calculate the remaining amount chemical left after including 
# leav-on/rins-off factor
amount_rinse <- function(x,y,z){
  set.seed(123)
  MC <- 1000
  for (u in 1:MC) {
    for (i in 1:nrow(x)) {
      y[,(i+2),u] <- z[,(i+2),u] * rtriangle(nrow(z), 
                                             a=x$LB_D[i], 
                                             b=x$UB_D[i],
                                             c=x$mean_D[i])
    }
  } 
  return(y)
}

# function to calculate the total chemical per MC per ID from all PCPs
 sum_amount_PCP <- function(x,y) {
   for (i in 1:1000){
     x [,i] <- rowSums(y[,3:16,i]) # Need to spsify the coolumns to sumarise to exclude adding the BW and ID
   }
   return(x)
 }
  


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### Calculation for males ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MC <- 1000
# Amounts in gram per day per person for each PCP product 

use_amount_m <- array(unlist(ProductUsePCP_male),dim = c(44,16,MC))

# x= AmountsMale, y=use_amount_m, z = ProductUsePCP_male

use_amount_m <- pcp_amount(AmountsMale, use_amount_m, ProductUsePCP_male)

print(use_amount_m[2,,5])
print(use_amount_m[10,,5])

####################PCP EXPOSURE MALE LB ####################################


# Amount of PFOA (total gram used * ng/g conc = ng)


amount_PFOA_m_LB <- use_amount_m
amount_PFOA_m_LB <- compound_amount(PFOAconcPCP_LB_male, amount_PFOA_m_LB, use_amount_m)

amount_PFOA_m_MB <- use_amount_m
amount_PFOA_m_MB <- compound_amount(PFOAconcPCP_MB_male, amount_PFOA_m_MB, use_amount_m)

amount_PFOA_m_UB <- use_amount_m
amount_PFOA_m_UB <- compound_amount(PFOAconcPCP_UB_male, amount_PFOA_m_UB, use_amount_m)


print(amount_PFOA_m_LB[3,,5])
print(amount_PFOA_m_LB[24,,5])

print(amount_PFOA_m_MB[3,,5])
print(amount_PFOA_m_MB[24,,5])

print(amount_PFOA_m_UB[3,,5])
print(amount_PFOA_m_UB[24,,5])

################################### Adding rinse factors ###########


amount_PFOA_m_LB_eef <- amount_PFOA_m_LB
amount_PFOA_m_LB_eef <- amount_rinse(EEFsMale,amount_PFOA_m_LB_eef, amount_PFOA_m_LB )

amount_PFOA_m_MB_eef <- amount_PFOA_m_MB
amount_PFOA_m_MB_eef <- amount_rinse(EEFsMale,amount_PFOA_m_MB_eef, amount_PFOA_m_MB )

amount_PFOA_m_UB_eef <- amount_PFOA_m_UB
amount_PFOA_m_UB_eef <- amount_rinse(EEFsMale,amount_PFOA_m_UB_eef, amount_PFOA_m_UB )


print(amount_PFOA_m_LB_eef[3,,5]) # still ug
print(amount_PFOA_m_LB_eef[24,,5])

print(amount_PFOA_m_MB_eef[3,,5]) # still ug
print(amount_PFOA_m_MB_eef[24,,5])

print(amount_PFOA_m_UB_eef[3,,5]) # still ug
print(amount_PFOA_m_UB_eef[24,,5])



## sumarise the external exposure over all PCP products for all MC #####

total_m_LB <- as.data.frame(matrix(NA,nrow = 44,ncol = 1000))
total_m_LB <- sum_amount_PCP(total_m_LB, amount_PFOA_m_LB_eef)

print(total_m_LB [3,5]) # still 

total_m_MB <- as.data.frame(matrix(NA,nrow = 44,ncol = 1000))
total_m_MB <- sum_amount_PCP(total_m_MB, amount_PFOA_m_MB_eef)

print(total_m_MB [25,5]) # still ng

total_m_UB <- as.data.frame(matrix(NA,nrow = 44,ncol = 1000))
total_m_UB <- sum_amount_PCP(total_m_UB, amount_PFOA_m_UB_eef)

print(total_m_UB [25,5]) # still ng



###### calculate the exposure in ug/day ########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

total_m_ID_LB <- as.data.frame(matrix(NA, nrow = 44, ncol = 1001))
total_m_ID_LB[,1] <- ProductUsePCP_male[,1] #get ID 
total_m_ID_LB[,2:1001] <- total_m_LB[,1:1000] # in ng/day

names(total_m_ID_LB)[1] <- "IDkode"

total_m_ID_MB <- as.data.frame(matrix(NA, nrow = 44, ncol = 1001))
total_m_ID_MB[,1] <- ProductUsePCP_male[,1] #get ID 
total_m_ID_MB[,2:1001] <- total_m_MB[,1:1000] # in ng/day

names(total_m_ID_MB)[1] <- "IDkode"

total_m_ID_UB <- as.data.frame(matrix(NA, nrow = 44, ncol = 1001))
total_m_ID_UB[,1] <- ProductUsePCP_male[,1] #get ID 
total_m_ID_UB[,2:1001] <- total_m_UB[,1:1000] # in ng/day

names(total_m_ID_UB)[1] <- "IDkode"


### make summary data for input to the PBPK model ug/day ##########################################

# LB
total_m_ID_LB_long <- total_m_ID_LB %>% pivot_longer(
  !IDkode, names_to = "MC", values_to = "value")  # make a long table

Sumtotal_m_ID_LB <-Sum(total_m_ID_LB_long) # Make summary data for males LB

head(Sumtotal_m_ID_LB)


write.xlsx(Sumtotal_m_ID_LB,
           file = file.path(newday,"Sumtotal_m_ID_LB.xlsx"),
           colNames = TRUE, borders = "rows"
)

# MB
total_m_ID_MB_long <- total_m_ID_MB %>% pivot_longer(
  !IDkode, names_to = "MC", values_to = "value")  # make a long table

Sumtotal_m_ID_MB <-Sum(total_m_ID_MB_long) # Make summary data for males MB

head(Sumtotal_m_ID_MB)


write.xlsx(Sumtotal_m_ID_MB,
           file = file.path(newday,"Sumtotal_m_ID_MB.xlsx"),
           colNames = TRUE, borders = "rows"
)

# UB
total_m_ID_UB_long <- total_m_ID_UB %>% pivot_longer(
  !IDkode, names_to = "MC", values_to = "value")  # make a long table

Sumtotal_m_ID_UB <-Sum(total_m_ID_UB_long) # Make summary data for males UB

head(Sumtotal_m_ID_UB)


write.xlsx(Sumtotal_m_ID_UB,
           file = file.path(newday,"Sumtotal_m_ID_UB.xlsx"),
           colNames = TRUE, borders = "rows"
)
### make summary data as ug/kg bw/day ##########################################

# LB
total_m_ID_LB_kg <- total_m_ID_LB
total_m_ID_LB_kg$BW <- ProductUsePCP_male[,2]
total_m_ID_LB_kg[,2:1001] <- total_m_ID_LB_kg[,2:1001]/total_m_ID_LB_kg$BW


total_m_ID_LB_kg  <- total_m_ID_LB_kg %>% dplyr::select(-BW) # remove BW
total_m_ID_LB_kg_long <- total_m_ID_LB_kg %>%  pivot_longer(
  !IDkode, names_to = "MC", values_to = "value") # make a long table

Sumtotal_m_ID_LB_kg <- Sum(total_m_ID_LB_kg_long)

head(Sumtotal_m_ID_LB_kg)

write.xlsx(Sumtotal_m_ID_LB_kg,
           file = file.path(newday,"Sumtotal_m_ID_LB_kg.xlsx"),
           colNames = TRUE, borders = "rows"
)

# MB
total_m_ID_MB_kg <- total_m_ID_MB
total_m_ID_MB_kg$BW <- ProductUsePCP_male[,2]
total_m_ID_MB_kg[,2:1001] <- total_m_ID_MB_kg[,2:1001]/total_m_ID_MB_kg$BW


total_m_ID_MB_kg  <- total_m_ID_MB_kg %>% dplyr::select(-BW) # remove BW
total_m_ID_MB_kg_long <- total_m_ID_MB_kg %>%  pivot_longer(
  !IDkode, names_to = "MC", values_to = "value") # make a long table

Sumtotal_m_ID_MB_kg <- Sum(total_m_ID_MB_kg_long)

head(Sumtotal_m_ID_MB_kg)

write.xlsx(Sumtotal_m_ID_MB_kg,
           file = file.path(newday,"Sumtotal_m_ID_MB_kg.xlsx"),
           colNames = TRUE, borders = "rows"
)

# UB
total_m_ID_UB_kg <- total_m_ID_UB
total_m_ID_UB_kg$BW <- ProductUsePCP_male[,2]
total_m_ID_UB_kg[,2:1001] <- total_m_ID_UB_kg[,2:1001]/total_m_ID_UB_kg$BW


total_m_ID_UB_kg  <- total_m_ID_UB_kg %>% dplyr::select(-BW) # remove BW
total_m_ID_UB_kg_long <- total_m_ID_UB_kg %>%  pivot_longer(
  !IDkode, names_to = "MC", values_to = "value") # make a long table

Sumtotal_m_ID_UB_kg <- Sum(total_m_ID_UB_kg_long)

head(Sumtotal_m_ID_UB_kg)

write.xlsx(Sumtotal_m_ID_UB_kg,
           file = file.path(newday,"Sumtotal_m_ID_UB_kg.xlsx"),
           colNames = TRUE, borders = "rows"
)


## making vectors of the melted MC's as basis for plotting and summary ##

# LB
total_m_ID_LB_vector <- as.vector(total_m_ID_LB_long$value) # ug/day will be included in the main figure
total_m_ID_LB_kg_vector <- as.vector(total_m_ID_LB_kg_long$value) 

plot(density(total_m_ID_LB_vector))

# MB
total_m_ID_MB_vector <- as.vector(total_m_ID_MB_long$value) # ug/day will be included in the main figure
total_m_ID_MB_kg_vector <- as.vector(total_m_ID_MB_kg_long$value) 

plot(density(total_m_ID_MB_vector))

# UB
total_m_ID_UB_vector <- as.vector(total_m_ID_UB_long$value) # ug/day will be included in the main figure
total_m_ID_UB_kg_vector <- as.vector(total_m_ID_UB_kg_long$value) 

plot(density(total_m_ID_UB_vector))

## make summary data for all males based on all MC interations  #######

# LB
describe(total_m_ID_LB_vector) #(ug/day)
describe(total_m_ID_LB_kg_vector) 

# MB
describe(total_m_ID_MB_vector) #(ug/day)
describe(total_m_ID_MB_kg_vector) 

# UB
describe(total_m_ID_UB_vector) #(ug/day)
describe(total_m_ID_UB_kg_vector) 


### PLOTT EXPOSURE TO PFOA FROM PCP FOR MALES LB, MB and UB ####

## make a data frame of the three vectors with summary iterations for LB, MB and UB (ug/day) ##

A = c(total_m_ID_LB_vector, total_m_ID_MB_vector,total_m_ID_UB_vector)
B = replicate(44000,"LB")
C = replicate(44000,"MB")
D = replicate(44000,"UB")
E = c(B,C,D)

pfoa_m_pcp <- as.data.frame(matrix(NA,nrow = 132000, ncol = 2))
pfoa_m_pcp[1] <- A
pfoa_m_pcp[2] <- E

names(pfoa_m_pcp)[1] <- "PFOAconc"
names(pfoa_m_pcp)[2] <- "Exposure"

## Make a cumulative density plot with the LB, MB and UB exposure ##

Plot_pfoa_m_pcp <- ggplot(data = pfoa_m_pcp, aes(PFOAconc, colour = Exposure)) +
  geom_line(stat = "ecdf", size = 0.5)+
  scale_colour_hue()+
  theme_minimal()+
  scale_x_continuous(limits=c(0.000001,0.05))+
  xlab("PFOA (ng/day)")+
  ylab("Cumulative probability")+
  annotate("text", x=0.1, y=0.98, label="Males", size=6)

Plot_pfoa_m_pcp


ggsave(file = file.path(newday,"Plot_pfoa_m_pcp.jpeg"),
       device = NULL,
       width=NA,
       height=NA,
       units="mm")

##### Calculate the contribution from the different PCP groups for the LB #########
###################################################################################

# try to melt the array across the MC simulatons, and then plot each PCP category with all MC for all ID

dimnames(amount_PFOA_m_LB_eef)[[2]] <- colnames(ProductUsePCP_male) # give names to the second dimension PCP categories

provideDimnames(amount_PFOA_m_LB_eef) # check if I succeeded 


# melting the array to a two dimentional matrix with all the MC simulation per food categiry

pfoa_cat_m_LB<- melt(amount_PFOA_m_LB_eef, varnames = names(dimnames(amount_PFOA_m_LB_eef)),
                           na.rm = FALSE, value.name = "value")


names(pfoa_cat_m_LB)[2] <- "PCP" # give the second coloumn the name "Foods"
pfoa_cat_m_LB <- pfoa_cat_m_LB %>% dplyr::select(-Var1)# remove Var1, the IDs. Not needed for plotting the PFOA from each food cat
pfoa_cat_m_LB <- pfoa_cat_m_LB %>% dplyr::select(-Var3) # remove the number of the MC simulations
pfoa_cat_m_LB <- pfoa_cat_m_LB[!grepl("IDkode", pfoa_cat_m_LB$PCP),] # delee all rows with IDkode
pfoa_cat_m_LB <- pfoa_cat_m_LB[!grepl("con_weight", pfoa_cat_m_LB$PCP),] # delee all rows with BW
pfoa_cat_m_LB<- pfoa_cat_m_LB [which(pfoa_cat_m_LB$value > 0), ] # remove all PCP categories that have no exposure


head(pfoa_cat_m_LB)

### Box-plot of all MC simulations for all IDs per food cat ###



Plot_pfoa_cat_m_LB <- ggplot(data = pfoa_cat_m_LB, aes(y=value, x=PCP, colour = PCP)) +
  geom_boxplot()+
  scale_colour_hue()+
  theme_minimal()+
  scale_y_log10(limits=c(0.000001,10), labels = scales::comma)+
  theme(axis.text.x = element_text(angle = 90, size = 15), axis.title = element_text(size = 24))+
  xlab("PCP category")+
  ylab("PFOA exposure (ng/day)")+
  annotate("text",x=3, y=1000, label="Males LB", size=6)

Plot_pfoa_cat_m_LB


ggsave(file = file.path(newday,"Plot_pfoa_cat_m_LB.jpeg"),
       device = NULL,
       width=NA,
       height=NA,
       units="mm")


#
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ##### Calculation for females ####
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Amounts in gram per day per person for each PCP product 

use_amount_f <- array(unlist(ProductUsePCP_female),dim = c(100,27,MC))

# x= AmountsFemale, y=use_amount_f, z = ProductUsePCP_male

use_amount_f <- pcp_amount(AmountsFemale, use_amount_f, ProductUsePCP_female)

print(use_amount_f[14,,5])


####################PCP EXPOSURE Female ####################################


# Amount of PFOA (total gram used * ng/g conc = ng)


amount_PFOA_f_LB <- use_amount_f
amount_PFOA_f_LB <- compound_amount(PFOAconcPCP_LB_female, amount_PFOA_f_LB, use_amount_f)

amount_PFOA_f_MB <- use_amount_f
amount_PFOA_f_MB <- compound_amount(PFOAconcPCP_MB_female, amount_PFOA_f_MB, use_amount_f)

amount_PFOA_f_UB <- use_amount_f
amount_PFOA_f_UB <- compound_amount(PFOAconcPCP_UB_female, amount_PFOA_f_UB, use_amount_f)


print(amount_PFOA_f_LB[14,,5])


print(amount_PFOA_f_MB[3,,5])
print(amount_PFOA_f_MB[24,,5])

print(amount_PFOA_f_UB[3,,5])
print(amount_PFOA_f_UB[24,,5])

################################### Adding rinse factors ###########


amount_PFOA_f_LB_eef <- amount_PFOA_f_LB
amount_PFOA_f_LB_eef <- amount_rinse(EEFsFemale,amount_PFOA_f_LB_eef, amount_PFOA_f_LB )

amount_PFOA_f_MB_eef <- amount_PFOA_f_MB
amount_PFOA_f_MB_eef <- amount_rinse(EEFsFemale,amount_PFOA_f_MB_eef, amount_PFOA_f_MB )

amount_PFOA_f_UB_eef <- amount_PFOA_f_UB
amount_PFOA_f_UB_eef <- amount_rinse(EEFsFemale,amount_PFOA_f_UB_eef, amount_PFOA_f_UB )


print(amount_PFOA_f_LB_eef[14,,5]) # still ng


print(amount_PFOA_f_MB_eef[3,,5]) # still ng
print(amount_PFOA_f_MB_eef[24,,5])

print(amount_PFOA_f_UB_eef[3,,5]) # still ng
print(amount_PFOA_f_UB_eef[24,,5])



## sumarise the external exposure over all PCP products for all MC #####

total_f_LB <- as.data.frame(matrix(NA,nrow = 100,ncol = 1000))
total_f_LB <- sum_amount_PCP(total_f_LB, amount_PFOA_f_LB_eef)

print(total_f_LB [3,5]) # still 

total_f_MB <- as.data.frame(matrix(NA,nrow = 100,ncol = 1000))
total_f_MB <- sum_amount_PCP(total_f_MB, amount_PFOA_f_MB_eef)

print(total_f_MB [25,5]) # still ng

total_f_UB <- as.data.frame(matrix(NA,nrow = 100,ncol = 1000))
total_f_UB <- sum_amount_PCP(total_f_UB, amount_PFOA_f_UB_eef)

print(total_f_UB [25,5]) # still ng



###### calculate the exposure in ng/day # This can be the input for the PBPK modelling ########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

total_f_ID_LB <- as.data.frame(matrix(NA, nrow = 100, ncol = 1001))
total_f_ID_LB[,1] <- ProductUsePCP_female[,1] #get ID 
total_f_ID_LB[,2:1001] <- total_f_LB[,1:1000] # in ng/day

names(total_f_ID_LB)[1] <- "IDkode"

total_f_ID_MB <- as.data.frame(matrix(NA, nrow = 100, ncol = 1001))
total_f_ID_MB[,1] <- ProductUsePCP_female[,1] #get ID 
total_f_ID_MB[,2:1001] <- total_f_MB[,1:1000] # in ng/day

names(total_f_ID_MB)[1] <- "IDkode"

total_f_ID_UB <- as.data.frame(matrix(NA, nrow = 100, ncol = 1001))
total_f_ID_UB[,1] <- ProductUsePCP_female[,1] #get ID 
total_f_ID_UB[,2:1001] <- total_f_UB[,1:1000] # in ng/day

names(total_f_ID_UB)[1] <- "IDkode"


### make summary data for input to the PBPK model ng/day ##########################################

# LB
total_f_ID_LB_long <- total_f_ID_LB %>% pivot_longer(
  !IDkode, names_to = "MC", values_to = "value")  # make a long table

Sumtotal_f_ID_LB <-Sum(total_f_ID_LB_long) # Make summary data for females LB

head(Sumtotal_f_ID_LB)


write.xlsx(Sumtotal_f_ID_LB,
           file = file.path(newday,"Sumtotal_f_ID_LB.xlsx"),
           colNames = TRUE, borders = "rows"
)

# MB
total_f_ID_MB_long <- total_f_ID_MB %>% pivot_longer(
  !IDkode, names_to = "MC", values_to = "value")  # make a long table

Sumtotal_f_ID_MB <-Sum(total_f_ID_MB_long) # Make summary data for females MB

head(Sumtotal_f_ID_MB)


write.xlsx(Sumtotal_f_ID_MB,
           file = file.path(newday,"Sumtotal_f_ID_MB.xlsx"),
           colNames = TRUE, borders = "rows"
)

# UB
total_f_ID_UB_long <- total_f_ID_UB %>% pivot_longer(
  !IDkode, names_to = "MC", values_to = "value")  # make a long table

Sumtotal_f_ID_UB <-Sum(total_f_ID_UB_long) # Make summary data for females UB

head(Sumtotal_f_ID_UB)


write.xlsx(Sumtotal_f_ID_UB,
           file = file.path(newday,"Sumtotal_f_ID_UB.xlsx"),
           colNames = TRUE, borders = "rows"
)
### make summary data as ng/kg bw/day ##########################################

# LB
total_f_ID_LB_kg <- total_f_ID_LB
total_f_ID_LB_kg$BW <- ProductUsePCP_female[,2]
total_f_ID_LB_kg[,2:1001] <- total_f_ID_LB_kg[,2:1001]/total_f_ID_LB_kg$BW


total_f_ID_LB_kg  <- total_f_ID_LB_kg %>% dplyr::select(-BW) # remove BW
total_f_ID_LB_kg_long <- total_f_ID_LB_kg %>%  pivot_longer(
  !IDkode, names_to = "MC", values_to = "value") # make a long table

Sumtotal_f_ID_LB_kg <- Sum(total_f_ID_LB_kg_long)

head(Sumtotal_f_ID_LB_kg)

write.xlsx(Sumtotal_f_ID_LB_kg,
           file = file.path(newday,"Sumtotal_f_ID_LB_kg.xlsx"),
           colNames = TRUE, borders = "rows"
)

# MB
total_f_ID_MB_kg <- total_f_ID_MB
total_f_ID_MB_kg$BW <- ProductUsePCP_female[,2]
total_f_ID_MB_kg[,2:1001] <- total_f_ID_MB_kg[,2:1001]/total_f_ID_MB_kg$BW


total_f_ID_MB_kg  <- total_f_ID_MB_kg %>% dplyr::select(-BW) # remove BW
total_f_ID_MB_kg_long <- total_f_ID_MB_kg %>%  pivot_longer(
  !IDkode, names_to = "MC", values_to = "value") # make a long table

Sumtotal_f_ID_MB_kg <- Sum(total_f_ID_MB_kg_long)

head(Sumtotal_f_ID_MB_kg)

write.xlsx(Sumtotal_f_ID_MB_kg,
           file = file.path(newday,"Sumtotal_f_ID_MB_kg.xlsx"),
           colNames = TRUE, borders = "rows"
)

# UB
total_f_ID_UB_kg <- total_f_ID_UB
total_f_ID_UB_kg$BW <- ProductUsePCP_female[,2]
total_f_ID_UB_kg[,2:1001] <- total_f_ID_UB_kg[,2:1001]/total_f_ID_UB_kg$BW


total_f_ID_UB_kg  <- total_f_ID_UB_kg %>% dplyr::select(-BW) # remove BW
total_f_ID_UB_kg_long <- total_f_ID_UB_kg %>%  pivot_longer(
  !IDkode, names_to = "MC", values_to = "value") # make a long table

Sumtotal_f_ID_UB_kg <- Sum(total_f_ID_UB_kg_long)

head(Sumtotal_f_ID_UB_kg)

write.xlsx(Sumtotal_f_ID_UB_kg,
           file = file.path(newday,"Sumtotal_f_ID_UB_kg.xlsx"),
           colNames = TRUE, borders = "rows"
)


## making vectors of the melted MC's as basis fr plotting and summary ##

# LB
total_f_ID_LB_vector <- as.vector(total_f_ID_LB_long$value) # ng/day will be included in the main figure
total_f_ID_LB_kg_vector_kg_vector <- as.vector(total_f_ID_LB_kg_long$value) # pg/kg bw/day

# MB
total_f_ID_MB_vector <- as.vector(total_f_ID_MB_long$value) # ng/day will be included in the main figure
total_f_ID_MB_kg_vector_kg_vector <- as.vector(total_f_ID_MB_kg_long$value) # pg/kg bw/day

# UB
total_f_ID_UB_vector <- as.vector(total_f_ID_UB_long$value) # ng/day will be included in the main figure
total_f_ID_UB_kg_vector_kg_vector <- as.vector(total_f_ID_UB_kg_long$value) # pg/kg bw/day


## make summary data for all females based on all MC interations  #######

# LB
describe(total_f_ID_LB_vector) #(ng/day)
describe(total_f_ID_LB_kg_vector_kg_vector) #(pg/kg bw/day)

# MB
describe(total_f_ID_MB_vector) #(ng/day)
describe(total_f_ID_MB_kg_vector_kg_vector) #(pg/kg bw/day)

# UB
describe(total_f_ID_UB_vector) #(ng/day)
describe(total_f_ID_UB_kg_vector_kg_vector) #(pg/kg bw/day)



### PLOTT EXPOSURE TO PFOA FROM PCP FOR females LB, MB and UB ####

## make a data frame of the three vectors with summary iterations for LB, MB and UB (ng/day) ##

A = c(total_f_ID_LB_vector, total_f_ID_MB_vector,total_f_ID_UB_vector)
B = replicate(100000,"LB")
C = replicate(100000,"MB")
D = replicate(100000,"UB")
E = c(B,C,D)

pfoa_f_pcp <- as.data.frame(matrix(NA,nrow = 300000, ncol = 2))
pfoa_f_pcp[1] <- A
pfoa_f_pcp[2] <- E

names(pfoa_f_pcp)[1] <- "PFOAconc"
names(pfoa_f_pcp)[2] <- "Exposure"

## Make a cumulative density plot with the LB, MB and UB exposure ##

Plot_pfoa_f_pcp <- ggplot(data = pfoa_f_pcp, aes(PFOAconc, colour = Exposure)) +
  geom_line(stat = "ecdf", size = 0.5)+
  scale_colour_hue()+
  theme_minimal()+
  scale_x_continuous(limits=c(0.001,1.5))+
  xlab("PFOA (ng/day)")+
  ylab("Cumulative probability")+
  annotate("text", x=0.1, y=0.98, label="Females", size=6)

Plot_pfoa_f_pcp


ggsave(file = file.path(newday,"Plot_pfoa_f_pcp.jpeg"),
       device = NULL,
       width=NA,
       height=NA,
       units="mm")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  Combining plot for males and females and LB, MB and UB
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Combined_PFOA_Exp_External_PCP <- ggarrange(Plot_pfoa_m_pcp, Plot_pfoa_f_pcp+ rremove("ylab"), 
                                                   ncol = 2, nrow = 1, common.legend = TRUE) 
Combined_PFOA_Exp_External_PCP


ggsave(filename=file.path(newday,"Combined_PFOA_Exp_External_PCP_kg.jpeg"),
       device = NULL,
       width=NA,
       height=NA,
       units="mm")


