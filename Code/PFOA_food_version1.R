#########################################################################################
## Exposusre assessment for PFOA from diet using the EuroMix stuy
## Trine Husøy
## Date: 090523
##########################################################################################

# Sett work directory and organise results

HOME <- "C:/Users/TRHU/Documents/R/PFAS_exposure"

setwd(HOME)


# Create a folder with current date in the Result folder
newday <- file.path('C:/Users/TRHU/Documents/R/PFAS_exposure/Results', Sys.Date())
dir.create(newday)



# Load required packages and libraries

library(ggplot2)
library(writexl)
library(ggpubr)
library(Hmisc)
library(openxlsx)
library(tidyverse)
library(flextable)



# Calculate the location parameter (Loc or logmean) of the lognormal distribution


Lognorm_Loc <- function(x){
  log(
    x$mean^2/sqrt(x$sd^2+x$mean^2)
  )
  
}


# Calculate the shape (logsd) parameter of the lognormal distribution. 


Lognorm_shape <- function(x){
  sqrt(
    log(
      1 + x$sd^2/x$mean^2
    )
  )
}



# A function for summarising the data


Sum <-function(x){
  x %>%
    summarise(
      N = n(),
      mean = mean(Conc, na.rm=TRUE),
      sd=sd(Conc, na.rm=TRUE),
      min=min(Conc, na.rm=TRUE),
      P05=quantile(Conc, .05, na.rm=TRUE),
      P50=quantile(Conc, .50, na.rm=TRUE),
      P95=quantile(Conc, .95, na.rm=TRUE),
      max=max(Conc, na.rm=TRUE)
    )
}  

Sum_1 <-function(x){
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



# The MC_Sim function performs the Monte Carlo (MC) simulation. 


MC_sim <- function(x, y ,z){
  set.seed(123)
  MC <- 1000
  x <- array(dim = c(nrow(y),nrow(z),MC))
  y[,1:2] <- NULL
  for (u in 1:MC){
    for(i in 1:nrow(z)){
      x[,i,u] <- y[,i]*
        rlnorm(n=1,z$Loc[i], z$Shape[i])
    }
  }
  return(x)
}


# Read data 

FoodIntakeDiaryDay1 <-  read.csv2("./Data/foodintake_dummy_day1.csv")
FoodIntakeDiaryDay1 <- data.frame(sapply(FoodIntakeDiaryDay1, function(x) as.numeric(as.character(x)))) # convert to nummeric
FoodIntakeDiaryDay2 <- read.csv2("./Data/foodintake_dummy_day2.csv")
FoodIntakeDiaryDay2 <- data.frame(sapply(FoodIntakeDiaryDay2, function(x) as.numeric(as.character(x))))  #convert to numeric
PFAS_LB <- read.csv2("./Data/3-SumPFAS food conc_LB.csv")
PFAS_MB <- read.csv2("./Data/3-SumPFAS food conc_MB.csv")
PFAS_UB <- read.csv2("./Data/3-SumPFAS food conc_UB.csv")
SexWeight <- read.csv2("./Data/EuroMix_dummy_sex_weight.csv")

############################################
## Cleaning of the data
#################################################

FoodIntakeDiaryBothDays <- aggregate(.~ IDkode, rbind(FoodIntakeDiaryDay1,FoodIntakeDiaryDay2), sum)/2
FoodIntakeDiaryBothDays[,1] <- FoodIntakeDiaryBothDays[,1]*2

# Merge FoodIntake with sex and weight

FoodIntakeDiaryBothDays <- merge(SexWeight, FoodIntakeDiaryBothDays, by = "IDkode", all = TRUE)


# Rename columns

PFAS_LB <- PFAS_LB %>% rename(food = X)
PFAS_MB <- PFAS_MB %>% rename(food = X)
PFAS_UB <- PFAS_UB %>% rename(food = X)

# Delete the food category flatfish, since this was not reported eaten in the diaries

FoodIntakeDiaryBothDays <- FoodIntakeDiaryBothDays %>%  dplyr::select(-FlatFish)
PFAS_LB <- PFAS_LB %>% filter(food != "FlatFish")
PFAS_MB <- PFAS_MB %>% filter(food != "FlatFish")
PFAS_UB <- PFAS_UB %>% filter(food != "FlatFish")


# Delete the food category "fish", since we will build up this later from the other fish categories

FoodIntakeDiaryBothDays <- FoodIntakeDiaryBothDays %>%  dplyr::select(-Fish)
PFAS_LB <- PFAS_LB %>% filter(food != "Fish")
PFAS_MB <- PFAS_MB %>% filter(food != "Fish")
PFAS_UB <- PFAS_UB %>% filter(food != "Fish")

# Extract food eaten for males and females separately

FoodIntakeDiaryBothDays_male <- FoodIntakeDiaryBothDays %>% filter(cat_male == 2)
FoodIntakeDiaryBothDays_female <- FoodIntakeDiaryBothDays %>%  filter(cat_male == 1)

# Delete the categories not needed any longer 

FoodIntakeDiaryBothDays_male <- FoodIntakeDiaryBothDays_male %>%  dplyr::select(-cat_male)
FoodIntakeDiaryBothDays_female <- FoodIntakeDiaryBothDays_female %>%  dplyr::select(-cat_male)


# Extract the PFOA data from the concentration data


PFOA_LB <- PFAS_LB %>% dplyr::select(food, N, PFOA_lb_mean, PFOA_lb_SD)
colnames(PFOA_LB) <- c("Foods", "N", "mean", "sd")
PFOA_MB <- PFAS_MB %>% dplyr::select(food, N, PFOA_mb_mean, PFOA_mb_SD)
colnames(PFOA_MB) <- c("Foods", "N", "mean", "sd")
PFOA_UB <- PFAS_UB %>% dplyr::select(food, N, PFOA_ub_mean, PFOA_ub_SD)
colnames(PFOA_UB) <- c("Foods", "N", "mean", "sd")


# Calculate the location parameter of the log normal the distribution

PFOA_LB$Loc <- Lognorm_Loc(PFOA_LB)
PFOA_LB$Shape <- Lognorm_shape(PFOA_LB)
colnames(PFOA_LB) <- c("Foods", "N", "mean", "sd", "Loc", "Shape")

PFOA_MB$Loc <- Lognorm_Loc(PFOA_MB)
PFOA_MB$Shape <- Lognorm_shape(PFOA_MB)
colnames(PFOA_MB) <- c("Foods", "N", "mean", "sd", "Loc", "Shape")

PFOA_UB$Loc <- Lognorm_Loc(PFOA_UB)
PFOA_UB$Shape <- Lognorm_shape(PFOA_UB)
colnames(PFOA_UB) <- c("Foods", "N", "mean", "sd", "Loc", "Shape")


# Replace NAs by zero

FoodIntakeDiaryBothDays_male[is.na(FoodIntakeDiaryBothDays_male)] <- 0
FoodIntakeDiaryBothDays_female[is.na(FoodIntakeDiaryBothDays_female)] <- 0
PFOA_LB[is.na(PFOA_LB)] <- 0
PFOA_MB[is.na(PFOA_MB)] <- 0
PFOA_UB[is.na(PFOA_UB)] <- 0

############################################################
# Probabilistic exposure estimate for males
###############################################################

FoodBothDays_male_PFOA_LB <- MC_sim(FoodBothDays_male_PFOA_LB, FoodIntakeDiaryBothDays_male, PFOA_LB)/1000000 # divide on 1000000 to get in ug/day
FoodBothDays_male_PFOA_MB <- MC_sim(FoodBothDays_male_PFOA_MB, FoodIntakeDiaryBothDays_male, PFOA_MB)/1000000
FoodBothDays_male_PFOA_UB <- MC_sim(FoodBothDays_male_PFOA_UB, FoodIntakeDiaryBothDays_male, PFOA_UB)/1000000

print(FoodBothDays_male_PFOA_LB[1,,7])
print(FoodBothDays_male_PFOA_MB[1,,7])
print(FoodBothDays_male_PFOA_UB[1,,7])


## Calculate the total PFOA exposure in ug/day per individual for each MC iteration


# LB
PFOA_BothDays_LB_male<-as.data.frame(matrix(NA, nrow = 44, ncol = 1000))
 for (i in 1:1000){
  PFOA_BothDays_LB_male[,i] <- rowSums(FoodBothDays_male_PFOA_LB[,,i]) 
  }
PFOA_BothDays_LB_male$IDkode <- FoodIntakeDiaryBothDays_male[,1]
rownames(PFOA_BothDays_LB_male) <- PFOA_BothDays_LB_male$IDkode
 
#MB
PFOA_BothDays_MB_male<-as.data.frame(matrix(NA, nrow = 44, ncol = 1000))
 for (i in 1:1000){
  PFOA_BothDays_MB_male[,i] <- rowSums(FoodBothDays_male_PFOA_MB[,,i]) 
  }
PFOA_BothDays_MB_male$IDkode <- FoodIntakeDiaryBothDays_male[,1]
rownames(PFOA_BothDays_MB_male) <- PFOA_BothDays_MB_male$IDkode

#UB
PFOA_BothDays_UB_male<-as.data.frame(matrix(NA, nrow = 44, ncol = 1000))
 for (i in 1:1000){
  PFOA_BothDays_UB_male[,i] <- rowSums(FoodBothDays_male_PFOA_UB[,,i]) 
  }
PFOA_BothDays_UB_male$IDkode <- FoodIntakeDiaryBothDays_male[,1]
rownames(PFOA_BothDays_UB_male) <- PFOA_BothDays_UB_male$IDkode



## Make summary data from total exposure for each individual for input to the PBPK model ug/day


PFOA_BothDays_LB_male_long <- PFOA_BothDays_LB_male %>%  
  pivot_longer(cols = -IDkode, values_to = "value") %>%
  as_tibble()
PFOA_BothDays_LB_male_long <- PFOA_BothDays_LB_male_long %>% dplyr::select(-name)

PFOA_BothDays_MB_male_long <- PFOA_BothDays_MB_male %>%  
  pivot_longer(cols = -IDkode, values_to = "value") %>% 
  as_tibble()
PFOA_BothDays_MB_male_long <- PFOA_BothDays_MB_male_long %>% dplyr::select(-name)

PFOA_BothDays_UB_male_long <- PFOA_BothDays_UB_male %>% 
  pivot_longer(cols = -IDkode, values_to = "value") %>% 
  as_tibble()
PFOA_BothDays_UB_male_long <- PFOA_BothDays_UB_male_long %>% dplyr::select(-name)


SumPFOA_BothDays_LB_male <- Sum_1(PFOA_BothDays_LB_male_long)
SumPFOA_BothDays_MB_male <- Sum_1(PFOA_BothDays_MB_male_long)
SumPFOA_BothDays_UB_male <- Sum_1(PFOA_BothDays_UB_male_long)




### Make summary data across all exposure assessmets for all males (ug/day)


PFOA_LB_males_BothDays_vector <- as.vector(PFOA_BothDays_LB_male_long$value)
plot(ecdf(PFOA_LB_males_BothDays_vector))


SumPFOA_LB_male_BothDays <- describe(PFOA_LB_males_BothDays_vector)
SumPFOA_LB_male_BothDays


## Make summary data for all males based on all MC iterations expressed as ug/kg bw/day

# Divide the ug PFOA on individual BW

#LB
PFOA_BothDays_LB_male_kg <- PFOA_BothDays_LB_male_long %>% 
  left_join(SexWeight, by = c("IDkode" = "IDkode")) %>% 
  dplyr::select(IDkode, value, con_weight) 
PFOA_BothDays_LB_male_kg$value <- PFOA_BothDays_LB_male_kg$value/PFOA_BothDays_LB_male_kg$con_weight
PFOA_BothDays_LB_male_kg <- PFOA_BothDays_LB_male_kg %>%  dplyr::select(-con_weight)

SumPFOA_BothDays_LB_male_kg <- Sum_1(PFOA_BothDays_LB_male_kg)

#MB
PFOA_BothDays_MB_male_kg <- PFOA_BothDays_MB_male_long %>% 
  left_join(SexWeight, by = c("IDkode" = "IDkode")) %>% 
  dplyr::select(IDkode, value, con_weight) 
PFOA_BothDays_MB_male_kg$value <- PFOA_BothDays_MB_male_kg$value/PFOA_BothDays_MB_male_kg$con_weight
PFOA_BothDays_MB_male_kg <- PFOA_BothDays_MB_male_kg %>%  dplyr::select(-con_weight)

SumPFOA_BothDays_MB_male_kg <- Sum_1(PFOA_BothDays_MB_male_kg)


#UB
PFOA_BothDays_UB_male_kg <- PFOA_BothDays_UB_male_long %>% 
  left_join(SexWeight, by = c("IDkode" = "IDkode")) %>% 
  dplyr::select(IDkode, value, con_weight) 
PFOA_BothDays_UB_male_kg$value <- PFOA_BothDays_UB_male_kg$value/PFOA_BothDays_UB_male_kg$con_weight
PFOA_BothDays_UB_male_kg <- PFOA_BothDays_UB_male_kg %>%  dplyr::select(-con_weight)

SumPFOA_BothDays_UB_male_kg <- Sum_1(PFOA_BothDays_UB_male_kg)


### Save the table as excel file


write.xlsx(SumPFOA_BothDays_LB_male_kg, file = file.path(newday,"SumPFOA_BothDays_LB_male_kg.xlsx"),
           colNames = TRUE, borders = "rows")

write.xlsx(SumPFOA_BothDays_MB_male_kg, file = file.path(newday,"SumPFOA_BothDays_MB_male-kg.xlsx"),
           colNames = TRUE, borders = "rows")

write.xlsx(SumPFOA_BothDays_UB_male_kg, file = file.path(newday,"SumPFOA_BothDays_UB_male_kg.xlsx"),
           colNames = TRUE, borders = "rows")


## Plott external exposure to PFOA for males for LB, MB, and UB


A = c(PFOA_BothDays_LB_male_long$value, PFOA_BothDays_MB_male_long$value, PFOA_BothDays_UB_male_long$value)
B = replicate(44000,"LB")
C = replicate(44000,"MB")
D = replicate(44000,"UB")
E = c(B,C,D)

PFOA_males_BothDays_LB_MB_UB <- as.data.frame(matrix(NA,nrow = 132000, ncol = 2))
PFOA_males_BothDays_LB_MB_UB[1] <- A
PFOA_males_BothDays_LB_MB_UB[2] <- E

names(PFOA_males_BothDays_LB_MB_UB)[1] <- "PFOAconc"
names(PFOA_males_BothDays_LB_MB_UB)[2] <- "Exposure"


### Make a cumulative density plot with the LB, MB and UB exposure

Plot_PFOA_males_diary_LB_MB_UB <- ggplot(data = PFOA_males_BothDays_LB_MB_UB, aes(PFOAconc, colour = Exposure)) +
  geom_line(stat = "ecdf", size = 0.5)+
  scale_colour_hue()+
  theme_minimal()+
  scale_x_continuous(limits=c(0.01,0.5))+
  xlab("PFOA (ng/kg bw/day)")+
  ylab("Cumulative probability")+
  annotate("text", x=0.45, y=0.95, label="Males", size=4)

Plot_PFOA_males_diary_LB_MB_UB

ggsave(filename=file.path(newday,"Plot_PFOA_males_diary_LB_MB_UB.jpeg"),
       device = NULL,
       width=NA,
       height=NA,
       units="mm")

###################################################
# Probabilistic exposure for females
###################################################


FoodBothDays_female_PFOA_LB <- MC_sim(FoodBothDays_female_PFOA_LB, FoodIntakeDiaryBothDays_female, PFOA_LB)/1000000 # divide on 1000000 to get in ug/day
FoodBothDays_female_PFOA_MB <- MC_sim(FoodBothDays_female_PFOA_MB, FoodIntakeDiaryBothDays_female, PFOA_MB)/1000000
FoodBothDays_female_PFOA_UB <- MC_sim(FoodBothDays_female_PFOA_UB, FoodIntakeDiaryBothDays_female, PFOA_UB)/1000000

print(FoodBothDays_female_PFOA_LB[1,,7])
print(FoodBothDays_female_PFOA_MB[1,,7])
print(FoodBothDays_female_PFOA_UB[1,,7])


## Calculate the total PFOA exposure in ug/day per individual for each MC iteration


# LB
PFOA_BothDays_LB_female<-as.data.frame(matrix(NA, nrow = 100, ncol = 1000))
 for (i in 1:1000){
  PFOA_BothDays_LB_female[,i] <- rowSums(FoodBothDays_female_PFOA_LB[,,i]) 
  }
PFOA_BothDays_LB_female$IDkode <- FoodIntakeDiaryBothDays_female[,1]
rownames(PFOA_BothDays_LB_female) <- PFOA_BothDays_LB_female$IDkode
 
#MB
PFOA_BothDays_MB_female<-as.data.frame(matrix(NA, nrow = 100, ncol = 1000))
 for (i in 1:1000){
  PFOA_BothDays_MB_female[,i] <- rowSums(FoodBothDays_female_PFOA_MB[,,i]) 
  }
PFOA_BothDays_MB_female$IDkode <- FoodIntakeDiaryBothDays_female[,1]
rownames(PFOA_BothDays_MB_female) <- PFOA_BothDays_MB_female$IDkode

#UB
PFOA_BothDays_UB_female<-as.data.frame(matrix(NA, nrow = 100, ncol = 1000))
 for (i in 1:1000){
  PFOA_BothDays_UB_female[,i] <- rowSums(FoodBothDays_female_PFOA_UB[,,i]) 
  }
PFOA_BothDays_UB_female$IDkode <- FoodIntakeDiaryBothDays_female[,1]
rownames(PFOA_BothDays_UB_female) <- PFOA_BothDays_UB_female$IDkode


## Make summary data from total exposure for each individual for input to the PBPK model ug/day


PFOA_BothDays_LB_female_long <- PFOA_BothDays_LB_female %>%  
  pivot_longer(cols = -IDkode, values_to = "value") %>%
  as_tibble()
PFOA_BothDays_LB_female_long <- PFOA_BothDays_LB_female_long %>% dplyr::select(-name)

PFOA_BothDays_MB_female_long <- PFOA_BothDays_MB_female %>%  
  pivot_longer(cols = -IDkode, values_to = "value") %>% 
  as_tibble()
PFOA_BothDays_MB_female_long <- PFOA_BothDays_MB_female_long %>% dplyr::select(-name)

PFOA_BothDays_UB_female_long <- PFOA_BothDays_UB_female %>% 
  pivot_longer(cols = -IDkode, values_to = "value") %>% 
  as_tibble()
PFOA_BothDays_UB_female_long <- PFOA_BothDays_UB_female_long %>% dplyr::select(-name)


SumPFOA_BothDays_LB_female <- Sum_1(PFOA_BothDays_LB_female_long)
SumPFOA_BothDays_MB_female <- Sum_1(PFOA_BothDays_MB_female_long)
SumPFOA_BothDays_UB_female <- Sum_1(PFOA_BothDays_UB_female_long)



### Make summary data across all exposure assessments for all females

PFOA_LB_females_BothDays_vector <- as.vector(PFOA_BothDays_LB_female_long$value)
plot(ecdf(PFOA_LB_females_BothDays_vector))

## make summary data for all males based on all MC interations (ng/day) #######

SumPFOA_LB_female_BothDays <- describe(PFOA_LB_females_BothDays_vector)
SumPFOA_LB_female_BothDays
SumPFOA_LB_female_BothDays <- summary(PFOA_LB_females_BothDays_vector)
SumPFOA_LB_female_BothDays
quantile(PFOA_LB_females_BothDays_vector, probs = c(.05, .5, .95))


## Make summary data for all females based on all MC iterations expressed as ug/kg bw/day

# Divide the ug PFOA on individual BW

#LB
PFOA_BothDays_LB_female_kg <- PFOA_BothDays_LB_female_long %>% 
  left_join(SexWeight, by = c("IDkode" = "IDkode")) %>% 
  dplyr::select(IDkode, value, con_weight) 
PFOA_BothDays_LB_female_kg$value <- PFOA_BothDays_LB_female_kg$value/PFOA_BothDays_LB_female_kg$con_weight
PFOA_BothDays_LB_female_kg <- PFOA_BothDays_LB_female_kg %>%  dplyr::select(-con_weight)

SumPFOA_BothDays_LB_female_kg <- Sum_1(PFOA_BothDays_LB_female_kg)

#MB
PFOA_BothDays_MB_female_kg <- PFOA_BothDays_MB_female_long %>% 
  left_join(SexWeight, by = c("IDkode" = "IDkode")) %>% 
  dplyr::select(IDkode, value, con_weight) 
PFOA_BothDays_MB_female_kg$value <- PFOA_BothDays_MB_female_kg$value/PFOA_BothDays_MB_female_kg$con_weight
PFOA_BothDays_MB_female_kg <- PFOA_BothDays_MB_female_kg %>%  dplyr::select(-con_weight)

SumPFOA_BothDays_MB_female_kg <- Sum_1(PFOA_BothDays_MB_female_kg)


#UB
PFOA_BothDays_UB_female_kg <- PFOA_BothDays_UB_female_long %>% 
  left_join(SexWeight, by = c("IDkode" = "IDkode")) %>% 
  dplyr::select(IDkode, value, con_weight) 
PFOA_BothDays_UB_female_kg$value <- PFOA_BothDays_UB_female_kg$value/PFOA_BothDays_UB_female_kg$con_weight
PFOA_BothDays_UB_female_kg <- PFOA_BothDays_UB_female_kg %>%  dplyr::select(-con_weight)

SumPFOA_BothDays_UB_female_kg <- Sum_1(PFOA_BothDays_UB_female_kg)



### Save the table as excel file


write.xlsx(SumPFOA_BothDays_LB_female_kg, file = file.path(newday,"SumPFOA_BothDays_LB_female_kg.xlsx"),
           colNames = TRUE, borders = "rows")

write.xlsx(SumPFOA_BothDays_MB_female_kg, file = file.path(newday,"SumPFOA_BothDays_MB_female-kg.xlsx"),
           colNames = TRUE, borders = "rows")

write.xlsx(SumPFOA_BothDays_UB_female_kg, file = file.path(newday,"SumPFOA_BothDays_UB_female_kg.xlsx"),
           colNames = TRUE, borders = "rows")


## Plott external exposure to PFOA for females for LB, MB, and UB

A = c(PFOA_BothDays_LB_female_long$value, PFOA_BothDays_MB_female_long$value, PFOA_BothDays_UB_female_long$value)
B = replicate(100000,"LB")
C = replicate(100000,"MB")
D = replicate(100000,"UB")
E = c(B,C,D)

PFOA_females_BothDays_LB_MB_UB <- as.data.frame(matrix(NA,nrow = 300000, ncol = 2))
PFOA_females_BothDays_LB_MB_UB[1] <- A
PFOA_females_BothDays_LB_MB_UB[2] <- E

names(PFOA_females_BothDays_LB_MB_UB)[1] <- "PFOAconc"
names(PFOA_females_BothDays_LB_MB_UB)[2] <- "Exposure"


### Make a cumulative density plot with the LB, MB and UB exposure

Plot_PFOA_females_diary_LB_MB_UB <- ggplot(data = PFOA_females_BothDays_LB_MB_UB, aes(PFOAconc, colour = Exposure)) +
  geom_line(stat = "ecdf", size = 0.5)+
  scale_colour_hue()+
  theme_minimal()+
  scale_x_continuous(limits=c(0.01,0.5))+
  xlab("PFOA (ng/kg bw/day)")+
  ylab("Cumulative probability")+
  annotate("text", x=0.45, y=0.95, label="Females", size=4)

Plot_PFOA_females_diary_LB_MB_UB

ggsave(filename=file.path(newday,"Plot_PFOA_females_diary_LB_MB_UB.jpeg"),
       device = NULL,
       width=NA,
       height=NA,
       units="mm")


######################################
# Combine plots for males and females
######################################


PFOA_males_females <- ggarrange(Plot_PFOA_males_diary_LB_MB_UB, Plot_PFOA_females_diary_LB_MB_UB + rremove("ylab"),
                                              ncol = 2, nrow = 1, common.legend = TRUE) 
PFOA_males_females

ggsave(filename=file.path(newday,"PFOA_males_females.jpeg"),
       device = NULL,
       width=NA,
       height=NA,
       units="mm")


###################################
# Merge data for males and females
###################################

SumPFOA_LB_food <- rbind(SumPFOA_BothDays_LB_female_kg, SumPFOA_BothDays_LB_male_kg)
SumPFOA_LB_food <- arrange(SumPFOA_LB_food, IDkode)
SumPFOA_LB_food <- SumPFOA_LB_food %>% dplyr::select(-N)
SumPFOA_LB_food <- SumPFOA_LB_food %>%  rename_with(.fn = function(.x){paste0( .x, "_oral")},
                                                    .cols = c(mean, sd, min, P05, P50, P95, max))


SumPFOA_LB_food_PCP <- merge(SumPFOA_LB_food, SumPFOA_LB_PCP, by = "IDkode", all = TRUE)
SumPFOA_LB_food_PCP <- merge(SumPFOA_LB_food_PCP, SexWeight, by = "IDkode", all = TRUE)
SumPFOA_LB_food_PCP$Tmc_Sex <-ifelse(SumPFOA_LB_food_PCP$cat_male == 1, 5550, 6500)
SumPFOA_LB_food_PCP <- rename(SumPFOA_LB_food_PCP, Sex=cat_male, BW=con_weight)
SumPFOA_LB_food_PCP <- SumPFOA_LB_food_PCP %>%  relocate(c("Sex", "BW", "Tmc_Sex"), .before = mean_oral)

write.xlsx(SumPFOA_LB_food_PCP, file = file.path(newday,"SumPFOA_LB_food_PCP.xlsx"),
           colNames = TRUE, borders = "rows")


