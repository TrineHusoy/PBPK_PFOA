####################################################################################################
# Sensitivity analysis of the PBPK model of PFOA using pksensi
# Trine Husøy
# Date: 090523
#####################################################################################################

HOME <- "your work directory"

setwd(HOME)


newday <- file.path('your work directory/Results', Sys.Date())
dir.create(newday)


library(data.table)
library(ggplot2)
library(openxlsx)
library(deSolve)
library(writexl)
library(sensitivity)
library(pksensi)
library(PKNCA)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### PBPK model PFOA ####
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### START PBPK model ###

## Time parameters ##

starttime = 0.1 # in hours
stoptime = (50*365*24)+0.1  # end of simmulation (years) - (50 years) 
dtout <-8760.1  # resolution of the output time (h)
times <- seq(starttime, stoptime, by=dtout) # not in the original code 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Physiological parameters (from Brown, et al 1997)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QCC = 12.5  # Cardiac blood output (L/h/kg^0.75)
QFC = 0.052 # Fraction cardiac output going to fat
QLC = 0.069 # Fraction cardiac output going to liver suggested by EFSA is 0.069, since they corrected for output to the gut
# Luccisano used 0.25, but Brown actually report 0.23
QKC = 0.175 # Fraction cardiac output going to kidney
QSkC = 0.058 # Fraction cardiac output going to skin
QGC = 0.181 # Fraction of cardiac output going to gut and the liver via portal arthery

BW = 59      # Body weight from the EuroMix study

## fractional tissue columes ##

VLC = 0.026  # Fraction liver volume
VFC = 0.214  # Fraction fat volume
VKC = 0.004  # Fraction kidney volume
VfilC = 0.0004  # Fraction filtrate compartment volume (10% of kidney volume)
VGC = 0.0171  # Fraction gut volume
VPlasC = 0.0428  # Fraction plasma volume (58% of blood)
Htc = 0.44  # hematocrit

Skinarea = 972  # Exposed area on skin (cm^2); Changes from 5 since 2018 EFSA
SkinTarea = 9.1*(BW*1000)^0.666  # Total area of the skin (cm^2)
Skinthickness = 0.1 # Skin thickness (cm)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Scaled cardiac output and blood flows##
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QC = QCC*BW^0.75  # Cardiac output adjusted for BW (L/h)
QCP = QC*(1-Htc)  # Cardiac output adjusted for plasma flow (L/h)

QTC=QLC+QFC+QKC+QGC+QSkC # Preparation for scaling

QL = QLC*QCP/QTC  # Scaled plasma flow to liver (L/h)
QF = QFC*QCP/QTC  # Scaled plasma flow to fat (L/h)
QK = QKC*QCP/QTC  # Scaled plasma flow to kidney (L/h)

QG = QGC*QCP/QTC  # Scaled plasma flow to gut (L/h)

QSk <- QSkC*QCP*(Skinarea/SkinTarea)/QTC # scaled plasma flow to the skin

QR = QCP-QL-QF-QK-QG-QSk # Plasma flow to the rest of the body.


## Scaled tissue volumes ##
VTC = VLC+VFC+VKC+VfilC+VGC+VPlasC +(Skinarea*Skinthickness)/1000 # Preparation for scaling 
VL = VLC*BW/VTC  # Scaled liver volume (L)
VF = VFC*BW/VTC  # Scaled fat volume (L)
VK = VKC*BW/VTC  # Scaled kidney volume (L)
Vfil = VfilC*BW/VTC  # Scaled filtrate compartment volume
VG = VGC*BW/VTC  # Gut volume (L)
VPlas = VPlasC*BW/VTC  # Scaled plasma volume(L)

VSk = (Skinarea*Skinthickness)/1000/VTC  # Skin volume (L)
VR = 0.84*BW-VL-VF-VK-Vfil-VG-VPlas-VSk  # Rest of the body volume (L). Need to know where the number 0.84 comes from???



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Chemical specific parameters (PFOA)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Free = 0.02 # Free fraction of PFOA in plasma
PL =2.2  # Plasma/liver partition coefficient
PF = 0.04  # Plasma/fat partition coefficient
PK = 1.05  # Plasma/kidney partition coefficient
PSk = 0.1  # Plasma/skin partition coefficient
PR = 0.12  # Plasma/rest of the body partition coefficient
PG = 0.05  # Plasma/gut partition coefficient

AbsPFOA = 0.016  # Dermal absorption fraction take from Fasano et al 2008. This is in line with the absorption reported from Franko et al 2012 divided on 1000 since PFOA is ionized 

# kidney and urine #

Tmc = 6250  # ug/h/kg^0.75 Maximum resorption rate, changed from 6 in the original Loccisano 2011 model (ug)
# representing a half-life of 2.3 years

Tm = Tmc*BW^0.75 # transporter maximum

Kt = 55  # Resorption affinity, changed from 0.055 in the original Loccisano 2011 model (ug)
# Expressed in ug to be consistent with the other parameters
kurinec = 0.00000183  # 0.044/24/1000 urinary clearance L/h/kg calculated from the urinary clearance of 0.044 ml/h/kg  taken from Fujii et al 2015 
kurine = kurinec*BW^(-0.25) # clearance urine (L/h), 

# Clearance parameters# 

kbiliaryc = 0.000109 # 2.62/24/1000 biliary clearance L/h/kg calculated from the biliary clearance of 2.62 ml/day  taken from Fujii et al 2015 
kfaecesc = 0.00000217 # 0.052/24/1000 faeces clearance L/h/kg clearance in faeces taken from Fujii et al 2015, calculated from 0.052 ml/day/kg

kbiliary = kbiliaryc*BW^0.1# L/h biliary clearence, BW adjusted to the volume from liver
kfaeces = kfaecesc*BW^0.001  # L/h faeces clearance, BW adjusted to the volume of GI tract

kfil = 0.2*QK  # Clearance from the kidney to the filtrate compartment (L/h); 20% of bloodstream to QK is cleared for 

## Free fraction of chemical in tissues ##

FreeL = Free/PL  # liver
FreeF = Free/PF  # fat
FreeK = Free/PK  # kidney
FreeSk = Free/PSk  # skin
FreeR = Free/PR  # rest of the body
FreeG = Free/PG  # gut


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Dosing parameters
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tchng = 50*365*24  # Duration of exposure (h); 50 years; turn dose on and off

## Dermal exposure ##


Dermconc = 4.328542e-03  # Dermal concentration (ug/kg bw/day) for one individual
Dermdose = Dermconc*BW*AbsPFOA     # Internal doseby multiplying with body weight and skin absorption (ug/day)


## Oral exposure ##

Oralconc = 9.129959e-05  # Oral concentration (ug/kg bw/day. Mean from ID1 from the EuroMix study
Oraldose = Oralconc*BW  # multiply with body weight ug/day

Tinput = 24  # duration of dose (h)

tinterval = 24 # the interval the dose should be repeated (h). 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compile parameters 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

para <- unlist(c(data.frame(
  Htc,
  Tmc,
  Kt,
  Free,
  BW,
  kurinec,
  PL,
  PF,
  PK,
  PSk,
  PR,
  PG,
  kurine,
  kbiliaryc,
  kfaecesc,
  kbiliary, 
  kfaeces,
  FreeL,
  FreeF,
  FreeK,
  FreeSk,
  FreeR,
  FreeG,
  tchng,
  QC,
  QCP,
  QL,
  QF,
  QK,
  kfil,
  QG,
  QSk,
  QR,
  VL,
  VF,
  VK,
  Vfil,
  VG,
  VPlas,
  VSk,
  VR,
  Tm,
  SkinTarea,
  AbsPFOA
  
)))

para

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Initial conditions
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

yini <- unlist(c(data.frame(
  Input1 = 0.0,  # amount of PFOA entering the body via oral dosing
  Input2 = 0.0, # amount of PFOA entering the body via dermal exposure
  APlas = 0.0,  # amount of PFOA in Plasma
  AG = 0.0,  # amount of PFOA in gut
  AL = 0.0,  # amount of PFOA in liver
  AF = 0.0,  # amount of PFOA in fat
  AK = 0.0,  # amount of PFOA in kidney
  Afil = 0.0,  # amount of PFOA in kidney filtrate
  Adelay = 0.0,  # amount of PFOA in storage compartment of urine
  Aurine = 0.0,  # amount of PFOA in urine
  Afaeces = 0.0, # amount of PFOA in faeces
  ASk = 0.0,  # amount of PFOA in skin
  AR = 0.0  # amount of PFOA in the rest of the body
)))


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Model for PFOA
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PBPKmodPFOA <- function(t,state,parameters){
  with(as.list(c(state,parameters)), {
    
    if(t<tchng){DoseOn=1} else{DoseOn=0}
    
    Input1 <- Oraldose/Tinput*(t %% tinterval<Tinput) 
    Input2 <- Dermdose/Tinput*(t %% tinterval<Tinput)
    
    # concentrations in the organs #
    CAFree <- APlas/VPlas  # free concentration of PFOA in plasma in (ug/L)
    CA <- CAFree/Free  # total concentration in plasma
    CG <- AG/VG  # concentration of PFOA in gut (ug/L)
    CVG <- CG/PG # concentration of PFOA leaving gut (ug/L)
    CL <- AL/VL  # concentration of PFOA in liver (ug/L)
    CVL <- CL/PL  # concentration of PFOA leaving liver (ug/L)
    CF <- AF/VF  # concentration of PFOA in fat (ug/L)
    CVF <- CF/PF # concentration of PFOA leaving fat (ug/L)
    CK <- AK/VK  # concentration of PFOA in kidney (ug/L)
    CVK <- CK/PK  # concentration of PFOA leaving the kidney (ug/L)
    CSk <- ASk/VSk  # concentration of PFOA in skin (ug/L)
    CVSk <- CSk/PSk  # Concentration of PFOA leaving the skin (ug/L)
    CR <- AR/VR  # concentration of PFOA in rest of the body (ug/L)
    CVR <- CR/PR  # concentration of PFOA leaving the rest of the body (ug/L)
    Cfil <- Afil/Vfil # concentration of PFOA in filtrate compartment
    
    
    ## Plasma compartment
    RPlas <- QF*CF*FreeF+(QL+QG)*CL*FreeL+QR*CR*FreeR+QSk*CSk*FreeSk+
      QK*CK*FreeK-QCP*CA*Free #-Qfil*CA*Free 
    # Rate of PFOA amount change in plasma (ug/h)
    
    
    ## Gut compartment
    RG <- QG*(CA*Free-CG*FreeG)+CL*FreeL*kbiliary+Input1*DoseOn-CG*FreeG*kfaeces
    # Rate of PFOA amount change in gut (ug/h)
    
    ## Liver compartment
    RL <- QL*CA*Free+QG*CG*FreeG-(QL+QG)*CL*FreeL-CL*FreeL*kbiliary
    # Rate of PFOA amount change in the liver (ug/h)
    
    ## Fat compartment
    RF <- QF*(CA*Free-CF*FreeF) # Rate of PFOA amount change in fat (ug/h)
    
    ## Kidney compartment
    RK <- QK*(CA*Free-CK*FreeK)+Tm*Cfil/(Kt+Cfil)-kfil*CK*Free # Qfil*CK*Free was introduced to reflect clearance to filtrate compartment
    # Rate of PFOA amount change in kidney (ug/h)
    
    ## Filtrate compartment
    Rfil <- kfil*(CK*Free-Cfil)-Tm*Cfil/(Kt+Cfil) # changed from Qfil*CA*Free 
    # Rate of PFOA amount change in filtrate compartment (ug/h)
    
    ## Storage compartment for urine
    Rdelay <- kfil*Cfil-kurine*Adelay # ug/h
    
    ## Urine
    Rurine <- kurine*Adelay # ug/h
    
    ## Feaces
    
    Rfaeces <- CG*FreeG*kfaeces #B5g/h
    
    ## Skin compartment
    RSk <- QSk*(CA*Free-CSk*FreeSk)+Input2*DoseOn # Rate of PFOA amount change in skin
    
    ## Rest of the body
    RR <- QR*(CA*Free-CR*FreeR)
    
    PFOAPlasTot <- APlas/(VPlasC*BW)/Free #microg/litre = ng/ml
    PFOAUrine <- Aurine
    PFOALiver <- AL/VL
    PFOAKidney <- AK/VK
    
    
    list(c(Input1, Input2, RPlas, RG, RL, RF, RK, Rfil, Rdelay, Rurine,  Rfaeces, RSk, RR),
         "PFOAPlasTot" = PFOAPlasTot) # This will be changed to measure the uncertainty/sensitivity in other organs
    
    
  })
}

outputs = c("PFOAPlasTot") # This will be changed to measure the uncertainty/sensitivity in other organs


# Solve the system of differential equations


results <- ode(y=yini, times=times, func=PBPKmodPFOA, parms=para, method="lsoda")

head(results)


PFOAamount <- as.data.frame(results)

### END PBPK ###


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Check mass balance
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Qbal = QCP-(QR+QL+QF+QK+QG+QSk) # Mass balance check for the cardiac output 
print(Qbal)

Vbal = (0.84*BW)-(VL+VF+VK+Vfil+VG+VPlas+VSk+VR)  # Mass balance check for the volumes
print(Vbal)

PFOA_bal <- sum(PFOAamount[,"Input1"]+ PFOAamount[,"Input2"]- PFOAamount[,"APlas"]- # Mass balance for PFOA
                  PFOAamount[,"AG"]-PFOAamount[,"AL"]-PFOAamount[,"AF"]-PFOAamount[,"AK"]-PFOAamount[,"AR"]-
                  PFOAamount[,"ASk"]-PFOAamount[,"Afil"]-PFOAamount[,"Aurine"]-PFOAamount[,"Adelay"]-PFOAamount[,"Afaeces"])

print(PFOA_bal)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     Sensitivity analysis using the pksensi package in R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Define the distribution of the parameters that you will analyse in the sensitivity test 

q <- c( "qunif" , "qunif" , "qunif" , "qunif", "qunif", "qunif", "qunif" , "qunif" , "qunif" , "qunif", 
        "qunif" , "qunif" , "qunif" , "qunif", "qunif", "qunif", "qunif") # , "qunif" , "qunif" , "qunif", 
        #"qunif", "qunif" , "qunif")



## Set parameter distribution ##
# we use 10% change in all parameters

LL <- 0.9 # 10% lower limit
UL <- 1.1 # 10% upper limit

q.arg <- list(list(min = para["Htc"]*LL, max= para["Htc"]*UL),
              list(min = para["Tmc"]*LL, max= para["Tmc"]*UL),
              list(min = para["Kt"]*LL, max =  para["Kt"]*UL),
              list(min = para["Free"]*LL, max =  para["Free"]*UL),
              list(min = para["BW"]*LL, max =  para["BW"]*UL),
              list(min = para["kurinec"]*LL, max = para["kurinec"]*UL),
              list(min = para["kbiliaryc"]*LL, max = para["kbiliaryc"]*UL),
              list(min = para["kfaecesc"]*LL, max = para["kfaecesc"]*UL),
              list(min = para["kfil"]*LL, max = para["kfil"]*UL),
              list(min = para["PL"]*LL, max = para["PL"]*UL),
              list(min = para["PF"]*LL, max = para["PF"]*UL),
              list(min = para["PK"]*LL, max = para["PK"]*UL),
              list(min = para["PSk"]*LL, max = para["PSk"]*UL),
              list(min = para["PR"]*LL, max = para["PR"]*UL),
              list(min = para["PG"]*LL, max = para["PG"]*UL),
              # list(min = para["VLC"]*LL, max = para["VLC"]*UL),
              # list(min = para["VFC"]*LL, max = para["VFC"]*UL),
              # list(min = para["VKC"]*LL, max = para["VKC"]*UL),
              # list(min = para["VfilC"]*LL, max = para["VfilC"]*UL),
              list(min = para["VGC"]*LL, max = para["VGC"]*UL),
              # list(min = para["VPlasC"]*LL, max = para["VPlasC"]*UL),
              # list(min = para["Skinarea"]*LL, max = para["Skinarea"]*UL),
              # list(min = para["Skinthickness"]*LL, max = para["Skinthickness"]*UL),
              list(min = para["AbsPFOA"]*LL, max = para["AbsPFOA"]*UL)
)



## Create parameter matrix ##  
set.seed(1234)
params <- c("Htc", "Tmc", "Kt", "Free", "BW", "kurinec", "kbiliaryc", "kfaecesc", "kfil", "PL", "PF", "PK", "PSk", "PR", "PG", "VGC", "AbsPFOA") #"VLC", "VFC", "VKC","VfilC", "VGC", "VPlasC", "Skinarea", "Skinthickness",
length(params)==length(q)
x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, rep = 10)

dim(x$a) # the array of c(model evaluation, replication, parameters)

## Conduct simulation ##
out <- solve_fun(x, time=times, func = PBPKmodPFOA, initState = yini, outnames = outputs)

saveRDS(object=out, file="out_scaled.rds")
#out <- readRDS("out_scaled.rds")

## Output of the Uncertainty analysis ##
pdf("out_scaled.pdf")
pksim(out)
dev.off()

## Output from the sensitivity analysis ##
pdf("out_scaled.pdf")
plot(out)
dev.off()

ResultsSI <- as.data.frame(print(out["tSI"]))
ResultsSI$Times <- rownames(ResultsSI)
ResultsSI$Times <- as.numeric(as.character(ResultsSI$Times))
ResultsSI$Year <- ResultsSI$Times/(365*24)

write.xlsx(ResultsSI,
           file =file.path(newday,"ResultsSI.xlsx"),
           colNames = TRUE, borders = "rows"
)

write.xlsx(ResultsSI,
           file = "Results/2023-09-24/ResultsSI.xlsx",
           colNames = TRUE
)

check(out)
pdf("heat_check_CI_scaled.pdf")
heat_check(out, index = "CI")
dev.off()

pdf("heat_check_all.pdf")
heat_check(out, show.all = TRUE)
dev.off()
