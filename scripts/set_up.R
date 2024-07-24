
## load packages
library(INLA)
library(ggplot2)
library(cowplot)
library(ISOweek)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(forcats)


## read in data
data<-readRDS("data/Rt_variables.RDS")

## load shapefile
malaysia.1 <- readRDS("data/gadm36_MYS_1_sp.rds") 


## read in functions
source("R/functions.R")
