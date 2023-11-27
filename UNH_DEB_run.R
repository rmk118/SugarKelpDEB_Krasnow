#Running DEB model on UNH outdoor growth expt
#Ruby Krasnow
#Last updated: Nov 27, 2023

#Import libraries
library(deSolve)
library(tidyverse)
library(lubridate)
library(patchwork)

#Required for model runs
source("SolveR_R.R")
source("KelpDEB_model.R")
source("./outdoorExpt/outdoorHOBO/outdoor_HOBO.R")
source("./outdoorExpt/outdoor_expt.R")