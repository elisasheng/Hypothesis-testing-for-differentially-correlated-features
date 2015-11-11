##################################################################################
### This script shows how to use functions defined in exampleFunctions.r to obtain
### the Minimally Dysregulated (MD) set for an example dataset with 34 features. 
##################################################################################

source("exampleFunctions.r")

# Load datasets: two data.frames that represent two groups.
# Rows correspond to observations, columns correspond to features.
load("exampleData.Rdata")

# Estimate the MD Set for a fixed value of alpha, two ways:

# 1. this function returns the features in the MD set, or zero if the MD set is 
# empty. Computations are done serially.
MDset(data1,data2,alpha=0.001)  

# 2. this function is identical to MDset() but performs computations in parallel.
MDset2(data1,data2,alpha=0.001)  

