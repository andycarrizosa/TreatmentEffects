# TreatmentEffects
This R Library re-codes what our team used to call "OneClarifyComb". The first big change 
is that the title of this function "TreatmentEffects" is now meant to be more intuitive and
less jargon-y with the title specifying exactly what the function will do.  That is, this function
allows us to easily estimate treatment effects of survey experiments. In order to produce these
estimates we rely on ["Clarify"](https://gking.harvard.edu/clarify) a parametric bootstrapping method
developed by Gary King. 

## Dependencies

**Treatment Effects** relies on several packages.  If this is your first time using TreatmentEffects
please run the following code in your R console before proceeding:

```
install.packages(c("stringr", "Hmisc", "ggplot2", "clarify", "tidyverse", "RColorBrewer","gtools",
"doParallel", "foreach", "remotes", "shiny"))
```

## Installing TreatmentEffects

TreatmentEffects can be installed as an R package using the following code:
```
library(remotes)
install_github("https://github.com/andycarrizosa/TreatmentEffects")
```
After that first installation it can be loaded as a regular R package for subsequent uses:

```
library(TreatmentEffects)
```

## Function Arguments

TreatmentEffects has the following arguments:

```
TreatmentEffects(data, treats = "", DV = "", model = c("linear", "logit"), 
controls = "", subgroups = "", sims = 10000, comb = "perm", clust=1) 
```

- **data:** is where you supply the dataframe to be used by the function
- **treats:** is the name of the variable that contains your treatment conditions.  
This variable must be categorical, and it can be either a string or a factor variable.
- **DV:**is the variable that you would like to treat as your dependent variable.  
This variable must be numeric and it must be either a continous or a binary variable (see **model**).
- **model:** Specifies the model you wish to run.  If the **DV** is continuous then specify "linear" in this argument,
and if the **DV** is binary then specify "logit" in this option
- **controls:** This argument recieves the variables you would like to include in your model as control variables.  This argument only accepts string vectors, 
and it can be a single item (for example controls="var1"), or it can be multiple items through concatenation (for example controls=c("var1", "var2", "var3")).
- **subgroups:** If left empty ("") the function finds the treatment effect for your full sample. But subgroups can be specified to find treatment effects within
different subgroups.  For example, if we specify subgroups="gender", then the function will estimate treatment effects for every category of that subgroup, such as "male",
"female", or "other".  Only one subgroup can be specified at a time.
- **sims:** This argument accepts the number of simulations you would like to perform in clarify's parametric bootstrapping. By default this value is set to 10000 iterations.
- **comb:** TreatmentEffects estimates effects comparing all treatment effects to the control condition, and comparing treatment effects against each other.  The **comb**
argument can accept two values. "perm" can be set if you want all permutations of treatment comparisons.  "comb" can be set if you only want combination comparisons
rather than permutation comparisons.  By default this argument is set to "perm".
- **clust:** clust specifies the number of clusters to use for estimation.  By default it is set at 1, which 
means estimations will take place by a single cluster. Raising this number can quicken estimation, but make sure
to not rase this number to higher than the total number of threads of your computer.





