For the WRI project, a repo that shows how to select a minimum number of stations
that don't compromise model performance.

## File list

#### `00_create_simulated_dataset.R`

Creates the simulated data. The key part is where you can define `beta_daymet`, 
`beta_green`, and `beta_albedo` which essentially 
changes the relationship between predictor and the response. Otherwise this
isn't really important -- we have these data already, just needed something.

However, we do make several assumptions about the dataset:

1. Every monitor has a value for every predictor every day. So the predictor matrix
has `N_monitors * N_days` of rows.

1. There is a monitor ID column, which is an integer, that is `1:N_monitors`

1. There is a day ID column which is an integer and is `1:N_days`


#### `01_model_airtemp.R`

Creates the linear model for the predictors. uses bootstrapping to calculate
wider confidence intervals from using fewer than all stations. 

This isn't used to influence the simulated annealing but I suppose serves as the
Benchmark for the individual components since the feature coverage should
essentially be 0 error as k approaches N

#### `02_annealing_prep.R`

Prepares the dataset for annealing. This also contains an R-version of the 
scoring function


#### `03_find_subset.R`

Finally run the simulated annealing, and many times per monitor to build the 
curve.

#### `simann.f90`

Simulated annealing in FORTRAN
the SCORE function is defined by model performance 
Remember that for FORTAN the inputs needs to be perfect for it to work
so if its an integer it needs to have as.integer() in the argument

