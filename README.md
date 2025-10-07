For the WRI project, a repo that shows how to select a minimum number of stations
that don't compromise model performance.

## File list

* `00_create_simulated_dataset.R`

Creates the simulated data. The key part of this is at the bottom where you 
can define `beta_green` and `beta_albedo` which essentially changes the relationship 
between either and the response. 

* `01_model_airtemp.R`

Creates the linear model for the predictors, and then does the extra step of
created a smoothed surface on a nearby grid using 2-d thin plate splines. Probably
don't need this part since we are just worried about monitors but essentially 
this is what Ian's code is doing.

* `02_annealing_prep.R`

Prepares the dataset for annealing, by creating monitor networks

* `03_find_subset.R`


* `simann.f90`

Simulated annealing in FORTRAN
the SCORE function is defined by model performance 


## Notes

* The big picture here is to run this for each city, because the number of 
monitors that are "representative" may change city by city. The 2nd task then
is using the GoogleAlpha to determine the characteristics of the _networks_ that
are chosen as representative. So essentially there are 2 tasks: tasks 1 - find 
representative stations using simulated annealing; task 2 - assess trends using 
Google things and maybe pick out network trends.

* The un-answered question here is whether the station is "representative", given
that we know that stations are more likely in some places than others (e.g., 
rich suburbs). I don't think that we can assess this without some external measure
of the model output, but we can assess similarity in predictor values and see
if the range etc is capture. But we can't assess where we don't have measured
air temperature.