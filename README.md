For the WRI project, a repo that shows how to select a minimum number of stations
that don't compromise model performance.

## File list

* `00_create_simulated_dataset.R`

Creates the simulated data. The key part is where you can define `beta_daymet`, 
`beta_green`, and `beta_albedo` which essentially 
changes the relationship between predictor and the response. Otherwise this
isn't really important -- we have these data already, just needed something.

* `01_model_airtemp.R`

Creates the linear model for the predictors. uses bootstrapping to calculate
wider confidence intervals from using fewer than all stations

* `02_annealing_prep.R`

Prepares the dataset for annealing

* `03_find_subset.R`



* `simann.f90`

Simulated annealing in FORTRAN
the SCORE function is defined by model performance 
Remember that for FORTAN the inputs needs to be perfect for it to work
so if its an integer it needs to have as.integer() in the argument

## Chad ToDo

* Right now you are just comparing the values directly, you aren't using the model
  I guess eventually you will want to see how the modeled values perform?
  Like does the subset affect your ability to MODEL, not just cover the datapoints themselves
  So you need to change this up slightly
  And do you assess this based on the model coefficients?
  
* Ok there are two additional thoughts here - you can't really do the RMSE alone, because you
  might just get lucky and get two points where it predicts really well. probably what you need
  is estimate the coefficients using the subset and then predict at all places you have data. 
  Again this doesn't solve the issue of over-representation in some areas, but again, we can't solve that.
  
  Another way to think about this is the bias-variance trade-off -- in some ways if you 
  leave the part about the distribution of the inputs in, you'll also account for this.
  
  So i think the ultimate solution is a scoring function that includes both the
  distribution of the predictors and how well the beta coefficients do at predicting at 
  every monitor, with gamma tuning params in front of each that I can turn on or off as 
  necessary
  
  Hmm another thing to consider is that the current model doesn't contain an 
  intercept for station -- but it probably doesn't need it, since all of the
  geographic co-variates probably uniquely identify it anyway. But a good
  thing to confirm
  
  Another issue you need to worry about is the units -- the units of each of 
  these will have an impact on how the penalty looks. the units of the predictors
  and the units of the MSE for both z1 and z2. Not sure exactly how to solve this
  but the best idea would probably be to (a) z-score the predictor matrix and then
  (b) test a range of z1 and z2 and see what answers are robust to that. 
  

## Notes

* The big picture here is to run this for each city, because the number of 
monitors that are "representative" may change city by city. The 2nd task then
is using the GoogleAlpha to determine the characteristics of the _networks_ that
are chosen as representative. So essentially there are 2 tasks: 
  
  * task 1 - find representative stations using simulated annealing
  
  * task 2 - assess trends using Google things and maybe pick out network trends.

* The un-answered question here is whether the station is "representative", given
that we know that stations are more likely in some places than others (e.g., 
rich suburbs). I don't think that we can assess this without some external measure
of the model output, but we can assess similarity in predictor values and see
if the range etc is captured. But we can't assess where we don't have measured
air temperature.

* I also have code to create a smoothed surface on a nearby grid using 2-d thin plate splines.
Probably don't need this part since we are just worried about monitors but essentially 
this is what Ian's code is doing.

* This annealing problem is hard because you are actually optimizing on 2 things;
  the number of stations, and which ones are being chosen. 
  Would this be the same as just running them independently by the number of stations
  and then comparing the minimum by station?

  * answer: seems like "joint annealing" is a thing. so the "moves" can be either
  swap around, add stations, or remove stations, and the probability of this can be adaptive
  OR you can just do this for 1 N at a time. 
  
  * Is it useful to know what the error is by station #? I think so. You can add the 
above as a discussion point.  

  * The other reason  you don't want to do joint is the choice of penalty function
  will be have a huge impact on the results. so actually doing it k by k allows us 
  to observe the surface, and users can choose more what they want (80%, 90%) etc.

  * so solution is 1 N at a time, in parallel: Declaring all dummy variables as 
  THREADPRIVATE seems to help. see [here](https://stackoverflow.com/questions/39196532/calling-subroutine-in-parallel-environment)
