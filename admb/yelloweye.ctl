## -------------------------------------------------------------------------- ##
##
## THIS IS THE CONTROL FILE FOR YELLOWEYE ROCKFISH COMBINED MODEL
## 
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
#  DEBUG FLAG (0 = FALSE, 1 = TRUE)
#
	0
#
## -------------------------------------------------------------------------- ##

## -------------------------------------------------------------------------- ##
#  NBOOT for number of bootstraps to run (0 = NOT RUN)
#
	0
#
## -------------------------------------------------------------------------- ##


## —————————————————————————————————————————————————————————————————————————— ##
##                  DESIGN MATRIX FOR PARAMETER CONTROLS                      ##
##  Prior descriptions   Parameter values                                     ##                                                                            
##  -0 uniform           (0,0)                                                ##
##  -1 normal            (p1=mu,p2=sig)                                       ##
##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
##  -3 beta              (p1=alpha,p2=beta)                                   ##
##  -4 gamma             (p1=alpha,p2=beta)                                   ##
## —————————————————————————————————————————————————————————————————————————— ##
##  init   lower    upper    est   prior
## value   bound    bound    phz    type    p1      p2    # PARAMETER         ##
## —————————————————————————————————————————————————————————————————————————— ##
## - init_int n_theta
   10
#
## - initial population
   -1.05  -15.0     0.00      2      1    -3.64     1   # log_natural_mortality
    4.0     2.00    8.00      1      1     4.0      1   # log_mean_rec
    4.0     2.00    8.00      1      1     4.0      1   # log_mean_y1
    0.5     0.00    1.5       3      0      0       0   # sig1
#
# - fishing mortality
   -5.0   -10.0     0.0       2      1    -3.64     1   # log_avg_F   
   -5.0   -10.0     0.0       2      1    -3.64     1   # log_avg_Fs 
#
## - selectivity 
   -1.5    -3.00    0.0       3      0      0     0   # slope
    2.0     1.0     3.0       3      0      0     0   # 50_sel
#
## - catchability
   -2.0    -7.0     0.00      1      0      0     0   # comm_q
   -2.0    -7.0     0.00      1      0      0     0   # iphc_q  
#

## —————————————————————————————————————————————————————————————————————————— ##


## —————————————————————————————————————————————————————————————————————————— ##
##                        OTHER MISCELLANEOUS CONTROLS                        ##
## —————————————————————————————————————————————————————————————————————————— ##
## number of controls to read in.
   8
## Value    # - Description
## - fixed value
 1.0      # sigr        -fixed
 0.1      # sigma_sport -fixed
 0.05     # sigma_catch -fixed
 4        # ph_rec - recruitment deviations vector phase
 3        # ph_init - Year 1 deviations vector  phase
 2        # ph_Fdev - commercial fisheries deviations vector phase
 3        # ph_FdevS - sport fish deviations vector phase
 5        # ph_spr - spawner-recruit estimates phase

#EOF
42
