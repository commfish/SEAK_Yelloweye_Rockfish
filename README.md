# This is the age-structured yelloweye rockfish model for Southeast Alaska Outside Waters  
# DIRECTORY CONTENTS
<hr>  
1.  **master**   
      The ADMB stock assessment model for Southeast Alaska Waters in a single file containing:  
      * a switch for estimating or fixing natural mortality  
      * retrospective analysis code
      * parametric bootstrap
      * mcmc output
2.  **est_M**   
      The ADMB stock assessment model for Southeast Alaska Waters with estimated natural mortality  
      * Static structure to make it easier to follow
      * parametric bootstrap 
      * mcmc output
3.  **fixed_M**   
      The ADMB stock assessment model for Southeast Alaska Waters with fixed natural mortality  
      * Static structure to make it easier to follow
      * parametric bootstrap 
      * mcmc output 
4.  **model_output**  
      Output files from the est_M and fixed_M structures for comparison. Also includes output from MCMC runs for selected derived quantities and parametric bootstraps of the main free parameters.  
5.  **assessment_graphics.R**  
      A set of graphic code to examine model outputs. Include read-in of the MCMC derived quantities, but only displays uncertainty as a function of the .std file values. Requires  'globals.R' file (included here)

#SUMMARY
The ADMB stock assessment model in this directory was recommended for formal adoption at the November 2016 Plan Team meeting in Seattle. Accordingly, the current work is to prepare the model and output for formal review by the Plan Team and Science and Statistical Committee.

The model itself is a region-wide model. Older versions contained distinct models for each region. The PT has requested that results from a model structure in which natural mortality is estimated and from a model structure in which nartural mortality is fixed at Tier 4 assumptions (*M* = 0.026) be presented and discussed.  

# DATA SOURCES:  
1. Commercial fishery and bycatch in the Pacific halibut directed fishery (merged):  
    * Total annual catch
    * Age-composition 
    * CPUE

2. Recreational fishery:  
    * Total annual catch;  

3. ADF&G Submersible and ROV density survey:  
    * Density (adults and sub-adults per square kilometer)  

4. International Pacific Halibut longline survey:  
    * CPUE

# OBJECTIVE FUNCTION  
1. Likelihoods:  
    a. Total annual commercial catch  
    b. Total annual recreational catch  
    c. Density  
    d. Fishery CPUE  
    e. IPHC survey CPUE  

2. Penalties:  
    a. Annual recruitment deviations  
    b. Year 1 abundance deviations  
    c. Mean *F*  
    d. Mean sport fishing *F*
    e. Mean recruitment  
    f. Mean year 1 abundance  

For a complete description of model structure, data, and function, please review the 'Assessment of yelloweye rockfish in SE AK Outside Waters_Nov' document.**


#SPECIFIC MODEL NOTES:  

1. The '*master*' directory has the full model - it includes a switch in the control file for fixing or estimating natural mortality, and it also includes the retrospective analysis code and the modifications to the objective functions and other loops required for the '-retro' flag to run. That can be a bit complicated when trying to get a handle on the code and model initially, so the '*master*' directory has been branched into two sub-models: '*est_M*', and '*fixed_M*'. Beginning with the '*est_M*' code, it is a bit easier to wade through all of it. The '*fixed_M*' model simply sets natural mortality to the Tier 4 assumption that *M* = 0.026. Distinct directories also makes the initial examination and graphics calls for the assessment document a bit easier.  

2. The current model structure only has the region-wide data  in the .dat file. 

3. Much of the code, especially in the objective function, is written out long-hand; it could be streamlined by using 'dnorm' calls, removing explicit loops, etc. I left it this way to faciliate my own error checking when I was wrestling with some ROV likelihood errors as well as implementation of the retrospective analysis code. Although that is all now resolved, I don't have the time to go through and clean it all up. By all means, feel free to clean this up once you are comfortable with the code workings.  


4. The .dat files produced by the R scripts, especially the RUNFILE and admb_dat_file_creation.R script, write the square root of the sample size for the age-composition matrix. Running the model, however, implements a MANUAL sample size correction: the variances of the residuals in the age-composition data are compared to the assumed input sample sizes by assessing the standard deviation of normalized residuals (Breen at el. 2003). Under the assumption that the normalized residuals are normally distributed, the acceptable limit for SDNR values, following Francis (2011) is given as  

$max(sdnr)<[\chi_{0.95}^2/(m-1)]^{0.5}$

for which m = the number of years in the age-composition data set. If the SNDR for any given year exceeds this limit, the input age composition sample size vector is divided by the SDNR vector, and the model re-run with the revised sample sizes. The process is iteratively repeated until the target maximum SDNR value is reached. **This means that the model is initially run using the square root of the sample size as output by the R script, the resulting SDNR examined (in the model.rep file), and if any single element of the vector exceeds the limit, the sample sizes are divided by the SDNR values, AND THE RESULTING REVISED SAMPLE SIZE ARE MANUALLY INSERTED IN THE MODEL.DAT FILE, REPLACING THE SAMPLE SIZE VECTOR (COMMENT THAT LINE OUT).**   
Note that the original sample size has been commented out and replaced by the SDNR modified sample size on lines 97 and 98 in admb\yelloweye.dat, which is the current model template.  

5. Retrospective analyses, in which a single year of data is removed for ten successive years one at a time, is a standard stock assessment inclusion. The retrospective model forms are in the master directory. The retro_yrs calls are at the very head of the .tpl file, under 'COMMAND LINE ARGUMENT FOR -RETRO'. The function removes one year of model data for each retro_yrs count (i.e. retro_yrs = 1 removes one year of data, retro_yrs = 10 removes ten years of data), but continues to estimate the full time series of parameters and derived quantities. The goal is to determine whether any specific data point(s) are exerting signficant forcing on model outputs and estimates. Under normal circumstances, implenetation of the retrospective analysis would entail simply running the model with the commend line flag '-retro X', where X is 1 to 10. In the yelloweye model, however, the submarine/ROV survey data that serve to scale absolute abundance are not annually updated. As those data are not present for some years, we need a separate counter for stripping them relative to the main 'retro_yrs' counter.
I have set up a secondary counter 'retro-mod' in the retro code that is called by a second flag (-step) following the '-retro' flag in the command line. The proper arguments for each removed year are given in the tpl files as well as here:  
 

    -RETRO -STEP  flag values                   
                                                                        
      For the current model ending in 2015::                            
                                                                        
      Retrospective 1  (remove 1 year of data)   ->   -retro 1 -step 1  
      Retrospective 2  (remove 2 year of data)   ->   -retro 2 -step 1  
      Retrospective 3  (remove 3 year of data)   ->   -retro 3 -step 2  
      Retrospective 4  (remove 4 year of data)   ->   -retro 4 -step 3  
      Retrospective 5  (remove 5 year of data)   ->   -retro 5 -step 3  
      Retrospective 6  (remove 6 year of data)   ->   -retro 6 -step 3  
      Retrospective 7  (remove 7 year of data)   ->   -retro 7 -step 4  
      Retrospective 8  (remove 8 year of data)   ->   -retro 8 -step 4  
      Retrospective 9  (remove 9 year of data)   ->   -retro 9 -step 5  
      Retrospective 10 (remove 10 years of data) ->   -retro 10 -step 5 


Note that saving output from this will be done manually. It would not be difficult to write an R loop for this either. Note also that not all of these will converge to positive Hessian definites, and each year will be different, so the relevant graphics code will have to account for this each year. 
  
