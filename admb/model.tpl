// --------------------------------------------------------------------------- //
//               Age-structured model for Yelloweye rockfish                   //
//                    Outside Southeast Alaskan Waters                         //
//                                                                             //
//                                 VERSION C.1                                 //
//                                  March 2016                                 //
//                                                                             //
//                                  AUTHORS                                    //
//                               Kray Van Kirk                                 //
//                           kray.vankirk@alaska.gov                           //
//                                                                             //
//                 Includes code from Pacific Ocean Perch model                //
//                               Dana Hanselman                                //
//                           dana.hanselman@noaa.gov                           //
//                                                                             //
//                   Scripts in FINAL and GLOBAL sections                      //
//                         for run time statistics                             //
//                              Steven Martell                                 //
//                          martell.steve@gmail.com                            //
//                                                                             //
//                                                                             //
//                                                                             //
// --------------------------------------------------------------------------- //


DATA_SECTION

 
  // |--------------------------------------------------------------------------|
  // |MCMC OUTPUT FILE
  // |--------------------------------------------------------------------------|

     !!CLASS ofstream evalout("evalout.prj");
     !!time(&start);



  // |--------------------------------------------------------------------------|
  // | STRINGS FOR INPUT FILES                                                  
  // |--------------------------------------------------------------------------|
  // |
  // | DataFile               : data to condition the assessment model    
  // | ControlFile            : controls for years, phases, and block options 

     init_adstring DataFile;      
     init_adstring ControlFile;    

  // | BaseFileName           : file prefix used for all  model output
  // | ReportFileName         : file name to which report file is printed

     !! BaseFileName = stripExtension(DataFile);  
     !! ReportFileName = BaseFileName + adstring(".rep");

     !! cout<<"You are modeling the "<<BaseFileName<<" stock of yelloweye rockfish"<<endl;
     !! cout<<""<<endl;



  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM CONTROL FILE                                             
  // |--------------------------------------------------------------------------|
  // | 
  // | DEBUG_FLAG  : Boolean Flag used for manual debugging
  // | nboot       : number of bootstrap draws to run. Refer to
  // |                PARAMETRIC_BOOTSTRAP section in FINAL_CALCS

     !! ad_comm::change_datafile_name(ControlFile);

     init_int DEBUG_FLAG;
     init_int nboot;

 // |-------------------------------------------------------------------------|
 // | DESIGN MATRIX FOR PARAMETER CONTROLS                                    |
 // |-------------------------------------------------------------------------|
 // | - theta_DM -> theta is a vector of estimated parameters.

    init_int    n_theta;
    init_matrix theta_DM(1,n_theta,1,7);
    vector      theta_ival(1,n_theta);
    vector      theta_lb(1,n_theta);
    vector      theta_ub(1,n_theta);
    ivector     theta_phz(1,n_theta);
    ivector     theta_iprior(1,n_theta);
    vector      theta_p1(1,n_theta);
    vector      theta_p2(1,n_theta);
    
    !! theta_ival = column(theta_DM,1);
    !! theta_lb  = column(theta_DM,2);
    !! theta_ub  = column(theta_DM,3);
    !! theta_phz = ivector(column(theta_DM,4));
    !! theta_iprior = ivector(column(theta_DM,5));
    !! theta_p1 = column(theta_DM,6);
    !! theta_p2 = column(theta_DM,7);


// |---------------------------------------------------------------------------|
// | Miscellaneous Controls
// |---------------------------------------------------------------------------|
// | nMiscCont  » Number of controls to read in.
// | dMiscCont  » Vector of miscellaneous controls,

   init_int nMiscCont;
   init_vector dMiscCont(1,nMiscCont);

   number sigr;
   number sigma_sport;
   number sigma_catch;
   int     ph_rec;
   int     ph_init;
   int     ph_Fdev;
   int     ph_FdevS;
   int     ph_spr;
   

  // |--------------------------------------------------------------------------|
  // | END OF FILE MARKER                             
  // |--------------------------------------------------------------------------|

     init_number eof_ctl
     
     LOCAL_CALCS

          sigr        = dMiscCont(1);
          sigma_sport = dMiscCont(2);
          sigma_catch = dMiscCont(3);
          ph_rec      = dMiscCont(4);
          ph_init     = dMiscCont(5);
          ph_Fdev     = dMiscCont(6);
          ph_FdevS    = dMiscCont(7);
          ph_spr      = dMiscCont(8);
   

    if(eof_ctl==42) cout << BaseFileName<<".ctl has been read correctly!"<<endl;
      else 
        {    
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|      Red alert! Captain to bridge! The .ctl file is compromised!     |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof_ctl<<", but the file *should* end with 42       |"<<endl;
         cout <<"| Please check the .ctl file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
         exit(1); 
    }

 END_CALCS



  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM DATA FILE                                                
  // |--------------------------------------------------------------------------|
  // | This calls the data file

     !! ad_comm::change_datafile_name(DataFile);



  // |--------------------------------------------------------------------------|
  // | MODEL STRUCTURAL DIMENSIONS  - CONSTANT OVER ALL REGIONS                                            
  // |--------------------------------------------------------------------------|
  // | 
  // | nages          -> number of ages
  // | styr           -> data start year
  // | endyr          -> data end year
  // | recage         -> age at recruitment into the model
  // | spawn_fract    -> spawning month (set to May [5])
  // | srv_fract      -> survey month (set to August [8])
  // | fsh_fract      -> month when fishery takes place (set to January [1])

     init_int      nages
     init_int      styr
     init_int      endyr
     init_int      recage
     init_number   spawn_fract;
     init_number   srv_fract;
     init_number   fshy_fract;



  // |--------------------------------------------------------------------------|
  // | CALCULATED CONTAINERS                                                    
  // |--------------------------------------------------------------------------|
  // |
  // | 1) styr_rec       -> recruitment start year
  // | 2) endyr_rec      -> recruitment end year
  // |
  // | 3) endyr_sp       -> last year of estimated spawn
  // | 4) nyrs           -> number of model years
  // |
  // | 5) yy             -> vector of specific model years
  // | 6) aa             -> vector of specific model ages
  // | 7) myy            -> vector for natural mortality cubic spline
  // |
  // | 8) nrecs_est      -> number of recruitment estimates; not currently used
  // | 9) wt_mature      -> weight of mature females (for spawning biomass)
  // |                      for each region and then in toto

     int       styr_rec;
     int       styr_sp;
     
     int       endyr_sp;
     int       nyrs;
     
     ivector   yy(styr,endyr);
     ivector   aa(1,nages);
     ivector   myy(styr+1,endyr);
     
     int       nrecs_est;
     vector    wt_mature(1,nages);




  // |--------------------------------------------------------------------------|
  // | PHYSIOLOGY                                                               
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) p_mature       -> proportion of females mature-at-age
  // | 5) wt             -> mean weight-at-age, regions (C,S,E,Global)
  // | 6) morphology     -> appearance for counting in the ROV survey

     init_vector   p_mature(1,nages)
     init_vector   wt(1,nages)
     init_vector   morphology(1,nages)


  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL LONGLINE COMMERCIAL CATCH                                   
  // |--------------------------------------------------------------------------|
  // | 
  // | obs_catch       -> total annual commercial catch and halibut bycatch
  // |                  
  // | nyrs_sport      -> number of years for the sport fishery
  // | yrs_sport       -> specific years of sport fishery
  // | obs_sportcatch  -> total annual sportfish catch

     init_vector   obs_catch(styr,endyr)

     init_int      nyrs_sport
     init_ivector  yrs_sport(1,nyrs_sport)
     init_vector   obs_sportcatch(1,nyrs_sport)


  // |--------------------------------------------------------------------------|
  // | SUBMERSIBLE / ROV SURVEYS - for each region and Global            
  // |--------------------------------------------------------------------------|
  // | 
  // | nyrs_srv        -> number of years of submersible / ROV surveys
  // | yrs_srv         -> specific survey years
  // | obs_srv_biom    -> mean biomass estimate
  // | obs_srv_se      -> standard error of DISTANCE
  // | area_skm        -> area of region in thousands of square kilometers
  // |
  // | l_var           -> vector for log transform of variance for -L

     init_int      nyrs_srv
     init_ivector  yrs_srv(1,nyrs_srv)
     init_vector   obs_srv_biom(1,nyrs_srv)
     init_vector   obs_srv_se(1,nyrs_srv)  
     init_number   area_skm

     vector l_var(1,nyrs_srv);


  // |--------------------------------------------------------------------------|
  // | AGE COMPOSITION LONGLINE FISHERY AND HALIBUT BYCATCH regional and Global                     
  // |--------------------------------------------------------------------------|
  // | 
  // | nyrs_fish_age       -> number of years of fishery age comps
  // | yrs_fish_age        -> specific years age comps
  // | nsamples_fish_age   -> some measure of relative sample size
  // | oac_fish            -> age comp data as %
  // | nmulti_fish_age     -> multinomial sample size
  // |
  // | Relative sample size = square root of number of fish aged to scale data impact

     init_int      nyrs_fish_age 
     init_ivector  yrs_fish_age(1,nyrs_fish_age)
     init_vector   nsamples_fish_age(1,nyrs_fish_age)  
     init_matrix   oac_fish(1,nyrs_fish_age,1,nages)
     
     vector        nmulti_fish_age(1,nyrs_fish_age) 



  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL LONGLINE COMMERCIAL CPUE - regional and Global                                 
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_cpue     -> number of years of CPUE estimates
  // | 2) yrs_cpue      -> specific years CPUE estimated
  // | 3) obs_cpue      -> mean CPUE by year
  // | 4) var_cpue      -> variance of CPUE
  // |
  // | 5) l_var_cpue    -> vector for variance log-transform

     init_int      nyrs_cpue    
     init_ivector  yrs_cpue(1,nyrs_cpue)  
     init_vector   obs_cpue(1,nyrs_cpue) 
     init_vector   var_cpue(1,nyrs_cpue)
     
     vector        l_var_cpue(1,nyrs_cpue)


  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL IPHC SURVEY CPUE                                 
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_cpue_iphc     -> number of years of CPUE estimates
  // | 2) yrs_cpue_iphc      -> specific years CPUE estimated 
  // | 3) obs_cpue_iphc      -> mean CPUE by year
  // | 4) var_cpue_iphc      -> variance of CPUE
  // |
  // | 5) l_var_iphc         -> vector for variance log-transform
  
     init_int      nyrs_cpue_iphc
     init_ivector  yrs_cpue_iphc(1,nyrs_cpue_iphc)
     init_vector   obs_cpue_iphc(1,nyrs_cpue_iphc)
     init_vector   var_cpue_iphc(1,nyrs_cpue_iphc)
     
     vector        l_var_iphc(1,nyrs_cpue_iphc)


  // |--------------------------------------------------------------------------|
  // | TRANSITION MATRICES                               
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) ageage                    -> aging error matrix

     init_matrix   ageage(1,nages,1,nages)



  // |--------------------------------------------------------------------------|
  // | END OF FILE MARKER                             
  // |--------------------------------------------------------------------------|

     init_int eof_dat;
     


  // |--------------------------------------------------------------------------|
  // | MINUTIAE                              
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) integer definitions 
  // | 2) 'small number' for adding to log functions
  // | 3) offsets for multinomials to allow perfect fit = 0
  
     int iyr
     int i
     int j
     int k
  
     number oo;

     number  offset;



 LOCAL_CALCS

    if(eof_dat==42) cout << BaseFileName<<".dat has been read correctly!"<<endl;
    else 
    {       
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|   ** Red alert! Captain to bridge! The .dat file is compromised! **  |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof_dat<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .dat file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }


    nyrs        = endyr - styr + 1;              // Number of years for model run
    spawn_fract = (spawn_fract - 1) / 12;        // Fraction of year at which spawning occurs
    srv_fract   = (srv_fract - 1)   / 12;        // Fraction of year at which survey occurs
    fshy_fract  = (fshy_fract - 1)  / 12;        // Fraction of year at which fishery occurs

    for(j=1; j<=nages; j++)
      {
         wt_mature(j)   = (wt(j) * p_mature(j))/2;      // Assumption of equal sex division
      }


    yy.fill_seqadd(styr,1);
    aa.fill_seqadd(recage,1);
    myy.fill_seqadd(styr+1,1);

    oo = 0.000001;

    nmulti_fish_age = nsamples_fish_age;


       // |---------------------------------------------------------------------|
       // | OFFSETS
       // |---------------------------------------------------------------------|
       // | Calculate "offset" for multinomials: age-comp for fishery and 
       // | halibut fishery bycatch COMBINED
       // | "Offset" value lets the multinomial likelihood equal zero when the 
       // | observed and predicted are equal as in Fournier (1990) "robustifies"

       offset = 0.0;


       for (i=1; i<=nyrs_fish_age; i++)
        {
          offset -= nmulti_fish_age(i) *((oac_fish(i) + oo)
                       *log(oac_fish(i) + oo)); 
        }


 END_CALCS


INITIALIZATION_SECTION

	theta theta_ival;



PARAMETER_SECTION

  // |--------------------------------------------------------------------------|
  // | POPULATION PARAMETERS
  // |--------------------------------------------------------------------------|
  // | - theta(1)  -> log natural mortality
  // | - theta(2)  -> log initial mean recruitment age-8
  // | - theta(3)  -> log average age-8 recruitment from styr to endyr
  // | - theta(4)  -> variance of log_mean_y1
  // | - theta(5)  -> log of mean F commercial + halibut bycatch
  // | - theta(6)  -> log of mean F sport fish
  // | - theta(7)  -> alpha (selectivity)
  // | - theta(8)  -> beta (selectivity)
  // | - theta(9)  -> delta (selectivity)
  // | - theta(10) -> gamma (selectivity)
  // | - theta(11) -> q catchability (commercial)
  // | - theta(12) -> q catchability (iphc)

     init_bounded_number_vector theta(1,n_theta,theta_lb,theta_ub,theta_phz);
     number log_natural_mortality;
     number log_mean_rec;
     number log_mean_y1;
     number sig1;
     number log_avg_F;
     number log_avg_Fs;
     number fish_sel_slope;
     number fish_sel_a50;
     number log_q1;
     number log_q2;


  // |-------------------------------------------------------------------------|
  // | STOCK-RECRUITMENT                                                       
  // |-------------------------------------------------------------------------|
  // | log_rec_dev     -> Vector of annual deviations from mean recruitment
  // | init_pop        -> Vector of Year 1 deviations from mean abundance
  // |
  // | recruitment     -> SD report vector for recruitment

     init_bounded_vector      log_rec_dev(styr,endyr,-5.,5.,ph_rec);
     init_bounded_vector      init_pop(1,nages-1,-5,5,ph_init);

     sdreport_vector          recruitment(styr,endyr);
    

  // |-------------------------------------------------------------------------|
  // | GEAR SELECTIVITY AND CATCHABILITY                                       
  // |-------------------------------------------------------------------------|
  // | fish_sel           -> commercial fishery selectivity 
  // | q1                 -> Exponentiated
  // | q2                 -> Exponentiated
  // |

     vector     fish_sel(1,nages);					

     number     q1
     number     q2
				

  // |-------------------------------------------------------------------------|
  // | FISHING MORTALITY                                                       
  // |-------------------------------------------------------------------------|
  // | log_F_devs         -> deviation vector about mean F
  // | log_F_devs_sport   -> deviation vector about mean sport F
  // |
  // | Fmort              -> commercial fishing mortality by year
  // | Fmort2             -> sport fishing mortality by year
  // |
  // | F                  -> commercial fishing mortality by year and age
  // | F2                 -> sport fishing mortality by year and age

     init_bounded_dev_vector  log_F_devs(styr,endyr,-15.0,15.0,ph_Fdev); 
     init_bounded_dev_vector  log_F_devs_sport(1,nyrs_sport,-15.0,15.0,ph_FdevS);
     
     vector                   Fmort(styr,endyr); 
     vector                   Fmort2(styr,endyr);
     
     matrix                   F(styr,endyr,1,nages); 
     matrix                   F2(styr,endyr,1,nages);


  // |-------------------------------------------------------------------------|
  // | MORTALITY                                                               
  // |-------------------------------------------------------------------------|
  // | M                 -> expnentiated log natural mortality
  // |
  // | Z                 -> Matrix of total mortality by year and age
  // | S                 -> Matrix of total survival by year and age

     number    M;
 
     matrix    Z(styr,endyr,1,nages);
     matrix    S(styr,endyr,1,nages);


  // |-------------------------------------------------------------------------|
  // | POPULATION VECTORS AND MATRICES                                         
  // |-------------------------------------------------------------------------|
  // | natage           -> numbers at age
  // | batage           -> biomass at age
  // | spbatage         -> spawning biomass at age
  // |
  // | spawn_biomass    -> total annual spawning biomass
  // | tot_biomass      -> total annual biomass
  // | tot_N            -> total annual abundance

     matrix                natage(styr,endyr,1,nages);
     matrix                batage(styr,endyr,1,nages);

     sdreport_vector       spawn_biom(styr,endyr);
     sdreport_vector       tot_biomass(styr,endyr);
     sdreport_vector       tot_N(styr,endyr);


  // |-------------------------------------------------------------------------|
  // | CATCH AND SURVEY VECTORS AND MATRICES                                   
  // |-------------------------------------------------------------------------|
  // | catage           -> commercal catch at age
  // | sportcatage      -> sport catch at age (yah, rite! Don't even tell me...)
  // |              
  // | pred_catch       -> total annual commercial catch
  // | pred_sportcatch  -> total annual sport catch
  // |
  // | pred_srv         -> total survey density (ROV/sub years only)
  // | pred_dens        -> total annual density
  // |
  // | pred_cpue        -> commercial fishery CPUE
  // | pred_cpue_iphc   -> IPHC survey CPUE
  // |
  // | age_comp         -> commercal fishery proportions at age
  // |
  // | expl_rate        -> exploitation rate
  // |
  // | res_age          -> catch-age residuals for graphics

     matrix                catage(styr,endyr,1,nages) 
     matrix                sportcatage(1,nyrs_sport,1,nages)
     
     vector                pred_catch(styr,endyr)
     vector                pred_sportcatch(1,nyrs_sport)

     sdreport_vector       pred_srv(1,nyrs_srv)
     sdreport_vector       pred_dens(styr,endyr)

     sdreport_vector       pred_cpue(1,nyrs_cpue)
     sdreport_vector       pred_cpue_iphc(1,nyrs_cpue_iphc)

     matrix                age_comp(1,nyrs_fish_age,1,nages)

     vector                expl_rate(styr,endyr)

     matrix                res_age(1,nyrs_fish_age,1,nages);
     

  // |-------------------------------------------------------------------------|
  // | DERIVED STATISTICAL CATCH VECTORS                                       |
  // |-------------------------------------------------------------------------|
  // | effn_fish_age      -> effective sample size, fishery removals          
  // | sdnr_fish_age      -> SDNR, commercial fishery
  // | 
  // | FROM HANSELMAN ET AL. 2013 GOA SABLEFISH STOCK ASSESSMEN
  // | 'SDNR' = standardized deviation of normalized results'
  // | * related to root mean squared error (RMSE)
  // | * compares variance of model residuals to variance of data set
  // | * SDNR = 1 >> model fit to data matches assumed data variance

     vector      effn_fish_age(1,nyrs_fish_age)	
     vector      sdnr_fish_age(1,nyrs_fish_age)


  // |-------------------------------------------------------------------------|
  // | PROJECTED POPULATIONS, SPR AND F LEVELS                                 |
  // |-------------------------------------------------------------------------|
  // | mF                 -> vector of potential F levels
  // | SBx                -> spawning biomass from mF            
  // | Bx                 -> biomass from mF
  // |
  // | Nspr               -> Number of projected spawners
  // | N0                 -> Number projected with F = 0

     init_bounded_vector mF(1,9,0.01,1,ph_spr)
     vector SBx(1,10)
     vector Bx(1,10)
					       
     matrix Nspr(1,10,1,nages)			
     matrix N0(1,10,1,nages)


  // |---------------------------------------------------------------------|
  // | PROJECTIONS @ F_ABC AND F_OFL
  // |---------------------------------------------------------------------|
  // | N_ABC              -> Projected N @ ABC
  // | N_OFL              -> Projected N @ OFL
  // |
  // | FABC               -> F vector @ ABC
  // | FOFL               -> F vector @ OFL
  // | ZABC               -> Z vevtor @ ABC
  // | ZOFL               -> Z vector @ OFL
  // |
  // | B40                -> B40             
  // | ABC                -> Allowable Biological Catch for next year
  // | OFL                -> Overfishing level for next year
  // | F55
  // | F50
  // |
  // | stdev_rec          -> st. dev. of estimated recruitment for
  // |                       variation in MCMC draws
  // | 
  // | catage_ABC         -> Catch at age ABC
  // | catage_OFL         -> Catch at age OFL
  // | catch_ABC          -> Total ABC catch (wt)
  // | catch_OFL          -> Total OFL catch (wt)
  // |
  // | spawn_biom_ABC     -> ABC spawning biomass
  // | spawn_biom_OFL     -> OFL spawning biomass
  // | tot_biom_ABC       -> ABC total biomass
          
     matrix N_ABC(endyr+1,endyr+15,1,nages);
     matrix N_OFL(endyr+1,endyr+15,1,nages);
     
     vector FABC(1,nages);
     vector FOFL(1,nages);
     vector ZABC(1,nages);
     vector ZOFL(1,nages);

     number stdev_rec;

     sdreport_number B40;
     sdreport_number ABC;
     sdreport_number OFL;
     number F65;
     number F60;

     matrix catage_ABC(endyr+1,endyr+15,1,nages);
     matrix catage_OFL(endyr+1,endyr+15,1,nages);
     vector catch_ABC(endyr+1,endyr+15);
     vector catch_OFL(endyr+1,endyr+15);

     sdreport_vector spawn_biom_ABC(endyr+1,endyr+15);
     sdreport_vector spawn_biom_OFL(endyr+1,endyr+15);
     sdreport_vector tot_biom_ABC(endyr+1,endyr+15);


  // |-------------------------------------------------------------------------|
  // | OBJECTIVE FUNCTION COMPONENS                                           
  // |-------------------------------------------------------------------------|
  // | c_catch_like         -> commercial catch
  // | s_catch_like         -> sport catch
  // |
  // | surv_like            -> survey biomass & CPUEs
  // | age_like             -> age composition          
  // | rec_like             -> recruitment stability
  // |              
  // | penalties            -> Mean rec, mean year 1, mean F, mean F2, M
  // |                         and deviation vectors for rec, Y1
  // | sprpen               -> spawner-recruit F values
  // | 
  // | Like                 -> Summed likelihoods
  // | obj_fun              -> Yo... I really need to say it?

     number c_catch_like;
     number s_catch_like;
     vector surv_like(1,3);
     number age_like;

     vector penalties(1,7);
     number sprpen
     
     number Like;	
     objective_function_value obj_fun;

        
PROCEDURE_SECTION
     initializeModelParameters();
 
     Selectivity();
            if(DEBUG_FLAG == 1) cout<<"**Selectivity**"<<endl;
            
     Mortality();
            if(DEBUG_FLAG == 1) cout<<"**Mortality**"<<endl;
            
     Abundance();
            if(DEBUG_FLAG == 1) cout<<"**Abundance**"<<endl;
            
     Catch();
            if(DEBUG_FLAG == 1) cout<<"**Catch**"<<endl;
            
     Predicted();
            if(DEBUG_FLAG == 1) cout<<"**Predicted**"<<endl;
            
     Population_Summaries();
            if(DEBUG_FLAG == 1) cout<<"**Population**"<<endl; 
     
     if (last_phase())
       {
         spr_rates();
     	    if(DEBUG_FLAG == 1) cout<<"**spr**"<<endl;
                
           Population_Projection();
            if(DEBUG_FLAG == 1) cout<<"**Projection**"<<endl; 
       }
       
     Objective_Function();
     


     if (mceval_phase())
       {
         writePosteriorSamples();

         evalout<<theta<<" "<<endl;
       }


FUNCTION void writePosteriorSamples()
	/**
	- This function is only envoked when the -mceval
		command line option is implemented.
	*/
	//if(nf==1){
		//ofstream ofs("ssb.ps");
	//}
	ofstream ofs0("bm.ps",ios::app);
	ofs0<<tot_biomass<<endl;
	ofstream ofs1("sbm.ps",ios::app);
	ofs1<<spawn_biom<<endl;
	ofstream ofs2("rec.ps",ios::app);
	ofs2<<recruitment<<endl;
	ofstream ofs3("dens.ps",ios::app);
	ofs3<<pred_dens<<endl;
	ofstream ofs4("catch_ABC.ps",ios::app);
	ofs4<<catch_ABC<<endl;
	ofstream ofs5("catch_OFL.ps",ios::app);
	ofs5<<catch_OFL<<endl;
	ofstream ofs8("spawn_biom_ABC.ps",ios::app);
	ofs8<<spawn_biom_ABC<<endl;
	ofstream ofs9("spawn_biom_OFL.ps",ios::app);
	ofs9<<spawn_biom_OFL<<endl;

       

FUNCTION void initializeModelParameters()
	//fpen = 0;

	log_natural_mortality = theta(1);
	log_mean_rec          = theta(2);
	log_mean_y1           = theta(3);
	sig1                  = theta(4);
	log_avg_F             = theta(5);
	log_avg_Fs            = theta(6);
    	fish_sel_slope        = theta(7);
    	fish_sel_a50          = theta(8);
        log_q1                = theta(9);
        log_q2                = theta(10);



FUNCTION Selectivity
  fish_sel.initialize();

  // |--------------------------------------------------------------------------|
  // | SELECTIVITY SWITCH
  // |--------------------------------------------------------------------------|
  // | It is assumed that selectivities for the commercial longline fishery
  // | and the bycatch in the halibut longline fishery are equivalent.

     for (j=1;j<=nages;j++)
       {
         fish_sel(j)  = 1/(1+mfexp(-mfexp(fish_sel_slope)  *
                          (j-mfexp(fish_sel_a50))));
       }
               
     fish_sel   =  fish_sel  / max(fish_sel);   //Scale to 1
               

FUNCTION Mortality
  Fmort.initialize();
  Fmort2.initialize();
  F.initialize();
  F2.initialize();
  Z.initialize();
  S.initialize();


  // |--------------------------------------------------------------------------|
  // | MORTALITY 
  // |--------------------------------------------------------------------------|
  // | Note here that we assume the same selectivity for the commercial
  // | fishery and the sport fishery. This is likely wrong, but we have
  // | no sport fish age data and only a few years of catch data

     M = mfexp(log_natural_mortality);


  // |----------------------------|
  // | ANNUAL FISHERY DEVIATIONS
  // |----------------------------|

     Fmort   =  mfexp(log_avg_F +  log_F_devs);
            
     if(DEBUG_FLAG == 1) cout<<"Fmort"<<endl;


  // |------------------------------|
  // | ANNUAL SPORTFISH DEVS (2006+)
  // |------------------------------|
        
     for(j=styr; j<=endyr; j++)
       {
         if(j<2006)
           {
             Fmort2(j) = 0;
           }
         else
           {
              Fmort2(j)  =  mfexp(log_avg_Fs + log_F_devs_sport(j-2005));
           }
       }
          

     if(DEBUG_FLAG == 1) cout<<"Sport"<<endl;


  // |-------------------|
  // | Fishing Mortality
  // |-------------------|

     for (j=styr; j<=endyr; j++)
       {
         F(j)  = Fmort(j)  * fish_sel;
         F2(j) = Fmort2(j) * fish_sel;  
       }
          
     if(DEBUG_FLAG == 1) cout<<"F"<<endl;
     
      
  // |-----------------|
  // | Total Mortality
  // |-----------------|

     Z = F + F2 + M;
     S = mfexp(-1.0*Z);

     if(DEBUG_FLAG == 1) cout<<"M"<<endl;
     
         

FUNCTION Abundance
  natage.initialize();
  
  // |---------------------------|
  // | Year 1 recruitment
  // |---------------------------|
  
     natage(styr,1)= mfexp(log_rec_dev(styr) + log_mean_rec + (sigr*sigr)/2)/2;

  // |---------------------------|
  // | Year 1 ages 9 - 75+
  // |---------------------------|

     for(j=1;j<=nages-1;j++)
       {
         natage(styr,j+1)=mfexp(log_mean_y1 + init_pop(j) + (sig1*sig1)/2)/2;
       }
       
  // |---------------------------|
  // | Years 2+
  // |---------------------------|
  
     for(i=styr+1; i<=endyr; i++) // 'i' loop
       {
         natage(i,1) = mfexp(log_rec_dev(i) + log_mean_rec + (sigr*sigr)/2)/2;            

           for(j=2; j<=nages; j++)
             {
               natage(i,j)  = natage(i-1,j-1)*S(i-1,j-1); 
             }                                                                         

  // |---------------------------|
  // | Years 2+ plus class
  // |---------------------------|

         natage(i,nages)    += natage(i,nages)*S(i,nages);

       } // close 'i' loop


     if(DEBUG_FLAG == 1) cout<<"Years plus"<<endl;
        


FUNCTION Catch
  catage.initialize();
  sportcatage.initialize();
  pred_catch.initialize();
  pred_sportcatch.initialize();

  // |------------------|
  // | COMMERCIAL CATCH
  // |------------------|

     for (i=styr; i<=endyr; i++)
       {
         catage(i) = elem_div(elem_prod(elem_prod(natage(i),F(i)),
                               (1.-S(i))),Z(i));
         pred_catch(i) =  catage(i)*wt;
       }

  // |-------------|
  // | SPORT CATCH
  // |-------------|
    
     for (i=1; i<=nyrs_sport; i++)
       {
         sportcatage(i) = elem_div(elem_prod(elem_prod(natage(yrs_sport(i)),
                           F2(yrs_sport(i))),(1.-S(yrs_sport(i)))),Z(yrs_sport(i)));
         pred_sportcatch(i) =  (sportcatage(i)*wt);
       }


     if(DEBUG_FLAG == 1) cout<<"Sport catch"<<endl;


FUNCTION Predicted
  pred_cpue.initialize();
  pred_cpue_iphc.initialize();
  pred_srv.initialize();
  pred_dens.initialize();
  age_comp.initialize();
  effn_fish_age.initialize();
  sdnr_fish_age.initialize();

  q1 = mfexp(log_q1);
  q2 = mfexp(log_q2);

     for (i=1;i<=nyrs_srv;i++)
       {
         pred_srv(i) =  natage(yrs_srv(i)) * morphology;
         pred_srv(i) /= area_skm;
       }


     // all years
     
     for(i=styr; i<=endyr; i++)
       {
         pred_dens(i) = natage(i) * morphology;
         pred_dens(i) /= area_skm;
       }


  // |-----------------------------------------------------------------------|
  // | FISHERY AGE COMPS, N, EFFN, SDNR
  // |-----------------------------------------------------------------------|

     for (i=1;i<=nyrs_fish_age;i++) 
       {
         age_comp(i)      = catage(yrs_fish_age(i))/
                            sum(catage(yrs_fish_age(i)))*ageage;

         effn_fish_age(i) = ((1-age_comp(i))*age_comp(i))/
                            norm2(oac_fish(i)-age_comp(i));
                            
         sdnr_fish_age(i) = sdnr(age_comp(i),oac_fish(i),
                            double(nmulti_fish_age(i)));
       }
 

  // |----------------------------------------------------------------------|
  // | COMMERCIAL FISHERY CPUE
  // |----------------------------------------------------------------------|

     for (i=1;i<=nyrs_cpue;i++)
       {
         pred_cpue(i) = q1* sum(elem_prod(natage(yrs_cpue(i)),wt)) / 1000; 
       }


  // |----------------------------------------------------------------------|
  // | IPHC SURVEY CPUE
  // |----------------------------------------------------------------------|  
  
     for (i=1;i<=nyrs_cpue_iphc;i++)
       {
         pred_cpue_iphc(i)  = q2* sum(natage(yrs_cpue_iphc(i))) / 1000;
       }


FUNCTION Population_Summaries
  tot_biomass.initialize();
  tot_N.initialize();
  expl_rate.initialize();
  spawn_biom.initialize();
  recruitment.initialize();


     for (i=styr;i<=endyr;i++)
       {
         tot_biomass(i) = (natage(i)*wt);
         tot_N(i)       = sum(natage(i));

         expl_rate(i)   = pred_catch(i)/tot_biomass(i);
                  
         spawn_biom(i)  = natage(i) * mfexp(-spawn_fract * M) * wt_mature;

         recruitment(i) = natage(i,1);

      }

             
            // spawn_biom_ABC(i) = sum((N_ABC(i) * mfexp(-spawn_fract * M)) *
                                             //wt_mature);
FUNCTION Penalties
  penalties.initialize();

    // |-----------------------------------------|
    // | Year 1 deviations penalty
    // |-----------------------------------------|
       
       penalties(1)  =norm2(init_pop+sig1*sig1/2.)/(2.*square(sig1)) +
                       size_count(init_pop)*log(sig1);



    // |-----------------------------------------|
    // | Recruitment deviations penalty
    // |-----------------------------------------|

       penalties(2)  = norm2(log_rec_dev+sigr*sigr/2.)/(2.*square(sigr)) +
                       size_count(log_rec_dev)*log(sigr);



    // |-----------------------------------------|
    // | Natural mortality penalty
    // |-----------------------------------------|
       int mpen = 0.2;
       if (last_phase())
         {
           mpen = 2;
         }
      
       penalties(3)  = 0.5*log(2*M_PI) + log(mpen) +
                       0.5*(square(log_natural_mortality - log(0.026)))
                       / (2*square(mpen));



    // |-----------------------------------------|
    // | Stabilize F and Fs estimates
    // |-----------------------------------------|
       int fpen = 1;
       if (last_phase())
         {
           fpen = 2;
         }
      
       penalties(4)  = dnorm(log_avg_F,log(0.02),fpen);
       penalties(5)  = dnorm(log_avg_Fs,log(0.01),fpen);
       penalties(6)  = dnorm(log_mean_rec,4,fpen);
       penalties(7)  = dnorm(log_mean_y1,4,fpen);


FUNCTION Catch_Like
  c_catch_like.initialize();
  s_catch_like.initialize();

  // |----------------------------------------------------------------------|
  // | CATCH LIKELIHOODS
  // |----------------------------------------------------------------------|

     c_catch_like  +=  0.5*log(2*M_PI) + log(sigma_catch) +
                       0.5*(norm2(log(obs_catch+oo) -
                       log(pred_catch+oo))) / (2.*square(sigma_catch));


    s_catch_like  +=  0.5*log(2*M_PI)  +log(sigma_sport) +
                      0.5*(norm2(log(obs_sportcatch+oo) -
                      log(pred_sportcatch+oo))) / (2.*square(sigma_sport));



FUNCTION Surv_Like
  surv_like.initialize();
  
  // |----------------------------------------------------------------------|
  // | ROV numbers per square kilometer - log-normal
  // |----------------------------------------------------------------------|

     for (i=1; i<=nyrs_srv; i++)
       {
         l_var(i) = log(1. + (square(obs_srv_se(i))/square(obs_srv_biom(i))));


         surv_like(1) += 0.5*log(2*M_PI) + log(sqrt(l_var(i))) +
                         0.5*(square(log(obs_srv_biom(i)) - log(pred_srv(i))) /
                         (2*(l_var(i))));
       }
  // |----------------------------------------------------------------------|
  // | COMMERCIAL CPUE - normal
  // |----------------------------------------------------------------------|

     for (i=1; i<=nyrs_cpue; i++)
       {
     
         surv_like(2) += 0.5*log(2*M_PI) + log(var_cpue(i))+
                         0.5*(square(obs_cpue(i)-pred_cpue(i)))
                         / (2*var_cpue(i)+oo);
       }


  // |----------------------------------------------------------------------|
  // | IPHC SURVEY CPUE - log-normal
  // |----------------------------------------------------------------------|

     for (i=1; i<=nyrs_cpue_iphc; i++)
       {                         

         l_var_iphc(i) = log(1 + ((var_cpue_iphc(i)+oo)/square(obs_cpue_iphc(i)+oo)));

         surv_like(3)  +=  0.5*log(2*M_PI) + log(sqrt(l_var_iphc(i))+oo) +
                           0.5*(square(log(obs_cpue_iphc(i))-log(pred_cpue_iphc(i))))
                           / (2*l_var_iphc(i)+oo);
       }


FUNCTION Age_Like
  age_like.initialize();

  // |----------------------------------------------------------------------|
  // | MULTINOMIAL LIKELIHOODS FOR AGE COMPOSITIONS
  // |----------------------------------------------------------------------|

     for (i=1; i <= nyrs_fish_age; i++)
       {
         age_like -= nmulti_fish_age(i)*
                     ((oac_fish(i) + oo) * log(age_comp(i) + oo)) ;
       }

     age_like -= offset;



FUNCTION Objective_Function
  Like.initialize();

  // |----------------------------------------------------------------------|
  // | CALL LIKELIHOOD FUNCTIONS
  // |----------------------------------------------------------------------|

     Catch_Like();
          if(DEBUG_FLAG == 1) cout<<"Catch_Like"<<endl;
     Surv_Like();
          if(DEBUG_FLAG == 1) cout<<"Surv_Like"<<endl;
     Age_Like();
          if(DEBUG_FLAG == 1) cout<<"Age_Like"<<endl;
     Penalties();
          if(DEBUG_FLAG == 1) cout<<"Penalties"<<endl;

  // |----------------------------------------------------------------------|
  // | SUM DATA LIKELIHOODS
  // |----------------------------------------------------------------------|

     Like       = c_catch_like;
     Like      += s_catch_like;
     Like      += sum(surv_like);
     Like      += age_like;

     obj_fun   = Like;         // Put here to capture the data likelihood

     obj_fun   += sum(penalties);


     if (last_phase())
       obj_fun += sprpen; 

  // |--------------------------------------|
  // | Just for some graphics
  // |--------------------------------------|

   for (i=1; i <= nyrs_fish_age; i++)
     {
       res_age(i) = oac_fish(i) - age_comp(i);
     }


FUNCTION spr_rates
  SBx=0.;
  Bx.initialize();
  Nspr.initialize();
  N0.initialize();
  sprpen.initialize();

  // |----------------------------------------------------------------------|
  // | SPR RATES - ADDED TO LIKELIHOOD AS 'sprpen'
  // |----------------------------------------------------------------------|
           

  // |-----------------------------------------|
  // | RECRUITMENT (AGE 8)
  // |-----------------------------------------|

     for (i=1;i<=10;i++)
       {
         Nspr(i,1) = 1.;
         N0(i,1) = mean(recruitment(1987,endyr-recage));
       }

  // |-----------------------------------------|
  // | AGES 9 - 75+
  // |-----------------------------------------|

     for (j=2;j<=nages;j++)  // j loop
       {
         Nspr(1,j) = Nspr(1,j-1)*mfexp(-1.*M);
         N0(1,j)   = N0(1,j-1)*mfexp(-1.*M);

       for (i=2; i<=10; i++)
         {
           Nspr(i,j) = Nspr(i,j-1)*mfexp(-1.*(M+mF(i-1)*fish_sel(j-1)));
           N0(i,j)   = N0(i,j-1)*mfexp(-1.*(M+mF(i-1)*fish_sel(j-1)));
         }
       } // close 'j' loop


  // |-----------------------------------------|
  // | PLUS CLASS
  // |-----------------------------------------|

      Nspr(1,nages) += Nspr(1,nages)*mfexp(-1.*M);
      N0(1,nages)   += N0(1,nages)*mfexp(-1.*M);


     for(i=2; i<=10; i++)
      {
         Nspr(i,nages) += Nspr(i,nages)*mfexp(-1.*(M+mF(i-1)*fish_sel(nages)));
         N0(i,nages)   += N0(i,nages)*mfexp(-1.*(M+mF(i-1)*fish_sel(nages)));
       }


  // |----------------------------------------------------------------------|
  // | PROJECTED BIOMASS FOR F LEVELS
  // |----------------------------------------------------------------------|

     SBx(1) +=  Nspr(1)*wt_mature*mfexp(-spawn_fract*M);
     

     for (j=1;j<=nages;j++)
       {
         for(i=2; i<=10; i++)
           {
             SBx(i) += Nspr(i,j)*wt_mature(j)*
                       mfexp(-spawn_fract*(M+mF(i-1)*fish_sel(j)));
           }
       } // close 'j' loop


     for(i=1; i<=10; i++)
       {
         Bx(i) +=  N0(i)*wt;
       }


  // |----------------------------------------------------------------------|
  // | CALCULATE F LEVELS
  // |----------------------------------------------------------------------|
  
     sprpen    = 100.*square(SBx(10)/SBx(1)-0.7);
     sprpen   += 100.*square(SBx(9)/SBx(1)-0.65);
     sprpen   += 100.*square(SBx(8)/SBx(1)-0.6);
     sprpen   += 100.*square(SBx(7)/SBx(1)-0.55);
     sprpen   += 100.*square(SBx(6)/SBx(1)-0.5);
     sprpen   += 100.*square(SBx(5)/SBx(1)-0.45);
     sprpen   += 100.*square(SBx(4)/SBx(1)-0.4);
     sprpen   += 100.*square(SBx(3)/SBx(1)-0.35);
     sprpen   += 100.*square(SBx(2)/SBx(1)-0.3);
  
     B40 = SBx(4)*mean(recruitment(1987,endyr-recage));
     F65 = mF(8);
     F60 = mF(7);

     if(DEBUG_FLAG == 1) cout<<"SPR"<<endl;


FUNCTION Population_Projection

  // |------------------------------------------------------------------------|
  // | ABUNDANCE FOR FIRST PROJECTION YEAR (endyr+1)
  // |------------------------------------------------------------------------|


  // |------------------------------------------------------------------------|
  // | AGE ONE (8)
  // |------------------------------------------------------------------------|

     int k;

     if(mceval_phase()) 
       {
         // random_number_generator r(1000);
         // standard deviation of mean recruitment
         stdev_rec   = sqrt(norm2(value(log_rec_dev(1987,endyr-recage))-
                       mean(value(log_rec_dev(1987,endyr-recage))))/
                       (size_count(value(log_rec_dev))-1));
                          
         k=k+1;
            
         // randomized projected recruitment
         N_ABC(endyr+1,1) = mfexp(value(log(mean(recruitment(1987,endyr-recage)))+
                       stdev_rec*stdev_rec/2+stdev_rec*randn(k)));

         N_OFL(endyr+1,1) = mfexp(value(log(mean(recruitment(1987,endyr-recage)))+
                       stdev_rec*stdev_rec/2+stdev_rec*randn(k)));

       }
          
     else 
       {
         N_ABC(endyr+1,1) = value(mean(recruitment(1987,endyr-recage)));
         N_OFL(endyr+1,1) = value(mean(recruitment(1987,endyr-recage)));
       }
     
  // |------------------------------------------------------------------------|
  // | AGES 9 - 75
  // |------------------------------------------------------------------------|

     for(j=2; j<=nages;j++) 
       {
          N_ABC(endyr+1,j) = natage(endyr,j-1)*S(endyr,j-1);
          N_OFL(endyr+1,j) = natage(endyr,j-1)*S(endyr,j-1);

       }

  // |------------------------------------------------------------------------|
  // | PLUS CLASS
  // |------------------------------------------------------------------------|

     N_ABC(endyr+1,nages) += natage(endyr,nages)*S(endyr,nages);
     N_OFL(endyr+1,nages) += natage(endyr,nages)*S(endyr,nages);

    
  // |------------------------------------------------------------------------|
  // | PROPAGATE FORWARD WITH ABC / OFL MORTALITIES
  // |------------------------------------------------------------------------|


     // |---------------------------------------------------------------------|
     // | TIER 4 NOAA REGULATIONS
     // |---------------------------------------------------------------------|
     // | FABC = F65 - author's recommendation
     // | FOFL = F60 - author's recommendation
     // | You can change these easily above in CALCULATE F LEVELS by changing
     // |   the vector elements of mF
     
        for (j=1;j<=nages;j++)
          {  
            FOFL(j) = fish_sel(j) * F60;
            FABC(j) = fish_sel(j) * F65;

            ZABC(j) = FABC(j) + M;
            ZOFL(j) = FOFL(j) + M;

          }


  // |------------------------------------------------------------------------|
  // | AGES NINE THROUGH +
  // |------------------------------------------------------------------------|

     for (i=endyr+2;i<=endyr+15;i++)
       {
         N_ABC(i,1) = value(mean(recruitment(1987,endyr-recage)));
         N_OFL(i,1) = value(mean(recruitment(1987,endyr-recage)));
         
           for(j=2; j<=nages;j++) 
             {
               N_ABC(i,j) = N_ABC(i-1,j-1)* mfexp(-FABC(j-1)-M);
               N_OFL(i,j) = N_OFL(i-1,j-1)* mfexp(-FOFL(j-1)-M);
             }

         N_ABC(i,nages) += N_ABC(i,nages) * mfexp(-FABC(nages)-M);
         N_OFL(i,nages) += N_OFL(i,nages) * mfexp(-FOFL(nages)-M);
       }


  // |------------------------------------------------------------------------|
  // | CATCH AND BIOMASS
  // |------------------------------------------------------------------------|

     for (i=endyr+1; i<=endyr+15; i++)
       {
         for (j=1; j<=nages; j++)
           {
             catage_ABC(i) = elem_div(elem_prod(elem_prod(N_ABC(i),FABC),
                                                (1.-mfexp(-ZABC))),ZABC);
             catage_OFL(i) = elem_div(elem_prod(elem_prod(N_OFL(i),FOFL),
                                                (1.-mfexp(-ZOFL))),ZOFL);
           
             catch_ABC(i) = catage_ABC(i)*wt;
             catch_OFL(i) = catage_OFL(i)*wt;
  
             tot_biom_ABC(i)   = N_ABC(i)*wt;
             
             spawn_biom_ABC(i) = (N_ABC(i) * mfexp(-spawn_fract * M)) *
                                             wt_mature;
             spawn_biom_OFL(i) = (N_OFL(i) * mfexp(-spawn_fract * M)) *
                                             wt_mature;
           }
       }
       

     OFL=catch_OFL(endyr+1);
     ABC=catch_ABC(endyr+1);

     if(DEBUG_FLAG == 1) cout<<"Projections"<<endl;
          



TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  arrmblsize=39000000;

  

RUNTIME_SECTION
  maximum_function_evaluations 5000 5000  5000  5000 5000;
  convergence_criteria 0.0001;


FINAL_SECTION

   // |----------------------------------------------------------------------|
   // | Print run time statistics to the screen.
   // |----------------------------------------------------------------------|
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
        cout<<""<<endl;
	cout<<"--Objective function value: "<<obj_fun<<endl;
        cout<<""<<endl;
	cout<<"--Maximum gradient component: "<<objective_function_value::gmax<<endl;
        cout<<""<<endl;
	cout<<"*******************************************"<<endl;
	cout<<""<<endl;
 

  // |--------------------------------------------------------------------------|
  // |                                                                          |
                                "PARAMETRIC BOOTSTRAP"                          ;
  // |                                                                          |
  // |--------------------------------------------------------------------------|
  // | This is a good check for parameter distributions and whether they are
  // |  being artificially truncated by the bounds defined in the
  // |  PARAMETER_SECTION. It is also a good method for
  // |  estimating uncertainty of derived quantities. MCMC is also readily
  // |  available for that. If the parameters are not exceptionally constrained,
  // |  the results between the parametric boostrap and the MCMC should be fairly
  // |  similar, and the MCMC is a bit easier to use.

     "  |-----------------------------------------------------------------  |  ";
     "  |                        BOOTSTRAP CONTROLS                         |  ";
     "  |                                                                   |  ";
     "  | 'const int stuff' should be set to the number of entries in the   |  ";
     "  |   .cor file. If it barfs, make it +1.                             |  ";
     "  |                                                                   |  ";
     "  | 'nboot', defining the number of bootstrap draws, is set in the    |  ";
     "  |   .ctl file and referenced in line 70 above.                      |  ";
     "  |                                                                   |  ";
     "  | If you make it through all the output messages and reach the      |  ";
     "  | 'M5 tie-in' message and then it simply sits there, it is because  |  ";
     "  |  one or more of the random draw bounds are halting it. Comment    |  ";
     "  |  out each parameter in the section below 'LIMITS ON RANDOM DRAWS' |  ";
     "  |  and re-run until it suddenly returns 'Got vector x'              |  ";
     "  |  and you'll know which parameter(s) is(are) misbehaving.          |  ";

        const int stuff = 412;

     "  |-----------------------------------------------------------------  |  ";


  // |------------------------------------------------------------------------|
  // | CALL IN/OUT FILE AND NAMESPACE
  // |------------------------------------------------------------------------|

     using namespace std;

     ifstream infile;
     ofstream outfile;


  // |------------------------------------------------------------------------|
  // | STRUCTURES
  // |------------------------------------------------------------------------|
  // |
  // |

     const int maxlen=46;
     const int maxline=5000;
     const int maxpar=initial_params::nvarcalc();
     char header_std[maxlen];
     char header_cor[maxline];
     char pname[maxlen];
     int icount;
     int ii;
     int jj;

     int nwork;
     int ndump;
     int marker;
     int nparameter;
     int index;
     double estimate;
     double stdev;
     double tmp_cor[stuff];
     double tmp_shit[stuff];
     //int count;

     //int good_count;

     cout<<" "<<endl;
     cout<<"What's the story, mother?"<<endl;


  // |------------------------------------------------------------------------|
  // | RANDOM NUMBER GENERATOR
  // |------------------------------------------------------------------------|
  // | - uses current system time as random seed

     random_number_generator rng(time(NULL));


  // |------------------------------------------------------------------------|
  // | CONTAINERS
  // |------------------------------------------------------------------------|

     dvector Boot_index_init(1,stuff);
     dvector Boot_stdev_init(1,stuff);

     dvector Boot_index(1,maxpar);
     dvector Boot_stdev(1,maxpar);

     dmatrix Boot_corr(1,stuff,1,stuff);
     dmatrix Boot_correlation(1,maxpar,1,maxpar);
     dmatrix Boot_covariance(1,maxpar,1,maxpar);
     dmatrix Boot_choleski(1,maxpar,1,maxpar);

     dvector Boot_ran(1,maxpar);
     dvector Bootshit(1,maxpar);

     dvector draw(1,maxpar);
     dvector bootN(1,maxpar);
     dmatrix boot(1,nboot,1,maxpar);

     Boot_index_init.initialize();
     Boot_index.initialize();
     Boot_stdev_init.initialize();
     Boot_stdev.initialize();
     Boot_corr.initialize();
     Boot_correlation.initialize();
     Boot_covariance.initialize();
     Boot_choleski.initialize();
     Boot_ran.initialize();
     draw.initialize();
     bootN.initialize();

     cout<<" "<<endl;
     cout<<"Containers loaded"<<endl;


  // |------------------------------------------------------------------------|
  // | OPEN THE 'MODEL.COR' FILE
  // |------------------------------------------------------------------------|

     infile.open("model.cor");


  // |------------------------------------------------------------------------|
  // | EXISTENCE CHECK
  // |------------------------------------------------------------------------|

     if (!infile.is_open())
       {
         cout<<" "<<endl;
         cout << "I'm sorry Dave... I can't do that..."<<endl;;
         exit(1);
       }


  // |------------------------------------------------------------------------|
  // | REMOVE FIRST TWO LINES OF ADMB CORRELATION FILE
  // |------------------------------------------------------------------------|

     infile.getline(header_cor,maxline,'\n'); //removes first line
     infile.getline(header_cor,maxline,'\n'); //removes second line


  // |------------------------------------------------------------------------|
  // | READ IN PARAMETER VALUE, ST. DEV., AND CORR VALUES
  // |------------------------------------------------------------------------|


     while (!infile.eof())
      {
        infile >> index;
        infile >> pname;
        infile >> estimate;
        infile >> stdev;

        // cout<<"COR file index "<<index<<", parameter name "<<pname<<endl;

        for (icount=1; icount<=index; icount++)
        {
            infile >> tmp_cor[icount];

            ii=index;  //-(Boot_index-1);
            jj=icount; //-(Boot_index-1);

            Boot_stdev_init(ii) = stdev;
            Boot_index_init(ii) = estimate;

            Boot_corr(ii,jj) = tmp_cor[icount];

        }
      }

     infile.close();

     cout<<" "<<endl;
     cout<<"I say take off and nuke it from orbit - only way to be sure"<<endl;


  // |------------------------------------------------------------------------|
  // | TRIM AND MIRROR CORRELATION MATRIX
  // |------------------------------------------------------------------------|

    for(ii=1; ii<=maxpar; ii++)
      {
         Boot_stdev(ii) = Boot_stdev_init(ii);
         Boot_index(ii) = Boot_index_init(ii);

       for(jj=1; jj<=maxpar; jj++)
         {
           Boot_correlation(ii,jj) = Boot_corr(ii,jj);
         }
      }

     for (ii=1; ii<=(maxpar-1); ii++)
       {
        for (jj=(ii+1); jj<=maxpar; jj++)
          {
            Boot_correlation(ii,jj) = Boot_correlation(jj,ii);
          }
       }

     for (ii=1; ii<=maxpar; ii++)
       {
        for (jj=1; jj<=maxpar; jj++)
          {
            Boot_covariance(ii,jj) = Boot_correlation(ii,jj)*Boot_stdev(ii)*Boot_stdev(jj);
          }
       }

     cout<<" "<<endl;
     cout<<"Trim, correlation and covariance implemented."<<endl;


  // |------------------------------------------------------------------------|
  // | CHOLESKI DECOMPOSITION FOR SQUARE ROOT OF CORRELATION
  // |------------------------------------------------------------------------|

     Boot_choleski = choleski_decomp(Boot_covariance);

     cout<<" "<<endl;
     cout<<"Choleski decomp successful; launch sequence initiated"<<endl;


  // |------------------------------------------------------------------------|
  // | REMOVE OLD BOOTSTRAP RUNS
  // |------------------------------------------------------------------------|

     remove("boot.txt");
     remove("SB.txt");
     remove("B.txt");
     remove("R.txt");
     remove("D.txt");


  // |------------------------------------------------------------------------|
  // | MESSAGE OUTPUT
  // |------------------------------------------------------------------------|

     cout<<""<<endl;
     cout<<"Old files terminated"<<endl;


  // |------------------------------------------------------------------------|
  // | MAIN COUNTER
  // |------------------------------------------------------------------------|

      nwork = 0;             // number of reasonable bootstrap vectors
      ndump = 0;             // number of unreasonable bootstrap vectors


  // |------------------------------------------------------------------------|
  // | MESSAGE OUTPUT
  // |------------------------------------------------------------------------|

     cout<<""<<endl;
     cout<<"M-5 tie in... M-5... working... calling Bootstrap"<<endl;


  // |------------------------------------------------------------------------|
  // | BOOTSTRAP LOOP
  // |------------------------------------------------------------------------|
  // |
  // | - Note that the random draws are implemented INSIDE the
  // |   bootstrap loop, and that each call produces a vector
  // |   MAXPAR long which is replicated over each iteration
  // |   of 'nboot'


      while(nwork < nboot)   //Primary loop -
        {                           //Doesn't close
         if(nwork<=nboot)              //until the very end
        {

      draw.initialize();
      draw.fill_randn(rng);         //Random draws ~N(0,1)


  // |------------------------------------------------------------------------|
  // | ADD THE RANDOM STANDARD DEVIATION TO THE MLE PARAMETER
  // |------------------------------------------------------------------------|
  // |
  // | This implements the multivariate normal distribution
  // | by including the correlation between parameters
  // | contained in the choleski matrix. If there were no correlation,
  // | the choleski would contain only zeroes in the off-diagonal,
  // | and condense to a univariate normal
  // |


     bootN.initialize();
     bootN = Boot_index +  Boot_choleski * draw;


  // |------------------------------------------------------------------------|
  // | PARAMETER VALUES
  // |------------------------------------------------------------------------|
  // |
  // | The vector 'bootN' contains all estimated MLE parameters
  // | from the original model, now with random variability added.
  // | This code assigns the elements of the 'bootN' vector their
  // | proper names as per the PARAMETER_SECTION above.
  // |
  // | That way we can simply re-call the primary model functions
  // | already defined above to calculate the new set of bootstrap
  // | derived quantities for each bootstrap iteration

     log_natural_mortality.initialize();
     log_mean_rec.initialize();
     log_mean_y1.initialize();
     sig1.initialize();
     log_avg_F.initialize();
     log_avg_Fs.initialize();
     fish_sel_slope.initialize();
     fish_sel_a50.initialize();
     log_q1.initialize();
     log_q2.initialize();
     log_rec_dev.initialize();
     init_pop.initialize();        
     log_F_devs.initialize();
     log_F_devs_sport.initialize();


     log_natural_mortality      =  bootN(1);      // Natural mortality
     log_mean_rec               =  bootN(2);      // Mean recruitment;
     log_mean_y1                =  bootN(3);      // Mean Year 1 naa;
     sig1                       =  bootN(4);      // Year 1 sigma
     log_avg_F                  =  bootN(5);      // Mean F
     log_avg_Fs                 =  bootN(6);      // Mean F
     fish_sel_slope             =  bootN(7);      // Fishery selectivity slope
     fish_sel_a50               =  bootN(8);      // Fishery selectivity 50%
     log_q1                     =  bootN(9);      // Catchability (all fisheries CPUE)
     log_q2                     =  bootN(10);     // Catchability (IPHC survey CPUE)

        
     for(ii = 11; ii<=41; ii++)
       {
         log_rec_dev(ii+styr-11) = bootN(ii);     // Recruitment deviations
       }


     for(ii = 42; ii<=108; ii++)
       {
         init_pop(ii-41) = bootN(ii);             // Year 1 abundance deviations
       }

      
     for(ii = 109; ii<=139; ii++)
       {
         log_F_devs(ii - 109 + styr) = bootN(ii); // F deviations, comm. fishery
       }


     for(ii = 140; ii<=149; ii++)
       {
         log_F_devs_sport(ii-139) = bootN(ii);    // F deviations, sport fishery
       }



     // |---------------------------------------------------------------------|
     // | LIMITS ON RANDOM DRAWS AS PER ORIGINAL PARAMETER BOUNDS
     // |---------------------------------------------------------------------|
     // |
     // | TEST DISTRIBUTION OF DRAWS AGAINST MLE DISTRIBUTION FROM MLE MEAN
     // | AND ST. DEV. FROM THE .STD FILE - BOUNDS SHOULD SERVE **ONLY** TO
     // | PROVIDE GENERAL LIMITS TO THE ADMB MINIMIZATION FUNCTION, **NOT**
     // | TO COMPENSATE FOR BAD CODE, ERRORS IN STATISTICAL ASSUMPTIONS, OR
     // | UNINFORMATIVE DATA.
     // |
     // | EXPERIMENT WITH PENALIZED DEVIATIONS FROM PRIOR VALUES AND DIFFERENT
     // | ERROR DISTRIBUTIONS IN THE OBJECTIVE FUNCION IF HARD BOUNDS RESULT
     // | IN A TRUNCATED EXPLORATION OF THE PARAMETER SPACE, AND THEN ALTER THE
     // | BOUNDS ACCORDINGLY

    
     if (log_natural_mortality    >   0.00 ||
         log_natural_mortality    < -15.0  ||
         log_mean_rec             >   8.0  ||
         log_mean_rec             <   2.0  ||
         log_mean_y1              >   8.0  ||
         log_mean_y1              <   2.0  ||
         sig1                     >   1.5  ||
         sig1                     <   0.0  ||
         log_avg_F                < -10.0  ||
         log_avg_F                >   0.0  ||
         log_avg_Fs               < -10.0  ||
         log_avg_Fs               >   0.0  ||
         fish_sel_slope           <  -3.0  ||
         fish_sel_slope           >   0.0  ||
         fish_sel_a50             >   3.0  ||
         fish_sel_a50             <   1.0  ||
         log_q1                   <  -7.0  ||
         log_q1                   >   0.0  ||
         log_q2                   <  -7.0  ||
         log_q2                   >   0.0  ||

         min(log_rec_dev)         <  -5.0  ||
         max(log_rec_dev)         >   5.0  ||
         min(init_pop)            <  -5.0  ||
         max(init_pop)            >   5.0  ||
         min(log_F_devs)          < -15.0  ||
         max(log_F_devs)          >  15.0  ||
         min(log_F_devs_sport)    < -15.0  ||
         max(log_F_devs_sport)    >  15.0   )

     {
         //cout<<"Dumping vector"<<endl;
          //go on to the next vector
         ndump++;
         continue;
     }
     else
     {
         nwork++;
         cout<<"Got vector "<<nwork<<endl;
     }


  // |-----------------------------------------------------------|
  // | CALL MODEL FUNCTIONS FROM MLE STRUCTURE ABOVE
  // |-----------------------------------------------------------|
  // |
  // | As we are not estimating anything in the bootstrap, we simply
  // | call the population dynamics functions from the MLE model
  // | but without the objective function. The randomized MLE
  // | parameters are used to calculate the derived quantities
  // | (abundance, spawning biomass, etc.) which are written to
  // | external files for assessment of variance.
  
     Selectivity();	
     Mortality();
     Abundance();
     Catch();
     Predicted();
     Population_Summaries();


  // |-----------------------------------------------------------|
  // | SAID EXTERNAL WRITING
  // |-----------------------------------------------------------|


    ofstream ofs0("boot.txt",ios::app);
    {
        ofs0<<bootN<<endl;
    }

    ofstream ofs1("SB.txt",ios::app);
    {
        ofs1<<spawn_biom<<endl;
    }

    ofstream ofs2("B.txt",ios::app);
    {
        ofs2<<tot_biomass<<endl;
    }

    ofstream ofs3("R.txt",ios::app);
    {
        ofs3<<recruitment<<endl;
    }

    ofstream ofs4("D.txt",ios::app);
    {
        ofs4<<pred_dens<<endl;
    }


   
     //} //if-else loop for 'good_count = 1'
     } // close conditional 'if' g++ loop
     } // close primary g++ loop




        cout<< "O frabjous day!"<<endl;
        cout<< "The sheep frolic!"<<endl;
        cout<<""<<endl;
        cout<<"        ...moo..."<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"             _.%%%%%%%%%%%%%             " <<endl;
        cout<<"            //-_%%%%%%%%%%%%%            " <<endl;
        cout<<"           (_ %\\%%%%%%%%%%%%%%~             "<<endl;
        cout<<"               %%%%%%%%%%%%%%             "<<endl;
        cout<<"                 %%%%%*%%%%              "<<endl;
        cout<<"            ,,,,,,||,,,,||,,,,,          "<<endl;
        cout<<""<<endl;

FUNCTION double sdnr(const dvar_vector& pred,const dvector& obs,double m)
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred);
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;


GLOBALS_SECTION

	#include <admodel.h>
        #include <string.h>
        #include <string>
        #include <time.h>
        #include <sstream>
        #include <fstream>
        adstring model_name;
        adstring data_file;

	#undef REPORT
	#define REPORT(object) report << #object "\n" << setw(8) \
	<< setprecision(4) << setfixed() << object << endl;

	#undef COUT
	#define COUT(object) cout << #object "\n" << setw(6) \
	<< setprecision(3) << setfixed() << object << endl;

        time_t start,finish;
        long hour,minute,second;
        double elapsed_time;

        adstring BaseFileName;
        adstring ReportFileName;
        adstring NewFileName;

        adstring stripExtension(adstring fileName)
	 {
           /*
	     This function strips the file extension
	     from the fileName argument and returns
	     the file name without the extension.
           */

           const int length = fileName.size();
           for (int i=length; i>=0; --i)
             {
               if (fileName(i)=='.')
                 {
                   return fileName(1,i-1);
                 }
             }

           return fileName;
         }
         

 void function_minimizer::mcmc_eval(void)
        {
                // |---------------------------------------------------------------------------|
                // | Added DIC calculation.  Martell, Jan 29, 2013                             |
                // |---------------------------------------------------------------------------|
                // | DIC = pd + dbar
                // | pd  = dbar - dtheta  (Effective number of parameters)
                // | dbar   = expectation of the likelihood function (average f)
                // | dtheta = expectation of the parameter sample (average y)

          gradient_structure::set_NO_DERIVATIVES();
          initial_params::current_phase=initial_params::max_number_phases;
          uistream * pifs_psave = NULL;

        #if defined(USE_LAPLACE)
        #endif

        #if defined(USE_LAPLACE)
            initial_params::set_active_random_effects();
            int nvar1=initial_params::nvarcalc();
        #else
          int nvar1=initial_params::nvarcalc(); // get the number of active parameters
        #endif
          int nvar;

          pifs_psave= new
            uistream((char*)(ad_comm::adprogram_name + adstring(".psv")));
          if (!pifs_psave || !(*pifs_psave))
          {
            cerr << "Error opening file "
                    << (char*)(ad_comm::adprogram_name + adstring(".psv"))
               << endl;
            if (pifs_psave)
            {
              delete pifs_psave;
              pifs_psave=NULL;
              return;
            }
          }
          else
          {
            (*pifs_psave) >> nvar;
            if (nvar!=nvar1)
            {
              cout << "Incorrect value for nvar in file "
                   << "should be " << nvar1 << " but read " << nvar << endl;
              if (pifs_psave)
              {
                delete pifs_psave;
                pifs_psave=NULL;
              }
              return;
            }
          }

          int nsamp = 0;
          double sumll = 0;
          independent_variables y(1,nvar);
          independent_variables sumy(1,nvar);

          do
          {
            if (pifs_psave->eof())
            {
              break;
            }
            else
            {
              (*pifs_psave) >> y;
              sumy = sumy + y;
              if (pifs_psave->eof())
              {
                double dbar = sumll/nsamp;
                int ii=1;
                y = sumy/nsamp;
                initial_params::restore_all_values(y,ii);
                initial_params::xinit(y);
                double dtheta = 2.0 * get_monte_carlo_value(nvar,y);
                double pd     = dbar - dtheta;
                double dic    = pd + dbar;
                double dicValue      = dic;
                double dicNoPar      = pd;

                cout<<"Number of posterior samples    = "<<nsamp    <<endl;
                cout<<"Expectation of log-likelihood  = "<<dbar     <<endl;
                cout<<"Expectation of theta           = "<<dtheta   <<endl;
                cout<<"Number of estimated parameters = "<<nvar1    <<endl;
                    cout<<"Effective number of parameters = "<<dicNoPar <<endl;
                    cout<<"DIC                            = "<<dicValue <<endl;
                                        cout<<"y"<<y<<endl;
                    cout<<"sumy"<<sumy<<endl;
                    cout<<"nsamp"<<nsamp<<endl;
                    cout<<"nvar"<<nvar<<endl;
                    cout<<"gmcv_nvar_y"<<get_monte_carlo_value(nvar,y)<<endl;
                break;
              }
              int ii=1;
              initial_params::restore_all_values(y,ii);
              initial_params::xinit(y);
              double ll = 2.0 * get_monte_carlo_value(nvar,y);
              sumll    += ll;
              nsamp++;
              // cout<<sumy(1,3)/nsamp<<" "<<get_monte_carlo_value(nvar,y)<<endl;
            }
          }
          while(1);
          if (pifs_psave)
          {
            delete pifs_psave;
            pifs_psave=NULL;
          }
          return;
        }
        


REPORT_SECTION
  REPORT(nages);
  REPORT(styr);
  REPORT(endyr);
  REPORT(recage);
  REPORT(p_mature);
  REPORT(wt);
  REPORT(morphology);
  REPORT(obs_catch);
  REPORT(nyrs_sport);
  REPORT(yrs_sport);
  REPORT(obs_sportcatch);
  REPORT(nyrs_srv);
  REPORT(yrs_srv);
  REPORT(obs_srv_biom);
  REPORT(area_skm);
  REPORT(nyrs_fish_age);
  REPORT(yrs_fish_age);
  REPORT(nsamples_fish_age);
  REPORT(oac_fish);
  REPORT(nyrs_cpue);
  REPORT(yrs_cpue);
  REPORT(obs_cpue);
  REPORT(var_cpue);
  REPORT(nyrs_cpue_iphc);
  REPORT(yrs_cpue_iphc);
  REPORT(obs_cpue_iphc);
  REPORT(var_cpue_iphc);
  //end DATA

  //likelihoods
  REPORT(yy);
  REPORT(aa);
  REPORT(c_catch_like);
  REPORT(s_catch_like);
  REPORT(age_like);
  REPORT(surv_like);
  REPORT(penalties);

  //END

  //M
  REPORT(mfexp(log_natural_mortality));
  //END

  //Predicted catch
  REPORT(catage);
  REPORT(pred_catch);
  REPORT(pred_sportcatch);
  //END


  //Catch age residuals
  REPORT(res_age);
  
  //Density
  REPORT(pred_srv);
  //END

  //CPUE
  REPORT(pred_cpue);
  REPORT(pred_cpue_iphc);
  //END

  //Selectivity & F
  REPORT(fish_sel);
  REPORT(Fmort);
  //END

  //Numbers at age
  REPORT(natage);
  //END

  //Spawning biomass
  REPORT(tot_biomass);
  REPORT(spawn_biom);
  //END

  //Recruitment
  REPORT(recruitment);
  //END

  //Forecast
  REPORT(catch_ABC);
  REPORT(catch_OFL);
  
  REPORT(ABC);
  REPORT(OFL);
  //END

  //F_40
  REPORT(mF);
  REPORT(B40);

  //SDNR
  REPORT(effn_fish_age);
  REPORT(sdnr_fish_age);
  REPORT(N_ABC)


