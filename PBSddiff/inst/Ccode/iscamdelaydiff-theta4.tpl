// ------------------------------------------------------------------------- //
//         integrated Statistical Catch Age Model (iSCAM)                    //
//                                                                           //
//                           VERSION 1.1                                     //
//               Tue Jul 19 22:23:58 PDT 2011                                //
//                                                                           //
//                                                                           //
//           Created by Steven Martell on 2010-04-09                         //
//           Copyright (c) 2010. All rights reserved.                        //
//                                                                           //
// AUTHORS: SJDM Steven Martell , RF Robyn Forrest (Delay Difference additions)                                           //
//                                                                           //
// CONVENTIONS: Formatting conventions are based on the The                  //
//               Elements of C++ Style (Misfeldt et al. 2004)                //
//                                                                           //
// NAMING CONVENTIONS:                                                       //
//             Macros       -> UPPERCASE                                     //
//             Constants    -> UpperCamelCase                                //
//             Functions    -> lowerCamelCase                                //
//             Variables    -> lowercase                                     //
//                                                                           //
// CHANGED add option for using empirical weight-at-age data                 //
// TODO:    add gtg options for length based fisheries                      //
// CHANGED add time varying natural mortality rate with splines              //
// TODO:    add cubic spline interpolation for time varying M               //
// CHANGED  Fix the type 6 selectivity implementation. not working.          //
// TODO:  fix cubic spline selectivity for only years when data avail        //
// CHANGED: fixed a bug in the simulation model log_ft_pars goes out         //
//        of bounds.                                                         //
// TODO: write a projection routine and verify equilibrium calcs             //
// TODO: add DIC calculation for MCMC routines (in -mcveal phase)            //
// CHANGED: add SOK fishery a) egg fishing mort 2) bycatch for closed ponds  //
//                                                                           //
//                                                                           //
//                                                                           //
// ------------------------------------------------------------------------- //
//-- CHANGE LOG:                                                           --//
//--  Nov 30, 2010 -modified survey biomass by the fraction of total       --//
//--                mortality that occurred during the time of the         --//
//--                survey. User specifies this fraction (0-1) in the      --//
//--                data file as the last column of the relative           --//
//--                abundance index.                                       --//
//--                                                                       --//
//--  Dec 6, 2010 -modified the code to allow for empiracle weight-        --//
//--               at-age data to be used.                                 --//
//--              -rescaled catch and relative abundance /1000, this       --//
//--               should be done in the data file and not here.           --//
//--                                                                       --//
//--  Dec 20, 2010-added prior to survey q's in control file               --//
//--                                                                       --//
//--  Dec 24, 2010-added random walk for natural mortality.                --//
//--                                                                       --//
//--  Jan 23, 2011-in Penticton Hospital with my mom in ICU, adopting      --//
//--               the naming conventions in The Elements of C++           --//
//--               style to keep my mind busy.                             --//
//--                                                                       --//
//-- May 5, 2011- added logistic selectcitivty as a fucntion of            --//
//--              mean body weight.  3 parameter logistic.                 --//
//--              NOT WORKING YET                                          --//
//--                                                                       --//
//-- May 6, 2011- added pre-processor commands to determin PLATFORM        --//
//--              either "Windows" or "Linux"                              --//
//--                                                                       --//
//-- use -mcmult 1.5 for MCMC with log_m_nodes with SOG herrning           --//
//--                                                                       --//
//--                                                                       --//
//-- Dec 11, 2011- added halibut branch to local git repository aim is to  --//
//--               add gender dimension and stock dimension.               --//
//--               This was created on the "twosex" branch in git merged   --//
//--                                                                       --//
//-- Dec 30, 2011- working on length-based selectivity for halibut.        --//
//--                                                                       --//
//-- Jan 5, 2012 - adding spawn on kelp fishery as catch_type ivector      --//
//--             - modified the following routines:                        --//
//--             - calcFisheryObservations                                 --//
//--             - calcTotalMortality                                      --//
//-- TODO: add catch_type to equilibrium calculations for reference points --//
//--                                                                       --//
//-- August 2015 - Delaydifference version: forced log(R0)=log(Rbar)       --//
//--                                                                       --//
// ------------------------------------------------------------------------- //


DATA_SECTION
	!! cout<<"iSCAM has detected that you are on a "<<PLATFORM<<" box"<<endl;
	init_adstring DataFile;
	init_adstring ControlFile;
	init_adstring ProjectFileControl; 
	
	!! ad_comm::change_datafile_name(ProjectFileControl);
	init_int n_tac;
	init_vector tac(1,n_tac);
	
	
	!! BaseFileName=stripExtension(ControlFile);
	!! cout<<BaseFileName<<endl;
	!! ReportFileName = BaseFileName + adstring(".rep");
	
	!! ad_comm::change_datafile_name(DataFile);
	
	int SimFlag;
	int rseed;
	int retro_yrs;
	int delaydiff;
	LOC_CALCS
		SimFlag=0;
		rseed=999;
		int on,opt;
		//the following line checks for the "-SimFlag" command line option
		//if it exists the if statement retreives the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			SimFlag=1;
			rseed=atoi(ad_comm::argv[on+1]);
			//if(SimFlag)exit(1);
		}
		
		// command line option for retrospective analysis. "-retro retro_yrs"
		retro_yrs=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1)
		{
			retro_yrs=atoi(ad_comm::argv[on+1]);
			cout<<"______________________________________________________\n"<<endl;
			cout<<"    **Implementing Retrospective analysis** "<<endl;
			cout<<"    **Number of retrospective years = "<<retro_yrs<<endl;
			cout<<"______________________________________________________"<<endl;
		}
		
		//RF May 22 2013
		// command line option for implementing delay difference model "-delaydiff"
		delaydiff=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-delaydiff",opt))>-1)
		{
			delaydiff=1;
			cout<<"______________________________________________________\n"<<endl;
			cout<<"    **Implementing Delay Difference model** "<<endl;
			cout<<"______________________________________________________"<<endl;
		}
	END_CALCS
	
	//Read in objects from data file using the init_ prefix
	init_int syr;
	init_int nyr;
	!! cout<<"syr\t"<<syr<<endl;
	!! cout<<"nyr\t"<<nyr<<endl;
	
	init_int sage;
	init_int nage;
	!! cout<<"sage\t"<<sage<<endl;
	!! cout<<"nage\t"<<nage<<endl;
	vector age(sage,nage);
	!! age.fill_seqadd(sage,1);
	
	init_int ngear;				//number of gear types with unique selectivities
	!! cout<<"ngear\t"<<ngear<<endl;
	init_vector allocation(1,ngear);
	init_ivector catch_type(1,ngear);  //RF NOTE: NEW
	ivector fsh_flag(1,ngear);
	LOC_CALCS
		//If allocation >0 then set fish flag =1 else 0
		int k;
		allocation = allocation/sum(allocation);
		for(k=1;k<=ngear;k++)
		{
			if(allocation(k)>0)
				fsh_flag(k)=1;
			else
				fsh_flag(k)=0;
		}
	END_CALCS
	
	//The following code has been deprecated
	//ivector ft_phz(1,ngear);
	//LOC_CALCS
	//	int k;
	//	ft_phz=1;
	//	for(k=1;k<=ngear;k++)
	//	{
	//		if(!fsh_flag(k))
	//			ft_phz(k)=-1;
	//	}
	//END_CALCS
	
	init_number fixed_m;		//FIXME: depricate this from data files
	init_number linf;
	init_number vonbk;
	init_number to;
	init_number a;
	init_number b;
	init_number ah;
	init_number gh;
	
	//DELAY DIFFERENCE MODEL PARAMETERS
	init_int kage;  //age at knife-edge recruitment 
        //DD growth parameters :: RF (02-Apr-2013)
	init_number alpha_g;  //growth alpha (intercept of Ford-Walford plot; derived from wk and wk-1, H&W 1992, p 334)
	init_number rho_g;  //growth rho (slope of Ford-Walford plot; H&W 1992, p 332)
  	init_number wk;
	
	vector la(sage,nage);		//length-at-age
	vector wa(sage,nage);		//weight-at-age
	LOC_CALCS
	  cout<<"linf\t"<<linf<<endl;
	  la=linf*(1.-exp(-vonbk*(age-to)));
	  wa=a*pow(la,b);
	  cout<<setprecision(2);		//2 decimal places for output
	  cout<<"la\n"<<la<<endl;
	  cout<<"wa\n"<<wa<<endl;
	  cout<<setprecision(5);
	END_CALCS
	
	//Time series data
	init_matrix catch_data(syr,nyr,1,ngear+1);
	matrix obs_ct(1,ngear,syr,nyr);
	int ft_count;
	int i;
	LOC_CALCS
		ft_count=0;
		for(k=1;k<=ngear;k++)
			obs_ct(k)=column(catch_data,k+1);
		
		for(k=1;k<=ngear;k++)
		{	
			for(i=syr;i<=nyr;i++)
				if( obs_ct(k,i)>0 ) ft_count++;
		}
		cout<<"ft_count\n"<<ft_count<<endl;
		cout<<"last row of catch \n"<<catch_data(nyr)<<endl;
		cout<<"Ok after catch extraction"<<endl;
	END_CALCS
	
	init_int nit;
	!! cout<<"Number of surveys "<<nit<<endl;
	init_ivector nit_nobs(1,nit);
	//#survey type 
	//## 1 = survey is proportional to vulnerable numbers
	//## 2 = survey is proportional to vulnerable biomass
	//## 3 = survey is proportional to spawning biomass (e.g., herring spawn survey)
	init_ivector survey_type(1,nit);
	init_3darray survey_data(1,nit,1,nit_nobs,1,5);
	//init_matrix survey_data(1,nit,1,4);
	imatrix iyr(1,nit,1,nit_nobs);
	imatrix igr(1,nit,1,nit_nobs);
	matrix it(1,nit,1,nit_nobs);
	matrix it_wt(1,nit,1,nit_nobs);		//relative weight
	matrix it_timing(1,nit,1,nit_nobs);	//timing of the survey (0-1)
	LOC_CALCS
		for(i=1;i<=nit;i++)
		{
			iyr(i)=ivector(column(survey_data(i),1));
			igr(i)=ivector(column(survey_data(i),3));
			it(i)=column(survey_data(i),2);
			it_wt(i)=column(survey_data(i),4)+1.e-10;//add a small constant to allow 0 weight
			it_timing(i)=column(survey_data(i),5);
		}
		cout<<"Last row of the relative abundance data\n"<<survey_data(nit)(nit_nobs(nit))<<endl;
		cout<<"OK after relative abundance index"<<endl;
		/*Normalize survey weights so estimated rho parameter is consistent with sd(epsilon)/sd(rec_devs)*/
		double mean_it_wt;
		mean_it_wt = mean(it_wt);
		for(i=1;i<=nit;i++)
		{
			it_wt(i) = it_wt(i)/mean_it_wt;
		}
	END_CALCS
	
	//Age comps for all gears.
	init_int na_gears;	//total number of aging observations
	init_ivector na_nobs(1,na_gears);
	init_ivector a_sage(1,na_gears);	//youngest age in the ageing matrix
	init_ivector a_nage(1,na_gears);	//oldest age in the ageing matrix

	init_3darray A(1,na_gears,1,na_nobs,a_sage-2,a_nage);
	!! cout<<"Last row of the age comp data\n"<<A<<endl;
	!! cout<<"OK after age data"<<endl;
	
	//Mean weight-at-age data (units are kg) (if exists)
	init_int n_wt_nobs;
	init_matrix tmp_wt_obs(1,n_wt_nobs,sage-1,nage);
	
	//@@@~ New data added by RF July 2012 for comparison with 2005 P cod assessment ~@@@
	init_number weight_sig;
	!! cout<<"weight_sig = "<<weight_sig<<endl;
	//PJS code to read non-contiguous mean weight array
	init_int n_weight;
	init_matrix obs_annual_mean_wt(1,n_weight,1,2);
	//RH code to create wyrs for report section
	ivector wyrs(1,n_weight);
	LOC_CALCS
		for(i=1;i<=n_weight;i++)
		{
			wyrs(i) = obs_annual_mean_wt(i,1);
		}
	END_CALCS
	!!cout<<"OK after non-contiguous mean weight array"<<endl;

	//Extra parameters for comparing to 2005 assessment - delete for 2013	
	init_vector anom_obs(syr,nyr);//observed sea level anomalies
	init_int eff_nyr; //final year of effort data
	init_vector effort(syr,eff_nyr); //observed effort data from 2005 assessment (derived from C/CPUE)
	!!cout<<"OK after extra P. cod input data"<<endl;
	
	matrix wt_obs(syr,nyr+1,sage,nage);		//weight-at-age
	matrix wt_dev(syr,nyr+1,sage,nage);		//standardized deviations in weight-at-age
	matrix fec(syr,nyr+1,sage,nage);		//fecundity-at-age
	vector avg_fec(sage,nage);				//average fecundity-at-age
	vector avg_wt(sage,nage);				//average weight-at-age
		
	LOC_CALCS
		int iyr;
		fec.initialize();
		avg_fec.initialize();
		cout<<endl<<"Maturity"<<endl;
		cout<<plogis(age,ah,gh)<<endl<<endl; //test
		
		for(i=syr;i<=nyr+1;i++)
		{
			wt_obs(i)=wa;			
			fec(i)=elem_prod(plogis(age,ah,gh),wt_obs(i));
		}
		//if empirical weight-at-age data exist, the overwrite wt_obs & fec - only if NOT delay difference model
		for(i=1;i<=n_wt_nobs;i++)
		{
			iyr=tmp_wt_obs(i,sage-1);  //index for year
			wt_obs(iyr)=tmp_wt_obs(i)(sage,nage);
			fec(iyr)=elem_prod(plogis(age,ah,gh),wt_obs(iyr));
		}
		
		//CHANGED average fecundity
		//int nfec = fec.rowmax()-fec.rowmin()+1;
		//avg_fec=colsum(fec)/nfec;
		
		//from Jake Schweigert: use mean-weight-at-age data
		//from the last 5 years for the projected mean wt.
		//dvector tmp=colsum(wt_obs.sub(nyr-5,nyr))/6.;
		//wt_obs(nyr+1) = tmp;
		
		/*June 8, 2012, average wt at age for all years*/
		avg_wt = colsum(wt_obs.sub(syr,nyr))/(nyr-syr+1.);
		wt_obs(nyr+1) = avg_wt;
		
		avg_fec   = colsum(fec.sub(syr,nyr))/(nyr-syr+1);
		fec(nyr+1) = avg_fec; 
		
		cout<<"n_wt_nobs\t"<<n_wt_nobs<<endl;
		cout<<"Ok after empirical weight-at-age data"<<endl;
		
		//July 14, 2011 Error handler to ensure non-zero wt_obs
		//in the data. Otherwise causes and error with weight-based
		//selectivities.
		if(min(wt_obs)==0)
		{
			cout<<"Cannot have a observed 0 mean weight at age\n";
			cout<<"in the data file.  Please fix.\n Aborting program!"<<endl;
			exit(2);
		}
		
		//May 5, 2011 SJDM: Calculating standardized deviates
		//in mean weights-at-age from empiracle data to use
		//with selectivity as a function of mean weight at age.
		//Idea borrowed from Vivian Haist in HCAM model.
		/*
			CHANGED: selectivity function based on average weight.
			Based on the logistic function;
			1) calculate a matrix of standardized deviates
			wt_dev = (wt_obs-mean(wt_obs))/sd(wt_obs);
			
			Not implemented, using the version provided by Vivian.
			
			Aug 5, 2011, noticed that Vivian's method was implemented
			in an awkward way where GN selectivities were a random walk
			sel(t+1)=sel(t) + tmp(2);
			
			Trying a compromize where estimating an additional parameter
			that attemps to explain residual variation in age-comps
			via changes in standardized mean weights at age. This is 
			implemented as:
			
			log_sel = log(plogis(age,s1,s2)) + s3*delta
			where delta is a matrix of standardized deviations
			(mu=0, std=1) of weights at age.  delta is calculated
			based on the wt_dev matrix above.
		*/
		wt_dev.initialize();
		dmatrix mtmp = trans(wt_obs);
		for(i=sage;i<=nage;i++)
		{
			dvector wa_dev = (mtmp(i)-mean(mtmp(i)))/sqrt(var(mtmp(i)));
			mtmp(i) = wa_dev;
		}
		wt_dev = trans(mtmp);	//each column has mean=0 sd=1
		
	END_CALCS
	
	//End of data file
	init_int eof;	
	LOC_CALCS
	  cout<<"eof = "<<eof<<endl;
	  if(eof==999){
		cout<<"\n -- END OF DATA SECTION -- \n"<<endl;
	  }else{
		cout<<"\n *** ERROR READING DATA *** \n"<<endl; exit(1);
	  }
	END_CALCS

	number fmsy;					//Fishing mortality rate at Fmsy
	number msy;						//Maximum sustainable yield
	number bmsy;					//Spawning biomass at MSY
	number Umsy;					//Exploitation rate at MSY
	number Vmsy;					//Vulnerable biomass at MSY
	vector age_tau2(1,na_gears);	//MLE estimate of the variance for the age comps
	//catch-age for simulation model (could be declared locally 3d_array)
	3darray d3C(1,ngear,syr,nyr,sage,nage);		
	
	
	// ***************************************************
	// ** Read parameter controls from ControlFile
	// ***************************************************
	!! ad_comm::change_datafile_name(ControlFile);
	
	init_int npar;
	init_matrix theta_control(1,npar,1,7);
	!! cout<<theta_control<<endl;
	
	vector theta_ival(1,npar);
	vector theta_lb(1,npar);
	vector theta_ub(1,npar);
	ivector theta_phz(1,npar);
	ivector theta_prior(1,npar);
	LOC_CALCS
		theta_ival = column(theta_control,1);
		theta_lb = column(theta_control,2);
		theta_ub = column(theta_control,3);
		theta_phz = ivector(column(theta_control,4));
		theta_prior=ivector(column(theta_control,5));
	END_CALCS
	
	
	// ***************************************************
	// ** Read selectivity parameter options
	// ***************************************************
	// type 1 = logistic (2pars)
	// type 2 = selcoffs (A-1 pars)
	// type 3 = cubic spline (age_nodes)
	// type 4 = time varying cubic spline (age_nodes for each year)
	// type 5 = bicubic spline with age_nodes adn yr_nodes
	// type 6 = fixed logistic by turning sel_phz to (-ve)
	// type 7 = logistic (3pars) as a function of body weight.
	init_ivector isel_type(1,ngear);  	//Switch for selectivity
	ivector isel_npar(1,ngear);			//ivector for # of parameters for each gear.
	ivector jsel_npar(1,ngear);			//ivector for the number of rows for time-varying selectivity.
	init_vector ahat(1,ngear);			//age-at-50% vulnerbality for logistic function
	init_vector ghat(1,ngear);			//std at 50% age of vulnerability for logistic funciton
	init_vector age_nodes(1,ngear);		//No. of age-nodes for bicubic spline.
	init_vector yr_nodes(1,ngear);		//No. of year-nodes for bicubic spline.
	init_ivector sel_phz(1,ngear);		//Phase for estimating selectivity parameters.
	init_vector sel_2nd_diff_wt(1,ngear);	//Penalty weight for 2nd difference in selectivity.
	init_vector sel_dome_wt(1,ngear);		//Penalty weight for dome-shaped selectivity.
	!! cout<<isel_type<<endl;
	
	LOC_CALCS
		//cout up the number of selectivity parameters
		//depending on the value of isel_type
		isel_npar.initialize();
		for(i=1;i<=ngear;i++)
		{	
			jsel_npar(i)=1;
			switch(isel_type(i))
			{
				case 1:
					// logistic selectivity
					isel_npar(i) = 2; 
					break;
					
				case 2:
					// age-specific coefficients
					isel_npar(i) = (nage-sage); 
					break;
					
				case 3:
				 	// cubic spline 
					isel_npar(i) = age_nodes(i);
					break;
					
				case 4:	 
					// annual cubic splines
					isel_npar(i) = age_nodes(i);
					jsel_npar(i) = (nyr-syr-retro_yrs)+1;
					break;
					
				case 5:  
					// bicubic spline
					jsel_npar(i) = age_nodes(i);
					isel_npar(i) = yr_nodes(i);
					break;
				
				case 6:
					// fixed logistic (no parameters estimated)
					// ensure sel_phz is set to negative value.
					isel_npar(i) = 2;
					if(sel_phz(i)>0) sel_phz(i) = -1;
					break;
					
				case 7:
					// CHANGED: Now working, Vivian Haist fixed it.
					// logistic (3 parameters) with mean body 
					// weight deviations. 
					isel_npar(i) = 2;
					break;
					
				case 8:
					// Alternative logistic selectivity with wt_dev coefficients.
					isel_npar(i) = 3;
					break;
					
				case 11:
					// Logistic length-based selectivity.
					isel_npar(i) = 2;
					break;
					
				case 12:
					// Length-based selectivity coeffs with cubic spline interpolation
					isel_npar(i) = (nage-sage);
					break;
					
				default: break;
			}
		}
		//cout<<"Number of estimated selectivity parameters\n"<<isel_npar<<endl;
	END_CALCS
	
	//Controls for prior on survey q.
	init_int nits;					//FIXME (redundant with nit, could be deprecated)
	init_ivector q_prior(1,nits);
	init_vector q_mu(1,nits);
	init_vector q_sd(1,nits);
	!! cout<<"nits\n"<<nits<<endl;
	!! cout<<"q Prior\n"<<q_mu<<endl<<q_sd<<endl;
	

	//Miscellaneous controls
	// 1 -> verbose
	// 2 -> recruitment model (1=beverton-holt, 2=rickers)
	// 3 -> std in catch first phase
	// 4 -> std in catch in last phase
	// 5 -> assumed unfished in first year (0=FALSE, 1=TRUE, 2 = AT EQUILIBRIUM WITH FISHING MORTALITY IN SYR - IMPLEMENTED ONLY IN DELAY DIFF MODEL))
	// 6 -> minimum proportion at age to consider in the dmvlogistic likelihood
	// 7 -> mean fishing mortality rate to regularize the solution
	// 8 -> standard deviation of mean F penalty in first phases
	// 9 -> standard deviation of mean F penalty in last phase.
	// 10-> phase for estimating deviations in natural mortality.
	// 11-> std in natural mortality deviations.
	// 12-> number of estimated nodes for deviations in natural mortality
	// 13-> fraction of total mortality that takes place prior to spawning - NOT IMPLEMENTED IN DELAY DIFFERENCE MODEL
	// 14-> switch for age-composition likelihood (1=dmvlogistic,2=dmultinom, 3 = Ageless (RF))
	// 15-> 1= fit to annual mean weights for commercial catch; 0=do not do this
	//extra parameters for comparing with 2005 assessment
	// 16 1= use environmental driver for recruitment, 0 = do not do this (default)
	// 17 1= fix observed effort for GEAR 1 commercial catch 0 = do not do this RF TEMP FOR DD MODEL TESTING
		
	init_vector cntrl(1,17);
	
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	//RF temporary code for comparing to 2005 assessment
	!!if(cntrl(16)==0) theta_phz(9) = -1; //turn off estimation of dt parameter
	!!cout<<"Phase for estimation of dt = "<<theta_phz(9)<<endl;
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	int verbose;
	
	init_int eofc;
	LOC_CALCS
		verbose = cntrl(1);
		cout<<"cntrl\n"<<cntrl<<endl;
		cout<<"eofc\t"<<eofc<<endl;
		if(eofc==999){
			cout<<"\n -- END OF CONTROL FILE -- \n"<<endl;
		}else{
			cout<<"\n ***** ERROR CONTROL FILE ***** \n"<<endl; exit(1);
		}
	END_CALCS
	
	int nf;
	ivector ilvec(1,7); //RF added 7th element (for mean weight likelihood)  
	!!ilvec=ngear;			//number of fisheries
	!!ilvec(2)=nit;			//number of surveys
	!!ilvec(3)=na_gears;	        //number of age-comps
	!!ilvec(4)=1;
	!!ilvec(7)=1;
	
	// SM Oct 31, 2010.  Implementing retrospective analysis.
	//Dont read in any more data below the retrospective reset of nyr
	!! nyr = nyr - retro_yrs;
	
INITIALIZATION_SECTION
 theta theta_ival;
	
PARAMETER_SECTION
	//Leading parameters
	//theta[1]		log_ro, or log_msy
	//theta[2]		steepness(h), or log_fmsy
	//theta[3]		log_m
	//theta[4]		log_avgrec
	//theta[5]		log_recinit
	//theta[6]		rho
	//theta[7]		vartheta
	//theta[8]		qc (RF TEMPORARY PARAMETER scalar scaling fishing mortality to effort)
	//theta[9]		dt (RF TEMPORARY PARAMETER scalar for observed recruitment anomalies - added by RF July 2012 for P. cod)
	

	init_bounded_number_vector theta(1,npar,theta_lb,theta_ub,theta_phz);
	!! for(int i=1;i<=npar;i++) theta(i)=theta_ival(i);

	//Selectivity parameters (A very complicated ragged array)
	init_bounded_matrix_vector sel_par(1,ngear,1,jsel_npar,1,isel_npar,-15.,15.,sel_phz);
			
	LOC_CALCS
			//initial values for logistic selectivity parameters
			//set phase to -1 for fixed selectivity.
			for(int k=1;k<=ngear;k++)
			{
				if( isel_type(k)==1 || isel_type(k)==6 || isel_type(k)>=7 )
				{
					sel_par(k,1,1) = log(ahat(k));
					sel_par(k,1,2) = log(ghat(k));
				}
			}
	END_CALCS
	
	!! int log_ft_phz = 1;
	//Fishing mortality rate parameters
	init_bounded_vector log_ft_pars(1,ft_count,-30.,3.0,log_ft_phz);
	
	LOC_CALCS
		if(!SimFlag) log_ft_pars = log(0.1);
	END_CALCS
	
	
	//CHANGED Vivian Haist: possible bias in Bo due to dev_vector for log_rec_devs
	//	-Try and compare estimates using init_bounded_vector version.
	//	-Now estimating Rinit, Rbar and Ro.
	
	// May 16, Just looked at Ianelli's code, he used bounded_vector, and 
	// separates the problem into init_log_rec_devs and log_rec_devs
	
	//Annual recruitment deviations
	//!! int ii;
	//!! ii=syr-nage+sage;
	//!! if(cntrl(5)) ii = syr;  //if initializing with ro
	//init_bounded_vector log_rec_devs(ii,nyr,-15.,15.,2);
	
	!! int init_dev_phz = 2;
	!! if(cntrl(5)==1) init_dev_phz = -1; //initialised at eqm no fishing mortality
	!! if(cntrl(5)==2) init_dev_phz = -1; //initialised at eqm with fishing mortality (delay diff only)
	init_bounded_vector init_log_rec_devs(sage+1,nage,-15.,15.,init_dev_phz);
	init_bounded_vector log_rec_devs(syr,nyr,-5.,5.,2);
	
	//Deviations for natural mortality
	!! int m_dev_phz = -1;
	!! m_dev_phz = cntrl(10);
	!! int n_m_devs = cntrl(12);
	init_bounded_vector log_m_nodes(1,n_m_devs,-5.0,5.0,m_dev_phz);
	//init_bounded_vector log_m_devs(syr+1,nyr,-5.0,5.0,m_dev_phz);
	
	objective_function_value f;
    
	number ro;					//unfished age-1 recruits
	number bo;					//unfished spawning stock biomass
	number kappa;				//Goodyear compensation ratio
	number h;					//steepness
	number m;					//initial natural mortality rate
	number m_bar;				//average natural mortality rate
	number log_avgrec;			//log of average recruitment.
	number log_recinit;			//log of initial recruitment in syr.
	//number log_avg_f;			//log of average fishing mortality DEPRICATED
	number rho;					//proportion of the observation error
	number varphi				//total precision in the CPUE & Rec anomalies.
	number sig;					//std of the observation errors in CPUE.
	number tau; 				//std of the process errors.
	number no; //unfished numbers
	number wbar; //mean weight delay diff model
	//@@@@@@@@@@@@@@ RF Parameter for comparing to 2005 assessment
	number qc;
	number dt;				//scalar for observed recruitment anomalies - added by RF July 2012 for P. cod
	//@@@@@@@@@@@@@@
	
	//recruitment parameters
	number so;
	number beta;
	
	vector log_rt(syr-nage+sage,nyr);
	vector vax(sage,nage);		//survey selectivity coefficients
	vector q(1,nit);			//survey catchability coefficients
	
	vector sbt(syr,nyr+1);		//spawning stock biomass
	vector tbt(syr,nyr+1); 			//total biomass
	vector rt(syr+sage,nyr); 	//predicted sage recruits from S-R curve	
	vector delta(syr+sage,nyr);	//residuals for stock recruitment
	vector avg_log_sel(1,ngear);//conditional penalty for objective function
	vector log_m_devs(syr+1,nyr);// log deviations in natural mortality
		
	matrix nlvec(1,7,1,ilvec);	//matrix for negative loglikelihoods RF added one component to the likelihood (mean weight)
	//!!cout<<"nlvec = "<<nlvec<<endl; exit(1);
	//matrix jlog_sel(1,ngear,sage,nage);		//selectivity coefficients for each gear type.
	//matrix log_sur_sel(syr,nyr,sage,nage);	//selectivity coefficients for survey.
	 
	matrix N(syr,nyr+1,sage,nage);			//Numbers at age
	matrix F(syr,nyr,sage,nage);			//Age-specific fishing mortality
	matrix M_tot(syr,nyr,sage,nage);		//Age-specific natural mortality
	matrix ft(1,ngear,syr,nyr);				//Gear specific fishing mortality rates
	matrix log_ft(1,ngear,syr,nyr);			//Gear specific log fishing mortlity rates
	matrix Z(syr,nyr,sage,nage);
	matrix S(syr,nyr,sage,nage);
	matrix ct(1,ngear,syr,nyr);			//predicted catch biomass
	matrix eta(1,ngear,syr,nyr);			//residuals for catch
	matrix epsilon(1,nit,1,nit_nobs);		//residuals for survey abundance index
	matrix pit(1,nit,1,nit_nobs);			//predicted relative abundance index
	matrix qt(1,nit,1,nit_nobs);			//catchability coefficients (time-varying)
	
	3darray Ahat(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//predicted age proportions by gear & year
	3darray A_nu(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//residuals for age proportions by gear & year
	3darray log_sel(1,ngear,syr,nyr,sage,nage);		//selectivity coefficients for each gear type.
	3darray Chat(1,ngear,syr,nyr,sage,nage);		//predicted catch-at-age
	
	//Additional objects for delay difference model
	vector numbers(syr,nyr+1); //RF added numbers in the delay difference model - in the ASM this is set to vulnerable numbers in com fishery
	vector biomass(syr,nyr+1); //RF added biomass in the delay difference model - in the ASM this is set to spawning biomass
	vector surv(syr,nyr);	     //delaydiff only
	vector vbcom(syr,nyr); //RF added vulnerable biomass in gear 1 - commercial fishery
	vector vncom(syr,nyr); //RF added vulnerable numbers in gear 1 - commercial fishery
	vector annual_mean_wt(syr,nyr);  //RF addition for P cod 
	
	sdreport_number sd_depletion;
	
PRELIMINARY_CALCS_SECTION
  //Run the model with input parameters to simulate real data.
  nf=0;
  if(SimFlag) 
  {
    initParameters();
    calcSelectivities();
    calcTotalMortality();
    simulation_model(rseed);
  }

RUNTIME_SECTION
    maximum_function_evaluations 100,200,500,25000,25000
    convergence_criteria 0.01,0.01,1.e-5,1.e-5

PROCEDURE_SECTION
	
	if(!delaydiff){	
		if(cntrl(5)==2) cntrl(5)=0; //This control determines whether population is unfished in syr (0=false). The delay diff model also has option 2 where the population is at equilibrium with fishing mortality - not implemented in ASM.
	
		initParameters();
		calcSelectivities();
		calcTotalMortality();
		calcNumbersAtAge();
		calcFisheryObservations();
		calcAgeProportions();
		calcSurveyObservations();
		calc_stock_recruitment();
		calc_annual_mean_weight(); //RF added this for P cod - only gets added to objective function if cntrl(15)==1
	}	

	if(delaydiff){
		cntrl(14)=3; //If using the delay difference model, switch age comp likelihood type to 3 - i.e. set the likelihood of the age comps to zero in the objective function
		initParameters_deldiff();
		calcTotalMortality();
		calcNumbersBiomass_deldiff();
		calcFisheryObservations_deldiff();
		calcSurveyObservations_deldiff();
		calc_stock_recruitment_deldiff();
		calc_annual_mean_weight(); //RF added this for P cod - only gets added to objective function if cntrl(15)==1
	}
	
	calc_objective_function();

	sd_depletion=sbt(nyr)/bo;
	
	if(mc_phase())
	{
		mcmcPhase=1;
	}
	
	if(mceval_phase())
	{
		mcmcEvalPhase=1;
		mcmc_output();
	}
	
	//The following causes a linker error
	//duplicate symbol in libdf1b2o.a
	//dvariable a=3.0;
	//cout<<"testing gammln(dvariable)"<<gammln(a)<<endl;

FUNCTION initParameters
  {
	/*
	This function is used to extract the specific parameter values
	from the init_bounded_number_vector to the specific variables
	used in the code.
	
	Note that you must call this routine before runnning the 
	simulation model to generate fake data.
	*/
	ro          = mfexp(theta(1));
	h           = theta(2);
	m           = mfexp(theta(3));
	log_avgrec  = theta(4);
	log_recinit = theta(5);
	  
	/*
	Variance partitioning:		SJDM new code - RF discovered this Sept 13
	-Estimating total variance as = 1/precision
	 and partition variance by rho = sig^2/(sig^2+tau^2).

	 E.g. if sig = 0.2 and tau =1.12 then
	 rho = 0.2^2/(0.2^2+1.12^2) = 0.03090235
        the total variance is kappa^2 = sig^2 + tau^2 = 1.2944
	*/
	rho         = theta(6);
	varphi      = sqrt(1.0/theta(7));
	sig         = sqrt(rho) * varphi;       
	tau         = sqrt(1-rho) * varphi;
	   
	//@@@@@@@@@@@@@RF Temporary parameter for comparing to 2005 assessment
	qc = theta(8);
	dt = theta(9); //scalar for observed recruitment anomalies - added by RF May 2013 for P. cod
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	switch(int(cntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			kappa = (4.*h/(1.-h));
			break;
		case 2:
			//Ricker model
			kappa = pow((5.*h),1.25);
		break;
	}
	
	
	//TODO Alternative parameterization using MSY and FMSY as leading parameters
	
	

	
	if(verbose)cout<<"**** Ok after initParameters ****"<<endl;
	
  }
	
FUNCTION initParameters_deldiff
  {
	/*
	This function is used to extract the specific parameter values
	from the init_bounded_number_vector to the specific variables
	used in the code.
	
	Note that you must call this routine before runnning the 
	simulation model to generate fake data.
	*/
	ro          = mfexp(theta(4));
	h           = theta(2);
	m           = mfexp(theta(3));
	log_avgrec  = theta(4);            // PJS changed this to set ln(R0)=ln(Rbar)
	log_recinit = theta(5);
	/*
	Variance partitioning:		SJDM new code - RF discovered this Sept 13
	-Estimating total variance as = 1/precision
	 and partition variance by rho = sig^2/(sig^2+tau^2).

	 E.g. if sig = 0.2 and tau =1.12 then
	 rho = 0.2^2/(0.2^2+1.12^2) = 0.03090235
	the total variance is kappa^2 = sig^2 + tau^2 = 1.2944
	*/
	rho         = theta(6);
	varphi      = sqrt(1.0/theta(7));
	sig         = sqrt(rho) * varphi;
	tau         = sqrt(1-rho) * varphi;
	
	
	//@@@@@@@@@@@@@RF Temporary parameter for comparing to 2005 assessment
	qc = theta(8);
	dt = theta(9); //scalar for observed recruitment anomalies - added by RF May 2013 for P. cod
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	switch(int(cntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			kappa = (4.*h/(1.-h)); //compensation ratio CR
			break;
		case 2:
			//Ricker model
			kappa = pow((5.*h),1.25);
		break;
	}
	
	if(verbose)cout<<"**** Ok after initParameters ****"<<endl;
	
  }
		
	
FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
  {
	RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	dvector fa(sage,nage);
	ia.fill_seqadd(0,1./(nodes-1));
	fa.fill_seqadd(0,1./(nage-sage));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	
	//some testing here
	/*dvar_vector spline_nodes(1,nodes);
		spline_nodes.fill_seqadd(-0.5,1./(nodes-1));
		cout<<spline_nodes<<endl;
		vcubic_spline_function test_ffa(ia,spline_nodes);
		cout<<test_ffa(fa)<<endl;
		exit(1);*/
	return(ffa(fa));
  }

FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs, const dvector& la)
  {
	/*interplolation for length-based selectivity coefficeients*/
	RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	return(ffa(fa));
  }

FUNCTION dvar_matrix cubic_spline_matrix(const dvar_matrix& spline_coffs)
  {
	RETURN_ARRAYS_INCREMENT();
	int nodes= spline_coffs.colmax()-spline_coffs.colmin()+1;
	int rmin = spline_coffs.rowmin();
	int rmax = spline_coffs.rowmax();
	
	dvector ia(1,nodes);
	dvector fa(sage,nage);
	ia.fill_seqadd(0,1./(nodes-1));
	//fa.fill_seqadd(sage,1);
	fa.fill_seqadd(0,1./(nage-sage));
	vcubic_spline_function_array fna(rmin,rmax,ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	return(fna(fa));
	
  }

FUNCTION calcSelectivities
  {
	/*
		This function loops over each ngear and calculates the corresponding
		selectivity for that gear type. It first uses a switch statement 
		to calculate selectivity curves based on isel_type where:
		1) logistic selectivity with 2 parameters
		2) age-specific selectivity coefficients with (nage-sage) parameters
		   and the last two age-classes are assumed to have the same selectivity.
		3) a reduced age-specific parameter set based on a bicubic spline.
		4) Time varying cubic spline.
		5) Time varying bicubic spline (2d version)
		6) Fixed logistic
		7) logistic selectivity based on relative changes in mean weight at age
		8) Time varying selectivity based on logistic with deviations in 
		   weights at age (3 estimated parameters).
		11) logistic selectivity with 2 parameters based on mean length
		12) length-based selectivity using cubic spline interpolation
		
		Following the initialization of the selectivity curves, time-varying 
		considerations are implemented.
		
		CHANGED: Add penality (10.*square(avg_log_sel)) to objective function 
		in cases where estimating sel_coffs to mimic an init_bounded_dev vector.
		
		CHANGED: Problem with case 7: turns out to be a random walk, so there
		is changes in selectivity even with constant growth.  Remove the
		random walk outside the switch statement.
		
		TODO: add an option for length-based selectivity.  Use inverse of
		allometric relationship w_a = a*l_a^b; to get mean length-at-age from
		empirical weight-at-age data, then calculate selectivity based on 
		mean length. 
	
	*/
	int i,j,k;
	double tiny=1.e-10;
	dvariable p1,p2,p3;
	dvar_vector age_dev=age;
	dvar_matrix t1;
	dvar_matrix tmp(syr,nyr-1,sage,nage);
	dvar_matrix tmp2(syr,nyr,sage,nage);
	dvar_matrix ttmp2(sage,nage,syr,nyr);
	//jlog_sel.initialize();
	log_sel.initialize();
	avg_log_sel.initialize();
	
	for(j=1;j<=ngear;j++)
	{
		tmp.initialize(); tmp2.initialize();
		dvector iy(1,yr_nodes(j));
		dvector ia(1,age_nodes(j));
		
		switch(isel_type(j))
		{
			case 1:
				// logistic selectivity for case 1 or 6
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				for(i=syr; i<=nyr; i++)
				{
					log_sel(j)(i) = log( plogis<dvar_vector>(age,p1,p2)+tiny );    // <dvar_vector> New in 2014 to allow compile with new statslib.h
				}
				break;
			
			case 6:
				// logistic selectivity for case 1 or 6
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				for(i=syr; i<=nyr; i++)
				{
					log_sel(j)(i) = log( plogis<dvar_vector>(age,p1,p2) );	   //<dvar_vector> New in 2014 to allow compile with new statslib.h
				}
				break;
				
			case 2:		
				// age-specific selectivity coefficients
				for(i=syr; i<=nyr; i++)
				{
					for(k=sage;k<=nage-1;k++)
					log_sel(j)(i)(k) = sel_par(j)(1)(k-sage+1);
					log_sel(j)(i,nage) = log_sel(j)(i,nage-1);
				}
				break;
				
			case 3:		
				// cubic spline
				log_sel(j)(syr)=cubic_spline( sel_par(j)(1) );
				for(i=syr; i<nyr; i++)
				{
					log_sel(j)(i+1) = log_sel(j)(i);
				}
				break;
				
			case 4:		
				// time-varying cubic spline every year
				for(i=syr; i<=nyr; i++)
				{
					log_sel(j)(i) = cubic_spline(sel_par(j)(i-syr+1));
				}
				
				//jlog_sel(j) = cubic_spline(sel_par(j)(1));
				//t1 = cubic_spline_matrix(sel_par(j).sub(2,jsel_npar(j)));
				//for(i = t1.indexmin(); i <= t1.indexmax(); i++)
				//{
				//	tmp( syr+(i-t1.indexmin()) ) = t1(i);
				//}
				break;
				
			case 5:		
				// time-varying bicubic spline
				ia.fill_seqadd( 0,1./(age_nodes(j)-1) );
				iy.fill_seqadd( 0,1./(yr_nodes(j)-1) );
				// bicubic_spline function is located in stats.cxx library
				bicubic_spline( iy,ia,sel_par(j),tmp2 );
				log_sel(j) = tmp2; 
				break;
				
			case 7:
				// time-varying selectivity based on deviations in weight-at-age
				// CHANGED This is not working and should not be used. (May 5, 2011)
				// SJDM:  I was not able to get this to run very well.
				// AUG 5, CHANGED so it no longer has the random walk component.
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				
				for(i = syr; i<=nyr; i++)
				{
					dvar_vector tmpwt=log(wt_obs(i)*1000)/mean(log(wt_obs*1000.));
					log_sel(j)(i) = log( plogis(tmpwt,p1,p2)+tiny );
				}	 
				break;
				
			case 8:
				//Alternative time-varying selectivity based on weight 
				//deviations (wt_dev) wt_dev is a matrix(syr,nyr+1,sage,nage)
				//p3 is the coefficient that describes variation in log_sel.
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				p3 = sel_par(j,1,3);
				
				for(i=syr; i<=nyr; i++)
				{
					tmp2(i) = p3*wt_dev(i);
					log_sel(j)(i) = log( plogis<dvar_vector>(age,p1,p2)+tiny ) + tmp2(i);	     // <dvar_vector> New in 2014 to allow compile with new statslib.h
				}
				break;
				
			case 11:
				//logistic selectivity based on mean length-at-age
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				
				for(i=syr; i<=nyr; i++)
				{
					dvector len = pow(wt_obs(i)/a,1./b);
					log_sel(j)(i) = log( plogis<dvar_vector>(len,p1,p2) );		  // <dvar_vector> New in 2014 to allow compile with new statslib.h
				}
				break;
				
			case 12:
				//length-specific selectivity coefficients
				//based on cubic spline interpolation
				for(i=syr; i<=nyr; i++)
				{
					dvector len = pow(wt_obs(i)/a,1./b);
					log_sel(j)(i)=cubic_spline( sel_par(j)(1), len );
				}
				break;
				
				
			default:
				log_sel(j)=0;
				break;
				
		}  // switch
		

		//log_sel(j)(syr) = jlog_sel(j)+tmp2(syr);
		//for(i=syr;i<nyr;i++)
		//{
		//	log_sel(j)(i+1)(sage,nage) = log_sel(j)(i)(sage,nage)+tmp(i)+tmp2(i+1);			
		//}
		
		//subtract mean to ensure mean(exp(log_sel))==1
		//substract max to ensure exp(log_sel) ranges from 0-1
		
		//RF MODIFICATION
		//Correction not applied for cases 1,6,7 - was returning values greater than one when exponentiated
		if(isel_type(j)==2 || isel_type(j)==3|| isel_type(j)==4|| isel_type(j)==5) {
				for(i=syr;i<=nyr;i++)
					log_sel(j)(i) -= log( mean(mfexp(log_sel(j)(i)))+tiny );
			//log_sel(j)(i) -= log(max(mfexp(log_sel(j)(i))));
		}	
		//cout<<"log_sel \t"<<j<<"\n"<<log_sel(j)<<"\n \n"<<endl;
		//testing bicubic spline  (SM Checked OCT 25,2010.  Works on the example below.)
		/*ia.fill_seqadd(0,1./(age_nodes(j)-1));
				iy.fill_seqadd(0,1./(yr_nodes(j)-1));
				dvar_matrix tn(1,age_nodes(j),1,yr_nodes(j));
				tn.colfill_seqadd(1,-.5,0.1);
				tn.colfill_seqadd(2,-.4,0.1);
				bicubic_spline(iy,ia,tn,tmp2);
				cout<<ia<<endl;
				cout<<iy<<endl;
				cout<<tn<<endl;
				cout<<tmp2<<endl;
				exit(1);*/
	
	}
	
	if(verbose)cout<<"**** Ok after calcSelectivities ****"<<endl;
	
  }	
	
FUNCTION calcTotalMortality
  {
	/*
	This routine calculates fishing mortality, total mortality
	and annual survival rates (exp(-Z)) for each age in each
	year.
	
	There is a complication in that if there is a survey gear
	then the fishing mortality rate is set to an extremely small
	value: exp(-70.)~3.975e-31, which is effectively 0.
	
	
	The above issue is no longer an issue b/c now estimating Fs
	for years where there is non-zero catch.
	
	SJDM.  Dec 24, 2010.  Adding time-varying natural mortality
	
	CHANGED May 20, 2011  Add cubic spline to the time-varying natural mortality
	
	Jan 5, 2012 Adding catch_type to allow for catch in numbers, weight or spawn.
	In the case of spawn on kelp (roe fisheries), the Fishing mortality does not
	occur on the adult component.  Added if(catch_type(k)!=3) //exclude roe fisheries
	*/
	int i,k,ki;
	dvariable ftmp;
	F.initialize();
	ft.initialize();
	log_ft.initialize();
	
	//Fishing mortality
	ki=1;
	for(k=1;k<=ngear;k++)
	{
		
		for(i=syr;i<=nyr;i++)
		{	
			ftmp=0;
			if( obs_ct(k,i)>0  )
				ftmp = mfexp(log_ft_pars(ki++));
			
			ft(k,i)=ftmp;

						
			//@@@@@@RF temporary code for comparing with 2005 assessment@@@@@@@@@@@@@@
			//~~~~~~~~~~~~~~~FOR P. cod~~~~~~~~~~~~~~~~~~~~~~
			//fix ft at q*obs_effort if cntrl 17 = 1 : FOR GEAR 1 ONLY
			//OVERWRITE ft
			if(cntrl(17)){
				if(i<=eff_nyr) {
					if(k==1){
						ftmp = effort(i)*qc;
						ft(k,i)=ftmp; 
					}
				}
			}
			
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			
			if(catch_type(k)!=3){	//exclude roe fisheries
				F(i)+=ftmp*mfexp(log_sel(k)(i));
			}
		}
	}
	
	//cout<<qc<<endl;
	//cout<<ft(1)<<endl<<endl; 
	
	//Natural mortality (year and age specific)
	//M_tot(syr,nyr,sage,nage);
	M_tot = m;

	// Cubic spline to interpolate log_m_devs (log_m_nodes)
	log_m_devs = 0.;
	if(active(log_m_nodes))
	{
		int nodes = size_count(log_m_nodes);
		dvector im(1,nodes);
		dvector fm(syr+1,nyr);
		im.fill_seqadd(0,1./(nodes-1));
		fm.fill_seqadd(0,1./(nyr-syr));
		vcubic_spline_function m_spline(im,log_m_nodes);
		//m_spline(fm);
		log_m_devs = m_spline(fm);
	}
	
	// Random walk in natural mortality.
	for(i=syr;i<=nyr;i++)
	{
		// if(active(log_m_devs)&&i>syr)
		if(active(log_m_nodes)&&i>syr)
		{
			M_tot(i)=M_tot(i-1)*exp(log_m_devs(i));
		}
	}
	m_bar = mean(M_tot);
	
	Z=M_tot+F;
	S=mfexp(-Z);
	if(verbose) cout<<"**** OK after calcTotalMortality ****"<<endl;
	
  }
	
	
FUNCTION calcNumbersAtAge
  {
	/*
		TODO Need to check the difference between the initialization 
		of the numbers at age here at the margins in comparison to the
		simulation model.
	*/
	
	int i,j;
	N.initialize();
	tbt.initialize(); //total biomass
	dvar_vector lx(sage,nage);
	
	if(cntrl(5)==1){	//If initializing in at unfished conditions
		log_rt(syr) = log(ro);
		
		for(j=sage;j<=nage;j++)
		{
			N(syr,j)=ro*exp(-m_bar*(j-sage)); //RF changed j-1 to j-sage to allow for recruits to be defined as fish of ansage
		}
	}
	if(cntrl(5)==0){
		
		log_rt(syr) = log_avgrec+log_rec_devs(syr);
		
		//@@@@@@@ RF temporary code for comparing to 2005 assessment @@@@@@@@
		if(cntrl(16)) log_rt(syr) = log_avgrec + log_rec_devs(syr) + dt*anom_obs(syr);
		//@@@@@@@@@@@@@@@@@@@@@@
		
		N(syr,sage)=mfexp(log_rt(syr));
		
		for(j=sage+1;j<=nage;j++)
		{
			N(syr,j)=mfexp(log_recinit+init_log_rec_devs(j))*exp(-m_bar*(j-sage));
		}
	}
	
	N(syr,nage)/=(1.-exp(-m_bar));
	//dvector ma = plogis(age,ah,gh); //maturity
	//cout<<"ma = "<<ma<<endl;
	//no=sum(elem_prod(N(syr),ma));//test - spawning numbers
	no=sum(N(syr));
	
	//initial number of sage recruits from year syr+1, nyr;
	for(i=syr+1;i<=nyr;i++){
		log_rt(i)=log_avgrec+log_rec_devs(i);
		
		//@@@@@@@ RF temporary code for comparing to 2005 assessment @@@@@@@@
		if(cntrl(16)) log_rt(i) = log_avgrec+log_rec_devs(i)+dt*anom_obs(i);
		//@@@@@@@@@@@@@@@@@@@@@@
		
		N(i,sage)=mfexp(log_rt(i));
	}
	N(nyr+1,sage)=mfexp(log_avgrec);    //RF doesn't like this projection step - prefers to stick to projection in projection model - this one calculates recruitment inconsistently with projection model
	
	/* Dynamic state variables */
	for(i=syr;i<=nyr;i++)
	{
		N(i+1)(sage+1,nage)=++elem_prod(N(i)(sage,nage-1),S(i)(sage,nage-1));//+1.e-10;
		N(i+1,nage)+=N(i,nage)*S(i,nage);
	}
	
	for(i=syr; i<=nyr; i++) numbers(i) = sum(elem_prod(N(i),mfexp(log_sel(1)(i)))); //vulnerable numbers commercial
	tbt = rowsum(elem_prod(N,wt_obs)); //total biomass
	if(verbose)cout<<"**** Ok after calcNumbersAtAge ****"<<endl;
	
  }
  
FUNCTION calcNumbersBiomass_deldiff
    {
  	int i;
  	int j;
  	numbers.initialize();
  	biomass.initialize();
  	tbt.initialize(); //total biomass same as biomass - needed for graphing comparison with ASM
  	dvariable snat;
  	dvariable sfished;
    	    	
  	 //DD initialisation
	 //Equlibrium mean weight - obtained from weq = Beq/Neq and solving for weq
	 //i.e. weq = [surv(alpha.Neq + rho.Beq + wk.R] / [surv.Neq + R]
	 // with substitutions Neq = Beq/weq and R = Neq(1 - surv)
	 //From SJDM, also used by Sinclair in 2005 p cod assessment
	 snat=mfexp(-m); //unfished survival rate
	
	 wbar=(snat*alpha_g+wk*(1.-snat))/(1.-rho_g*snat);
	  no=ro/(1.-snat); //H&W 1992 p339
          bo=no*wbar;      
         
         /*
      cout<<"m snat wbar ro no bo"<<endl;
       cout<<m<<endl;
        cout<<snat<<endl;
        cout<<wbar<<endl;
        cout<<ro<<endl;
        cout<<no<<endl;
       cout<<bo<<endl; exit(1);  
       */ 
       
          //stock-recruit parameters
          if(cntrl(2)==1){ //Beverton-Holt
          	so=kappa*(ro/bo);
          	beta=(kappa-1)/bo;}
	  
	  if(cntrl(2)==2){ //Ricker
	  	
	  	so=kappa*ro/bo;
		 beta=log(kappa)/bo;
	  }
	
	//Initialise dynamic equations
	surv(syr) = mfexp(-m-ft(1, syr)); 
	
	 if(cntrl(5)==0){
	   	//Unfished and not at equilibrium - Initialise as for ASM
		   log_rt(syr) = log_avgrec+log_rec_devs(syr);

		//@@@@@@@ RF temporary code for comparing to 2005 assessment @@@@@@@@
		if(cntrl(16)) log_rt(syr) = log_avgrec + log_rec_devs(syr) + dt*anom_obs(syr);
		//@@@@@@@@@@@@@@@@@@@@@@

	   	N(syr,sage)=mfexp(log_rt(syr));
		for(j=sage+1;j<=nage;j++)
		{
			N(syr,j)=mfexp(log_recinit+init_log_rec_devs(j))*exp(-m*(j-sage));
		}
		N(syr,nage)/=(1.-exp(-m));
		
		numbers(syr) = sum(N(syr));
		biomass(syr) = sum(elem_prod(N(syr),wt_obs(syr))); //total biomass
		annual_mean_wt(syr) = biomass(syr)/numbers(syr);
		surv(syr) = mfexp(-m-ft(1, syr)); 
   	   } //end if
			
	//check these two options are the same in the absence of fishing mortality	
	if(cntrl(5)==1){
		//start at equlibrium unfished
		numbers(syr)= no;
		biomass(syr) = bo;
		annual_mean_wt(syr) = wbar;
		log_rt(syr) = log(ro);
		
		//@@@@@@@ RF temporary code for comparing to 2005 assessment @@@@@@@@
		if(cntrl(16)) log_rt(syr) = log(ro) + dt*anom_obs(syr);
  		//@@@@@@@@@@@@@@@@@@@@@@
   	  }  //end if
   	 
   	  if(cntrl(5)==2){
   	  	//start at equlibrium with fishing mortality - different approach to ASM
   	  	sfished = surv(syr); //equilibrium survivorship at initial fishing mortality (gear 1 commercial fishery)
   	  	annual_mean_wt(syr) = (sfished*alpha_g + wk*(1-sfished))/(1-rho_g*sfished);
   	  	biomass(syr) = -(annual_mean_wt(syr)*(wk*so-1)+sfished*(alpha_g+rho_g* annual_mean_wt(syr)))/(beta*(sfished*alpha_g+sfished*rho_g* annual_mean_wt(syr)- annual_mean_wt(syr)));
   	  	numbers(syr) = biomass(syr)/annual_mean_wt(syr);
   	  	sbt(syr) = biomass(syr);
   	 	   	  	
   	  }//end if
   	  
    	for(i=syr+1;i<=nyr;i++){
    		log_rt(i)=log_avgrec+log_rec_devs(i); 
  		//@@@@@@@ RF temporary code for comparing to 2005 assessment @@@@@@@@
  		if(cntrl(16)) log_rt(i) = log_avgrec+log_rec_devs(i)+dt*anom_obs(i);
  		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  		
  		 //Update biomass and numbers	
	   	biomass(i) =(surv(i-1)*(rho_g*biomass(i-1)+alpha_g*numbers(i-1))+wk*mfexp(log_rt(i)));
	    	numbers(i)=surv(i-1)*numbers(i-1)+mfexp(log_rt(i));
	    	annual_mean_wt(i)=biomass(i)/numbers(i);		//calculate predicted weight in dynamics - possible option to fit to it
		surv(i) = mfexp(-m-ft(1, i));
   	}	
  	  //RF doesn't like this projection step - prefers to stick to projection in projection model - this one calculates recruitment inconsistently with projection model
  	dvariable rnplus=mfexp(log_avgrec); //assume recruits nyr+1 average - same as for ASM
  	biomass(nyr+1)=(surv(nyr)*(rho_g*biomass(nyr)+alpha_g*numbers(nyr))+wk*rnplus); 
  	numbers(nyr+1)=surv(nyr)*numbers(nyr)+rnplus;
  	
  	sbt = biomass; //set spawning biomass to biomass
  	tbt = biomass; //total biomass - for testing
  	//cout<<ft(1)<<endl<<endl;
  	if(verbose)cout<<"**** Ok after calcNumbersBiomass_deldiff ****"<<endl;
  	
  }

FUNCTION calcAgeProportions
  {
	/*This function loops over each gear and year
	and calculates the predicted proportions at age
	sampled based on the selectivity of that gear and
	the numbers-at-age in the population.*/
	
	int i,k,iyr,ig;
	for(k=1;k<=na_gears;k++)
	{
		for(i=1;i<=na_nobs(k);i++)
		{
			iyr=A(k,i,a_sage(k)-2);	//index for year
			ig=A(k,i,a_sage(k)-1);	//index for gear
			if(iyr>nyr)break;		//trap for retrospective analysis
			
			A_nu(k,i,a_sage(k)-2)=iyr;
			A_nu(k,i,a_sage(k)-1)=ig;
			Ahat(k,i,a_sage(k)-2)=iyr;
			Ahat(k,i,a_sage(k)-1)=ig;
			Ahat(k)(i)(a_sage(k),a_nage(k))=Chat(ig)(iyr)(a_sage(k),a_nage(k))
										/sum(Chat(ig)(iyr)(a_sage(k),a_nage(k)));
		}
		
	}
	
	if(verbose)cout<<"**** Ok after calcAgeProportions ****"<<endl;

  }	

FUNCTION calcFisheryObservations
  {
	/*
	Dec 6, 2010.  Modified ct calculations to include
				  empirical weight at age data (wt_obs);
				
	Jan 16, 2011. 	modified this code to get age-comps from surveys, rather than 
					computing the age-comps in calc_fisheries_observations
	
	Jan 6, 2012. 	modified code to allow for catch observations in numbers,
					biomass, and harvest of roe.  
	
	Jun 22, 2012.	added eta variable for catch residuals.
	
	*/
	
	/*
		FIXED Reconcile the difference between the predicted catch 
		here and in the simulation model.
	*/
	int i,k;
	ct.initialize();
	eta.initialize();
	
	for(i=syr;i<=nyr;i++)
	{
		for(k=1;k<=ngear;k++)
		{
			dvar_vector log_va=log_sel(k)(i);
			
			//SJDM Jan 06, 2012 Modification as noted above.
			if(obs_ct(k,i)>0)
			{
				dvar_vector fa=ft(k,i)*mfexp(log_va);
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
				switch(catch_type(k))
				{
					case 1:	//catch in weight
						ct(k,i) = Chat(k,i)*wt_obs(i);
					break;
					case 2:	//catch in numbers
						ct(k,i) = sum(Chat(k,i));
					break;
					case 3:	//catch in roe that does not contribute to SSB
						dvariable ssb = elem_prod(N(i),exp(-Z(i)*cntrl(13)))*fec(i);
						ct(k,i) = ( 1.-mfexp(-ft(k,i)) )*ssb;
					break;
				}
				eta(k,i) = log(obs_ct(k,i))-log(ct(k,i));
			}
			else
			{/*If there is no commercial fishery then set Chat equal to 
			   the expected proportions at age.*/
				dvar_vector fa = mfexp(log_va);
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
			}
			
			
			/*
			Changed this on Jan 16, 2011 as it was preventing
			convergence for simuation with zero error due to the tiny
			constant added to F.
			
			need to add a tiny constant to deal with catch-age data 
			for non-extractive survey age comps
			dvar_vector fa=ft(k,i)*mfexp(log_va)+1.e-30;
			
			//Catch-at-age by commercial gear
			if(fsh_flag(k))
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
			
			//Catch weight by gear
			//ct(k,i)=Chat(k,i)*wa;  
			ct(k,i)=Chat(k,i)*wt_obs(i);  //SM Dec 6, 2010*/
		}
	}
	if(verbose)cout<<"**** Ok after calcFisheryObservations ****"<<endl;

  }	
	
FUNCTION calcSurveyObservations
  {
	/*This code needs to be modified to accomodate
	multiple surveys or block changes in survey q.
	
	Oct 31, 2010, added retrospective counter.
	
	Nov 22, 2010, adding multiple surveys. Still need to check with retrospective option
	
	Nov 30, 2010, adjust the suvery biomass by the fraction of Z that has occurred 
	when the survey was conducted. For herring spawning biomass this would be after the 
	fishery has taken place.
	
	Dec 6, 2010, modified predicted survey biomass to accomodate empirical weight-at-age 
	data (wt_obs).
	
	May 11, 2011.  CHANGED Vivian Haist pointed out an error in survey biomass comparison.
	The spawning biomass was not properly calculated in this routine. I.e. its different 
	than the spawning biomass in the stock-recruitment routine. (Based on fecundity which
	changes with time when given empirical weight-at-age data.)
	
	Jan 6, 2012.  CHANGED corrected spawn survey observations to include a roe 
	fishery that would remove potential spawn that would not be surveyed.
	
	*/
	/*
		CHANGED add capability to accomodate priors for survey q's.
		DONE
	*/
	int i,j,ii,k,kk;
	
	//survey abudance index residuals
	epsilon.initialize();
	pit.initialize();
	
	for(i=1;i<=nit;i++)
	{	
		int nx=0;		//counter for retrospective analysis
		dvar_matrix V(1,nit_nobs(i),sage,nage);  //vulnerable units for survey comparison
		V.initialize(); 
		
		//RF testing - to delete
		//dvar_matrix V2(1,nit_nobs(i),sage,nage);  
		//V2.initialize();
		
		//dvar_matrix VB(1,nit_nobs(i),sage,nage); //vulnerable biomass
		for(j=1;j<=nit_nobs(i);j++)
		{
			ii=iyr(i,j);
			k=igr(i,j);
			if(ii>nyr) break;	//trap for retrospective analysis.
			dvar_vector log_va=log_sel(k)(ii);
			
			//Adjust survey biomass by the fraction of the mortality that 
			//occurred during the time of the survey.
			dvar_vector Np = elem_prod(N(ii),exp( -Z(ii)*it_timing(i,j))); // if it_timing is 0 Np - N(ii)
			//cout<<N(ii)<<endl;
			//cout<<Np<<endl<<endl; //exit(1);
			//get fishing mortality rate on spawn.
			dvariable ftmp = 0;
			for(kk=1;kk<=ngear;kk++)
				if(catch_type(kk)==3)
					ftmp += ft(kk,ii);
					
			switch(survey_type(i))
			{
				case 1:
					V(j)=elem_prod(Np,mfexp(log_va));
				break;
				case 2:
					V(j)=elem_prod(elem_prod(Np,mfexp(log_va)),wt_obs(ii));
					
				break;
				case 3:
					//SJDM Jan 6, 2012 Modified for roe fishery
					V(j)=elem_prod(Np,fec(ii))*exp(-ftmp);
				break;
			}
			
			//VB(j)=elem_prod(V(j),wt_obs(ii));		//SM Dec 6, 2010
			
			//If the survey is a spawn index, then need to account
			//for changes in fecundity.
			
			if(iyr(i,j)<=nyr) nx++;
		}
		dvar_vector t1 = rowsum(V);//V*wa;
		//cout<<"V\n"<<V<<endl<<endl;
		//cout<<"t1\n"<<t1<<endl<<endl;
		
		//RF testing - to delete
		//dvar_vector t2 = rowsum(V2);
		//cout<<"t2"<<endl<<t2<<endl;
		
		//See Ludwig & Walters 1994
		//Note this is incorrect if each survey has different weights.
		dvar_vector zt=log(it(i).sub(1,nx))-log(t1(1,nx));
		//cout<<"zt\n"<<t1(1,nx)<<endl;
		epsilon(i).sub(1,nx) = zt-mean(zt);
		q(i) = exp(mean(zt));
		pit(i).sub(1,nx)=t1(1,nx)*q(i);	//predicted index
		
		//TODO, this might not be working correctly, simulation test it.
		if(q_prior(i)==2)
		{
			//random walk in q
			epsilon(i)=0;
			dvar_vector fd_zt=first_difference(zt);
			epsilon(i).sub(1,nx-1) = fd_zt-mean(fd_zt);
			//dvar_vector qt(1,nx);
			qt(i,1) = exp(zt(1));
			for(j=2;j<=nx;j++)
				qt(i,j) = qt(i,j-1)*exp(fd_zt(j-1));
				
			pit(i).sub(1,nx)=elem_prod(t1(1,nx),qt(i)(1,nx));
			//cout<<sum(epsilon(i))<<endl;
			//exit(1);
		}
	}
   if(verbose)cout<<"**** Ok after calcSurveyObservations ****"<<endl;
	
  }
  
FUNCTION calc_stock_recruitment
  {
	/*
	The following code is used to derive unfished
	spawning stock biomass bo and the stock-
	recruitment parameters for the:
	Beverton-Holt Model
		-Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau)
	Ricker Model
		-Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau)
		
	Dec 6, 2010.  Modified code to allow empirical weight-at-age data
				This required a fecundity-at-age matrix.  Need to 
				project spawning biomass into the future.
	CHANGED bo should be a function of the average natural mortality
	TODO update phib calculation if age-specific M's are used.
	
	May 6, 2010.  Changed phib calculation based on average M 
	in cases where M varies over time. Otherwise bo is biased based
	on the initial value of M.
	
	CHANGED Need to adjust spawning biomass to post fishery numbers.
	CHANGED Need to adjust spawners per recruit (phib) to average fecundity.
	
	Jan 6, 2012.  Need to adjust stock-recruitment curvey for reductions 
	in fecundity associated with removal of roe from a spawn on kelp fishery.
	*/ 
	int i,j,k;
	dvariable   phib, phitb;
	dvariable   tau = sqrt(1.-rho)*varphi;
	dvar_vector     ma(sage,nage);
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector     lx(sage,nage); lx(sage) = 1.0;
	dvar_vector     lw(sage,nage); lw(sage) = 1.0;
	
	// -steps (1),(2),(3)
	dvar_matrix t_M_tot = trans(M_tot);
	for(j=sage; j<=nage;j++)
	{
		ma(j) = mean(t_M_tot(j));
		if(j>sage)
		{
			lx(j) = lx(j-1)*mfexp(-ma(j-1));
		} 
		lw(j) = lx(j) * mfexp(-ma(j)*cntrl(13));
	}
	lx(nage) /= 1.0 - mfexp(-ma(nage));
	lw(nage) /= 1.0 - mfexp(-ma(nage));
	phib      = lw * avg_fec;
	phitb      = lw * avg_wt;	    //RF this is total biomass per recruit 
	
	// step (4)
	for(i=syr;i<=nyr;i++)
	{
		sbt(i) = elem_prod(N(i),exp(-Z(i)*cntrl(13)))*fec(i);
		//Adjustment to spawning biomass for roe fisheries
		for(k=1;k<=ngear;k++)
		{
			if(catch_type(k)==3)
			{
				sbt(i) *= mfexp(-ft(k,i));
			}
		}
	}
	// sbt(nyr+1) = elem_prod(N(nyr+1),exp(-M_tot(nyr))) * fec(nyr+1);	
	// Changes as per Nathan Taylor's email.
	sbt(nyr+1) = elem_prod(N(nyr+1),exp(-M_tot(nyr))) * avg_fec;	
	dvar_vector tmp_st=sbt(syr,nyr-sage).shift(syr+sage);

	so = kappa/phib;		//max recruits per spawner
       	bo = ro*phib;  			//unfished spawning biomass
	
	switch(int(cntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			beta   = (kappa-1.)/bo;
			tmp_rt = elem_div(so*tmp_st,1.+beta*tmp_st);
			break;
		case 2:
			//Ricker model
			beta   = log(kappa)/bo;
			tmp_rt = elem_prod(so*tmp_st,exp(-beta*tmp_st));
		break;
	}
	
	//residuals in stock-recruitment curve
	rt    = mfexp(log_rt(syr+sage,nyr));
	delta = log(rt)-log(tmp_rt)+0.5*tau*tau;
	
	if(verbose)cout<<"**** Ok after calc_stock_recruitment ****"<<endl;
	
  }

FUNCTION calcFisheryObservations_deldiff
  {
	int i,k;
	ct.initialize();
	eta.initialize();
	for(i=syr;i<=nyr;i++)
	{
		for(k=1;k<=ngear;k++)
		{
			//SJDM Jan 06, 2012 Modification as noted above.
			if(obs_ct(k,i)>0)
			{
								
				switch(catch_type(k))
				{
					case 1:	//catch in weight
						ct(k,i) = biomass(i)*(1-mfexp(-m-ft(k,i)))*(ft(k,i)/(m + ft(k,i)));
					break;
					case 2:	//catch in numbers
						ct(k,i) = numbers(i)*(1-mfexp(-m-ft(k,i)))*(ft(k,i)/(m + ft(k,i)));
					break;
					case 3:	//roe - NOT IMPLEMENTED FOR DELAY DIFFERENCE MODEL
						cout<<"WARNING: CATCH TYPE 3 (ROE FISHERY) NOT IMPLEMENTED FOR DELAY DIFFERENCE MODEL:"<<endl;
						cout<<"USE THE AGE-STRUCTURED MODEL - TERMINATING PROGRAM"<<endl; exit(1);
					break;
				}
				eta(k,i) = log(obs_ct(k,i))-log(ct(k,i));
			}
						
		}
	}
    if(verbose)cout<<"**** Ok after calcFisheryObservations_deldiff ****"<<endl;
   
  }	
	
FUNCTION calcSurveyObservations_deldiff
  {
	int i,j,ii,k,kk;
	
	//survey abudance index residuals
	epsilon.initialize();
	pit.initialize();
		
	for(i=1;i<=nit;i++)
	{	
		int nx=0;		//counter for retrospective analysis
		dvar_vector V(1,nit_nobs(i));  //vulnerable units for survey comparison
		V.initialize();
		
		//dvar_matrix VB(1,nit_nobs(i),sage,nage); //vulnerable biomass
		for(j=1;j<=nit_nobs(i);j++)
		{
			ii=iyr(i,j);
			k=igr(i,j);
			if(ii>nyr) break;	//trap for retrospective analysis.
			
			//Adjust survey biomass by the fraction of the mortality that 
			//occurred during the time of the survey.
			dvariable z = ft(k,ii)+m;
			dvariable Np = numbers(ii) * exp( -z * it_timing(i,j));
			dvariable Bp = biomass(ii) *exp( -z * it_timing(i,j));
			
			switch(survey_type(i))
			{
				case 1:
					V(j)=Np;
				break;
				case 2:
					V(j)=Bp;
				break;
				case 3:
					V(j)=Bp;
				break;
			}
			
			if(iyr(i,j)<=nyr) nx++;
		}
		dvar_vector t1 = V;
		//cout<<"t1\n"<<t1<<endl; 
		
		//See Ludwig & Walters 1994
		//Note this is incorrect if each survey has different weights.
		dvar_vector zt=log(it(i).sub(1,nx))-log(t1(1,nx));
		//cout<<"zt\n"<<t1(1,nx)<<endl;
		epsilon(i).sub(1,nx) = zt-mean(zt);
		q(i) = exp(mean(zt));
		pit(i).sub(1,nx)=t1(1,nx)*q(i);	//predicted index
		
	}
   if(verbose)cout<<"**** Ok after calcSurveyObservations_deldiff ****"<<endl;
	
  }
  
FUNCTION calc_stock_recruitment_deldiff
  {
	/*
	The following code is used to derive unfished
	spawning stock biomass bo and the stock-
	recruitment parameters for the:
	Beverton-Holt Model
		-Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau)
	Ricker Model
		-Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau)
	*/
	
	//for the delay difference model so and beta are calculated in the function calcNumbersBiomass_deldiff
	int i,k;
	dvariable tau = sqrt(1.-rho)*varphi;
	dvar_vector tmp_rt(syr,nyr);	  //need to calc recruits in all years, then offset by kage yrs
	dvar_vector tmp_st(syr,nyr);
	int iicount=0;
	
	switch(int(cntrl(2)))
	{
		 case 1:
			//Beverton-Holt model
					
			for(i=syr; i<=nyr; i++){
				iicount++;
				tmp_st(i)=biomass(i);
				
				if(iicount <=kage){ 
					tmp_rt(i) =  so*tmp_st(syr)/(1.+beta*tmp_st(syr));  
				}else
					tmp_rt(i) = so*tmp_st(i-kage)/(1.+beta*tmp_st(i-kage));
				}
		break;
		case 2:
			//Ricker model
			for(i=syr; i<=nyr; i++){
				iicount++;
				tmp_st(i)=biomass(i);
				
				if(iicount <=kage){ 
					tmp_rt(i) =so*tmp_st(syr)*exp(-beta*tmp_st(syr));
				}else
					tmp_rt(i) =so*tmp_st(i-kage)*exp(-beta*tmp_st(i-kage));
			}
		break;
	}
	
	//residuals in stock-recruitment curve
	//    = mfexp(log_rt(syr+sage,nyr));  //rf changed bounds  log_rt(syr+sage,nyr
	//delta = log(rt)-log(tmp_rt(syr+sage,nyr))+0.5*tau*tau;
	rt    = mfexp(log_rt(syr+sage,nyr));
	delta = log(rt)-log(tmp_rt(syr+sage,nyr))+0.5*tau*tau;
	
	if(verbose)cout<<"**** Ok after calc_stock_recruitment_deldiff ****"<<endl;
	
  }
  
//RF added for comparison to delay difference model
FUNCTION calc_annual_mean_weight
  {
  	//******** FOR GEAR 1 ONLY (commercial catch) **********:
  	//1. calculate vulnerable biomass bat each year t = sum(Nat,Vat,Wat)
  	//2. calculate annual mean weight at age from bat/sum(Nat,Vat)
  	//This is only called in the Catch Age model - annual mean weight is calculated as part of calcNumbersBiomass_deldiff for the delay diff model
  	int i;
  	if(!delaydiff) {
		for(i=syr; i<=nyr; i++)
		{
			vbcom(i) = sum(elem_prod(N(i),elem_prod(mfexp(log_sel(1)(i)),wt_obs(i)))); //vulnerable biomass commercial
			vncom(i) = sum(elem_prod(N(i),mfexp(log_sel(1)(i)))); //vulnerable numbers commercial
			annual_mean_wt(i) = vbcom(i)/vncom(i);
		}
	}else {
		 vbcom(syr,nyr) = biomass(syr,nyr); //vulnerable biomass commercial - for calculating ut in report section
		 vncom(syr,nyr) = numbers(syr,nyr); //vulnerable numbers commercial
	}
  	
  	//cout<<"annual_mean_wt = "<<annual_mean_wt<<endl;
  	//exit(1);
  	
  	if(verbose)cout<<"**** Ok after calc_annual_mean_weight ****"<<endl;
  }
  
	
FUNCTION calc_objective_function
  {
	//Dec 20, 2010.  SJDM added prior to survey qs.
	/*q_prior is an ivector with current options of 0 & 1.
	0 is a uniform density (ignored) and 1 is a normal
	prior density applied to log(q).*/

	/*
	There are several components to the objective function
	Likelihoods (nlvec):
		-1) likelihood of the catch data
		-2) likelihood of the survey abundance index
		-3) likelihood of the survey age comps
		-4) likelihood of the fishery age comps
		-5) likelihood for stock-recruitment relationship
		-6) likelihood for fishery selectivities
		-7) likelihood for fitting to mean weight (RF added this for P cod)
	*/
	int i,j,k;
	double o=1.e-10;

	dvar_vector lvec(1,6); //RF comment - this seems to only have a selectivity likelihood in it
	lvec.initialize();
	nlvec.initialize();
	
	//1) likelihood of the catch data (retro)
	double sig_c =cntrl(3);
	if(last_phase())sig_c=cntrl(4);
	if(active(log_ft_pars)) 
	for(k=1;k<=ngear;k++){
		for(i=syr;i<=nyr;i++)
		{
			if(obs_ct(k,i)!=0)
				nlvec(1,k)+=dnorm(eta(k,i),0.0,sig_c);
				//nlvec(1,k)+=dnorm(log(ct(k,i)),log(obs_ct(k,i)),sig_c);
		}
		//if(active(log_ft_pars))
		//	nlvec(1,k)=dnorm(log(obs_ct(k).sub(syr,nyr)+o)-log(ct(k).sub(syr,nyr)+o),sig_c);
	}
	
	
	//2) likelihood of the survey abundance index (retro)
	for(k=1;k<=nit;k++)
	{
		dvar_vector sig = (sqrt(rho)*varphi)/it_wt(k);
		nlvec(2,k)=dnorm(epsilon(k),sig);
	}
	
	//3) likelihood for age-composition data
	for(k=1;k<=na_gears;k++)
	{	
		if(na_nobs(k)>0){
			int naa=0;
			int iyr;
			//retrospective counter
			for(i=1;i<=na_nobs(k);i++)
			{
				iyr=A(k,i,a_sage(k)-2);	//index for year
				if(iyr<=nyr)naa++;
			}
			
			dmatrix O=trans(trans(A(k)).sub(a_sage(k),a_nage(k))).sub(1,naa);
			dvar_matrix P=trans(trans(Ahat(k)).sub(a_sage(k),a_nage(k))).sub(1,naa);
			//dvar_matrix nu=trans(trans(Ahat(k)).sub(a_sage(k),a_nage(k))).sub(1,naa); //residuals
			dvar_matrix nu(O.rowmin(),O.rowmax(),O.colmin(),O.colmax()); //residuals
			nu.initialize();
			
			//CHANGED add a switch statement here to choose form of the likelihood
			//RF Added 3rd switch: 3=no likelihood (for "ageless" model) - currently still need a dummy matrix to avoid read data errors
			switch(int(cntrl(14)))
			{
				case 1:
					nlvec(3,k) = dmvlogistic(O,P,nu,age_tau2(k),cntrl(6));
				break;
				case 2:
					nlvec(3,k) = dmultinom(O,P,nu,age_tau2(k),cntrl(6));
				break;
				case 3:
					nlvec(3,k) = 0.; //RF added this May 22 2013 for 'ageless' model
				break;
			}
			
			for(i=1;i<=naa/*na_nobs(k)*/;i++)
			{
				iyr=A(k,i,a_sage(k)-2);	//index for year
				A_nu(k)(i)(a_sage(k),a_nage(k))=nu(i);
			}
		}
	}
	
	//4) likelihood for stock-recruitment relationship
	dvariable tau = sqrt(1.-rho)/varphi;
	if(active(theta(1))){
			nlvec(4,1)=dnorm(delta,tau);
	}
	
	//5-6) likelihood for selectivity paramters
	if(!delaydiff){
		for(k=1;k<=ngear;k++)
		{
			if(active(sel_par(k))){
				//if not using logistic selectivity then
				//CHANGED from || to &&  May 18, 2011 Vivian
				if( isel_type(k)!=1 && 
					isel_type(k)!=7 && 
					isel_type(k)!=8 &&
					isel_type(k)!=11 )  
				{
					for(i=syr;i<=nyr;i++)
					{
						//curvature in selectivity parameters
						dvar_vector df2=first_difference(first_difference(log_sel(k)(i)));
						nlvec(5,k)+=sel_2nd_diff_wt(k)/(nage-sage+1)*norm2(df2);

						//penalty for dome-shapeness
						for(j=sage;j<=nage-1;j++)
							if(log_sel(k,i,j)>log_sel(k,i,j+1))
								nlvec(6,k)+=sel_dome_wt(k)
											*square(log_sel(k,i,j)-log_sel(k,i,j+1));
					}
				}
			}
		}//end for(k=1;k<=ngear;k++) 
	
	
		// CONSTRAINT FOR SELECTIVITY DEV VECTORS
		// Ensure vector of sel_par sums to 0. (i.e., a dev_vector)
		// TODO for isel_type==2 ensure mean 0 as well (ie. a dev_vector)

		for(k=1;k<=ngear;k++)
		{
			if( active(sel_par(k)) &&
				isel_type(k)!=1 &&
				isel_type(k)!=7 &&
				isel_type(k)!=8 &&
				isel_type(k)!=11 )
			{
				dvariable s=0;
				if(isel_type(k)==5)  //bicubic spline version ensure column mean = 0
				{
					dvar_matrix tmp = trans(sel_par(k));
					for(j=1;j<=tmp.rowmax();j++)
					{
						s=mean(tmp(j));
						lvec(1)+=1000.*s*s;
					}
				}
				if( isel_type(k)==4 ||
					isel_type(k)==3 || 
					isel_type(k)==12 )
				{
					dvar_matrix tmp = sel_par(k);
					for(j=1;j<=tmp.rowmax();j++)
					{
						s=mean(tmp(j));
						lvec(1)+=1000.*s*s;
					}
				}
			}
		   }//end for(k=1;k<=ngear;k++) 
	}//end if(!delaydiff)
	
	//7) likelihood for mean weight - RF added for P cod + PJS additions
	dvar_vector pred_annual_mean_wt(1, n_weight);
	dvar_vector epsilon_wt(1, n_weight);
	epsilon_wt.initialize();
	int stepy;
 
	for(i=1; i<=n_weight; i++)
	{
     		stepy=obs_annual_mean_wt(i,1);
           	pred_annual_mean_wt(i) = annual_mean_wt(stepy);
           	epsilon_wt(i) = log(pred_annual_mean_wt(i)) - log(obs_annual_mean_wt(i,2));    
	}
	if(cntrl(15)) nlvec(7) = dnorm(epsilon_wt,weight_sig); //fit to annual mean weight if cntrl 15 is switched on  #This is the vector version of dnorm. Gets summed at L2232
 	//end of PJS additions

	/*
	PRIORS for estimated model parameters from the control file.
	*/
	dvar_vector priors(1,npar);
	dvariable ptmp; priors.initialize();
	for(i=1;i<=npar;i++)
	{
		if(active(theta(i)))
		{
			switch(theta_prior(i))
			{
			case 1:		//normal
				ptmp=dnorm(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
				ptmp=dlnorm(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			case 3:		//beta distribution (0-1 scale)
				double lb,ub;
				lb=theta_lb(i);
				ub=theta_ub(i);
				ptmp=dbeta((theta(i)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
				break;
				
			case 4:		//gamma distribution
				ptmp=dgamma(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			default:	//uniform density
				ptmp=1./(theta_control(i,3)-theta_control(i,2));
				break;
			}
			priors(i)=ptmp;	
		}
	}
	
	//Priors for survey q based on control file.
	dvar_vector qvec(1,nits);

	qvec.initialize();
	for(i=1;i<=nits;i++)
	{
		if(q_prior(i)==1) 
		{
			qvec(i)=dnorm(log(q(i)),q_mu(i),q_sd(i));
		}
	}
	
	//** Legacy **  By accident took Rick Methot's bag from Nantes.
	//301 787 0241  Richard Methot cell phone.
	//ibis charles du gaulle at
	//01 49 19 19 20
	
	/*
	The following are penalties that are invoked in early
	phases to ensure reasonable solutions are obtained,
	and to reduce the sensitivity of initial starting
	conditions.  Penalties include:
		-1) keep average fishing mortality rate near 
			0.2 and in the last phase relax this constraint.
		-3) normal prior for log rec devs with std=50.
	*/
	
	dvar_vector pvec(1,7);
	pvec.initialize();
	
	//Penalties to regularize the solution for fishing mortality rates
	dvariable log_fbar = mean(log_ft_pars);
	if(last_phase())
	{
		pvec(1) = dnorm(log_fbar,log(cntrl(7)),cntrl(9));
		//Penalty for log_rec_devs (large variance here)
		pvec(4) = dnorm(log_rec_devs,2.0);
		pvec(5) = dnorm(init_log_rec_devs,2.0);
	}
	else
	{
		pvec(1) = dnorm(log_fbar,log(cntrl(7)),cntrl(8));
		//Penalty for log_rec_devs (CV ~ 0.0707) in early phases
		pvec(4)=100.*norm2(log_rec_devs);
		pvec(5)=100.*norm2(init_log_rec_devs);
	}
	
	//Priors for deviations in natural mortality rates
	//if(active(log_m_devs))
	if(active(log_m_nodes))
	{
		double std_mdev = cntrl(11);
		dvar_vector fd_mdevs=first_difference(log_m_devs);
		pvec(2) = dnorm(fd_mdevs,std_mdev);
		pvec(2) += 0.5*norm2(log_m_nodes);
	}
	
	
	if(verbose)
	{
		cout<<"nlvec\t"<<nlvec<<endl;
		cout<<"lvec\t"<<lvec<<endl;
		cout<<"priors\t"<<priors<<endl;
		cout<<"penalties\t"<<pvec<<endl;
	}
	f=sum(nlvec)+sum(lvec)+sum(priors)+sum(pvec)+sum(qvec);
	//cout<<"Objective function value = "<<f<<endl;
	nf++;
	if(verbose)cout<<"**** Ok after calc_objective_function ****"<<endl;
	
  }

//RF - this one agrees with slow msy code	
FUNCTION void equilibrium(const double& fe,const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dvector& va,double& re,double& ye,double& be,double& phiq,double& dphiq_df, double& dre_df)
  {
	/*
	This is the equilibrium age-structured model that is 
	conditioned on fe (the steady state fishing mortality rate).
	
	In the case of multiple fisheries, fe is to be considered as the
	total fishing mortality rate and each fleet is given a specified
	allocation based on its selectivity curve.  The allocation to 
	each fleet must be specified a priori.
	
	args:
	fe	-steady state fishing mortality
	ro	-unfished sage recruits
	kap	-recruitment compensation ration
	m	-instantaneous natural mortality rate
	age	-vector of ages
	wa	-mean weight at age
	fa	-mean fecundity at age
	va	-mean vulnerablity at age for fe gear.
	ak	-allocation of total ye to gear k.
	
	Modified args:
	re	-steady state recruitment
	ye	-steady state yield
	be	-steady state spawning biomass
	phiq		-per recruit yield
	dre_df		-partial of recruitment wrt fe
	dphiq_df	-partial of per recruit yield wrt fe
	
	FIXME add Ricker model to reference points calculations.
	FIXME partial derivatives for dphif_df need to be fixed when cntrl(13)>0.
	*/
	int i;
	
	int nage=max(age);
	int sage=min(age);
	dvector lx=pow(exp(-m),age-double(sage));
	lx(nage)/=(1.-exp(-m));
	dvector lz=lx;
	dvector za=m+fe*va;
	dvector sa=1.-exp(-za);
	dvector qa=elem_prod(elem_div(va,za),sa);
	
	double phie = lx*fa;		//eggs per recruit
	double so = kap/phie;
	double beta = (kap-1.)/(ro*phie);
	
	
	double dlz_df = 0, dphif_df = 0;
	dphiq_df=0; dre_df=0;
	for(i=sage; i<=nage; i++)
	{
		if(i>sage) lz[i]=lz[i-1]*exp(-za[i-1]);
		if(i>sage) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*va[i-1]*exp(-za[i-1]);
		if(i==nage){ //6/11/2007 added plus group.
					lz[i]/=(1.-mfexp(-za[i]));
					
					dlz_df=dlz_df/(1.-mfexp(-za[i]))
							-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
					/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i])));
				}
		dphif_df=dphif_df+fa(i)*dlz_df;
		dphiq_df=dphiq_df+wa(i)*qa(i)*dlz_df+(lz(i)*wa(i)*va(i)*va(i))/za(i)*(exp(-za[i])-sa(i)/za(i));
	}
	//CHANGED need to account for fraction of mortality that occurs
	//before the spawning season in the recruitment calculation.
	//cout<<"lz\t"<<elem_prod(lz,exp(-za*cntrl(13)))<<endl;
	//exit(1);
	//double phif = lz*fa;
	double phif = elem_prod(lz,exp(-za*cntrl(13)))*fa;
	phiq=sum(elem_prod(elem_prod(lz,wa),qa));
	re=ro*(kap-phie/phif)/(kap-1.);
	//cout<<fe<<" spr ="<<phif/phie<<endl;
	if(re<=0) re=0;
	dre_df=(ro/(kap-1.))*phie/square(phif)*dphif_df;
	ye=fe*re*phiq;
	be=re*phif;	//spawning biomass
	
	//cout<<"Equilibrium\n"<<ro<<"\n"<<re<<"\n"<<ye<<endl;
		
  }
	
FUNCTION void calc_reference_points()
  {
	/**
	\file iscam.tpl
	\author Steven Martell
	Uses Newton_Raphson method to determine Fmsy and MSY 
	based reference points.  
	
	Code check: appears to find the correct value of MSY
	in terms of maximizing ye.  Check to ensure rec-devs
	need a bias correction term to get this right.
	
	Modification for multiple fleets:
		Need to pass a weighted average vector of selectivities
		to the equilibrium routine, where the weights for each
		selectivity is based on the allocation to each fleet.
		
		Perhaps as a default, assign an equal allocation to each
		fleet.  Eventually,user must specify allocation in 
		control file.  DONE
		
		Use selectivity in the terminal year to calculate reference
		points.
	
	June 8, 2012.  SJDM.  Made the following changes to this routine.
		1) changed reference points calculations to use the average
		   weight-at-age and fecundity-at-age.
		2) change equilibrium calculations to use the catch allocation
		   for multiple gear types. Not the average vulnerablity... this was wrong.
	*/
	int i,j;
	double re,ye,be,ve,phiq,dphiq_df,dre_df,fe;
	double dye_df,ddye_df,d2ye_df2,spr;
	
	/* Initial guess for fmsy */
	fe = 1.0*value(m_bar);
	
	/*Calculate average vulnerability*/
	dmatrix va(1,ngear,sage,nage);
	dvector va_bar(sage,nage);
	va_bar.initialize();
	
	
	/*CHANGED Allow for user to specify allocation among gear types.*/
	/*FIXME:  this allocation should be on the catch on the vulnerabilities*/
	
	for(j=1;j<=ngear;j++)
	{
		va_bar+=allocation(j)*value(exp(log_sel(j)(nyr)));
		va(j) = value(exp(log_sel(j)(nyr)));
	}
	
	/*CHANGED: SJDM June 8, 2012 fixed average weight-at-age for reference points
	           and average fecundity-at-age.
	*/
	
	//RF - DON'T USE THIS ONE - COMMENTED OUT IN GLOBALS SECTION
	#if defined(USE_NEW_EQUILIBRIUM)
		/* Newton-Raphson method to determine MSY-based reference points. */
		for(i=1;i<=15;i++)
		{
			equilibrium(fe,allocation,value(ro),value(kappa),value(m_bar),age,avg_wt,
					avg_fec,va,re,ye,be,ve,dye_df,d2ye_df2);
		
			fe = fe - dye_df/d2ye_df2;
			if(fabs(dye_df)<1e-6)break;
		}
		
		fmsy=fe;
		equilibrium(fe,allocation,value(ro),value(kappa),value(m_bar),age,avg_wt,
				avg_fec,va,re,ye,be,ve,dye_df,d2ye_df2);
	#endif
	
	#if !defined(USE_NEW_EQUILIBRIUM)
		for(i=1;i<=20;i++)
		{
			//equilibrium(fe,value(ro),value(kappa),value(m),age,wa,fa,value(exp(log_sel(1)(nyr))),re,ye,be,phiq,dphiq_df,dre_df);
			//equilibrium(fe,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),
			//			fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
			equilibrium(fe,value(ro),value(kappa),value(m_bar),age,avg_wt,
						avg_fec,va_bar,re,ye,be,phiq,dphiq_df,dre_df);
		
			dye_df = re*phiq+fe*phiq*dre_df+fe*re*dphiq_df;
			ddye_df = phiq*dre_df + re*dphiq_df;
			fe = fe - dye_df/ddye_df;
			if(verbose) cout<<"fe\t"<<fe<<"\t"<<dye_df<<"\t"<<ye<<endl;
			if(sfabs(dye_df)<1.e-5)break;
		}
		fmsy=fe;
		
		equilibrium(fmsy,value(ro),value(kappa),value(m_bar),age,avg_wt,
					avg_fec,va_bar,re,ye,be,phiq,dphiq_df,dre_df);
	#endif
	
	msy=ye;
	bmsy=be;
	Vmsy=be;
	Umsy=msy/Vmsy;
	
	/*TODO print this to the REPORT file for plotting.*/
	/*SM Loop over discrete value of fe and ensure above code is 
	finding the correct value of msy.*/
	
	if(!mceval_phase())
	{
		ofstream report_file("iscamdelaydiff.eql");
		if(report_file.is_open())
		{
			report_file<<"index\t fe \t ye \t be \t ve \t re \t spr\n";
		
			fe = 0; i=0;
			while(i < 1500)
			{
				#if !defined(USE_NEW_EQUILIBRIUM)
				equilibrium(fe,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),
							fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
				#endif
			
				#if defined(USE_NEW_EQUILIBRIUM)
				equilibrium(fe,allocation,value(ro),value(kappa),value(m_bar),age,avg_wt,
						avg_fec,va,re,ye,be,ve,dye_df,d2ye_df2);
				#endif
				if(re<=0)break;
			
				double spr = value(-ro/((kappa-1)*re-ro*kappa));
				report_file<<i++<<"\t"<<fe<<"\t"<<ye<<"\t"<<be<<"\t"<<ve<<"\t";
				report_file<<re<<"\t"<<spr<<endl;
				
				fe += 0.01;
			}
		}
	}//exit(1);
	
	if(verbose)cout<<"**** Ok after calc_reference_points ****"<<endl;
  }

FUNCTION void simulation_model(const long& seed)
  {
	/*
	Call this routine to simulate data for simulation testing.
	The random number seed can be used to repeat the same 
	sequence of random number for simulation testing.
	
	Implemented using the "-SimFlag 99" command line option where
	99 is the random number seed.
	
	-SimFlag 99 is a special case used for the manuscript for case 1.
	-SimFlag 000 is a special case with 0 error (exact data)
	
	-This routine will over-write the observations in memory
	with simulated data, where the true parameter values are
	the initial values.  Change the standard deviations of the 
	random number vectors epsilon (observation error) or 
	recruitment devs wt (process error).
	*/
	
	
	cout<<"___________________________________________________\n"<<endl;
	cout<<"  **Implementing Simulation--Estimation trial**    "<<endl;
	cout<<"___________________________________________________"<<endl;
	cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
	cout<<"___________________________________________________\n"<<endl;
	
	
	//Indexes:
	int i,j,k,ii,ki;

	//3darray Chat(1,ngear,syr,nyr,sage,nage);
	//C.initialize();
	
	
	
	/*----------------------------------*/
	/*	-- Generate random numbers --	*/
	/*----------------------------------*/
	random_number_generator rng(seed);
	dvector wt(syr-nage-1,nyr);			//recruitment anomalies
	dmatrix epsilon(1,nit,1,nit_nobs);  //observation errors in survey
	
	double sig = value(sqrt(rho)/varphi);
	double tau = value(sqrt(1.-rho)/varphi);
	
	if(seed==000)
	{
		cout<<"No Error\n";
		sig=0.;
		tau=0.;
	}
	wt.fill_randn(rng); wt *= tau;
	epsilon.fill_randn(rng); 
	
	//now loop over surveys and scale the observation errors
	for(k=1;k<=nit;k++)
	{
		for(j=1;j<=nit_nobs(k);j++)
			epsilon(k,j) *= sig/it_wt(k,j);
	}
	
	cout<<"	OK after random numbers\n";
	/*----------------------------------*/
	
	/*----------------------------------*/
    /*		--Initialize model--		*/
	/*CHANGED now calculating phie based on m_bar and avg_fec*/
	/*----------------------------------*/
	dvector lx=pow(exp(-value(m_bar)),age-min(age));
	lx(nage)/=(1.-exp(-value(m_bar)));
	double phie=(lx*exp(-value(m_bar)*cntrl(13)))*avg_fec;//fec(syr);
	so=kappa/phie;
	
	
	if(cntrl(2)==1) beta=(kappa-1.)/(ro*phie);
	if(cntrl(2)==2) beta=log(kappa)/(ro*phie);
	
	//Initial numbers-at-age with recruitment devs
	/*for(i=syr;i < syr+sage;i++)
			N(i,sage)=exp(log_avgrec+wt(i));
			
		for(j=sage+1;j<=nage;j++)
			N(syr,j)=exp(log_avgrec+wt(syr-j))*lx(j);
		*/
	
	N.initialize();
	if(cntrl(5)){	//If initializing in at unfished conditions
		log_rt(syr) = log(ro);
		for(j=sage;j<=nage;j++)
		{
			N(syr,j)=ro*exp(-m_bar*(j-1.));
		}
	}
	else{			//If starting at unfished conditions
		log_rt(syr) = log_avgrec;
		N(syr,sage)=mfexp(log_rt(syr));
		for(j=sage+1;j<=nage;j++)
		{
			N(syr,j)=mfexp(log_recinit+init_log_rec_devs(j))*exp(-m_bar*(j-sage));
		}
	}
	N(syr,nage)/=(1.-exp(-m_bar));
	
	//log_rt=log_avgrec+log_rec_devs;
	//log_rt(syr) = log(ro);

	for(i=syr+1;i<=nyr;i++){
		log_rt(i)=log_avgrec+log_rec_devs(i);
		N(i,sage)=mfexp(log_rt(i));
	}
	N(nyr+1,sage)=mfexp(log_avgrec);

	/*
	for(j=sage;j<=nage;j++) 
	{
		if(cntrl(5))  //if starting at unfished state
		{
			N(syr,j)=ro*exp(-m_bar*(j-1));
		}
		else{
			log_rt(syr-j+sage)=log_avgrec+log_rec_devs(syr-j+sage);
			N(syr,j)=mfexp(log_rt(syr-j+sage))*exp(-m_bar*(j-sage));
		}
	}

	N(syr,nage)/=(1.-exp(-m_bar));
	*/
	cout<<"	Ok after initialize model\n";
	/*----------------------------------*/
	
	
	
	/*----------------------------------*/
    /*		--    Selectivity   --		*/
	/*----------------------------------*/
	/*
		-Based on values in the control file.
		-Add autocorrelated random numbers
		for time varying or something to that
		effect.
		
		-If seed==99 then set up a special case
		for the cubic spline manunscript using
		the eplogistic function where g goes from
		strongly domed to asymptotic, e.g.,
		g = 0.2 * (nyr-i)/(nyr-syr);
		
	*/

	/*CHANGED May 15, 2011 calcSelectivities gets called from PRELIMINARY_CALCS*/
	dmatrix va(1,ngear,sage,nage);			//fishery selectivity
	d3_array dlog_sel(1,ngear,syr,nyr,sage,nage);
	dlog_sel=value(log_sel);
	/*
	for(k=1;k<=ngear;k++)
		for(i=syr;i<=nyr;i++)
		{
			//sel(k)(i)=plogis(age,ahat(k),ghat(k));
			log_sel(k)(i)=log(plogis(age,ahat(k),ghat(k)));
			log_sel(k)(i) -= log(mean(exp(log_sel(k)(i))));
		}
		//log_sel(j)(i) -= log(mean(mfexp(log_sel(j)(i))));
	*/
	cout<<"	Ok after selectivity\n";

	/*----------------------------------*/
	
	
	/*----------------------------------*/
    /*	--  Population dynamics  --		*/
	/*----------------------------------*/
	
	dmatrix zt(syr,nyr,sage,nage);			//total mortality
	zt.initialize();
	dmatrix ft(syr,nyr,1,ngear);
	ft.initialize();
	dvector sbt(syr,nyr+1);
	sbt.initialize();
	
	
	for(i=syr;i<=nyr;i++)
	{   
		
		//total biomass at age
		//dvector bt = elem_prod(value(N(i)),wa);
		dvector bt = elem_prod(value(N(i)),wt_obs(i));

		/*calculate instantaneous fishing mortalities
		based on Baranov's catch equation and the 
		observed catch from each fleet.*/
		dvector oct = trans(obs_ct)(i);
		
		
		for(k=1;k<=ngear;k++)
			va(k)=exp(dlog_sel(k)(i));
		
		//get_ft is defined in the Baranov.cxx file
		//CHANGED these ft are based on biomass at age, should be numbers at age
		//ft(i) = get_ft(oct,value(m),va,bt);
		ft(i) = get_ft(oct,value(m),va,value(N(i)),wt_obs(i));
		//cout<<trans(obs_ct)(i)<<"\t"<<oct<<endl;
		
		// overwrite observed catch incase it was modified by get_ft
		for(k=1;k<=ngear;k++)
			obs_ct(k,i)=oct(k);
		
		//total age-specific mortality
		//dvector zt(sage,nage);
		zt(i)=value(m);
		for(k=1;k<=ngear;k++){
			zt(i)+= ft(i,k)*exp(dlog_sel(k)(i));
		}
		
		
		//CHANGED definition of spawning biomass based on ctrl(12)
		sbt(i) = value(elem_prod(N(i),exp(-zt(i)*cntrl(13)))*fec(i));
		
		//Update numbers at age
		if(i>=syr+sage-1)
		{
			double rt;
			//double et=value(N(i-sage+1))*fec(i-sage+1);
			double et=sbt(i-sage+1);
			if(cntrl(2)==1)rt=value(so*et/(1.+beta*et));
			if(cntrl(2)==2)rt=value(so*et*exp(-beta*et));
			N(i+1,sage)=rt*exp(wt(i)-0.5*tau*tau);
			
			/*CHANGED The recruitment calculation above is incosistent
			  with the assessment model.  Below recruitment is based on
			  rt=exp(log_avgrec + wt + rt_dev), where the rt_dev calculation
			is based on the BH or Ricker model.*/
			//double rt_dev = log(rt)-value(log_avgrec);
			//N(i+1,sage)=exp(log_avgrec+wt(i));
			
		}


		N(i+1)(sage+1,nage)=++elem_prod(N(i)(sage,nage-1),exp(-zt(i)(sage,nage-1)));
		N(i+1,nage)+=N(i,nage)*exp(-zt(i,nage));
		
		
		//Catch & Catch-at-age
		for(k=1;k<=ngear;k++)
		{
			if(ft(i,k)>0)
			{
				dvector sel = exp(dlog_sel(k)(i));
				d3C(k)(i)=elem_prod(elem_div(ft(i,k)*sel,zt(i)),elem_prod(1.-exp(-zt(i)),value(N(i))));
				obs_ct(k,i)=d3C(k)(i)*wt_obs(i);
			}
			else	//if this is a survey
			{
				dvector sel = exp(dlog_sel(k)(i));
				d3C(k)(i)=elem_prod(elem_div(sel,zt(i)),elem_prod(1.-exp(-zt(i)),value(N(i))));
			}
		}
		
	}
	
	
	//initial values of log_ft_pars set to true values
	ki=1;
	for(k=1;k<=ngear;k++)
		for(i=syr;i<=nyr;i++)
			if(obs_ct(k,i)>0){
				log_ft_pars(ki++)=log(ft(i,k));
			}
	
	// Error handler to inform user population went extinct.
	if(min(sbt(syr,nyr))<=1.e-5)
	{
		cout<<"---------------------------------------\n";
		cout<<"Simulated population went extinct, try\n";
		cout<<"increasing steepness, Ro and Rbar\n";
		cout<<sbt<<endl;
		cout<<"Minimum spawning biomass="<<min(sbt(syr,nyr))<<endl;
		cout<<"---------------------------------------\n";
		exit(1);
	}
	
	//Average recruitment calculation
	cout<<"	log(mean(column(N,sage))) = "<<mean(log(column(N,sage)))<<endl;
	cout<<"	log_avgrec = "<<log_avgrec<<endl;
	cout<<"	Ok after population dynamics\n";
	/*----------------------------------*/
	
	/*----------------------------------*/
    /*	--  Observation models  --		*/
	/*----------------------------------*/
	
	//Simulated Age-compositions
	int ig;
	for(k=1;k<=na_gears;k++)
	{
		for(i=1;i<=na_nobs(k);i++)
		{
			ii=A(k,i,a_sage(k)-2);	//index for year
			ig=A(k,i,a_sage(k)-1);	//index for gear
			dvector pa = d3C(ig)(ii);	//
			pa/=sum(pa);
			
			
			dvector t1=pa(a_sage(k),a_nage(k));
			t1/=sum(t1);
			A(k)(i)(a_sage(k),a_nage(k))=rmvlogistic(t1,0.3,i+seed);
			if(seed==000)
			{
				A(k)(i)(a_sage(k),a_nage(k))=t1;
			}
			//cout<<iyr<<"\t"<<k<<endl;
		}
	}

	//cout<<Ahat<<endl;
	//Relative abundance indices
	//CHANGED fixed this to reflect survey timing etc & survey_type
	for(k=1;k<=nit;k++)
	{   
		for(i=1;i<=nit_nobs(k);i++)
		{
			ii=iyr(k,i);
			ig=igr(k,i);
			dvector sel = exp(dlog_sel(ig)(ii));
			dvector Np = value(elem_prod(N(ii),exp(-zt(ii)*it_timing(k,i))));
			switch(survey_type(k))
			{
				case 1: //survey based on numbers
					Np = elem_prod(Np,sel);
				break;
				case 2: //survey based on biomass
					Np = elem_prod(elem_prod(Np,sel),wt_obs(ii));
				break;
				case 3: //survey based on spawning biomass
					Np = elem_prod(Np,fec(ii));
				break;
			}
			it(k,i) = sum(Np) * exp(epsilon(k,i));
		}
	}
	cout<<"	OK after observation models\n";
	/*----------------------------------*/
	
	//CHANGED Fixed bug in reference points calc call from simulation model,
	//had to calculate m_bar before running this routine.
	
	calc_reference_points();
	//cout<<"	OK after reference points\n"<<fmsy<<endl;
	//exit(1);
	//	REPORT(fmsy);
	//	REPORT(msy);
	//	REPORT(bmsy);
	
	
	cout<<"___________________________________________________"<<endl;
	ofstream ofs("iscam.sim");
	ofs<<"fmsy\n"<<fmsy<<endl;
	ofs<<"msy\n"<<msy<<endl;
	ofs<<"bmsy\n"<<bmsy<<endl;
	ofs<<"va\n"<<va<<endl;
	ofs<<"sbt\n"<<sbt<<endl;//<<rowsum(elem_prod(N,fec))<<endl;
	ofs<<"rt\n"<<rt<<endl;
	ofs<<"ct\n"<<obs_ct<<endl;
	ofs<<"ft\n"<<trans(ft)<<endl;
	//ofs<<"ut\n"<<elem_div(colsum(obs_ct),N.sub(syr,nyr)*wa)<<endl;		//RF THIS ISN'T TRUE UT - NOT A FUNCTION OF SELECTIVITY
	//ofs<<"ut\n"<<1.-mfexp(-ft(1))<<endl;	   // RF added this but hasn't tested it - haven't fugured out whether it should be transposed or is needed
	ofs<<"iyr\n"<<iyr<<endl;
	ofs<<"it\n"<<it<<endl;
	ofs<<"N\n"<<N<<endl;
	ofs<<"A\n"<<A<<endl;
	ofs<<"dlog_sel\n"<<dlog_sel<<endl;
	cout<<"  -- Simuation results written to iscam.sim --\n";
	cout<<"___________________________________________________"<<endl;
	
	//cout<<N<<endl;
	//exit(1);
  }
	
FUNCTION dvector cis(const dvector& na)
  {
	//Cohort Influenced Selectivity
	//This function returns a vector of residuals from a
	//linear regression of log(pa)= a+b*age+res that can be
	//used to modify age-based selectivity according to relative
	//cohort strengths.
	
	//SM  Currently not used at all in iscam and should be deprecated.
	dvector y = log(na);
	dvector x = age;
	double b = sum(elem_prod(x-mean(x),y-mean(y)))/sum(square(x-mean(x)));
	double a = mean(y)-b*mean(x);
	dvector res = y - (a+b*x);
	return(res);
  }

REPORT_SECTION
  {
	//RF added CCAM's report section here so can use CCAM GUI
	 if(verbose)
  	cout<<"Start of Report Section..."<<endl;
  	//cout<<"Female sbo report =  "<<bo/2<<endl<<endl;
  	//CG RF change - a large number of changes to report section
  	int i,j,k,g;
  
   	report<<"npar\n"<<npar<<endl;
  	for(i=1;i<=npar;i++){
  		report<<"priorPars_theta"<<i<<endl<<theta_control(i)<<endl;
  	}
  	
  	REPORT(q_prior);
  	REPORT(q_mu);
  	REPORT(q_sd);
  
  	report<<"nf\n"<<nf<<endl;
  	report<<"runtime\n"<<(long(difftime(finish,start))%3600)%60<<endl;
  	report<<"ObjectiveFunction\n"<<f<<endl;
  
    	REPORT(ControlFile);
  	REPORT(f);
  	REPORT(nlvec);
  	REPORT(delaydiff) 	
  	if(delaydiff) report<<"#  Running delay difference model"<<endl;
  	if(!delaydiff) report<<"#  Running age-structured model"<<endl;
    
  	report<<"\n#______________________________\n"<<endl;
  	report<<"#______________________________"<<endl;
    	report<<"#  Estimated parameters"<<endl;
  	REPORT(theta);  // all 9 parameter estimates.
  	REPORT(ro);
  	double steepness=value(theta(2));
	REPORT(steepness);
  	REPORT(m);
  	REPORT(log_avgrec);
  	double rbar=value(exp(log_avgrec));
  	REPORT(rbar);
  	REPORT(log_recinit);
  	double rinit=value(exp(log_recinit));
  	REPORT(rinit);
  	REPORT(tau);
  	REPORT(sig);
    	REPORT(rho);
  	REPORT(varphi);
	REPORT(age_tau2);
  	//@@@@@@@@@Temporary parameters
  	REPORT(qc);
  	REPORT(dt); //RF: this is new parameter to scale effect of observed recruitment anomalies
  	
  	REPORT(q);
  	REPORT(qt);
  	REPORT(so); //SR alpha parameter
  	//REPORT(beta);//SR beta parameter
  	report<<"bta\n"<<beta<<endl; 
  	report<<"sbo\n"<<bo<<endl; //female spawning biomass
  	REPORT(no);
  	REPORT(kappa); //recruitment compensation ratio
  	//Delay difference model outputs
	if(delaydiff){
		REPORT(wbar);
		REPORT(alpha_g);
		REPORT(rho_g);
		REPORT(wk);
	}

  	if(!delaydiff){
		//Selectivity
		report<<"#______________________________"<<endl;
		report<<"#  Selectivity parameters"<<endl;

		report<<"ahat"<<endl;
		for(k=1;k<=ngear;k++) report<<mfexp(sel_par(k,1,1))<<" ";
		report<<endl<<"ghat"<<endl;
		for(k=1;k<=ngear;k++) report<<mfexp(sel_par(k,1,2))<<" ";
		report<<endl;

		//report<<"log_sel"<<endl;
		//for(k=1;k<=ngear;k++)
		//	for(i=syr;i<=nyr;i++)
		//		report<<k<<"\t"<<log_sel(k)(i)<<endl;
		report<<"selectivity"<<endl;
		for(g=1;g<=ngear;g++){
			for(i=syr;i<=nyr;i++) report<<g<<"\t"<<mfexp(log_sel(g)(i))<<endl;
		}

	}
  	ivector yr(syr,nyr);
  	ivector yrs(syr,nyr+1);
  	//ivector wyrs(annual_mean_wt_syr, annual_mean_wt_nyr);
  	yr.fill_seqadd(syr,1); 
  	yrs.fill_seqadd(syr,1); 
  	//wyrs.fill_seqadd(annual_mean_wt_syr,1); 

  	
  	REPORT(ngear);
  	REPORT(nits);
  	REPORT(yr);
  	REPORT(yrs);
  	REPORT(wyrs);
  	REPORT(iyr); //WHEN MULTIPLE SURVEYS - MAKE SURE THAT CODE IS ROBUST TO FIRST SURVEY BEING SHORTER THAN SUBSEQUENT SURVEYS
  	REPORT(age);
  	REPORT(sage);
  	REPORT(nage);
  	REPORT(la);
  	if(!delaydiff) REPORT(wa);
  	if(!delaydiff) REPORT(fec);
  	REPORT(obs_ct);
  	REPORT(ct);
  	REPORT(ft);
  	REPORT(effort);
  	REPORT(eta);
  	/*FIXED small problem here with array bounds if using -retro option*/
  	//if(!delaydiff) report<<"ut\n"<<elem_div(colsum(ct)(syr,nyr),N.sub(syr,nyr)*wa)<<endl; //RF Need to use predicted catch because model does not always fit all catch points
  	//if(delaydiff) report<<"ut\n"<<elem_div(colsum(ct)(syr,nyr),biomass(syr,nyr))<<endl;
  	//report<<"ut\n"<<elem_div(ct(1)(syr,nyr),vbcom(syr,nyr))<<endl; //RF Need to use predicted catch because model does not always fit all catch points
  	report<<"ut\n"<<1.-mfexp(-ft(1))<<endl;
  	int rectype=int(cntrl(2));
  	REPORT(rectype);
  	REPORT(rt);
	dvector ln_rt=value(log_rt(syr,nyr));
  	REPORT(ln_rt);
  	REPORT(delta);
  	REPORT(log_rec_devs);
  	REPORT(it); //WHEN MULTIPLE SURVEYS - MAKE SURE THAT CODE IS ROBUST TO FIRST SURVEY BEING SHORTER THAN SUBSEQUENT SURVEYS
  	REPORT(pit); //WHEN MULTIPLE SURVEYS - MAKE SURE THAT CODE IS ROBUST TO FIRST SURVEY BEING SHORTER THAN SUBSEQUENT SURVEYS
  	REPORT(epsilon);
	if(!delaydiff){
		REPORT(F);
		REPORT(M_tot);
		REPORT(na_gears);
		REPORT(na_nobs);
		REPORT(a_sage);
		REPORT(a_nage);
		REPORT(A); 
		REPORT(Ahat);
		REPORT(A_nu); //array of residuals
		REPORT(N);
		REPORT(wt_obs);
	}	
  	
  	//RF Added parameters
  	REPORT(weight_sig);
  	REPORT(annual_mean_wt);
	REPORT(obs_annual_mean_wt);
	REPORT(vbcom);
	REPORT(vncom);
	REPORT(numbers.sub(syr,nyr+1)); //RF in ASM this is the same as vulnerable numbers - 
	
	if(last_phase())
	{	
		if(!delaydiff){
			calc_reference_points();
			REPORT(fmsy);
			REPORT(msy);
			REPORT(bmsy);
		}
		
		if(delaydiff){
			//run_FRPdd();
			findmsy_dd(msy, fmsy, bmsy);
			REPORT(fmsy);
			REPORT(msy);
			REPORT(bmsy);
		}
		//run_FRP(); //USE THIS ONLY FOR TESTING REF POINTS - RF HAS DONE THIS FOR PCOD
	}
  	
  	//biomasses
  	cout<<"Report Section CHECK ***********************************..."<<endl;
  	cout<<"N min row index: "<<N.rowmin()<<endl;
  	cout<<"N max row index: "<<N.rowmax()<<endl;
  	cout<<"N min col index: "<<N.colmin()<<endl;
  	cout<<"N max col index: "<<N.colmax()<<endl;
  	cout<<"wt_obs min row index: "<<wt_obs.rowmin()<<endl;
  	cout<<"wt_obs max row index: "<<wt_obs.rowmax()<<endl;
  	cout<<"wt_obs min col index: "<<wt_obs.colmin()<<endl;
  	cout<<"wt_obs max col index: "<<wt_obs.colmax()<<endl<<endl;
  	dmatrix nt3=value(trans(N).sub(3,nage));
  	dmatrix wt3=value(trans(wt_obs).sub(3,nage));
  	dvector bt3=colsum(elem_prod(nt3,wt3));
  	
  	report<<"sbt\n"<<sbt<<endl; //SPAWNING BIOMASS
  	report<<"biomass\n"<<biomass.sub(syr, nyr)<<endl; //SPAWNING BIOMASS//from dd model in ASM this is the same as sbt	   //RF don't report projection year based on average recruits
  	report<<"bt3\n"<<bt3<<endl;
  	report<<"tbt\n"<<tbt.sub(syr, nyr)<<endl; //Total BIOMASS
  	dmatrix nt2=value(trans(N).sub(2,nage));
  	dmatrix wt2=value(trans(wt_obs).sub(2,nage));
  	dvector bt2=colsum(elem_prod(nt2,wt2));
  	report<<"bt2\n"<<bt2<<endl;
  	//dvector rt=value(column(N.sub(syr+1,nyr+1),sage)/exp(-m));
  	//report<<"rt_sageminus1\n"<<rt<<endl;
  	dvector rt3(1,3);
  	
  	dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),wt_obs(nyr+1)));
  	REPORT(future_bt);
  	  
  	//Parameter controls
  	dmatrix ctrl=theta_control;
  	REPORT(ctrl);
  	REPORT(cntrl);
    	
    	if(last_phase()) {
    		for(i=1;i<=n_tac;i++)
		{
			if (!delaydiff) projection_model(tac(i));
			if(delaydiff) projection_model_dd(tac(i));
		}
	}
  	
  	if(verbose)cout<<"END of Report Section..."<<endl;
  	
  	/*IN the following, I'm renaming the report file
	in the case where retrospective analysis is occurring*/
	if(retro_yrs && last_phase() && PLATFORM =="Linux")
	{
		//adstring rep="iscam.ret"+str(retro_yrs);
		//rename("iscam.rep",rep);
		adstring copyrep = "cp iscam.rep iscam.ret"+str(retro_yrs);
		system(copyrep);
	}
    }

	
FUNCTION decision_table
  {
	/*
	This function takes a vector of projected catches and computes the following
	Reference points: Bmsy, Bo, Fmsy, Umsy.
	
	Biomass Metrics for the decision table:
	1) P(SB_{t+1} < SB_{t})
	2) P(SB_{t+1} < 0.25 B_{0})
	3) P(SB_{t+1} < 0.75 B_{0})
	4) P(SB_{t+1} < 0.40 B_{MSY})
	5) P(SB_{t+1} < 0.80 B_{MSY})
	
	Harvest Metrics for the decision table:
	1) P(U_{t+1} > Target harvest rate)
	2) P(U_{t+1} > 1/2 Fmsy)
	3) P(U_{t+1} > 2/3 Fmsy)
	4) P(tac/2+  > 20%)
	
	Key to the harvest metric is the definition of Umsy and allocation to fleets.
	
	Pseudocode:
		1) Calculate reference points (Fmsy, Bmsy)
		2) Loop over vector of proposed catches
		3) Evaluate biomass metrics for each posterior sample
		4) Evaluate harvest metrics for each posterior sample
	
	*/
	int i;
	// 1) Calculate reference pionts.
	if (!delaydiff) calc_reference_points();
	//if (delaydiff) run_FRPdd();
	if (delaydiff) findmsy_dd(msy, fmsy, bmsy);
	
	// 2) Loop over vector of proposed catches
	//    This vector should is now read in from the projection file control (pfc).
	for(i=1;i<=n_tac;i++)
	{
		if (!delaydiff) projection_model(tac(i));
		if(delaydiff)   projection_model_dd(tac(i));
	}
  }
	
FUNCTION mcmc_output
  {
	if(nf==1){
		adstring str_q;
		str_q="q";
		ofstream ofs("iscamdelaydiff.mcmc");
		ofs<<"log.ro\t h\t log.m\t log.rbar\t log.rinit\t rho\t varphi\t qc\t dt\t";
		ofs<<"bo\t bmsy\t msy\t fmsy\t";
		ofs<<"SSB\t Age-4\t Poor\t Average\t Good\t";
		for(int i=1;i<=nit;i++)ofs<<str_q<<i<<"\t";
		ofs<<"f\t"<<endl;
		
		ofstream of1("sbt.mcmc");
		ofstream of2("rt.mcmc");
		//ofstream of3("1-spr.mcmc");
		ofstream of4("ft.mcmc");
		if(!delaydiff) ofstream of5("bt2.mcmc");
		ofstream of6("tbt.mcmc");
		ofstream of7("surveyres.mcmc");
		//ofstream of8("sprmsy_status.mcmc");
		//ofstream of9("spr40_status.mcmc");
		ofstream of10("sbtdepletion.mcmc");
		ofstream of12("recDevs.mcmc");
		ofstream of13("surveyres2.mcmc"); //only use if number of surveys (nits) > 1
		ofstream of14("surveyres3.mcmc"); //only use if number of surveys (nits) > 2
		ofstream of15("surveyres4.mcmc"); //only use if number of surveys (nits) == 4 
		
	}
	
	// leading parameters & reference points
	if (!delaydiff) calc_reference_points();
	//if (delaydiff) run_FRPdd();
	if (delaydiff) findmsy_dd(msy, fmsy, bmsy);
	
	// decision table output
	
		dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),wt_obs(nyr+1)));
		double future_bt4 = sum(future_bt(4,nage));
		dvector rt3 = age3_recruitment(value(column(N,3)),wt_obs(nyr+1,3),value(M_tot(nyr,3)));	
		//dvector bt = value(rowsum(elem_prod(N,wt_obs)));
		dmatrix nt2=value(trans(N).sub(2,nage));
		dmatrix wt2=value(trans(wt_obs).sub(2,nage));
		dvector bt2=colsum(elem_prod(nt2,wt2));
	
	
	ofstream ofs("iscamdelaydiff.mcmc",ios::app);
	ofs<<theta;
	ofs<<" "<<bo<<" "<<bmsy<<" "<<msy<<" "<<fmsy<<"\t\t";
	ofs<<sbt(nyr)<<" "<<future_bt4<<" "<<future_bt4+rt3<<"\t\t";
	ofs<<q<<" "<<f<<endl;
	
	// output spawning stock biomass
	ofstream of1("sbt.mcmc",ios::app);
	of1<<sbt<<endl;
	
	// output age-1 recruits
	ofstream of2("rt.mcmc",ios::app);
	of2<<rt<<endl;
	
	//ofstream of3("1-spr.mcmc",ios::app);
	//of3<<1-ispr<<endl;
	
	// output ft
	ofstream of4("ft.mcmc",ios::app);
	of4<<ft(1)<<endl;
	// output age-2 biomass
	if(!delaydiff){
		ofstream of5("bt2.mcmc",ios::app);
		of5<<bt2<<endl;
	}	
	// output total biomass
	ofstream of6("tbt.mcmc",ios::app);
	of6<<tbt<<endl;
	// output survey residuals
	ofstream of7("surveyres.mcmc",ios::app);
	of7<<epsilon(1)<<endl;
	
	// output spr status (fmsy policy)
	//ofstream of8("sprmsy_status.mcmc",ios::app);
	//of8<<(1.-ispr)/(1.-spr_msy)<<endl;
	// output spr status (f40 policy)
	//ofstream of9("spr40_status.mcmc",ios::app);
	//of9<<(1.-ispr)/(1.-0.4)<<endl;
	
	ofstream of10("sbtdepletion.mcmc",ios::app);
	of10<<sbt(syr,nyr+1)/bo<<endl;
	
	// output recruitment deviations 
	ofstream of12("recDevs.mcmc", ios::app);
	of12<<mfexp(log_rec_devs)<<endl;
	//Write extra files for survey residuals if more than  one survey
	//(can't have more than 4 right now - make dynamic)
	if(nits>1) {
		ofstream of13("surveyres2.mcmc", ios::app);
		of13<<epsilon(2)<<endl;
	}
	if(nits>2) {
		ofstream of14("surveyres3.mcmc", ios::app);
		of14<<epsilon(3)<<endl;
	} 
	if(nits==4) {
		ofstream of15("surveyres4.mcmc", ios::app);
		of15<<epsilon(4)<<endl;
	}
	
	/* June 12, 2012.  SJDM Call decision table. */
	decision_table();  
	
	
	// DEPRECATED	
	// Projection model.

	
	//for(int i=0;i<=10;i++)
	//{
	//	double tac = double(i)/10. * 1.5*msy;
	//	projection_model(tac);
	//}
	
	// Deviance Information Criterion
	/*
		DIC = pd + Dbar, where:
		pd is the effective number of parameters,
		Dbar is the expected (average deviance)
		
		Dbar = mean(D(theta))
		pd = Dbar - D(mean(theta))
		Agorithm:
			-for each sample compute Dtotal += 2*f;
			-compute a running Dbar = Dtotal/nf;
			
			
	
	//cout<<initial_params::nvarcalc()<<endl;
	int nvar = initial_params::nvarcalc();
	dvector y(1,nvar);
	static dvector sum_y(1,nvar);
	initial_params::xinit(y);
	sum_y = sum_y + y;
	cout<<y(1,3)<<endl;
	cout<<sum_y(1,3)<<endl;
	double fi = get_monte_carlo_value(nvar,y);
	if(nf==2)exit(1);
	*/
	
  }

FUNCTION dvector age3_recruitment(const dvector& rt, const double& wt,const double& M)
  {
	/*
	This routine returns the poor average and good age-3 recruits
	that is used in constructing the decision table for the pacific
	herring fisheries.
	
	-1) sort the rt vector from small to large
	-2) find the 33rd and 66th percentiles
	-3) take the averge of the 0-33, 34-66, 67-100
	-4) multiply by the average weight
	-5) return the age-3 recruitment biomass
	*/
	
	dvector s_rt = sort(rt);
	dvector rbar(1,3);
	
	double idx = floor((nyr-syr+1.)/3.);
	int ix1 = syr+int(idx);
	int ix2 = syr+int(2.*idx);
	rbar(1) = mean(s_rt(syr,ix1));
	rbar(2) = mean(s_rt(ix1+1,ix2));
	rbar(3) = mean(s_rt(ix2+1,nyr));
	rbar = rbar*wt*exp(-M);
	//cout<<rbar<<endl;
	return(rbar);
  }

FUNCTION void projection_model(const double& tac)
  {
	/*
	This routine conducts population projections based on 
	the estimated values of theta.  Note that all variables
	in this routine are data type variables.
	
	Arguments:
	tac is the total allowable catch that must be allocated 
	to each gear type based on allocation(k)
	
	theta(1) = log_ro
	theta(2) = h
	theta(3) = log_m
	theta(4) = log_avgrec
	theta(5) = log_recinit
	theta(6) = rho
	theta(7) = vartheta
	
	** NOTES **
	* Projections are based on average natural mortality and fecundity.
	* Selectivity is based on selectivity in terminal year.
	* Average weight-at-age is based on mean weight in the last 5 years.
	
	Aug 20, 2012 Found a bug, See issue 2 on github.
	*/
	static int runNo=0;
	runNo ++;
	int i,j,k;
	int pyr = nyr+1;	//projection year. 
	
	// --derive stock recruitment parameters
	// --survivorship of spawning biomass
	dvector lx(sage,nage);
	double   tau = value(sqrt(1.-rho)*varphi); 
	double   m_M = value(m_bar); 
	double m_rho = cntrl(13);
	lx(sage)     = 1;
	for(i=sage; i<=nage; i++)
	{
		lx(i) = exp( -m_M*(i-sage) -m_rho*m_M );
		if(i==nage) 
			lx(i) /= 1.0 - exp( -m_M );
	}	
	double phib = lx*avg_fec;
	double so   = value(kappa)/phib;
	double bo   = value(ro)*phib;
	
	double beta;
	switch(int(cntrl(2)))
	{
		case 1:  // Beverton-Holt
			beta = (value(kappa)-1.)/bo;
		break;
		case 2:  // Ricker
			beta = log(value(kappa))/bo;
		break;
	}
	
	/* Fill arrays with historical values */
	dvector p_sbt(syr,pyr);
	dvector  p_ct(1,ngear);
	dmatrix  p_ft(nyr+1,pyr,1,ngear);
	dmatrix   p_N(syr,pyr+1,sage,nage);
	dmatrix   p_Z(syr,pyr,sage,nage);
	p_N.initialize();
	p_sbt.initialize();
	p_Z.initialize();
	p_N.sub(syr,nyr)   = value( N.sub(syr,nyr) );
	p_sbt(syr,nyr)     = value( sbt(syr,nyr)   );
	p_Z.sub(syr,nyr)   = value( Z.sub(syr,nyr) );
	
	
	/* Selectivity and allocation to gears */
	dmatrix va_bar(1,ngear,sage,nage);
	for(k=1;k<=ngear;k++)
	{
		p_ct(k)   = allocation(k)*tac;
		va_bar(k) = exp(value(log_sel(k)(nyr)));
	}
		
	/* Simulate population into the future under constant tac policy. */
	for(i = nyr; i<=pyr; i++)
	{
		
		if(i > nyr)
		{
			// get_ft is defined in the Baranov.cxx file
			p_ft(i) = getFishingMortality(p_ct, value(m_bar), va_bar, p_N(i),avg_wt);
			// p_ft(i)(1,nfleet) = fmsy(1,nfleet);
			
			// calculate total mortality in future years
			p_Z(i) = value(m_bar);
			for(k=1;k<=ngear;k++)
			{
				p_Z(i)+=p_ft(i,k)*va_bar(k);
			}
		}
		
		
		// spawning biomass
		p_sbt(i) = elem_prod(p_N(i),exp(-p_Z(i)*cntrl(13))) * avg_fec;
		
		
		// sage recruits with random deviate xx
		// note the random number seed is repeated for each tac level.
		double  xx = randn(nf+i)*tau;
		if(i>=syr+sage-1)
		{
			double rt;
			double et = p_sbt(i-sage+1);		// lagged spawning biomass
			
			if(cntrl(2)==1)						// Beverton-Holt model
			{
				rt=(so*et/(1.+beta*et));
			}
			if(cntrl(2)==2)						// Ricker model
			{
				rt=(so*et*exp(-beta*et));
			}
			
			p_N(i+1,sage)=rt*exp(xx-0.5*tau*tau); 
		}
		
		/* Update numbers at age in future years */
		p_N(i+1)(sage+1,nage) =++ elem_prod(p_N(i)(sage,nage-1),exp(-p_Z(i)(sage,nage-1)));
		p_N(i+1,nage)        +=   p_N(i,nage)*exp(-p_Z(i,nage));
		
		//Predicted catch for checking calculations
		//for(k=1;k<=nfleet;k++)
		//{
		//	dvector ba = elem_prod(p_N(i),avg_wt);
		//	cout<<k<<" tac = "<<tac<<"\t ct = ";
		//	cout<<sum(elem_div(elem_prod(elem_prod(ba,p_ft(i,k)*va_bar(k)),1.-exp(-p_Z(i))),p_Z(i)));
		//	cout<<" fmsy = "<<fmsy<<" ft = "<<p_ft(i,k)<<endl;
		//}
	}	
	//cout<<"fmsy\n"<<fmsy<<endl;
	//cout<<"Spawning biomass\n"<<p_sbt<<endl;
	//exit(1);
	/* 
	  Write output to *.proj file for constructing decision tables. 
	
	  Biomass Metrics for the decision table:
	  1) P(SB_{t+1} < SB_{t})
	  2) P(SB_{t+1} < 0.25 B_{0})
	  3) P(SB_{t+1} < 0.75 B_{0})
	  4) P(SB_{t+1} < 0.40 B_{MSY})
	  5) P(SB_{t+1} < 0.80 B_{MSY})
	  
	  Harvest Metrics for the decision table:
	  1) P(U_{t+1} > Umsy)
	  2) P(U_{t+1} > 1/2 Umsy)
	  3) P(U_{t+1} > 2/3 Umsy)
	  4) P(tac/2+  > 20%)
	
	  Metric for spawning depletion and spawning trends:
	  1) P(5-year decline)
	  
	  Defn: Stock status is based on spawning biomass
	  Defn: Removal rate is based on removals/spawning biomass
	  Defn: Harvest rate    : U_{t+1}=TAC/SBio
	  Defn: MSY harvest rate: Umsy=MSY/SBmsy
	  
	
	*/
	if(nf==1 && runNo==1)
	{
		ofstream ofs(BaseFileName + ".proj");
		ofs<<" tac";
		ofs<<"      P(SB1)"; 
		ofs<<"      P(SB2)";
		ofs<<"      P(SB3)";
		ofs<<"      P(SB4)";
		ofs<<"      P(SB5)";
		ofs<<"       P(U1)";
		ofs<<"       P(U2)";
		ofs<<"       P(U3)";
		ofs<<"       P(U4)";
		ofs<<"       P(D5)";
		ofs<<endl;
		cout<<"Bo when nf==1 \t"<<bo<<endl;
	}
	
	double  ut  = tac / p_sbt(pyr);
	double u20  = tac / ( (p_N(pyr)(3,nage)*exp(-value(M_tot(nyr,3))))* avg_wt(3,nage) );
	
	/* Average rate of change in spawning biomass in last 5 years */
	double dSb5 = mean(log(p_sbt(pyr-5,pyr)) - log(p_sbt(pyr-6,pyr-1).shift(pyr-5)));
	
	ofstream ofs(BaseFileName + ".proj",ios::app);
	ofs<< setprecision(4)               <<setw(4) 
	   << tac                           <<setw(12)
	   << p_sbt(pyr-1)/p_sbt(pyr)       <<setw(12)
	   << 0.25*bo/p_sbt(pyr)            <<setw(12)
	   << 0.75*bo/p_sbt(pyr)            <<setw(12)
	   << 0.40*bmsy/p_sbt(pyr)          <<setw(12)
	   << 0.80*bmsy/p_sbt(pyr)          <<setw(12)
	   << ut/Umsy                       <<setw(12)
	   << ut/(0.5*Umsy)                 <<setw(12)
	   << ut/(2./3.*Umsy)               <<setw(12)
	   << u20/0.2                       <<setw(12)
	   << dSb5+1                        <<setw(12)
	   << endl;
	// cout<<"Finished projection model"<<endl;
  }

FUNCTION void findmsy_dd(double& msy, double& fmsy, double& bmsy)
 {
	/*
	Brute force method to determine msy, fmsy & bmsy (added by PJS and RH, Aug 2015)
	*/
	static int runNo=0;
	runNo ++;
	int i,j,k;
	j =0;
	int nproj = 200;
	int s = 1;
	int pyr = nyr+nproj; //projection years

	//get parameters - convert to data objects
	double pbo   = value(bo);
	double pso   = value (so);
	double pbeta = value(beta);

	//values for searching fishing mortality rates
	double stft	 = 0.01;
	double endft = 1.0; //0.4;
	double stepft= 0.01;
	int n = (endft - stft)/stepft;
	double ift;
	dvector ye(s,n);
	dvector be(s,n);
	dvector fe(s,n);
		  
	dvector p_bt(syr,pyr);
	dvector p_ft(syr,pyr);
	dvector p_N(syr,pyr);
	dvector p_S(syr,pyr);
	dvector p_rt(syr+sage,pyr);//dmatrix p_Ft(1,ngear,syr,pyr);

	p_bt.initialize();
	p_ft.initialize();
	p_N.initialize();
	p_S.initialize();
	p_rt.initialize();
	ye.initialize();
	be.initialize();
	fe.initialize();
	
	p_ft(syr,nyr)  = value(ft(1)(syr,nyr));
	p_N(syr,nyr)   = value(numbers(syr,nyr));
	p_bt(syr,nyr)  = value(biomass(syr,nyr)); //sbt and vul biomass all the same for delay diff
	p_S(syr,nyr)   = value(surv(syr,nyr));
	p_rt(syr+sage,nyr) = value(rt(syr+sage,nyr));
		
	/* Simulate population into the future under constant F policies */
	
	for(ift=stft;ift<=endft;ift=ift+stepft) {
		j=j+1;
		for(i = nyr+1; i<=pyr; i++)
		{
			//recruits
			double p_tau = value(tau); 
			//double xx = randn(nf+i)*p_tau;
			double xx = p_tau;
				
			double et=p_bt(i-kage); //delay diff
			if(cntrl(2)==1)p_rt(i)=value((so*et/(1.+beta*et))*exp(xx-0.5*p_tau*p_tau));
			if(cntrl(2)==2)p_rt(i)=value((so*et*exp(-beta*et))*exp(xx-0.5*p_tau*p_tau));

			//Update biomass and numbers	
			p_bt(i) = (p_S(i-1)*(rho_g*p_bt(i-1)+alpha_g*p_N(i-1))+wk*p_rt(i));
			p_N(i)  = p_S(i-1)*p_N(i-1)+p_rt(i);
			//cout<<i<<"  "<<p_rt(i)<<"  "<<p_bt(i)<<"  "<<p_N(i)<<endl;
			//cout<<i<<"  "<<p_bt<<endl;
		
			//Calculate mortality for next projection year
			p_S(i) = mfexp(-(value(m)+ift));
		}
		be(j) = p_bt(pyr);
		ye(j) = (1-mfexp(-ift))*p_bt(pyr);
		fe(j) = ift;
	}

	//get MSY and Fmsy
	msy=max(ye);
	double mtest;
	for(k=1; k<=n; k++)
	{
		mtest=ye(k);
		if(mtest==msy) fmsy=fe(k);
		if(mtest==msy) bmsy=be(k);
		
	}
	
	//cout<<"Slower"<<endl;
	//cout<<msy<<endl;
	//cout<<fmsy<<endl;
	//cout<<bmsy<<endl; 
 }

FUNCTION void projection_model_dd(const double& tac)
 {
	/*
	This routine conducts population projections based on 
	the estimated values of theta.  Note that all variables
	in this routine are data type variables.
	
	Arguments:
	tac is the total allowable catch that must be allocated 
	to each gear type based on allocation(k)
	
	theta(1) = log_ro
	theta(2) = h
	theta(3) = log_m
	theta(4) = log_avgrec
	theta(5) = log_recinit
	theta(6) = rho
	theta(7) = vartheta
	
	** NOTES **
	* Projections are based on estimated constant natural mortality 
	
	*/
	static int runNo=0;
	runNo ++;
	int i,j,k;
	//int pyr = nyr+2;	//projection year. THE YEAR OF INTEREST IS 2015. 2014 NUMBERS AND BIOMASS ARE DETERMINED BY 2013 CATCH WHICH HAS ALREADY OCCURRED (PAC 2013)
	// SST 2015: THE PROJECTED YEAR OF INTEREST STARTS IN 2017. 2016 NUMBERS AND BIOMASS ARE DETERMINED BY 2015 CATCH WHICH HAS ALREADY OCCURRED
	int pyr = nyr+6; //projection years

	//get parameters - convert to data objects
	double pbo   = value(bo);
	double pso   = value (so);
	double pbeta = value(beta);
		  
	dvector p_bt(syr,pyr);
	dvector p_ft(syr,pyr);
	dvector p_N(syr,pyr);
	dvector p_S(syr,pyr);
	dvector p_rt(syr+sage,pyr);//dmatrix p_Ft(1,ngear,syr,pyr);

	p_bt.initialize();
	p_ft.initialize();
	p_N.initialize();
	p_S.initialize();
	p_rt.initialize();
	
	p_ft(syr,nyr)  = value(ft(1)(syr,nyr));
	p_N(syr,nyr)   = value(numbers(syr,nyr));
	p_bt(syr,nyr)  = value(biomass(syr,nyr)); //sbt and vul biomass all the same for delay diff
	p_S(syr,nyr)   = value(surv(syr,nyr));
	p_rt(syr+sage,nyr) = value(rt(syr+sage,nyr));
		
	//control points    - these are "historical" control points based on biomass and F reconstruction (P.cod specific)
	int nshort=2004-syr+1;
	int nlong=2012-syr+1;
	double meanfshort;  // average F between 1956 and 2004
	double meanflong;   // average F between 1956 and 2012
	double meanbshort;  // average B between 1956 and 2004
	double meanblong;   // average B between 1956 and 2012
	double minb;	     // biomass in 1971 for 5CD or 1985 for 5AB

	dvector hist_ftshort(syr,2004);
	dvector hist_ftlong(syr,2012);
	dvector hist_btshort(syr,2004);
	dvector hist_btlong(syr,2012);
	hist_ftshort.initialize();  hist_ftlong.initialize();
	hist_btshort.initialize();  hist_btlong.initialize();

	hist_ftshort=value(ft(1)(syr,2004));
	hist_btshort=value(biomass.sub(syr,2004));
	if(nyr>=2012){
		hist_ftlong=value(ft(1)(syr,2012));
		hist_btlong=value(biomass.sub(syr,2012));
	}
	meanfshort=sum(hist_ftshort)/nshort;
	if(nyr>=2012) meanflong=sum(hist_ftlong)/nlong;
	meanbshort=sum(hist_btshort)/nshort;
	if(nyr>=2012) meanblong=sum(hist_btlong)/nlong;
	minb=hist_btshort(1996);

	/* Simulate population into the future under constant tac policy. */
	
	for(i = nyr+1; i<=pyr; i++)
	{
		//recruits
		//double p_tau = value(sqrt(1-rho)/varphi);
		double p_tau = value(tau); 
		double xx = randn(nf+i)*p_tau;
				
		double rt;
		double et=p_bt(i-kage); //delay diff
		if(cntrl(2)==1)p_rt(i)=value((so*et/(1.+beta*et))*exp(xx-0.5*p_tau*p_tau));
		if(cntrl(2)==2)p_rt(i)=value((so*et*exp(-beta*et))*exp(xx-0.5*p_tau*p_tau));

		//numbers and biomass
		//Update biomass and numbers	
		p_bt(i) = (p_S(i-1)*(rho_g*p_bt(i-1)+alpha_g*p_N(i-1))+wk*p_rt(i));
		p_N(i)  = p_S(i-1)*p_N(i-1)+p_rt(i);
		//cout<<i<<"  "<<p_rt(i)<<"  "<<p_bt(i)<<"  "<<p_N(i)<<endl;
		//cout<<i<<"  "<<p_bt<<endl;

		//get_ft is defined in the Baranov.cxx file
		p_ft(i) = get_ftdd(tac,value(m),p_bt(i));	    //hardwiring the catch to gear 1 for this assessment       m_bar is same as constant M

		/*
		//test get_ftdd with Baranov equation
		double testf =p_ft(i);
		double testc = p_bt(i)*(1-mfexp(-value(m)-testf))*(testf/(value(m) + testf));
		cout<<i<<"  "<<testf<<"  "<<tac<<"  "<<testc<<endl<<endl;
		*/

		//Calculate mortality for next projection year
		p_S(i) = mfexp(-(value(m)+p_ft(i)));
	}

	if(mceval_phase()){
		if(nf==1 && runNo==1){
			cout<<"Running MCMC projections"<<endl;
			ofstream ofsmcmc("iscammcmc.proj.csv");
			write_proj_headers(ofsmcmc, syr, nyr, pyr);
			ofsmcmc.flush();
		}
		ofstream ofsmcmc("iscammcmc.proj.csv", ios::app);
		write_proj_output(ofsmcmc, syr, nyr, tac, pyr, p_bt, p_ft, value(bo), fmsy, bmsy, msy);
		ofsmcmc.flush();
	}
	
	  //S= Short (1956-2004) 
	   //L-Long(1956-2012)
//	if(mceval_phase()){
//		if(nf==1 && runNo==1)
//		{
//			ofstream ofsP("iscammcmc.proj");
//			ofsP<<"tac" <<setw(6)     <<   "\t";
//			ofsP<<"B2014" <<setw(6)     <<   "\t";
//			ofsP<<"B2015" <<setw(6)     <<   "\t";
//			ofsP<<"B2015B2014" <<setw(6)     <<   "\t";		   //want probability B2015<B2014 - this will be < 1 if true
//			ofsP<<"F2013" <<setw(6)     <<   "\t";
//			ofsP<<"F2014" <<setw(6)     <<   "\t";
//			ofsP<<"F2014F2013" <<setw(6)     <<   "\t";		   //want probability F2014>F2013     - this will be > 1 if true
//			//MSY based ref points
//			ofsP<<"BMSY" <<setw(6)     <<   "\t";
//			ofsP<<"B2015BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<BMSY - this will be < 1 if true
//			ofsP<<"B201508BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<0.8BMSY - this will be< 1 if true
//			ofsP<<"B201504BMSY" <<setw(6)     <<   "\t";		   //want probability B2015<0.4BMSY - this will be < 1 if true
//			ofsP<<"FMSY" <<setw(6)     <<   "\t";
//			ofsP<<"F2014FMSY"<<setw(6)     <<   "\t";		   //want probability F2014>F2013 - this will be > 1 if true
//			//Historical ref points "short"	 1956-2004
//			ofsP<<"Bmin" <<setw(6)     <<   "\t";
//			ofsP<<"B2015Bmin" <<setw(6)     <<   "\t";		   //want probability B2015<Bmin 
//			ofsP<<"BAvg_S" <<setw(6)     <<   "\t";
//			ofsP<<"B2015BAvg_S" <<setw(6)     <<   "\t";		   //want probability B2015<Bavg 
//			ofsP<<"FAvg_S" <<setw(6)     <<   "\t";
//			ofsP<<"F2014FAvg_S"<<setw(6)     <<   "\t";	
//			//Historical ref points "long"	 1956-2012
//			ofsP<<"BAvg_L" <<setw(6)     <<   "\t";
//			ofsP<<"B2015BAvg_L" <<setw(6)     <<   "\t";		   //want probability B2015<Bavg - this will be < 1 if true
//			ofsP<<"FAvg_L" <<setw(6)     <<   "\t";
//			ofsP<<"F2014FAvg_L\n";		   //want probability F2014>F2013 - this will be > 1 if true
//		      
//			cout<<"Running MCMC evaluations"<<endl;
//			cout<<"In delay diff proj Bo when nf==1 \t"<<bo<<endl;
//		}

//		ofstream ofsP("iscammcmc.proj",ios::app);
//		ofsP <<tac <<setw(6)                            <<"\t"
//		  << p_bt(pyr-1) <<setw(6)       <<"\t"	      
//		  << p_bt(pyr) <<setw(6)       <<"\t"		 
//		  << p_bt(pyr)/p_bt(pyr-1) <<setw(6)      <<"\t"	     
//		 << p_ft(pyr-2) <<setw(6)      <<"\t"
//		 << p_ft(pyr-1)  <<setw(6)     <<"\t"
//		 << p_ft(pyr-1)/p_ft(pyr-2)  <<setw(6)     <<"\t"	 
//		//MSY based ref points
//		<<bmsy <<setw(6)     <<   "\t"
//		<<p_bt(pyr)/bmsy <<setw(6)     <<   "\t"		 
//		<<p_bt(pyr)/(0.8*bmsy) <<setw(6)     <<   "\t"		  
//		<<p_bt(pyr)/(0.4*bmsy) <<setw(6)     <<   "\t"		   
//		<<fmsy <<setw(6)     <<   "\t"
//		<<p_ft(pyr-1)/fmsy <<setw(6)     <<   "\t"		   
//		//Historical ref points "short"	 1956-2004
//		<<minb <<setw(6)     <<   "\t"
//		<<p_bt(pyr)/minb <<setw(6)     <<   "\t"		   
//		<<meanbshort <<setw(6)     <<   "\t"
//		<<p_bt(pyr)/meanbshort <<setw(6)     <<   "\t"		   
//		<<meanfshort <<setw(6)     <<   "\t"
//		<<p_ft(pyr-1)/meanfshort<<setw(6)     <<   "\t"		  
//		 //Historical ref points "long"	 1956-2012
//		<<meanblong <<setw(6)     <<   "\t"
//		<<p_bt(pyr)/meanblong <<setw(6)     <<   "\t"		   
//		<<meanflong <<setw(6)     <<   "\t"
//		<<p_ft(pyr-1)/meanflong<<   "\t"		   	   		   
//		 <<endl;
//	   }

	//S= Short (1956-2004) 
	//L-Long(1956-2012)
	if(last_phase() && !mceval_phase()){
		if(runNo==1) {
			cout<<"Running MPD projections"<<endl;
			ofstream ofsP("iscammpd.proj.csv");
			ofsP<<"tac"                   <<",";
			ofsP<<"B"<<nyr+1              <<",";
			ofsP<<"B"<<nyr+2              <<",";
			ofsP<<"B"<<nyr+2<<"B"<<nyr+1  <<",";   //want probability B2015<B2014 - this will be < 1 if true
			ofsP<<"F"<<nyr                <<",";
			ofsP<<"F"<<nyr+1              <<",";
			ofsP<<"F"<<nyr+1<<"F"<<nyr    <<",";   //want probability F2014>F2013     - this will be > 1 if true
			//MSY based ref points
			ofsP<<"MSY"                   <<",";
			ofsP<<"BMSY"                  <<",";
			ofsP<<"B"<<nyr+2<<"BMSY"      <<",";   //want probability B2015<BMSY - this will be < 1 if true
			ofsP<<"B"<<nyr+2<<"BMSY80"    <<",";   //want probability B2015<0.8BMSY - this will be< 1 if true
			ofsP<<"B"<<nyr+2<<"BMSY40"    <<",";   //want probability B2015<0.4BMSY - this will be < 1 if true
			ofsP<<"FMSY"                  <<",";
			ofsP<<"F"<<nyr+1<<"FMSY"      <<",";   //want probability F2014>F2013 - this will be > 1 if true
			//Historical ref points "short"	 1956-2004
			ofsP<<"Bmin"                  <<",";
			ofsP<<"B"<<nyr+2<<"Bmin"      <<",";   //want probability B2015<Bmin 
			ofsP<<"BAvg_S"      <<",";
			ofsP<<"B"<<nyr+2<<"BAvg_S"    <<",";   //want probability B2015<Bavg  
			ofsP<<"FAvg_S"                <<",";
			ofsP<<"F"<<nyr+1<<"FAvg_S"    <<",";
			//Historical ref points "long"	 1956-2012
			ofsP<<"BAvg_L"                <<",";
			ofsP<<"B"<<nyr+2<<"BAvg_L"    <<",";   //want probability B2015<Bavg - this will be < 1 if true
			ofsP<<"FAvg_L"                <<",";
			ofsP<<"F"<<nyr+1<<"FAvg_L"    <<endl;  //want probability F2014>F2013 - this will be > 1 if true
		}
		cout<<"tac = "<<tac<<endl;
		ofstream ofsP("iscammpd.proj.csv",ios::app);
		ofsP
		<<tac                            <<","
		<< p_bt(nyr+1)                   <<","
		<< p_bt(nyr+2)                   <<","
		<< p_bt(nyr+2)/p_bt(nyr+1)       <<","
		<< p_ft(nyr)                     <<","
		<< p_ft(nyr+1)                   <<","
		<< p_ft(nyr+1)/p_ft(nyr)         <<","
		//MSY based ref points
		<<msy                            <<","
		<<bmsy                           <<","
		<<p_bt(nyr+2)/bmsy               <<","
		<<p_bt(nyr+2)/(0.8*bmsy)         <<","
		<<p_bt(nyr+2)/(0.4*bmsy)         <<","
		<<fmsy                           <<","
		<<p_ft(nyr+1)/fmsy               <<","
		//Historical ref points "short"	 1956-2004
		<<minb                           <<","
		<<p_bt(nyr+2)/minb               <<","
		<<meanbshort                     <<","
		<<p_bt(nyr+2)/meanbshort         <<","
		<<meanfshort                     <<","
		<<p_ft(nyr+1)/meanfshort         <<","
		//Historical ref points "long"	 1956-2012
		<<meanblong                      <<","
		<<p_bt(nyr+2)/meanblong          <<","
		<<meanflong                      <<","
		<<p_ft(nyr+1)/meanflong          <<endl;
	}
 }
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
 //Extra test functions by RF to test ref points
FUNCTION void slow_msy(dvector& ftest, dvector& ye, dvector& be, double& msy, double& fmsy, double& bmsy )
	//THIS CODE VERIFIES THAT THE EQM CODE IS RETURNING CORRECT REF POINTS
	//DD CODE HAS ALSO BEEN CHECKED IN AN R MODEL AND IS CORRECT
	//THIS MEANS THAT THE DIFFRENCES BETWEEN DD AND ASM IS DUE TO GROWTH PARAMETERS
	 
	int i;
	int j;
	int k;
	int t;
	int NF=size_count(ftest);
	int Nyr=100; //number of years to run out the model
	ye.initialize();
	be.initialize();
	
	double sa;
	dvector za(sage,nage); za.initialize();
	dvector saf(sage,nage); saf.initialize();
	dvector lx(sage,nage); lx.initialize();
	dvector vd(sage,nage); vd.initialize();
	
	vd=value(mfexp(log_sel(1)(nyr)));
	
	double Ro=value(ro);
	double CR=value(kappa);
	double M=value(m);
	
	dmatrix Nn(1,Nyr+1,sage,nage);			//Numbers at age
	dmatrix Ff(1,Nyr+1,sage,nage);			//Age-specific fishing mortality
	dmatrix Zz(1,Nyr+1,sage,nage);
	dvector Ss(sage,nage);
	dmatrix Cc(1,Nyr,sage,nage);
	dvector Ssb (1,Nyr+1);
	dvector Bb(1,Nyr+1);
	dvector Y(1,Nyr);				//predicted catch biomass
	
	dvector finaly(1,NF);
	dvector finalb(1,NF);
	
	//unfished
	sa=mfexp(-M);
	lx(sage)=1.0;
	for(i=(sage+1); i<=nage; i++)
		lx(i)=lx(i-1)*sa;
	lx(nage)/=(1.-sa); 
		
	//Initialize model - same for all F scenarios
	
	for(j=sage;j<=nage;j++) Nn(1,j)=Ro*lx(j);
	Ssb(1)= sum(elem_prod(Nn(1),avg_fec));
	
	for(k=1;k<=NF;k++){ 
		
		za=(M+ftest(k)*vd);
		Ss=mfexp(-za);
		
		/*
		cout<<"Nn "<<endl<<Nn(1)<<endl;
		cout<<"SS "<<endl<<Ss<<endl;
		cout<<"SSb "<<endl<<Ssb(1)<<endl;
		cout<<"wt "<<avg_wt<<endl; 
		*/
			
		for(t=1;t<=Nyr;t++){
			
			Nn(t+1)(sage+1,nage)=++elem_prod(Nn(t)(sage,nage-1),Ss(sage,nage-1));
			Nn(t+1,nage)+=Nn(t,nage)*Ss(nage);
			if(t==1) Nn(t+1)(sage)=Ro;
			if(t>1) Nn(t+1)(sage)=value(so)*Ssb(t-1)/(1.+value(beta)*Ssb(t-1));
			Ssb(t+1)= sum(elem_prod(Nn(t+1),avg_fec));
			
			//catch
			for(j=sage;j<=nage;j++) Cc(t,j) = ((ftest(k)*vd(j))/za(j))*(1.-exp(-(za(j))))*Nn(t,j)*avg_wt(j);
			Y(t)=sum(Cc(t));
		}//end t
		
	 finaly(k)=Y(Nyr);
	 finalb(k)=Ssb(Nyr);
	} //end k
	
	//cout<<"finalY "<<finaly<<endl; 
	//cout<<"finalB "<<finalb<<endl; 
	
	//get MSY and Fmsy
	msy=max(finaly);
	double mtest;
	for(k=1; k<=NF; k++)
	{
		mtest=finaly(k);
		if(mtest==msy) fmsy=ftest(k);
		if(mtest==msy) bmsy=Ssb(k);
	}
	
	cout<<"Ref points from running out model"<<endl;
	cout<<msy<<endl;
	cout<<fmsy<<endl;
	cout<<bmsy<<endl;
	
	
FUNCTION void ddiff_msy(dvector& ftest, dvector& ye, dvector& be, double& msy, double& fmsy, double& bmsy )
	//for testing delay diff ref points - slow way RF too lazy to find derivatives!
	int i;
	int k;
	int NF=size_count(ftest);
	ye.initialize();
	be.initialize();
	double se;
	double we;
	double rec_a=value(so);
	double rec_b=value(beta);
	double M = value(m);
	
	// Calculate equilibrium survivorship as function of FMSY
	for(k=1; k<=NF; k++)
	{
		 se = exp(-M -ftest(k));
		  we = (se*alpha_g+wk *(1.-se))/(1.-rho_g*se);
		  be(k) = -1*((-we + se*alpha_g + se*rho_g*we + wk*rec_a*we)/(rec_b*(-we + se*alpha_g + se*rho_g*we))); //Martell
		 
		  // Calculate equilibrium yield
		  ye(k)   = be(k)*(1-mfexp(-ftest(k)-M))*(ftest(k)/(ftest(k)+M));
		  if(ye(k)<0) ye(k)=0.;
		 if(be(k)<0) be(k)=0.;
	}	
	
	//get MSY and Fmsy
	msy=max(ye);
	double mtest;
	for(k=1; k<=NF; k++)
	{
		mtest=ye(k);
		if(mtest==msy) fmsy=ftest(k);
		if(mtest==msy) bmsy=be(k);
		
	}
	
	/*
	cout<<"Slow"<<endl;
	cout<<msy<<endl;
	cout<<fmsy<<endl;
	cout<<bmsy<<endl; 
	*/
	

FUNCTION void run_FRP()
	//Reference points
	if(last_phase()){
	cout<<"                                                            "<<endl;
	cout<<"*********Getting reference points************"<<endl;
	cout<<"*******************************************"<<endl;
	}
	
	dvector ftest(1,1001);
	ftest.fill_seqadd(0,0.001);
	int Nf;
	Nf=size_count(ftest);
	double Fmsy;
	double MSY;
	double Bmsy;
	dvector Ye(1,Nf); //Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
	dvector Be(1,Nf); //Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
	double fmsy,msy,bmsy,msy2;
	dvector ye(1,Nf);
	dvector be(1,Nf);

	
	slow_msy(ftest, ye, be, msy, fmsy, bmsy);
	Fmsy=fmsy;
	MSY=msy;
	Bmsy=bmsy;
	Ye=ye;
	Be=be;

	ofstream ofsr("TEST_frp.rep");
	ofsr<<"Fmsy"<<endl<<Fmsy<<endl;
	ofsr<<"MSY"<<endl<<MSY<<endl;
	ofsr<<"Bmsy"<<endl<<Bmsy<<endl;
	ofsr<<"ftest"<<endl<<ftest<<endl;
	ofsr<<"Ye"<<endl<<Ye<<endl;	
	ofsr<<"Be"<<endl<<Be<<endl;	

	
FUNCTION void run_FRPdd()
	//Reference points
	dvector ftest(1,1001);
	ftest.fill_seqadd(0,0.001);
	int Nf;
	Nf=size_count(ftest);
	dvector Ye(1,Nf); //Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
	dvector Be(1,Nf); //Matrix for putting numerically derived equilibrium catches for calculating MSY and FMSY (in R)
	dvector ye(1,Nf);
	dvector be(1,Nf);

	ddiff_msy(ftest, ye, be, msy, fmsy, bmsy);
	
		
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

	

GLOBALS_SECTION
 //	#define USE_NEW_EQUILIBRIUM
	typedef double DP;
	//#include <dfridr.cpp>
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
    
	#undef COUT
	#define COUT(object) cout << #object "\n" << object <<endl;

	#if defined(_WIN32) && !defined(__linux__)
		const char* PLATFORM = "Windows";
	#else
		const char* PLATFORM = "Linux";
	#endif

	#include <admodel.h>
	#include <time.h>
	#include <string.h>
	#include <statsLib.h>
	//#include "stats.cxx"
	#include "baranov.cxx"

	#undef NA
	#define NA -99.0
	#include "utilities.cpp"

	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	bool mcmcPhase = 0;
	bool mcmcEvalPhase = 0;
	
	adstring BaseFileName;
	adstring ReportFileName;
	
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
	
	class Selex
	{
		
	};
	
	
FINAL_SECTION
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
	cout<<"--Number of function evaluations: "<<nf<<endl;
	cout<<"--Results are saved with the base name:\n"<<"\t"<<BaseFileName<<endl;
	cout<<"*******************************************"<<endl;

	//Make copies of the report file using the ReportFileName
	//to ensure the results are saved to the same directory 
	//that the data file is in. This should probably go in the 
	//FINAL_SECTION
	
	//CHANGED only copy over the mcmc files if in mceval_phase()
	
	if(last_phase() && PLATFORM =="Linux" && !retro_yrs)
	{
		adstring bscmd = "cp iscam.rep " +ReportFileName;
		system(bscmd);
		
		bscmd = "cp iscam.par " + BaseFileName + ".par";
		system(bscmd); 
		
		bscmd = "cp iscam.std " + BaseFileName + ".std";
		system(bscmd);
		
		bscmd = "cp iscam.cor " + BaseFileName + ".cor";
		system(bscmd);
		
		if( mcmcPhase )
		{
			bscmd = "cp iscam.psv " + BaseFileName + ".psv";
			system(bscmd);
			
			cout<<"Copied binary posterior sample values"<<endl;
		}
		
		if( mcmcEvalPhase )
		{		
			bscmd = "cp iscam.mcmc " + BaseFileName + ".mcmc";
			system(bscmd);
		
			bscmd = "cp sbt.mcmc " + BaseFileName + ".mcst";
			system(bscmd);
		
			bscmd = "cp rt.mcmc " + BaseFileName + ".mcrt";
			system(bscmd);
		
			cout<<"Copied MCMC Files"<<endl;
		}
	}

	if( last_phase() && PLATFORM =="Linux" && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "cp iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}


