/*******************************************************************************************
Created By:         Brandy Ringham
Date:               6/22/13
Description:        Create power curve for hypothetical oral cancer biomarkers study. Use 
                    inputs from Elashoff et al., 2012.
                    
				    Fixed predictor = disease status (normal, case)
                    Dependent variables = mRNA (molecule1, molecule2, molecule3);
                    Hypothesis test = group difference for any mRNA
*******************************************************************************************/
title1 h = 12pt "Oral Cancer Biomarkers Demonstration";

footnote1 "&sysdate";
footnote2 "OCBiomarkersDemonstration.sas";

/*initiate SAS/IML*/
proc iml;

	/*define planned per group sample size*/
    unadjustRepNt = 75;

	/*vector of anticipated proportion of missing data*/
    piVec = { 0, .05, .1 };

	/*correlation between repeated measures*/
	rho = .4;

	/*define 7 inputs of power analysis*/

	/*1.  alpha*/
	alpha = .05;

	/*2.  between subject contrast matrix*/
	/*dimensions 1x2*/
	c = { 1 -1 };

	/*3. within subject contrast matrix*/
	/*dimensions 3x2*/
	u = { 1 0 0,
          0 1 0,
          0 0 1 };
	
	/*4. sigma*/
	/*see below--values will vary so it resides in the do loop*/

	/*5. beta */
	/*see below--values will vary so it resides in the do loop*/

	/*6. null hypothesis*/
	/*set by POWERLIB as a 1x2 matrix of 0's by default*/

	/*7. essenceX (Muller and Stewart, 2006)*/
	/*dimensions 2x2*/
	essenceX = { 1 0,
	             0 1 };

	/*calculate total planned sample size - unadjusted*/
	unadjustNt = unadjustRepNt * nrow( essenceX );
	
	/*include NMK code*/
	/*NMK module calculate adjusted sample size, or E(Nmk)*/
	%include "NMK.SAS" / nosource2;

	/*include POWERLIB code*/
	/*POWERLIB module calculates power*/
	%include "POWERLIB22.IML" / nosource2;

	/*set options for POWERLIB*/
	opt_off = { HF UN BOX WARN }; 
	opt_on = { UNIFORCE DS FRACREPN NOPRINT };

	/*loop through rows of proportion missing vector*/
	do piID = 1 to nrow( piVec );

		/*define proportion missing as the value in the current row of the vector*/
		pi = piVec[ piID, ];

		/*loop through different scale factors for the beta matrix*/
		do k = 0 to 2 by .005;

		    /*variances are from Table 3, Elashoff et al., 2012*/
			var = 2.9**2;

			/*define sigma matrix*/
			oneMat = j( 3, 1, 1 );
			sigma = var * ( oneMat * oneMat` * rho + I( 3 ) * ( 1 - rho ) );

			/*cell means are from for 
			1) IL-1B, 
		    2) IL-8, and 
			3) SAT, Cohort 4 reported in Table 3, Elashoff et al., 2012*/

			/*first set means for normals*/
			mu1n = 20.1;
			mu2n = 19.8;
			mu3n = 21.3;

			/*define difference between normals and cases*/
			delta1 = -1.3;
			delta2 = -2.1;
			delta3 = -1.4;

			/*calculate mean for cases based on scale factor*/
			mu1c = mu1n + k * delta1;
			mu2c = mu2n + k * delta2;
			mu3c = mu3n + k * delta3;

			/*define beta matrix from means*/
			beta = ( mu1c || mu2c || mu3c ) //
                   ( mu1n || mu2n || mu3n );

			/*calculate sample sizes adjusted for anticipated amount of missing data*/
			/*use Nm1*/
			call NMK( unadjustNt, 2, pi, 1, 0, ENmk, error );

			/*calculated adjusted number of participants per group*/
			/*set repn equal to adjusted repn to pass into POWERLIB*/
			repn = ENmk / nrow( essenceX );

			/*save power results in a dataset*/
			dsname = { work adjustPower };

			/*calculate power adjusted for missing data using POWERLIB*/
			run power;

			use work.adjustPower;

			/*create an IML matrix from the POWERLIB output dataset*/
			read all var{ ALPHA TOTAL_N POWER_MULT EPSILON EXEPS_GG POWER_GG } into adjustPower[ colname = name ];

			/*delete dataset after use*/
			call delete( "work", "adjustPower" );

			/*row of output containing experimental conditions and power results*/
			powerRow = rho || k || pi || repn || 2 || adjustPower[ , 1:3 ];

			/*create matrix of results over all proportion missing*/
			power = power // powerRow;

			/*free variables to be reused*/
			free powerRow ENmk repn dsname adjustPower;

		end; /*end k loop*/

	end; /*end pi loop*/

	/*list of column names for output dataset*/
	expNames = { "rho" "k" "PropMiss" "unadjustRepNt" "NumDepVars" };
	allNames = expNames || name[ , 1:3 ];

	/*output matrix as SAS dataset*/
	create out01.appPow from power[ colname = allNames ];
	append from power;

quit;

/*set size of graph*/
goptions hsize = 4in vsize = 4in;

/*plot power curve*/
proc gplot data = out01.appPow;

	title2 h = 12pt "Power Curve";

	/*plot power by scale factor*/
	plot power_mult * k = propMiss / noframe vaxis = axis1 haxis = axis2 legend = legend1
                                     href = 1 vref = .9 cvref = red chref = red;
	
	/*define linetypes for different proportion missing*/
	symbol1 i = j l = 1 c = black w = 1;
	symbol2 i = j l = 41 c = black w = 1;
	symbol3 i = j l = 35 c = black w = 1;

	/*define axes options*/
	axis1 order = ( 0 to 1 by .2 ) minor = none major = none label = ( angle = 90 font = times h = 12pt "Power" )
          value = ( h = 12pt font = times );
	axis2 order = ( 0 to 2 by .5 ) minor = none major = none label = ( font = times h = 12pt "k" )
          value = ( h = 12pt font = times );

	/*define legend*/
	legend1 value = ( height = 12pt font = "times" "0%" "5%" "10%" ) position = ( inside bottom right )
            label = ( font = times h = 12pt position = top justify = center "Percent Missing" ) mode = protect 
            across = 1 down = 3 frame shape = line( .3in );

run;
quit;

title;
footnote;
