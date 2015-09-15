/****************************************************************************************
Date:					6/28/2013
Created By:				Brandy Ringham
Description:			This module calculates the expected value of N(m1), N(m2), and 
                        N(m9) (Catellier and Muller, 2000; Ringham et al., in 
                        submission).

						Copyright (C) 2010 Regents of the University of Colorado.

                        This program is free software; you can redistribute it and/or 
                        modify it under the terms of the GNU General Public License as 
                        published by the Free Software Foundation; either version 2 of 
                        the License, or (at your option) any later version. This program 
                        is distributed in the hope that it will be useful, but WITHOUT 
                        ANY WARRANTY; without even the implied warranty of 
                        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
                        General Public License for more details. You should have received 
                        a copy of the GNU General Public License along with this program; 
                        if not, write to the Free Software Foundation, Inc., 51 Franklin 
                        Street, Fifth Floor, Boston, MA 02110-1301, USA.
							
Inputs:					Nt		number of independent sampling units
						p		number of repeated measures
						pi		probability missing value
						k		index for missing data summary statistic (see index key)
						print	1 = print on, otherwise printing is turned off

Output:					ENmk	label for mean of the missing data summary statistic
						error	label for error code matrix

Error Code Key:			The error matrix contains error indicators for the inputs that
                        are checked in the beginning of the module. Error indicators are 
                        defined as follows:

						0		no error
						1		error

						The error indicator's position in the error matrix defines the 
                        referent input.

						Row		Input
						1		error state of Nt
						2		error state of p
						3		error state of pi
						4		error state of k

Usage:					call NMK( 12, 4, 0.10, 1, 1, ENm1, error );

****************************************************************************************/

/*begin module definition*/
start NMK( Nt, p, pi, k, print, ENmk, error );

	/*initialize error indicators to 0*/
	error = j( 4, 1, 0 );

	/*check that Nt is valid*/
	if Nt < 0 then do;

		error[ 1, 1 ] = 1;

		if print = 1 then do;

			print "invalid number of independent sampling units";
			print "Nt should be a positive integer";
			print "currently, Nt =" Nt;
			print "no expected values or variances were calculated";

		end;

	end;

	/*check that p is valid*/
	if p < 0 then do;

		error[ 2, 1 ] = 1;

	   if print = 1 then do;

			print "invalid number of repeated measures";
			print "p should be a positive integer";
			print "currently, p =" p;
			print "no expected values or variances were calculated";

		end;

	end;

	/*check that pi is valid*/
	if pi < 0 | pi > 1 then do;

		error[ 3, 1 ] = 1;

	   if print = 1 then do;

			print "invalid value for pi";
			print "pi should be a number between 0 and 1";
			print "currently, pi =" pi;
			print "no expected values or variances were calculated";

		end;

	end;

	/*check that k is valid*/
	if k ^= 1 &
	   k ^= 2 &
	   k ^= 9 then do;

	   error[ 4, 1 ] = 1;

	   if print = 1 then do;

			print "invalid index for the missing data summary statistic";
			print "k should be in the set {1, 2, 9}";
			print "currently, k = " k;
			print "No expected values or variances were calculated";

		end;

	end;

	/*only complete the next set of steps if there were no errors*/
	/*otherwise, program skips over the next block and ends with no calculations*/
	if sum( error ) = 0 then do; 

		/*calculate expected value of N(m1)*/
		if k = 1 then ENmk = Nt * ( 1 - pi )**p;

		/*calculate expected value of for N(m2)*/
		else if k = 2 then do;

			if pi = 0 then ENmk = Nt;

			else if pi = 1 then ENmk = 0;

			else

			/*regression model for expected value of N(m2)*/
			/*see Table 4, Ringham et al. (in submission)*/
			ENmk =     62.7318676
					-   5.1567678 * ( Nt / 10 )
					-   0.1324363 * ( Nt / 10 )**2
					-   1.6042196 * p
					+   0.0640387 * p**2
					- 147.7255861 * ( 1 - pi )
					+  87.3472243 * ( 1 - pi )**2
					-   0.4981166 * ( Nt / 10 ) * p 
					+   0.0019218 * ( ( Nt / 10 ) * p )**2
 					+   1.2421550 * p * ( 1 - pi )
					-   0.0423721 * ( p * ( 1 - pi ) )**2
					+  14.8545994 * ( Nt / 10 ) * ( 1 - pi )
					+   0.1812043 * ( ( Nt / 10 ) * ( 1 - pi ) )**2 
					+   0.4137540 * ( Nt / 10 ) * p * ( 1 - pi )
					-   0.0013272 * ( ( Nt / 10 ) * p * ( 1 - pi ) )**2; 

		end;

		/*calculate expected value of Nm9*/
		else if k = 9 then ENmk = Nt * ( 1 - pi );

	end;

/*end module*/
finish;
