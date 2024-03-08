The package RobustSharpe is a Matlab implementation of the robust Sharpe 
ratio testing based  on Ledoit & Wolf(2008, Journal of Empirical Finance).


Before running these commands, you need to store them in your Matlab working 
directory or change the working directory appropriately.


Of course, first of all, load your data into Matlab. 
For example, to load the data 
in the file rethedge.csv, write

data = load('...directory where the data file is...\rethedge.csv');



The main command is
				
                                          robustSharpe()

As an easy start, you can run
				
                                        robustSharpe(data)
				

to test the null hypothesis "H0: Difference of Sharpe ratios is zero". 

The matrix data MUST BE of size [T x 2], otherwise you get an error.
 
In this reduced form input of robustSharpe(data), all other parameters of 
robustSharpe() are set to the default. 

Note that the computation of the 
optimal block size as in optimalblrobustSharpe() is computed automatically 
as part of robustSharpe(data). 
The remaining time in seconds is printed in 
the console. On most machines, this should only take some minutes. 


The full input and help file can be found by typing

				
                                                      help robustSharpe

in Matlab.




The second main command is
				
                                         optimalblrobustSharpe()

Again, the minimal input is
 
    			               optimalblrobustSharpe(data)



As said before, if no prespecified block size is passed to robustSharpe(), 
the block size calibration is done automatically. 
If you have determined 
a prespecified block size by whatever means, you can pass it to robustSharpe() 
directly so that the computation time reduces to some seconds on most machines.
The command 
                       		       help optimalblrobustSharpe


gives details on possible inputs and outputs of the command. 




The main computational burden lies in the calibration of the block
size. 

The default number of bootstrap repetitions are M = 2,000 (for the
outer loop) and MBM = 200 (for the inner loop). 
The routine can be sped up
 by reducing these numbers (say to 1,000 and 100, respectively). 

Needless to say, the numbers should be picked as large as possible, 
given the computational resources, in view of accuracy.



All other routines are help routines. 
retagg.csv and rethedge.csv are the 
two data sets in the empirical applications of the paper. They correspond 
to ret.hedge and ret.agg in the R implementation.



In both robustSharpe() and optimalblrobustSharpe(), one can specify the 
kernel used to compute the prewhitened kernel variance estimator. 
The 
two options are Gallant/Parzen or Quadratic Spectral. 



The routine andmon6cvm2.m uses parts of the gmmopt package from Mike Cliff
version 1.1. 
However, you don't need to install the gmmopt package to run
 robustSharpe or any other of our routines, as the required parts of gmmopt
are is included in this package.



You are free to use this package as you see fit. Please give proper credit 
though. 



Bugs or comments may be reported to danwunderli@iew.uzh.ch.

Zurich, Institure for Empirical Research in Economics, U Zurich. March 2009.
