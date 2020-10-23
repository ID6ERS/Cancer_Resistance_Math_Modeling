------------------------ README ------------------------
Please place all R scripts and the excel file in the same directory before running

fit_all.R: 
------------
This R script performs fitting of cellular growth data to PSMSR using the GA algorithm.
Following are the inputs to this script:

1. datafile: An excel file containing mean cell population with time for different S:T ratios. 
Each sheet in the excel file represents different experiment conditions, e.g. mixed, mixed 3 weeks etc.

2. dataset_name: Name of the excel sheet from which the data should be read.

3. fit_cisplatin: True or False; Is the data to be fit for the cisplatin free 
(phenotype switch and stress removal turned on) or 'with cisplatin' model.

4. labels: labels of the different trends used in plotting

5. Initial guesses: parameters$K0, parameters$Kb0, etc.; Initial guesses to GA

6. fixed_param_list: If certain parameters should be kept fixed during GA optimization, 
replace NA with fixed values for those parameters only. All the other parameters (NAs)
will be optimized.

7. Kmax: Maximum K0 for the linear regime of phenotype switch - stress relationship.

8. upper, lower: Upper and lower limits of the search spaces used by GA.

9. ratios: Initial proportions of switched phenotypes in the population, 
set to 10% of the sensitive and tolerant populations.

Required libraries: readxl, parallel, foreach, GA, desolve


functions.R:
--------------
Miscellaneous functions used by the GA fitting.


parameters.R:
--------------
Sets the parameters corresponding to the fitting of different datasets


simulate_intermittent2_therapy.R
---------------------------------
Simulates intermittent therapy with the sensitive and tolerant cell populations using PSMSR
Following are the inputs to this script:

1. parameters_nodrug: parameters to be used to simulate cell growth without cisplatin
2. parameters_wdrug: parameters for simulating cell growth with cisplatin
Both these parameters are set to the parameter set obtained by fitting the mixed culture data in presence of cisplatin.

3. Ctot: Initial cell population
4. S_T_ratio: Initial sensitive to tolerant cell ratio
5. ratios: Initial proportions of switched phenotypes in the population, 
set to 10% of the sensitive and tolerant populations.
6. drug_conc: Cisplatin concentration used in the simulations
7. therapy_regimen: time points where drug concentration is changed or cell passaging is performed.
	times: time-points
	drug: drug concentration applied
	media: Stress level (set to 0.001 due to cell passaging)
	fraction: New cell population after passaging
	
Required libraries: readxl, parallel, foreach, plotrix