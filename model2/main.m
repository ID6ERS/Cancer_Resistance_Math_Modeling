%% Main matlab script
% this script should be used for plotting the relevant data and
% corresponding fitting achieved using "Density Games type" model. 
%  NOTE : the optimized parameters for model fitting have already been
%  stored in MAT files. we use these MAT files for plotting the fitting
%  obtained from model ( see modifiedDensityGames.m for definition)

%%
% Extract the data from CSV files. The CSV files have been created from
% original excel data itself.
[ncr,ncg,cr,cg,hours,ncx0,cx0]=extractData();

%%
% Plot with/without-cisplatin related data along with the corresponding "Density 
% games" inspired model fit. Also plot the errors in model fitting and payoffs 
% (changing with time in our model).

% for without-cisplatin plot, use the following
% plots_modelFit_payoffs_fitError(ncr,ncg,ncx0,hours,'withoutcis')
 
% for with-cisplatin plot, use the following 
 plots_modelFit_payoffs_fitError(cr,cg,cx0,hours,'withcis')
