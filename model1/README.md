# EGT_lotka-volterra
Codes used to fit the cancer cell growth to a lotka-volterra formalism

# Code files description
The codes attached here calculate the parameters resulting in the best fit of logistic and modified lotka-voterra model to the corresponding data. In each folder, the data used for fitting is stored in a file "data.mat". "model.m" contains the function that calculates the fitting error to be minimized. The mathematical model used for fittin in each case is included within the model function. The file "gr.m" has the code to iteratively estimate the parameters corresponding to best-fit corresponding to the data and the model. 

## Logistic_fits folder
In the folder Logistic_fits, the "model" calculates the fitting error between the data and a simple logistic model. In the gr.m file, notice two separate parameter estimates for red cells (sensitive cells, Bred) and green cells (resistant cells, Bgreen).

## Lotka-volterra_fits folder
In the folder lotka-volterra_fits, the "model" function calculates the fitting-error using the estimated model and 8 different experiments (each with different initial compositions of red and green cells). The experimental data is consolidated in data.m. The ODE function within the model, "DifEq" takes a parameter "opt", which corresponds to one of the four possible types of interactions tested between green and red cells. The equations are modified accordingly. Details of these options and the corresponding equations are given in the manuscript.

# Running the code
To run get the parameter estimates in either scenarios, one needs to run the gr.m file in MATLAB. The code has been tested in MATLAB versions R2018b and above, but should run fine in other versions as well. In the Logistic fits, the parameters are stored in Bred an Bgreen variables, whereas for Lotka-volterra fits, the parameters are stored in the table params, where each row corresponds to an "opt" and each column corresponds to a parameter. 
