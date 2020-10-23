library(readxl)
library(parallel)
library(foreach)
source("functions.R")

ncore = detectCores()
cl = parallel::makeForkCluster(ncore)
doParallel::registerDoParallel(cl)

# Excel file name
# The data is organized as: <Time, days> <sensitive cell count, ratio 1> <tolerant cell count, ratio 1> <sensitive cell count, ratio 2> <tolerant cell count, ratio 2> ...
datafile = "~/data/PostDoc Stuff/Core Projects/Salgia/misc/Atish_data/data_cleaned.xlsx"

# Excel sheet name containing the dataset to be used
dataset_name = "mixed 3 weeks cisplatin"

# Fit normal growth data or with cisplatin present
# Fitting with cisplatin requires setting community cooperation parameters to be set to zero (Kb0, Ksd)
fit_cisplatin = T
tdrug = numeric()

# Labels for graph legends
labels = c("1:1", "1:2", "1:4", "1:8")

parameters = list()
if(!fit_cisplatin) {
# Initial guesses
parameters$K0 = 0.0487
parameters$Kb0 = 3.54
parameters$Kgs0 = 0.713
parameters$Kgt0 = 0.682
parameters$Ks = 0.000714
parameters$Ksd = 0.00567
parameters$a = 0.045
parameters$b = 0.037
parameters$g = 0.015

parameters$Sdrug = 0

parameters$AUCs5 = 0
parameters$AUCs95 = 0
parameters$SCALEdrugs = 0

parameters$AUCt5 = 0
parameters$AUCt95 = 0
parameters$SCALEdrugt = 0

# Some parameters can be set to be fixed, while the others are optimized
fixed_param_list = list()
fixed_param_list$K0 = NA #0.01144903
fixed_param_list$Kb0 = NA #1.652723
fixed_param_list$Kgs0 = NA #0.3846179
fixed_param_list$Kgt0 = NA #0.1788549
fixed_param_list$Ks = NA #0.0007027514
fixed_param_list$Ksd = NA
fixed_param_list$a = NA #0.06679257
fixed_param_list$b = NA
fixed_param_list$g = NA #0.09161607

fixed_param_list$Sdrug = 0
fixed_param_list$AUCs5 = 0
fixed_param_list$AUCs95 = 0
fixed_param_list$SCALEdrugs = 0

fixed_param_list$AUCt5 = 0
fixed_param_list$AUCt95 = 0
fixed_param_list$SCALEdrugt = 0

sensitives_plot_ylim = c(-2000,40000)
tolerants_plot_ylim = c(-2000,65000)
} else {
  parameters$K0 = 0.03
  parameters$Kb0 = 0 #3.54
  parameters$Kgs0 = 0.713
  parameters$Kgt0 = 0.682
  parameters$Ks = 0.000714
  parameters$Ksd = 0
  parameters$a = 0.045
  parameters$b = 0.037
  parameters$g = 0.015
  
  parameters$Sdrug = 75.46712
  
  parameters$AUCs5 = 89.05143
  parameters$AUCs95 = 616.0868
  parameters$SCALEdrugs = 10
  
  parameters$AUCt5 = 180.1463
  parameters$AUCt95 = 376.5121
  parameters$SCALEdrugt = 4
  
  fixed_param_list = list()
  fixed_param_list$K0 = NA #0.01144903
  fixed_param_list$Kb0 = 0 #1.652723
  fixed_param_list$Kgs0 = NA #0.3846179
  fixed_param_list$Kgt0 = NA #0.1788549
  fixed_param_list$Ks = NA #0.0007027514
  fixed_param_list$Ksd = 0
  fixed_param_list$a = NA #0.06679257
  fixed_param_list$b = NA
  fixed_param_list$g = NA #0.09161607
  
  fixed_param_list$Sdrug = NA
  fixed_param_list$AUCs5 = NA
  fixed_param_list$AUCs95 = NA
  fixed_param_list$SCALEdrugs = NA
  
  fixed_param_list$AUCt5 = NA
  fixed_param_list$AUCt95 = NA
  fixed_param_list$SCALEdrugt = NA
  
  sensitives_plot_ylim = c(-2000,15000)
  tolerants_plot_ylim = c(-2000,40000)
}

Kmax = 0.04

# Upper and lower bounds for GA
#          K0     Kb0   Kgs0    Kgt0   Ks       Ksd     a     b     g   Sdrug   AUCs5   AUCs95  SCALEdrugs  AUCt5   AUCt95  SCALEdrugt
upper = c(0.1,     5,    2.0,    2.0, 0.001,   0.02,  0.1,  0.1,  0.1,  500,    500,     1500,      20,        500,      1500,      20)
lower = c(0,      0,     0,       0,     0,    0,     0,    0,    0,     1,      1,      10,      0.1,      1,        10,      0.1)

# Please ignore this
fit_all = T

param_names = names(fixed_param_list)
npar = length(parameters)

ratios = list(
  Cso0_ratio = 0.9,
  Css0_ratio = 0.1,
  Cto0_ratio = 0.9,
  Cts0_ratio = 0.1
)
initial_values = list( ratios = ratios, S0 = 0.001)

# Read exp. data
dataset_orig = read_excel(datafile, sheet = dataset_name)
t = as.numeric(dataset_orig[,1][[1]])
nset = (ncol(dataset_orig) - 1)/2
x_col = seq(2,2*nset,2)
ndata = length(t)
nrep = 1
x1_exp = matrix(nrow = ndata, ncol = nset)
x2_exp = matrix(nrow = ndata, ncol = nset)
for(i in 1:nset) {
  x1_exp[,i] = as.numeric(dataset_orig[,x_col[i]][[1]])
  x2_exp[,i] = as.numeric(dataset_orig[,x_col[i]+1][[1]])
}

colors = c("violet", "blue", "green", "orange", "red", "darkred", "black")
Ctot_list = matrix(nrow = ndata, ncol = nset * nrep)

parameters_preoptim = parameters

GA = NA
while(is.na(GA)) {
  GA <- tryCatch({
    ga(type = "real-valued",
       fitness =  function(x,param_names, t, x1_exp, x2_exp, fixed_param_list, initial_values, Kmax) -optimFunc_all(x,param_names, t, x1_exp, x2_exp, fixed_param_list, initial_values, Kmax),
       param_names, t, x1_exp, x2_exp, fixed_param_list, initial_values, Kmax,
       lower = lower, upper = upper,
       suggestions = as.vector(unlist(parameters_preoptim)),
       monitor = GAmonitor_new,
       popSize = 50, maxiter = 5000, run = 100, parallel = T)
  }, error = function(c) NA
  )
}

parameters_allatonce = ParamVecToList( c(GA@solution[1,]), param_names )

stopCluster(cl)

x1_pred = matrix(nrow = ndata, ncol = nset)
x2_pred = matrix(nrow = ndata, ncol = nset)
Cso = matrix(nrow = ndata, ncol = nset)
Css = matrix(nrow = ndata, ncol = nset)
Cto = matrix(nrow = ndata, ncol = nset)
Cts = matrix(nrow = ndata, ncol = nset)
S = matrix(nrow = ndata, ncol = nset)
x1_err = numeric(nset)
x2_err = numeric(nset)
for(i in 1:nset) {
  Cs0 = x1_exp[1,i]
  Ct0 = x2_exp[1,i]
  S0 = initial_values$S0
  initial_conditions = list(Cso0 = ratios$Cso0*Cs0, Css0 = ratios$Css0*Ct0, Cto0 = ratios$Cto0*Ct0, Cts0 = ratios$Cts0*Cs0, S0 = S0)
  if(fit_all) {
    parameters_vec = as.numeric(parameters_allatonce)
    if(sum(!is.na(fixed_param_list)) == 0) {
      parameters_allatonce = ParamVecToList(parameters_vec, param_names)
    } else {
      temp = as.numeric(fixed_param_list)
      variable_ind = which(is.na(temp))
      temp[variable_ind] = parameters_vec[variable_ind]
      parameters_optim = ParamVecToList(temp, param_names)
      parameters_allatonce = parameters_optim
    }
  } else {
    parameters_optim = ParamVecToList(as.numeric(all_parameters[i,2:10]), param_names)
  }
  output = ode(as.numeric(initial_conditions), t, calcDerivatives_new, list(parameters=parameters_optim, Kmax=Kmax),  method = "lsoda")
  x1_pred[,i] = output[,2] + output[,5]
  x2_pred[,i] = output[,4] + output[,3]
  Cso[,i] = output[,2]
  Css[,i] = output[,3]
  Cto[,i] = output[,4]
  Cts[,i] = output[,5]
  S[,i] = output[,6]
  x1_err[i] = mean(abs(x1_pred[,i] - x1_exp[,i])) / diff(range(c(x1_exp[,i])))
  x2_err[i] = mean(abs(x2_pred[,i] - x2_exp[,i])) / diff(range(c(x2_exp[,i])))
}

par.old = par()
par(mfcol = c(2,1), mar=c(4,4,1,0.8))
plot(0,0, xlim = c(0,max(t)+1), ylim = sensitives_plot_ylim, xlab = "Time, days", ylab = "Sensitive Cell count")
for(i in 1:nset) {
  points(t, x1_pred[,i], 'l', col = colors[i], lwd = 2)
  points(t, x1_exp[,i], 'p', col = colors[i])
}
legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
plot(0,0, xlim = c(0,max(t)+1), ylim = tolerants_plot_ylim, xlab = "Time, days", ylab = "Tolerant Cell count")
for(i in 1:nset) {
  points(t, x2_pred[,i], 'l', col = colors[i], lwd = 2)
  points(t, x2_exp[,i], 'p', col = colors[i])
}
legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
par(par.old)

