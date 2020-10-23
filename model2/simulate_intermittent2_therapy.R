library(readxl)
library(parallel)
library(foreach)
library(plotrix)
source("functions.R")
source("parameters.R")

#------------no drug-----------#
parameters_nodrug = parameters_mixed_cisplatin
#------------no drug-----------#

#---------- with drug----------#
parameters_wdrug = parameters_mixed_cisplatin
#---------- with drug----------#

Ctot = 50000
S_T_ratio = 4
Cs0 = S_T_ratio/(1+S_T_ratio) * Ctot
Ct0 = 1/(1+S_T_ratio) * Ctot

ratios = list(
  Cso0_ratio = 0.9,
  Css0_ratio = 0.1,
  Cto0_ratio = 0.9,
  Cts0_ratio = 0.1
)
initial_values = list( ratios = ratios, S0 = 0.001)
dt = 1
drug_conc = 50

tdrug = numeric()
Sdrug_therapy = numeric()

# #initial intermittent
# therapy_regimen = data.frame(times = c(3,         10),
#                              drug = c(drug_conc, 0),
#                              media =c(0.001,     0.001),
#                              fraction = c(1,         1/4))
# output_intermittent_initial = simulate_intermittent_therapy(Cs0, Ct0, initial_values, parameters_nodrug, parameters_wdrug, therapy_regimen, dt = dt, Kmax = 0.04)
# Cs0 = (output_intermittent_initial$Cso + output_intermittent_initial$Css)[length(output_intermittent_initial$t)]/4
# Ct0 = (output_intermittent_initial$Cto + output_intermittent_initial$Cts)[length(output_intermittent_initial$t)]/4

# intermittent 1
therapy_regimen = data.frame(times = c(4,     7,     7,     7),
                              drug = c(0,     0,     0,     0),
                             media = c(0.001, 0.001, 0.001, 0.001),
                         fraction = c(1,      1/4,   1/4,   1/4))
output_intermittent1 = simulate_intermittent_therapy(Cs0, Ct0, initial_values, parameters_nodrug, parameters_wdrug, therapy_regimen, dt = dt, Kmax = Kmax)

# intermittent 2
therapy_regimen = data.frame(times = c(4,     7,     7,     7),
                          drug = c(drug_conc, 0,     0,     0),
                             media = c(0.001, 0.001, 0.001, 0.001),
                             fraction = c(1,      1/4,   1/4,   1/4))
output_intermittent2 = simulate_intermittent_therapy(Cs0, Ct0, initial_values, parameters_nodrug, parameters_wdrug, therapy_regimen, dt = dt, Kmax = Kmax)

# continuous
therapy_regimen = data.frame(times = c(4,        7,         7,         7),
                             drug = c(drug_conc, drug_conc, drug_conc, drug_conc),
                             media =c(0.001,     0.001,     0.001,     0.001),
                         fraction = c(1,         1/4,       1/4,       1/4))
output_continuous = simulate_intermittent_therapy(Cs0, Ct0, initial_values, parameters_nodrug, parameters_wdrug, therapy_regimen, dt = dt, Kmax = Kmax)

Cs = output_intermittent1$Cso + output_intermittent1$Css
Ct = output_intermittent1$Cto + output_intermittent1$Cts
tolerant_sensitive_ratio_fc_intermittent1 = ((Ct/Cs)/(Ct0/Cs0))

Cs = output_intermittent2$Cso + output_intermittent2$Css
Ct = output_intermittent2$Cto + output_intermittent2$Cts
tolerant_sensitive_ratio_fc_intermittent2 = ((Ct/Cs)/(Ct0/Cs0))

Cs = output_continuous$Cso + output_continuous$Css
Ct = output_continuous$Cto + output_continuous$Cts
tolerant_sensitive_ratio_fc_continuous = ((Ct/Cs)/(Ct0/Cs0))

yrange = range(c(tolerant_sensitive_ratio_fc_intermittent1, tolerant_sensitive_ratio_fc_intermittent2)) #, tolerant_sensitive_ratio_fc_continuous))
ygap = c(max(tolerant_sensitive_ratio_fc_intermittent1)*1.2, max(tolerant_sensitive_ratio_fc_continuous)*0.8)

# plot intermittent1
t = output_intermittent1$t
#gap.plot(t,tolerant_sensitive_ratio_fc_intermittent, gap = ygap, type = "b", pch = 19, col = "black", xlab = "days", ylab = "Fold change in tolerant:sensitive ratio", ylim = yrange)
plot(t,tolerant_sensitive_ratio_fc_intermittent1, type = "b", pch = 19, col = "black", xlab = "days", ylab = "Fold change in tolerant:sensitive ratio", ylim = yrange)

# plot intermittent2
t = output_intermittent2$t
#gap.plot(t,tolerant_sensitive_ratio_fc_intermittent, gap = ygap, type = "b", pch = 19, col = "black", xlab = "days", ylab = "Fold change in tolerant:sensitive ratio", ylim = yrange)
points(t,tolerant_sensitive_ratio_fc_intermittent2, type = "b", pch = 19, col = "red", xlab = "days", ylab = "Fold change in tolerant:sensitive ratio", ylim = yrange)

# plot continuous
t = output_continuous$t
#gap.plot(t,tolerant_sensitive_ratio_fc_continuous, gap = ygap, type = "b", pch = 19, col = "red", add = T)
points(t[1:25],tolerant_sensitive_ratio_fc_continuous[1:25], type = "b", pch = 19, col = "blue")
legend("topleft", legend = c("Intermittent 1", "intermittent 2", "continuous"), col = c("black", "red", "blue"), lty = 1, pch = 19)
