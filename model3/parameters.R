#------------mixed------------#
parameters_mixed = list(
  K0 = 4.87E-02,
  Kb0 = 3.54,
  Kgs0 = 0.713,
  Kgt0 = 0.682,
  Ks = 7.14E-04,
  Ksd = 5.67E-03,
  a = 0.045,
  b = 0.037,
  g = 0.015,
  
  Sdrug = 90.73,
  
  AUCs5 = 258.68,
  AUCs95 = 918.86,
  SCALEdrugs = 0,
  
  AUCt5 = 221.24,
  AUCt95 = 947.96,
  SCALEdrugt = 0
)

#------------mixed, 3 weeks------------#
parameters_mixed_3weeks = list(
  K0 = 0.0498,
  Kb0 = 2.66,
  Kgs0 = 0.695,
  Kgt0 = 0.313,
  Ks = 0.000672,
  Ksd = 0.00538,
  a = 0.05,
  b = 0.038,
  g = 0.016,
  
  Sdrug = 93.09,
  
  AUCs5 = 242.58,
  AUCs95 = 1337.33,
  SCALEdrugs = 0,
  
  AUCt5 = 193.28,
  AUCt95 = 1304.72,
  SCALEdrugt = 0
)

#------------mixed, cisplatin------------#
parameters_mixed_cisplatin = list(
  K0 = 2.73E-02,
  Kb0 = 0,
  Kgs0 = 1.041,
  Kgt0 = 0.949,
  Ks = 6.16E-04,
  Ksd = 0,
  a = 0.06,
  b = 0.041,
  g = 0.054,
  
  Sdrug = 93.09,
  
  AUCs5 = 242.58,
  AUCs95 = 1337.33,
  SCALEdrugs = 8.48,
  
  AUCt5 = 193.28,
  AUCt95 = 1304.72,
  SCALEdrugt = 4.36
)

#------------mixed, 3 weeks, cisplatin------------#
parameters_mixed_3weeks_cisplatin = list(
  K0 = 0.02444,
  Kb0 = 0,
  Kgs0 = 0.816,
  Kgt0 = 0.964,
  Ks = 0.000726,
  Ksd = 0,
  a = 0.057,
  b = 0.055,
  g = 0.064,
  
  Sdrug = 90.73,
  
  AUCs5 = 258.68,
  AUCs95 = 918.86,
  SCALEdrugs = 5.64,
  
  AUCt5 = 221.24,
  AUCt95 = 947.96,
  SCALEdrugt = 3.53
)

Kmax = 0.04