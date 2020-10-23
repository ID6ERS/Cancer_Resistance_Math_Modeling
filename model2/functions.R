library(GA)
library(deSolve)
library(parallel)
library(memoise)

deriv1 = function(y,x) {
  n = length(x)
  h = x[2] - x[1]
  d = numeric(n)
  for(i in 1:n) {
    if(i == 1) {
      d[i] = (y[2] - y[1]) / h
    } else if(i == n) {
      d[i] = (y[n] - y[n-1]) / h
    } else {
      d[i] = (y[i+1] - y[i-1]) / 2/h
    }
  }
  return(d)
}

sigmoid = function(x,x5,x95,ymin,ymax,isplot = F) {
  ln19 = 2.944439
  y = ymin + (ymax-ymin) / (1 + exp(ln19*(1-2*(x-x5)/(x95-x5))))
  if(isplot) {
    xrange = x95-x5
    yrange = ymax-ymin
    plot(x,y,'l',lwd = 2, xlim = c(x5-0.1*xrange,x95+0.1*xrange), ylim = sort(c(ymin-0.1*yrange,ymax+0.1*yrange)))
    abline(v = c(x5,x95), lty = 2)
    abline(h = c(ymin,ymax), lty = 2, col = "red")
  }
  return(y)
}

AUCfunc = function(t, S) {
  if(length(tdrug) == 0) {
    #print("problem")
    return(t*S)
  } else {
    #print(t)
    return(integrate(f = function(x)(return(approx(tdrug,Sdrug_therapy,xout = x, rule = 2)$y)), lower = 0, upper = t)$value)
  }
}

calcRateConstants = function(t, Stress, parameters, Kmax) {
  K0 = parameters$K0  # Ka/Kb
  Kb0 = parameters$Kb0
  Kgs0 = parameters$Kgs0
  Kgt0 = parameters$Kgt0
  Ks = parameters$Ks
  Ksd = parameters$Ksd
  a = parameters$a
  b = parameters$b
  g = parameters$g
  
  Sdrug = parameters$Sdrug
  
  AUCs5 = parameters$AUCs5
  AUCs95 = parameters$AUCs95
  SCALEdrugs = parameters$SCALEdrugs
  
  AUCt5 = parameters$AUCt5
  AUCt95 = parameters$AUCt95
  SCALEdrugt = parameters$SCALEdrugt
  
  Si = Stress
  
  if(!is.na(Kmax)) {
    St = (Kmax/K0-1)/4/a
  } else {
    St = -1
  }
  if(St <= 0)
    St = 1e30
  a1 = 1/St
  
  if(Si > St) {
    K = K0*(1+a*St) -(Kmax-K0)/2 + (Kmax-K0) * a1*Si/(1+a1*Si)
  }
  else {
    K = K0 * (1+a*Si)
  }
  if(Si > St) {
    Kgs = Kgs0*(1 - 2*b*St + b*St*St/(0.001+Si))
    Kgt = Kgt0*(1 - 2*g*St + g*St*St/(0.001+Si))
  }
  else {
    Kgs = Kgs0 * (1 - b*Si)
    Kgt = Kgt0 * (1 - g*Si)
  }
  Kb = Kb0
  Ka = K * Kb
  AUC = AUCfunc(t, Sdrug)
  #AUC =  t * Sdrug
  Kdrugs = SCALEdrugs * sigmoid(AUC, AUCs5, AUCs95, 0.05, 0.95)
  Kdrugt = SCALEdrugt * sigmoid(AUC, AUCt5, AUCt95, 0.05, 0.95)
  return(list(Kgs = Kgs, Kgt = Kgt, Ka = Ka, Kb = Kb, Kdrugs = Kdrugs, Kdrugt = Kdrugt))
}

calcDerivatives_new = function(t,y,pars) {
  y1 = numeric(length(y))
  Csoi = y[1]
  Cssi = y[2]
  Ctoi = y[3]
  Ctsi = y[4]
  Stress = y[5]
  Csi = Csoi+Cssi
  Cti = Ctoi+Ctsi
  Ctot = sum(y[1:4])
  
  parameters = pars$parameters
  Kmax = pars$Kmax
  
  K0 = parameters$K0  # Ka/Kb
  Kb0 = parameters$Kb0
  Kgs0 = parameters$Kgs0
  Kgt0 = parameters$Kgt0
  Ks = parameters$Ks
  Ksd = parameters$Ksd
  a = parameters$a
  b = parameters$b
  g = parameters$g
  Sdrug = parameters$Sdrug
  
  Si = Stress + Sdrug
  rateConstants = calcRateConstants(t, Si, parameters, Kmax)
  Kgs = rateConstants$Kgs
  Kgt = rateConstants$Kgt
  Ka = rateConstants$Ka
  Kb = rateConstants$Kb
  Kdrugs = rateConstants$Kdrugs
  Kdrugt = rateConstants$Kdrugt
  
  if(is.nan(Kdrugs))
    Kdrugs = 0
  if(is.nan(Kdrugt))
    Kdrugt = 0
  
  Cmax = 1e7
  y1[1] = Kgs*Csoi*(1-Csoi/Cmax)-Ka*Csoi + Kb*Ctsi - Csoi * Kdrugs
  y1[2] = Kgs*Cssi*(1-Cssi/Cmax)-Ka*Cssi + Kb*Ctoi - Cssi * Kdrugs
  y1[3] = Ka*Cssi + Kgt*Ctoi*(1-Ctoi/Cmax)-Kb*Ctoi - Ctoi * Kdrugt #+ Kgt2*Ctoi*Ctoi
  y1[4] = Ka*Csoi + Kgt*Ctsi*(1-Ctsi/Cmax)-Kb*Ctsi - Ctsi * Kdrugt #+ Kgt2*Ctsi*Ctsi
  y1[5] = (Ks*Csi - Ksd*Cti)
  return(list(y1))
}

optimFunc_new = function(parameters_vec, parameter_names, dataset, initial_conditions, fixed_param_list, Kmax = 0.04) {
  t = as.vector(dataset[,1])
  if(sum(!is.na(fixed_param_list)) == 0) {
    parameters = ParamVecToList(parameters_vec, parameter_names)
  } else {
    temp = as.numeric(fixed_param_list)
    variable_ind = which(is.na(temp))
    temp[variable_ind] = parameters_vec[variable_ind]
    parameters = ParamVecToList(temp, parameter_names)
  }
  output = tryCatch({
    output = ode(as.numeric(initial_conditions), t, calcDerivatives_new, list(parameters=parameters, Kmax=Kmax), method = "lsoda")
  }, error = function(e) {
    output = 1e50
  })
  if(length(output) == 1) {
    return(output)
  }
  x1 = output[,2] + output[,5]
  x2 = output[,4] + output[,3]
  #return((sum((dataset[,2]-x1)^2))/(diff(range(dataset[,2]))^2) + (sum((dataset[,3]-x2)^2))/(diff(range(dataset[,3]))^2))
  return(sum((dataset[,2]-x1)^2) + sum((dataset[,3]-x2)^2))
}

optimFunc_all = function(parameters_vec, parameter_names, t, dataset_x1, dataset_x2, fixed_param_list, initial_values, Kmax = 0.04) {
  Cso0_ratio = initial_values$ratios$Cso0
  Css0_ratio = initial_values$ratios$Css0
  Cto0_ratio = initial_values$ratios$Cto0
  Cts0_ratio = initial_values$ratios$Cts0
  S0 = initial_values$S0
  
  nset = ncol(dataset_x1)
  rmsd = 0
  for(i in 1:nset) {
    dataset = cbind(t, dataset_x1[,i], dataset_x2[,i])
    Cs0 = dataset_x1[1,i]
    Ct0 = dataset_x2[1,i]
    S0 = 0.001
    initial_conditions = list(Cso0 = Cso0_ratio*Cs0, Css0 = Css0_ratio*Ct0, Cto0 = Cto0_ratio*Ct0, Cts0 = Cts0_ratio*Cs0, S0 = S0)
    rmsd = rmsd + optimFunc_new(parameters_vec, parameter_names, dataset, initial_conditions, fixed_param_list, Kmax)
  }
  return(rmsd)
}

ParamVecToList = function(parameters_vec, parameter_names) {
  temp = as.list(parameters_vec)
  names(temp) = parameter_names
  return(temp)
}

GAmonitor_new = function(GA) {
  gaMonitor(GA, digits = getOption("digits"))
  if((GA@iter %% 5) == 0) {
    x1_pred = matrix(nrow = ndata, ncol = nset)
    x2_pred = matrix(nrow = ndata, ncol = nset)
    Cso = matrix(nrow = ndata, ncol = nset)
    Css = matrix(nrow = ndata, ncol = nset)
    Cto = matrix(nrow = ndata, ncol = nset)
    Cts = matrix(nrow = ndata, ncol = nset)
    S = matrix(nrow = ndata, ncol = nset)
    x1_err = numeric(nset)
    x2_err = numeric(nset)
    parameters_vec = GA@population[which(GA@fitness == max(GA@fitness))[1],]
    if(sum(!is.na(fixed_param_list)) == 0) {
      parameters_allatonce = ParamVecToList(parameters_vec, param_names)
    } else {
      temp = as.numeric(fixed_param_list)
      variable_ind = which(is.na(temp))
      temp[variable_ind] = parameters_vec[variable_ind]
      parameters_allatonce = ParamVecToList(temp, param_names)
    }
    Cso0_ratio = initial_values$ratios$Cso0
    Css0_ratio = initial_values$ratios$Css0
    Cto0_ratio = initial_values$ratios$Cto0
    Cts0_ratio = initial_values$ratios$Cts0
    S0 = initial_values$S0
    for(i in 1:nset) {
      Cs0 = x1_exp[1,i]
      Ct0 = x2_exp[1,i]
      initial_conditions = list(Cso0 = Cso0_ratio*Cs0, Css0 = Css0_ratio*Ct0, Cto0 = Cto0_ratio*Ct0, Cts0 = Cts0_ratio*Cs0, S0 = S0)
      output = tryCatch({
        output = ode(as.numeric(initial_conditions), t, calcDerivatives_new, list(parameters=parameters_allatonce, Kmax=Kmax), method = "lsoda")
      }, error = function(e) {
        output = 1e50
      })
      if(length(output) > 1) {
        x1_pred[,i] = output[,2] + output[,5]
        x2_pred[,i] = output[,4] + output[,3]
        Cso[,i] = output[,2]
        Css[,i] = output[,3]
        Cto[,i] = output[,4]
        Cts[,i] = output[,5]
        S[,i] = output[,6]
      }
    }
    if(length(output) > 1) {
      par.old = par()
      par(mfcol = c(2,1), mar=c(4,4,1,0.8))
      plot(0,0, xlim = c(0,max(t)+1.5), ylim = sensitives_plot_ylim, xlab = "Time, days", ylab = "Sensitive Cell count")
      for(i in 1:nset) {
        points(t, x1_pred[,i], 'l', col = colors[i], lwd = 2)
        points(t, x1_exp[,i], 'p', col = colors[i])
      }
      legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
      plot(0,0, xlim = c(0,max(t)+1.5), ylim = tolerants_plot_ylim, xlab = "Time, days", ylab = "Tolerant Cell count")
      for(i in 1:nset) {
        points(t, x2_pred[,i], 'l', col = colors[i], lwd = 2)
        points(t, x2_exp[,i], 'p', col = colors[i])
      }
      legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
      par(par.old)
    }
  }
}

plot_properties = function(legend = F) {
  par.old = par()
  par(mfrow = c(3,2), mar=c(4,5,1,0.8))
  plot(0,0, xlim = c(0,max(t)+1), ylim = c(0,2000), xlab = "Time, days", ylab = "Tolerant cell population\n(switched from sensitives)")
  for(i in 1:nset) {
    points(t, Cts[,i], 'l', col = colors[i], lwd = 2)
  }
  if(legend)
    legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
  
  plot(0,0, xlim = c(0,max(t)+1), ylim = c(0,60000), xlab = "Time, days", ylab = "Sensitive cell population\n(switched from tolerants)")
  for(i in 1:nset) {
    points(t, Css[,i], 'l', col = colors[i], lwd = 2)
  }
  if(legend)
    legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
  
  plot(0,0, xlim = c(0,max(t)+1), ylim = c(0,0.15), xlab = "Time, days", ylab = "Fraction of phenotype-\nswitched cells, (from sensitives)")
  for(i in 1:nset) {
    points(t, Cts[,i]/(Cso[,i]+Cts[,i]), 'l', col = colors[i], lwd = 2)
  }
  if(legend)
    legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
  
  plot(0,0, xlim = c(0,max(t)+1), ylim = c(0,1.0), xlab = "Time, days", ylab = "Fraction of phenotype-\nswitched cells, (from tolerants)")
  for(i in 1:nset) {
    points(t, Css[,i]/(Cto[,i]+Css[,i]), 'l', col = colors[i], lwd = 2)
  }
  if(legend)
    legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
  
  plot(0,0, xlim = c(0,max(t)+1), ylim = c(-10,50), xlab = "Time, days", ylab = "Stress")
  for(i in 1:nset) {
    points(t, S[,i], 'l', col = colors[i], lwd = 2)
  }
  if(legend)
    legend("topright", legend = labels, lty = 1, col = colors[1:nset], inset = 0.01)
  par(par.old)
}

simulate_intermittent_therapy = function(Cs0, Ct0, initial_values, parameters_nodrug, parameters_wdrug, therapy_regimen, dt = 0.25, Kmax = 0.04, t_relax = 2) {
  tdrug <<- numeric()
  Sdrug_therapy <<- numeric()
  regime_rules = data.frame(from = numeric(), to = numeric(), parameters = list(), media = numeric(), fraction = numeric())
  nparam = length(parameters_nodrug)
  param_names = names(parameters_nodrug)
  told = 0
  for(i in 1:nrow(therapy_regimen)) {
    t_current = seq(0, therapy_regimen$times[i], dt) + told
    tdrug <<- c(tdrug, t_current)
    Sdrug_therapy <<- c(Sdrug_therapy, rep(therapy_regimen$drug[i], length(t_current)))
    told = max(t_current)
  }
  told = 0  
    # check for change in drug administration
  for(i in 1:nrow(therapy_regimen)) {
    if(i > 1) {
      if(therapy_regimen$drug[i] * therapy_regimen$drug[i-1] == 0) {
        current_param = if(therapy_regimen$drug[i-1] > 0) parameters_wdrug else parameters_nodrug
        current_param$Sdrug = therapy_regimen$drug[i-1]
        new_param = if(therapy_regimen$drug[i] > 0) parameters_wdrug else parameters_nodrug
        new_param$Sdrug = therapy_regimen$drug[i]
        t_transient = seq(0, t_relax, 0.1)
        temp = matrix(nrow = length(t_transient), ncol = nparam)
        for(j in 1:nparam) {
          temp[,j] = sigmoid(t_transient, max(t_transient)*0.05, max(t_transient)*0.95, current_param[[j]], new_param[[j]])
        }
        for(j in 1:(length(t_transient)-1)) {
          regime_rules = rbind(regime_rules, 
                               data.frame(from = t_transient[j]+told, to = t_transient[j+1]+told, 
                                          parameters = ParamVecToList(temp[j,], param_names), 
                                          media = if(j==1) therapy_regimen$media[i] else 0,
                                          fraction = therapy_regimen$fraction[i]))
        }
        told = t_transient[j+1]+told
        regime_rules = rbind(regime_rules, 
                             data.frame(from = told, to = sum(therapy_regimen$times[1:i]), 
                                        parameters = new_param, 
                                        media = therapy_regimen$media[i],
                                        fraction = therapy_regimen$fraction[i]))
        told = sum(therapy_regimen$times[1:i])
      } else {
        current_param = if(therapy_regimen$drug[i] > 0) parameters_wdrug else parameters_nodrug
        current_param$Sdrug = therapy_regimen$drug[i]
        regime_rules = rbind(regime_rules, 
                             data.frame(from = told, to = sum(therapy_regimen$times[1:i]), 
                                        parameters = current_param, 
                                        media = therapy_regimen$media[i],
                                        fraction = therapy_regimen$fraction[i]))
        told = sum(therapy_regimen$times[1:i])
      }
    } else {
      current_param = if(therapy_regimen$drug[i] > 0) parameters_wdrug else parameters_nodrug
      current_param$Sdrug = therapy_regimen$drug[i]
      regime_rules = rbind(regime_rules, 
                           data.frame(from = told, to = sum(therapy_regimen$times[1:i]), 
                                      parameters = current_param, 
                                      media = therapy_regimen$media[i],
                                      fraction = therapy_regimen$fraction[i]))
      told = sum(therapy_regimen$times[1:i])
    }
  }
  #plot(regime_rules$to, regime_rules$parameters.Kgs0, type = 'b')
  #print(regime_rules)
  ratios = initial_values$ratios
  Cso = ratios$Cso0*Cs0
  Css = ratios$Css0*Ct0
  Cto = ratios$Cto0*Ct0
  Cts = ratios$Cts0*Cs0
  S = initial_values$S0
  t = 0
  told = 0
  count = 1
  for(i in 1:nrow(regime_rules)) {
    # set initial values and parameters
    parameters = ParamVecToList(as.numeric(regime_rules[i,grep("parameters",names(regime_rules))]), param_names)
    #print(parameters)
    if(regime_rules$to[i] - regime_rules$from[i] < dt)
      t_current = c(regime_rules$from[i], regime_rules$to[i]) # - regime_rules$from[i]
    else
      t_current = seq(regime_rules$from[i], regime_rules$to[i], dt) # - regime_rules$from[i]
    if(regime_rules$media[i] == 0)
      S0 = S[length(S)]
    else
      S0 = regime_rules$media[i]
    # Check for cell passaging
    if(regime_rules$media[i] > 0) {
      Cs = (Cso[length(S)] + Css[length(S)]) * regime_rules$fraction[i]
      Ct = (Cto[length(S)] + Cts[length(S)]) * regime_rules$fraction[i]
      initial_conditions = list(Cso0 = ratios$Cso0*Cs, Css0 = ratios$Css0*Cs, Cto0 = ratios$Cto0*Ct, Cts0 = ratios$Cts0*Ct, S0 = S0)
    } else {
      initial_conditions = list(Cso0 = Cso[length(S)], Css0 = Css[length(S)], Cto0 = Cto[length(S)], Cts0 = Cts[length(S)], S0 = S0)
    }
    output = ode(as.numeric(initial_conditions), t_current, calcDerivatives_new, list(parameters=parameters, Kmax=Kmax),  method = "lsoda")
    t = c(t, output[,1])
    #t = c(t, t[length(t)] + output[,1])
    Cso = c(Cso, output[,2])
    Css = c(Css, output[,3])
    Cto = c(Cto, output[,4])
    Cts = c(Cts, output[,5])
    S = c(S, output[,6])
    #count = length(t)
    #told = t[count]
  }
  tdrug <<- numeric()
  Sdrug_therapy <<- numeric()
  return(list(t = t,
              Cso = Cso,
              Css = Css,
              Cto = Cto,
              Cts = Cts,
              S = S))
}

filled.contour_new <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)),
            z,
            xlim = range(x, finite=TRUE),
            ylim = range(y, finite=TRUE),
            zlim = range(z, finite=TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20,
            color.palette = cm.colors,
            col = color.palette(length(levels) - 1),
            plot.title, plot.axes, key.title, key.axes,
            asp = NA, xaxs = "i", yaxs = "i", las = 1, axes = TRUE,
            frame.plot = axes, ...)
  {
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    
    mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
    on.exit(par(par.orig))
    
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(las = las)
    
    ## Plot the 'plot key' (scale):
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 0.2
    mar[4L] <- mar[4L]+1.4
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0,1), ylim = range(levels), xaxs = "i", yaxs = "i")
    levels_colorbar = seq(min(levels), max(levels), diff(range(levels))/(length(levels)-1))
    rect(0, levels_colorbar[-length(levels)], 1, levels_colorbar[-1L], col = col)
    if (missing(key.axes)) {
      if (axes)
        axis(4,at=levels_colorbar,labels=levels) #axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title))
      key.title
    
    ## Plot contour-image::
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    
    .filled.contour(x, y, z, levels, col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) box()
    if (missing(plot.title))
      title(...)
    else
      plot.title
    invisible()
  }
