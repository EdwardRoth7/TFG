
#Steady state

initial_ss <- function(A, beta, theta){
  with(as.list(c(A, beta, theta)),{
    w = A
    b_ss = 1100
    taxn_ss = ((((1/beta)-1)*(theta^(1/2)))/A)*b_ss
    gama = 0.001
    N = ((1 - taxn_ss)/theta)^(1/2)
    C = A * N
    u = (log(C) - theta/2 * N^2)
    U = 1/(1 - beta) * u
    b_ss = 1100
    taxn_ss = ((((1/beta)-1)*(theta^(1/2)))/A)*b_ss
    gama = -0.001
    
    
    return(list('w'= w, 'C'= C, 'N'= N, 'ns'= N, 'ni'= N, 'nr'= N, 'U'= U, 'A'= A, 'beta'= beta, 'theta'= theta, 'Neff'= N, 'b_ss' = b_ss, 'taxn_ss' = taxn_ss ))
  }
  )
}

#SIR without Macro

td_sir <- function(pi3, H, eps, pid, pir, phi, theta, A){
  with(as.list(c(pi3, H, eps, pid, pir, phi, theta, A)),{
    S = integer(H)
    I = integer(H)
    R = integer(H)
    D = integer(H)
    
    S[1] = 1 - eps
    I[1] = eps
    R[1] = 0
    D[1] = 0
    
    Tr = integer(H)
    
    for(t in 2:H){
      Tr[t-1] = pi3 * S[t-1] * I[t-1]
      
      S[t] = S[t-1] - Tr[t-1]
      I[t] = I[t-1] + Tr[t-1] - (pir + pid) * I[t-1]
      R[t] = R[t-1] + pir * I[t-1]
      D[t] = D[t-1] + pid * I[t-1]
    }
    
    P = integer(H) + 1
    P = P - D
    
    n = theta^(-1/2)
    cs = A * n
    cr = A * n
    ci = phi * A * n
    
    N = n * (S + phi*I + R)
    C = S * cs + I * ci + R * cr
    
    return(list('S'= S, 'I'= I, 'R'= R, 'D'= D, 'P'= P, 'N'= N, 'C'= C))
  }
  )
}

#Evaluates the model (ERT(2020), Appendix A)

td_eval <- function(ns, ni, nr, ctax, he, U_ss, c_ss, n_ss, pr_treat, pr_vacc, pi1, pi2, pi3, H, eps, pidbar, pir, phi, theta, A, beta, kappa){
  with(as.list(c(ns, ni, nr, ctax, he, U_ss, c_ss, n_ss, pr_treat, pr_vacc, pi1, pi2, pi3, H, eps, pidbar, pir, phi, theta, A, beta, kappa)),{
    

    b_ss = 1100
    taxn_ss = ((((1/beta)-1)*(theta^(1/2)))/A)*b_ss
    gama = -0.001
    
    S = integer(H)
    I = integer(H)
    R = integer(H)
    D = integer(H)
    pid = integer(H)
    b = integer(H)
    pic = integer(H)
    pin = integer(H)
    taxn = integer(H)
    dp = integer(H)
    cr = integer(H)
    transfer = integer(H)
    ci = integer(H)
    cs = integer(H)
    comprobacion1 = integer(H)
    
    S[1] = 1 - eps
    I[1] = eps
    R[1] = 0
    D[1] = 0
    #pic[1] = pi1 - (1.5*10^(-9)) * he[1]
    #pin[1] = pi2 - (7.5*10^(-7)) * he[1]
    pic[1] = pi1 * (1.0045 + (0.1 - 0.9)/(1 + 0.5 * exp(- 4 * (he[1]/20 - 0.8)))^2)
    pin[1] = pi2 * (1.0045 + (0.1 - 0.9)/(1 + 0.5 * exp(- 4 * (he[1]/20 - 0.8)))^2)
    b[1] = b_ss + (he[1] + ((1/beta) - 1) * b_ss - taxn_ss * A * (S[1] * ns[1] + I[1] * ni[1] * phi + R[1] * nr[1]))/(1 + gama * A * (S[1] * ns[1] + I[1] * ni[1] * phi + R[1] * nr[1]))
    taxn[1] = taxn_ss + gama * (b[1] - b_ss)
    #dp[1] = he[1] + ((1/beta) - 1 ) * b_ss - taxn[1] * A * (S[1] * ns[1] + I[1] * ni[1] * phi + R[1] * nr[1])
    dp[1] = b[1] - b_ss
    cr[1] = (A * (1 - taxn[1]))/ ((1 + ctax[1]) * theta * nr[1])
    transfer[1] = (1 + ctax[1]) * cr[1] - A * (1 - taxn[1]) * nr[1] + b[1] - (1/beta) * b_ss 
    ci[1] = (A * phi * (1 - taxn[1]) * ni[1] + transfer[1] + (1/beta) * b_ss - b[1]) / (1 + ctax[1])
    cs[1] = (A * (1 - taxn[1]) * ns[1] + transfer[1] + (1/beta) * b_ss - b[1]) / (1 + ctax[1])
    comprobacion1[1] = -dp[1] + ((1/beta)-1) * b_ss - taxn[1] * A * (S[1] * ns[1] + I[1] * ni[1] * phi + R[1] * nr[1])
    
    Tr = integer(H)
    
    for (t in 2:H){
      
      Tr[t-1] = pic[t-1] * S[t-1] * (A * ns[t-1] - he[t-1]) * I[t-1] * (A * phi * ni[t-1] - he[t-1]) + pin[t-1] * S[t-1] * ns[t-1] * I[t-1] * ni[t-1] + pi3 * S[t-1] * I[t-1]
      
      pid[t-1] = pidbar + kappa * I[t-1]^2
      S[t] = S[t-1] - Tr[t-1]
      I[t] = I[t-1] + Tr[t-1] - (pir + pid[t-1]) * I[t-1]
      R[t] = R[t-1] + pir * I[t-1]
      D[t] = D[t-1] + pid[t-1] * I[t-1]
      #pic[t] = pi1 - (1.5*10^(-9)) * he[t]
      #pin[t] = pi2 - (7*10^(-7)) * he[t]
      pic[t] = pi1 * (1.0045 + (0.1 - 0.9)/(1 + 0.5 * exp(- 4 * (he[t]/20 - 0.8)))^2)
      pin[t] = pi2 * (1.0045 + (0.1 - 0.9)/(1 + 0.5 * exp(- 4 * (he[t]/20 - 0.8)))^2)
      b[t] = (he[t] + (1/beta) * b[t-1] + (gama * b_ss - taxn_ss) * A * (S[t] * ns[t] + I[t] * ni[t] * phi + R[t] * nr[t]))/(1 + gama * A * (S[t] * ns[t] + I[t] * ni[t] * phi + R[t] * nr[t]))
      taxn[t] = taxn_ss + gama * (b[t] - b_ss)
      #dp[t] = he[t] + ((1/beta) - 1 ) * b[t-1] - taxn[t] * A * (S[t] * ns[t] + I[t] * ni[t] * phi + R[t] * nr[t]) 
      dp[t] = b[t] - b[t-1]
      cr[t] = (A * (1 - taxn[t]))/ ((1 + ctax[t]) * theta * nr[t])
      transfer[t] = (1 + ctax[t]) * cr[t] - A * (1 - taxn[t]) * nr[t] + b[t] - (1/beta) * b[t-1] 
      ci[t] = (A * phi * (1 - taxn[t]) * ni[t] + transfer[t] + (1/beta) * b[t-1] - b[t]) / (1 + ctax[t])
      cs[t] = (A * (1 - taxn[t]) * ns[t] + transfer[t] + (1/beta) * b[t-1] - b[t]) / (1 + ctax[t])
      comprobacion1[t] = -dp[t] + ((1/beta)-1) * b[t-1] - taxn[t] * A * (S[t] * ns[t] + I[t] * ni[t] * phi + R[t] * nr[t])
    }
    
    #cr = (A * (1 - taxn)) / ((1 + ctax) * theta * nr)
    #transfer = ctax * cr 
    #ci = (A * phi * (1 - taxn)) / ((1 + ctax) * theta * ni)
    #cs = (A * (1 - taxn)) / ((1 + ctax) * theta * ns)
    
    
    P = integer(H) + 1
    P = P - D
    
    pid[length(pid)] = pidbar + kappa * I[length(I)]^2
    
    tau = pic * cs * I * ci + pin * ns * I * ni + pi3 * I
    
    Ur = integer(H+1)
    Ui = integer(H+1)
    Us = integer(H+1)
    
    Ur[length(Ur)] = U_ss
    Us[length(Us)] = U_ss
    Ui[length(Ui)] = (log(phi * c_ss) - theta / 2 * n_ss^2 + (1 - pr_treat[length(pr_treat)]) * beta * pir * Ur[length(Ur)] + beta * pr_treat[length(pr_treat)] * Ur[length(Ur)]) / (1 - beta * (1 - pir - pid[length(pid)]))
    
    for (t in H:1){
      
      Ur[t] = log(cr[t]) - theta / 2 * nr[t]^2 + beta * Ur[t+1]
      Ui[t] = log(ci[t]) - theta / 2 * ni[t]^2 + (1 - pr_treat[t]) * beta * ((1 - pir - pid[t]) * Ui[t+1] + pir * Ur[t+1]) + beta * pr_treat[t] * Ur[t+1]
      Us[t] = log(cs[t]) - theta / 2 * ns[t]^2 + (1 - pr_vacc[t]) * beta * ((1 - tau[t]) * Us[t+1] + tau[t] * Ui[t+1]) + pr_vacc[t] * beta * Ur[t+1]
      
    }
    
    mus = beta * (1 - pr_vacc) * (Us[2:length(Us)] - Ui[2:length(Ui)])    #mus = landa tau
    lams = (theta * ns + mus * pin * I * ni) / (A * (1 - taxn))           #lams = landa de R.P. de susceptibles
    lami = theta * ni / (phi * A * (1 - taxn))                            #lami = landa de R.P. de infectados
    
    R1 = ctax * (S * cs + I * ci + R * cr) - transfer * (S + I + R)
    R2 = lami * (1 + ctax) - 1 / ci
    R3 = lams * (1 + ctax) + mus * pic * I * ci - 1 / cs
    
    C = S * cs + I * ci + R * cr
    N = S * ns + I * ni * phi + R * nr
    Neff = S * ns + I * phi * ni + R * nr  # redundant
    walras = A * Neff - C - he
    mortality = pid / pir
    
    
    return(list('ns'= ns, 'ni'= ni, 'nr'= nr, 'cs'= cs, 'ci'= ci, 'cr'= cr, 'transfer'= transfer, 'Tr'= Tr, 'S'= S, 'I'= I,
                'R'= R, 'D'= D, 'tau'= tau, 'Ur'= Ur[1:(length(Ur)-1)], 'Ui'= Ui[1:(length(Ui)-1)], 'Us'= Us[1:(length(Us)-1)], 'mus'= mus, 'lami'= lami, 'P'= P,
                'lams'= lams, 'R1'= R1, 'R2'= R2, 'R3'= R3, 'C'= C, 'N'= N, 'Neff'= Neff, 'walras'= walras, 'pid'= pid,
                'pr_treat'= pr_treat, 'pr_vacc'= pr_vacc, 'mortality'= mortality, 'ctax'= ctax,
                'kappa'= kappa, 'pir'= pir, 'beta'= beta, 'pidbar'= pidbar, 'theta'= theta, 'A'= A, 'eps'= eps, 'pi1'= pi1,
                'pi2'= pi2, 'pi3'= pi3, 'pic'= pic, 'pin'= pin, 'b'= b, 'b_ss'= b_ss, 'he'= he,'dp'= dp, 'taxn'= taxn, 'taxn_ss'= taxn_ss, 'comprobacion1' = comprobacion1))  }
  )
}

#Compute the Jacobian matrix (gradient-based method of optimization)... in parallel

get_J <- function(ss, ctax, he, pr_treat, pr_vacc, pi1, pi2, pi3, H, eps, pidbar, pir, phi, theta, A, beta, kappa, h){
  
  td0 <- td_eval(ns=ss[['N']] * (integer(H)+1), ni=ss[['N']] * (integer(H)+1), nr=ss[['N']] * (integer(H)+1), pr_treat=pr_treat, pr_vacc=pr_vacc, U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, kappa=kappa, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']])
  
  n.cores <- parallel::detectCores() - 1
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  doParallel::registerDoParallel(cl = my.cluster)
  
  multiResultClass <- function(result1=NULL,result2=NULL,result3=NULL,result4=NULL,result5=NULL,result6=NULL,result7=NULL,result8=NULL,result9=NULL)
  {
    me <- list(
      result1 = result1,
      result2 = result2,
      result3 = result3,
      result4 = result4,
      result5 = result5,
      result6 = result6,
      result7 = result7,
      result8 = result8,
      result9 = result9
    )
    
    ## Set the name for the class
    class(me) <- append(class(me),"multiResultClass")
    return(me)
  }

  # Compute via direct method
  var <- foreach(t = 1:H, .combine= 'c', .export = "td_eval") %dopar% {
    
    result <- multiResultClass()
    
    result$result1 <- (td_eval(ns=ss[['N']] + h * (1:H == t), ni=ss[['N']] * (integer(H) + 1), nr=ss[['N']] * (integer(H) + 1), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R1']] - td0[['R1']])/h
    result$result2 <- (td_eval(ns=ss[['N']] + h * (1:H == t), ni=ss[['N']] * (integer(H) + 1), nr=ss[['N']] * (integer(H) + 1), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R2']] - td0[['R2']])/h
    result$result3 <- (td_eval(ns=ss[['N']] + h * (1:H == t), ni=ss[['N']] * (integer(H) + 1), nr=ss[['N']] * (integer(H) + 1), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R3']] - td0[['R3']])/h
    
    result$result4 <- (td_eval(ns=ss[['N']] * (integer(H) + 1), ni=ss[['N']] + h * (1:H == t), nr=ss[['N']] * (integer(H) + 1), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R1']] - td0[['R1']])/h
    result$result5 <- (td_eval(ns=ss[['N']] * (integer(H) + 1), ni=ss[['N']] + h * (1:H == t), nr=ss[['N']] * (integer(H) + 1), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R2']] - td0[['R2']])/h
    result$result6 <- (td_eval(ns=ss[['N']] * (integer(H) + 1), ni=ss[['N']] + h * (1:H == t), nr=ss[['N']] * (integer(H) + 1), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R3']] - td0[['R3']])/h
    
    result$result7 <- (td_eval(ns=ss[['N']] * (integer(H) + 1), ni=ss[['N']] * (integer(H) + 1), nr=ss[['N']] + h * (1:H == t), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R1']] - td0[['R1']])/h
    result$result8 <- (td_eval(ns=ss[['N']] * (integer(H) + 1), ni=ss[['N']] * (integer(H) + 1), nr=ss[['N']] + h * (1:H == t), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R2']] - td0[['R2']])/h
    result$result9 <- (td_eval(ns=ss[['N']] * (integer(H) + 1), ni=ss[['N']] * (integer(H) + 1), nr=ss[['N']] + h * (1:H == t), U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)[['R3']] - td0[['R3']])/h
    
    return(result)
  }
  
  x = rbind(matrix(unlist(var[seq(1, length(1:(H*H)), 9)]), H, H), matrix(unlist(var[seq(2, length(1:(H*H)), 9)]), H, H), matrix(unlist(var[seq(3, length(1:(H*H)), 9)]), H, H))
  y = rbind(matrix(unlist(var[seq(4, length(1:(H*H)), 9)]), H, H), matrix(unlist(var[seq(5, length(1:(H*H)), 9)]), H, H), matrix(unlist(var[seq(6, length(1:(H*H)), 9)]), H, H))
  z = rbind(matrix(unlist(var[seq(7, length(1:(H*H)), 9)]), H, H), matrix(unlist(var[seq(8, length(1:(H*H)), 9)]), H, H), matrix(unlist(var[seq(9, length(1:(H*H)), 9)]), H, H))
  
  J = cbind(x, y, z)
  
  stopCluster(my.cluster)
  
  return(J)
}


################################################################################
#UTILS
################################################################################

isEmpty <- function(x){
  if (length(x)==0){
    x= 0
  }
  return(x)
}


pack_vectors <- function(vs, names, Tr){
  
  v <- integer(length(names)*Tr)
  
  for (i in 1:length(names)){
    if (i == 1){
      v[1:Tr] <- vs[[ names[[1]] ]]
    }
    else{
      v[((i-1)*Tr+1):(i*Tr)] <- vs[[ names[[i]] ]]
    }
  }
  
  return(v)
}

unpack_vectors <- function(v, names, Tr){
  
  vs <- vector(mode = "list", length = length(names))
  names(vs) <- names
  for (name in names){
    vs[[name]] <- rep(NA, Tr)
  }
  
  i = 1
  for (name in names){
    if (i == 1){
      vs[[name]] = v[1:Tr]
    }else{
      vs[[name]] = v[((i-1)*Tr+1):(i*Tr)]
    }
    i <- i + 1
  }
  return(vs)
}


factor <- function(X){
  
  ex <- Matrix::expand(lu(X))
  P <- ex$P
  L <- ex$L
  LU <- U <- ex$U
  LU[lower.tri(U)] <- L[lower.tri(L)]
  
  
  Z = list(LU, P, L, U)
  
  return(Z)
}


factored_solve <- function(Z, b){
  
  (Y <- solve(Z[[3]], solve(Z[[2]]) %*% b))
  (X <- solve(Z[[4]], Y))
  
  return(X)
}


J_to_HU <- function(J, H, unknowns, targets){
  H_U_factored = factor(J)
  return(H_U_factored)
}

################################################################################
### SOLUTION
################################################################################

#Newton's method, iterate until convergence. System of eq. solved using LU factorization.

td_solve <- function(ctax, he, pr_treat, pr_vacc, pi1, pi2, pi3, eps, pidbar, pir, kappa, phi, 
                     theta, A, beta, maxit, h, tol, noisy, H_U){
  
  H = length(ctax)
  unknowns = list('ns', 'ni', 'nr')
  targets = list('R1', 'R2', 'R3')
  
  ss = initial_ss(A, beta, theta)
  
  print('Precomputing Jacobian...')
  J = get_J(ss=ss, ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, h=h, kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)
  H_U_factored = J_to_HU(J, H, unknowns, targets)
  print('Done!')
  
  Usinter = rep(list(integer(length = H)), length(unknowns))
  Us <- vector(mode = "list", length = length(Usinter))
  names(Us) <- unknowns
  
  for(i in unknowns){for (j in 1:length(Usinter)){Us[[j]] = Usinter[[j]] + ss[[i]]}}
  
  Uvec = matrix(pack_vectors(Us, unknowns, H), ncol = 1, byrow = T)
  
  for (it in 1:maxit){
    
    results = td_eval(ns=Us[['ns']],ni=Us[['ni']],nr=Us[['nr']], U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)
    
    for(k in targets){errors = max(abs(results[[k]]))}
    
    Hvec = pack_vectors(results, targets, H) 
    Uvec = Uvec - factored_solve(H_U_factored, Hvec)
    Us <- unpack_vectors(Uvec, unknowns, H)
    
    # Imprimir el progreso de la simulación
    cat("\rProgreso de la simulación: ", round(it/maxit*100, 2), "%", sep="")
    
  }
  
  return(list('results'=results, 'H_U_factored' = H_U_factored))
}


################################################################################

################################################################################


#he_arbit <- readRDS("he.arbit.rds")
#he.arbit = he_arbit[1:H]

#td3.opt <- td_solve(integer(H), integer(H), pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                  #eps=eps, pidbar=pid, pir=pir, kappa=0, phi=phi, theta=theta, 
                  #A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro, Best Containment Policy

#saveRDS(td3.opt[["H_U_factored"]], "H_U_factored.opt.RDS")

