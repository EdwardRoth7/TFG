#Same as td_eval in Funcs.R, modified for convenience.
H_U_factored.opt <- readRDS("H_U_factored.opt.RDS")

ram_solve <- function(ctax, he, pr_treat, pr_vacc, pi1, pi2, pi3, eps, pidbar, pir, kappa, phi, theta, A, beta, maxit, h, tol, noisy, H_U){
  
  H = length(ctax)
  unknowns = list('ns', 'ni', 'nr')
  targets = list('R1', 'R2', 'R3')
  
  ss = initial_ss(A, beta, theta)

  Usinter = rep(list(integer(length = H)), length(unknowns))
  Us <- vector(mode = "list", length = length(Usinter))
  names(Us) <- unknowns
  
  for(i in unknowns){for (j in 1:length(Usinter)){Us[[j]] = Usinter[[j]] + ss[[i]]}}
  
  Uvec = matrix(pack_vectors(Us, unknowns, H), ncol = 1, byrow = T)
  
  for (it in 1:maxit){

    results_ram = td_eval(ns=Us[['ns']],ni=Us[['ni']],nr=Us[['nr']], U_ss=ss[['U']], ctax=ctax, he=he, pi1=pi1, pi2=pi2, pi3=pi3, 
                          H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi, theta=theta, A=A, beta=beta, c_ss=ss[['C']], 
                          n_ss=ss[['N']], kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)
    
    for(k in targets){errors = max(abs(results_ram[[k]]))}
    
    Hvec = pack_vectors(results_ram, targets, H) 
    Uvec = Uvec - factored_solve(H_U, Hvec)
    Us <- unpack_vectors(Uvec, unknowns, H)
  }
  
  return(results_ram)
}

#Planner's social utility function

planner <- function(ctax, he, s0, i0, r0, pr_treat, pr_vacc, pi1, pi2, pi3, eps, pidbar, pir, kappa, phi, theta, 
                    A, beta, maxit, h, tol, noisy, H_U){
  
  #Objective function.
  
  # solve transition path for given guess
  out = ram_solve(ctax, he, pr_treat=pr_treat, pr_vacc=pr_vacc, pi1=pi1, pi2=pi2, pi3=pi3, eps=eps, pidbar=pidbar, pir=pir,
                  kappa=kappa, phi=phi, theta=theta, A=A, beta=beta, maxit=maxit, h=h, tol=tol, noisy=noisy, H_U=H_U)
  
  # welfare
  W = s0 * out[['S']][1] * out[['Us']][1] + i0 * out[['I']][1] * out[['Ui']][1] + r0 * out[['R']][1] * out[['Ur']][1]
  
  return(-W)
  
} 

#Planner's Jacobian matrix

planner_jac <- function(ctax, he, s0, i0, r0, pr_treat, pr_vacc, pi1, pi2, pi3, eps, pidbar, pir, kappa, phi, theta, 
                        A, beta, maxit, h, tol, noisy, H_U){
  
  
  H = length(ctax)
  H = length(ctax)
  dW_ctax = integer(H)
  dW_he = integer(H)
  
  W0 = planner(ctax=ctax, he=he, s0=s0, i0=i0, r0=r0, pr_treat=pr_treat, pr_vacc=pr_vacc, pi1=pi1, pi2=pi2, pi3=pi3,
               eps=eps, pidbar=pidbar, pir=pir, kappa=kappa, phi=phi, theta=theta, A=A, beta=beta, maxit=maxit, h=h,
               tol=tol, noisy=noisy, H_U=H_U_factored.opt)
  
  for (t in 1:H){
    W1_ctax = planner(ctax=ctax + h*(1:H == t), he=he, s0=s0, i0=i0, r0=r0, pr_treat=pr_treat,
                      pr_vacc=pr_vacc, pi1=pi1, pi2=pi2, pi3=pi3, eps=eps, pidbar=pidbar, pir=pir, 
                      kappa=kappa, phi=phi, theta=theta, A=A, beta=beta, maxit=maxit, h=h, tol=tol, 
                      noisy=noisy, H_U=H_U_factored.opt)
    dW_ctax[t] = (W1_ctax - W0) / h
    
    W1_he = planner(ctax=ctax, he=he + h*(1:H == t), s0=s0, i0=i0, r0=r0, pr_treat=pr_treat,
                    pr_vacc=pr_vacc, pi1=pi1, pi2=pi2, pi3=pi3, eps=eps, pidbar=pidbar, pir=pir, 
                    kappa=kappa, phi=phi, theta=theta, A=A, beta=beta, maxit=maxit, h=h, tol=tol, 
                    noisy=noisy, H_U=H_U_factored.opt)
    dW_he[t] = (W1_he - W0) / h
  }
  
  return(c(dW_ctax, dW_he))
}

#Solves Ramsey's problem, takes a while.

ramsey <- function(ctax0, he0, s0, i0, r0, pr_treat, pr_vacc, pi1, pi2, pi3, eps, pidbar, pir, kappa, phi, theta, A, 
                   beta, maxit, h, tol, tol_ramsey, noisy){ 
  
  
  H = length(ctax0)
  unknowns = list('ns', 'ni', 'nr')
  targets = list('R1', 'R2', 'R3')
  
  ss = initial_ss(A, beta, theta)
  
  print('Precomputing Jacobian...')
  J = get_J(ss=ss, ctax=ctax0, he=he0, pi1=pi1, pi2=pi2, pi3=pi3, H=H, eps=eps, pidbar=pidbar, pir=pir, phi=phi,
            theta=theta, A=A, beta=beta, h=h, kappa=kappa, pr_treat=pr_treat, pr_vacc=pr_vacc)
  H_U_factored = J_to_HU(J, H, unknowns, targets)
  print('Done!')
  
  
  obj = function(pars) {
    ctax = pars[1:H]
    he = pars[(H+1):(2*H)] 
    return(planner(ctax, he, s0=s0, i0=i0, r0=r0, pr_treat=pr_treat, pr_vacc=pr_vacc, pi1=pi1, pi2=pi2, pi3=pi3,
                   eps=eps, pidbar=pidbar, pir=pir, kappa=kappa, phi=phi, theta=theta, A=A, beta=beta,
                   maxit=maxit, h=h, tol=tol, noisy=noisy, H_U=H_U_factored.opt))
  }
  
  jac = function(pars) {
    ctax = pars[1:H]
    he = pars[(H+1):(2*H)] 
    return(planner_jac(ctax, he, s0=s0, i0=i0, r0=r0, pr_treat=pr_treat, pr_vacc=pr_vacc, pi1=pi1, pi2=pi2, pi3=pi3,
                       eps=eps, pidbar=pidbar, pir=pir, kappa=kappa, phi=phi, theta=theta, A=A, beta=beta,
                       maxit=maxit, h=h, tol=tol, noisy=noisy, H_U=H_U_factored.opt))
  }
  
  lower_bounds = c(rep(0, 2*H))  # Límite inferior de 0 para ctax y he
  upper_bounds = c(rep(1, H), rep(Inf, H))  # Límite superior de 1 para ctax y Inf para he
  
  res = optim(par = c(ctax0, he0), fn = obj, gr = jac, method='L-BFGS-B', lower = lower_bounds, upper = upper_bounds, control = list("ndeps" = tol_ramsey, "reltol" = tol))
  
  #res = optim(par = c(ctax0, he0), fn = obj, gr = jac, method='BFGS', control = list("ndeps" = tol_ramsey, "reltol" = tol))
  
  return(res)
  
}

# Resuelve el problema de Ramsey y guarda el resultado en resopt
#resopt <- ramsey(ctax0=integer(H), he0=integer(H), s0=1, i0=1, r0=1, pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM,
                 #eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
                 #A=A, beta=beta, maxit=100, h=1E-4, tol=1E-8, tol_ramsey=1E-3, noisy=FALSE) # Ramsey, SIR Macro

#ctax_opt = resopt$par[1:H]
#he_opt = resopt$par[(H+1):(2*H)]

