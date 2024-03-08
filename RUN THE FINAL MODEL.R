
library(tidyverse)
library(Matrix)
library(stats)
library(ggpubr)
library(foreach)
library(doParallel)

source("https://raw.githubusercontent.com/EdwardRoth7/TFG/main/Funcs%20he%20optimo.R?token=GHSAT0AAAAAACPITFZG4ITTHSBERHBH3AR4ZPK46QQ")
#source('~/Desktop/calibratepis (3).R')
#source('~/Desktop/Funcs (3).R')
#source('~/Desktop/Ramsey (3).R')
#source('~/Desktop/figures (3).R')

### KEY PARAMETERS #############################################################

H = 250

A = 39.8351648 
beta = 0.96^(1/52) 
theta = 0.00127551
eps = 0.001 
pid = 7*0.005/18
pir = 7/18 - 7*0.005/18 
phi = 0.8

################################################################################

### CALIBRATE PI'S #############################################################

scale1 = 1000000
scale2 = 1000

SIRMacpis <- calib(pis_guess= c(0.2,0.2,0.2),H = H,eps = eps,pir = pir, 
                   pid = pid,pis1_shr_target=1/6,pis2_shr_target=1/6,RplusD_target=0.6,
                   C = initial_ss(A=A,beta=beta,theta=theta)[["C"]], 
                   N = initial_ss(A=A,beta=beta,theta=theta)[["N"]],
                   scale1=scale1,scale2=scale2) #SIR Macro

pi1SM <- SIRMacpis$par[1]/scale1
pi2SM <- SIRMacpis$par[2]/scale2
pi3SM <- SIRMacpis$par[3]

SIRMacpis <- calib(pis_guess= c(0.2,0.2,0.2),H = H,eps = eps,pir = pir, 
                   pid = pid,pis1_shr_target=0,pis2_shr_target=0,RplusD_target=0.6,
                   C = initial_ss(A=A,beta=beta,theta=theta)[["C"]], 
                   N = initial_ss(A=A,beta=beta,theta=theta)[["N"]],
                   scale1=scale1,scale2=scale2) #SIR, no Macro

pi3S <- SIRMacpis$par[3]

################################################################################

### FIGURES ####################################################################

ssgraph <- initial_ss(A=A, beta=beta, theta=theta) #Steady State

td1 <- td_sir(pi3=pi3S, H, eps=eps, pid=pid, pir=pir, phi=phi, theta=theta, A=A) # SIR, no Macro

td3.opt <- td_solve(integer(H), integer(H), pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                    eps=eps, pidbar=pid, pir=pir, kappa=0, phi=phi, theta=theta, 
                    A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro, sin confinamiento y sin gasto en salud

saveRDS(td3.opt[["H_U_factored"]], "H_U_factored.opt.RDS")

# Resuelve el problema de Ramsey y guarda el resultado en resopt
resopt <- ramsey(ctax0=integer(H), he0=integer(H), s0=1, i0=1, r0=1, pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM,
                 eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
                 A=A, beta=beta, maxit=100, h=1E-4, tol=1E-8, tol_ramsey=1E-3, noisy=FALSE) # Ramsey, SIR Macro

ctax_opt = resopt$par[1:H]
he_opt = resopt$par[(H+1):(2*H)]


td3.opt <- td_solve(ctax_opt, he_opt, pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                    eps=eps, pidbar=pid, pir=pir, kappa=0, phi=phi, theta=theta, 
                    A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro, sin confinamiento y sin gasto en salud





