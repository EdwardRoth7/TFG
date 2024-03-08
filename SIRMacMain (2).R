
library(tidyverse)
library(Matrix)
library(stats)
library(ggpubr)
library(foreach)
library(doParallel)

source('~/Desktop/calibratepis (3).R')
source('~/Desktop/Funcs (3).R')
source('~/Desktop/Ramsey (3).R')
source('~/Desktop/figures (3).R')

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

td2 <- td_solve(integer(H), pr_treat=integer(H), pr_vacc=integer(H),pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, eps=eps, 
                pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, A=A, 
                beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro

taxload <- ramsey(ctax0=integer(H), s0=1, i0=1, r0=1, pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM,
                  eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
                  A=A, beta=beta, maxit=100, h=1E-4, tol=1E-8, tol_ramsey=1E-3, noisy=FALSE) # Ramsey, SIR Macro


#taxload = readRDS(file = "/Users/edward/Downloads/ctax.rds") #Solves Ramsey's problem. Any change in the parameters must be 
taxnormal = taxload$par[1:H] #followed by a re-running of ramsey() function (H + 2 Jacobians + optim()/base/ optimization)

td3 <- td_solve(tax_prueba, pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro, Best Containment Policy

#taxmedpreload <- ramsey(ctax0=taxload$par, s0=1, i0=1, r0=1, pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, 
#                        pi3=pi3SM, eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
#                        A=A, beta=beta, maxit=100, h=1E-4, tol=1E-8, tol_ramsey=1E-3, noisy=FALSE) # Ramsey, SIR Macro w/ endogenous mortality

taxmedpreload = readRDS(file = "/Users/edward/Downloads/ctaxmedpre (1).rds") #Solves Ramsey's problem. Any change in the parameters must be 
taxm = taxmedpreload[1:H] #followed by a re-running of ramsey() function (H + 2 Jacobians + optim()/base/ optimization)

td4 <- td_solve(integer(H), pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro w/ endogenous mortality

td5 <- td_solve(taxm, pr_treat=integer(H), pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # Best Containment Policy w/ endogenous mortality

td6 <- td_solve(integer(H), pr_treat=integer(H) + 1/52, pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro w/ treatment

#taxtreat <- ramsey(ctax0=taxnormal, s0=1, i0=1, r0=1, pr_treat=integer(H) + 1/52, pr_vacc=integer(H), pi1=pi1SM, pi2=pi2SM, 
#                   pi3=pi3SM, eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
#                   A=A, beta=beta, maxit=100, h=1E-4, tol=1E-8, tol_ramsey=1E-3, noisy=FALSE) # Ramsey, SIR Macro w/ treatment

taxtreatload = readRDS(file = "C:/.../taxtreat.rds") #Solves Ramsey's problem. Any change in the parameters must be
taxtreatment = taxtreatload[1:H] #followed by a re-running of ramsey() function (H + 2 Jacobians + optim()/base/ optimization)

#taxvacc <- ramsey(ctax0=taxtreatment, s0=1, i0=1, r0=1, pr_treat=integer(H), pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, 
#                   pi3=pi3SM, eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
#                   A=A, beta=beta, maxit=100, h=1E-4, tol=1E-8, tol_ramsey=1E-3, noisy=FALSE) # Ramsey, SIR Macro w/ vaccines

taxvaccload = readRDS(file = "C:/.../taxvacc.rds")
taxvacc = taxvaccload[1:H]

td7 <- td_solve(integer(H), pr_treat=integer(H), pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # SIR Macro w/ vaccines


td8 <- td_solve(taxvacc, pr_treat=integer(H), pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0.0, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # Best Containment Policy w/ vaccine

#taxbench <- ramsey(ctax0=taxvacc, s0=1, i0=1, r0=1, pr_treat=integer(H) + 1/52, pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, 
#                   pi3=pi3SM, eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
#                   A=A, beta=beta, maxit=100, h=1E-4, tol=1E-8, tol_ramsey=1E-3, noisy=FALSE)

taxbenchload = readRDS(file = "C:/.../taxbench.rds")
taxbench = taxbenchload[1:H]

td9 <- td_solve(integer(H), pr_treat=integer(H) + 1/52, pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # Benchmark SIR Macro Model


td10 <- td_solve(taxbench, pr_treat=integer(H) + 1/52, pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
                A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # Best Simple Containment Policy (Benchmark)

eexit12 <- c(taxbench[1:12], integer(H-12)) #Early exit
eexit44 <- c(taxbench[1:44], integer(H-44)) #Early exit

td11 <- td_solve(eexit12, pr_treat=integer(H) + 1/52, pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                 eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
                 A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # Early Exit from Best Simple Containment Policy (Benchmark), 12 weeks

td12 <- td_solve(eexit44, pr_treat=integer(H) + 1/52, pr_vacc=integer(H) + 1/52, pi1=pi1SM, pi2=pi2SM, pi3=pi3SM, 
                 eps=eps, pidbar=pid, pir=pir, kappa=0.9, phi=phi, theta=theta, 
                 A=A, beta=beta, maxit=50, h=1E-4, tol=1E-8, noisy=FALSE, H_U=NA) # Early Exit from Best Simple Containment Policy (Benchmark), 44 weeks


figure_one(H = H, ssgraph = ssgraph, t1 = td1, t2 = td2)

figure_two(H = H, ssgraph = ssgraph, t1 = td2)

figure_three(H = H, ssgraph = ssgraph, t1 = td2, t2 = td3, ctax = taxnormal)

figure_four(H = H, ssgraph = ssgraph, t1 = td2, t2 = td4, t3 = td5, ctax = taxm)

figure_five(H = H, ssgraph = ssgraph, t1 = td2, t2 = td6, t3 = td7, ctax = taxtreatment)

figure_six(H = H, ssgraph = ssgraph, t1 = td2, t2 = td7, t3 = td8, ctax = taxvacc)

figure_seven(H = H, ssgraph = ssgraph, t1 = td9, t2 = td10, ctax = taxbench)

figure_eight(H = H, ssgraph = ssgraph, t1 = td10, t2 = td11, t3 = td12, ctax = taxbench)


#Gráfico Deuda
plot(1:H,td3[["results"]][["b"]] , type="l", main="Evolución temporal de la deuda pública per cápita", xlab="Semana", ylab="b")


#Gráfico Déficit Público
plot(1:H, td3[["results"]][["dp"]], type="l", main="Evolución temporal del déficit público per cápita", xlab="Semana", ylab="dp")


#Gráfico Horas Trabajadas
plot(td3[["results"]][["N"]], type = "l", main="Comparación caida horas trabajadas sin déficit/con déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(min(td3N_normal), 28) )
lines(td3N_normal, type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


#Gráfico Consumo
plot(td3[["results"]][["C"]], type = "l", main="Comparación consumo con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(min(td3C_normal), 1115))
lines(td3C_normal, type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


#Gráfico Confinamiento
plot(td3[["results"]][["ctax"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max( taxnormal   )) )
lines(taxnormal, type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


#Gráfico Mortalidad
plot(td3[["results"]][["D"]], type = "l", main="Comparación mortalidad con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max( td3D_normal   )) )
lines(td3D_normal, type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")




############################################################################################################################################

#Comparación logistica con Normal

############################################################################################################################################

#Gráfico Deuda
plot(1:H, td3.1[["b"]], type="l", main="Evolución temporal de la deuda pública per cápita", xlab="Semana", ylab="b", col = "blue")

plot(1:H, td3.1[["walras"]], type="l", main="Evolución temporal del déficit público per cápita", xlab="Semana", ylab="dp", col = "blue")

#Gráfico Déficit Público
plot(1:H, td3.1[["dp"]], type="l", main="Evolución temporal del déficit público per cápita", xlab="Semana", ylab="dp", col = "blue")

plot(1:H, td3.1[["he"]], type="l", main="Evolución temporal del déficit público per cápita", xlab="Semana", ylab="dp", col = "blue")

#Gráfico Horas Trabajadas
plot(td3.1[["N"]], type = "l", main="Comparación caida horas trabajadas sin déficit/con déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(25, 28) )
lines(td2[["N"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


plot(td3.1[["C"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(min(td2[["C"]]), max(td3.1[["C"]])) )
lines(td2[["C"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")



plot(td3.1[["C"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(min(td2[["C"]]), max(td3.1[["C"]])) )
lines(A * td3.1[["N"]] - 0 * td3.1[["he"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")



#Gráfico Mortalidad
plot(td3.1[["pic"]], type = "l", main="Comparación mortalidad con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td3.1[["pic"]])) )
lines(td2[["pi1"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


plot(td3.1[["D"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(Dtd2_normal)) )
lines(td2[["D"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


plot(td3.1[["I"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td2[["I"]])) )
lines(td2[["I"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")

plot(td3.1[["he"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td3.1[["he"]])) )
lines(-td3.1[["comprobacion1"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")

plot(td3.1[["I"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes")


############################################################################################################################################

#Comparación logistica con anterior

############################################################################################################################################


#Gráfico Deuda
plot(1:H, td3.0[["b"]], type="l", main="Evolución temporal de la deuda pública per cápita", xlab="Semana", ylab="b", col = "blue")

plot(1:H, td3.0[["walras"]], type="l", main="Evolución temporal del déficit público per cápita", xlab="Semana", ylab="dp", col = "blue")

#Gráfico Déficit Público
plot(1:H, td3.0[["dp"]], type="l", main="Evolución temporal del déficit público per cápita", xlab="Semana", ylab="dp", col = "blue")

plot(1:H, td3.1[["he"]], type="l", main="Evolución temporal del déficit público per cápita", xlab="Semana", ylab="dp", col = "blue", ylim = c(0, max(td3.0[["he"]])))
lines(td3.0[["he"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


#Gráfico Horas Trabajadas
plot(td3.1[["N"]], type = "l", main="Comparación caida horas trabajadas sin déficit/con déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(27, 28) )
lines(td3.0[["N"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")

plot(td3.1[["C"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(min(td3.0[["C"]]), max(td3.1[["C"]])) )
lines(td3.0[["C"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


plot(td3.0[["C"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(min(td2[["C"]]), max(td3.1[["C"]])) )
lines(A * td3.0[["N"]] - 0 * td3.0[["he"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")



#Gráfico Mortalidad
plot(td3.1[["pic"]], type = "l", main="Comparación mortalidad con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td3.1[["pic"]])) )
lines(td3.0[["pic"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


plot(td3.1[["pin"]], type = "l", main="Comparación mortalidad con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td3.1[["pin"]])) )
lines(td3.0[["pin"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


plot(td3.1[["D"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td3.0[["D"]])) )
lines(td3.0[["D"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")


plot(td3.1[["I"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td3.0[["I"]])) )
lines(td3.0[["I"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")

plot(td3.0[["he"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes", ylim = c(0, max(td3.0[["he"]])) )
lines(-td3.0[["comprobacion1"]], type = "l", col = "red")
legend(x = 0.5, y = 1.1, legend = c("Muertes sin déficit", "Muertes con déficit"), col = c("blue", "red"), lty = 1, xpd = TRUE, horiz = TRUE, bty = "n")

plot(td3.1[["I"]], type = "l", main="Comparación confinamiento con déficit/sin déficit", col = "blue", xlab = "Semanas", ylab = "Muertes")




