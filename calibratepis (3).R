
#Similar to td_eval (Funcs.R)

calibrate_pis <- function(pis_guess,H,eps,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,C,N,scale1,scale2){

    pis1=pis_guess[1]/scale1
    pis2=pis_guess[2]/scale2
    pis3=pis_guess[3]
    
    S = integer(H)
    I = integer(H)
    R = integer(H)
    D = integer(H)
    Tr = integer(H)
    
    S[1] = 1 - eps
    I[1] = eps
    R[1] = 0
    D[1] = 0
    
    for (t in 2:H){
      Tr[t-1] = pis1*S[t-1]*C^2*I[t-1]+pis2*S[t-1]*N^2*I[t-1]+pis3*S[t-1]*I[t-1]
      S[t] = S[t-1]-Tr[t-1]
      I[t] = I[t-1]+Tr[t-1]-(pir+pid)*I[t-1]
      R[t] = R[t-1]+pir*I[t-1]
      D[t] = D[t-1]+pid*I[t-1]
    }
    
    err1 =pis1_shr_target-(pis1*C^2)/(pis1*C^2+pis2*N^2+pis3)
    err2 =pis2_shr_target-(pis2*N^2)/(pis1*C^2+pis2*N^2+pis3)
    err3 =RplusD_target-(R[length(R)]+D[length(D)])
    
    results  = 100*(err1)^2 + 100*(err2)^2 + 100*(err3)^2

    RnotSIR=T[1]/I[1]/(pir+pid)
    
    return(results)
  
}

#Minimizes calibrate_pis' error terms

calib <- function(pis_guess,H,eps,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,C,N,scale1,scale2){

  objec = function(pve) calibrate_pis(pve,H=H,eps=eps,pir=pir,pid=pid,pis1_shr_target=pis1_shr_target,
                                      pis2_shr_target=pis2_shr_target,RplusD_target=RplusD_target,C=C,
                                      N=N,scale1=scale1,scale2=scale2)

  respi = optim(par = pis_guess, fn = objec, gr= NULL, method='BFGS', control = 1E-6)
  
  return(respi)
  
}


