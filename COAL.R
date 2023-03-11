library(MASS) 
library(glmnet)
library(lqa) 
library(jackknifeKME)
###using lqa
coal <- function(Y, delta, A, W){
  srt <- order(Y)
  sy <- as.double(Y[srt]) #sorted Y
  sdelta <- as.integer(delta[srt])# sorted delta
  sA <- as.integer(A[srt])# sorted A
  sW=W[srt,]#sorted W
  ####
 # kmweights <- numeric(n)
  #kmweights[1] <- 1/n
  #for (i in 2:n) {
 #   kmweights[i] <- kmweights[i - 1] * (n - i + 2)/(n - i + 
 #                                                     1) * (((n - i + 1)/(n - i + 2))^sdelta[i - 1])
 # }
 # kmwts <- kmweights * sdelta
  #if (sdelta[n] == 0) 
    
   # kmwts[n] <- 1 - sum(kmwts)
  #####
  wei=kmweight(sy,sdelta)
  Data <- as.data.frame(cbind(Y = sy, A = sA, sW), row.names = NULL)
  p = dim(W)[2]
  var.list = colnames(Data)[(1:p)+2]
  y.form = formula(paste("Y~."))
  lm.Y = lm(y.form,data=Data, weight=wei)
  #lm.Y = glm(y.form,data=Data, family = "gaussian")
  betaXY = coef(lm.Y)[3:(dim(W)[2]+2)] 
  #----------------------------------------------------
  # set vector of possible lambda's to try
  lambda_vec = c(0.49, 0.1, 0.05, seq(0, -10, length.out = 11))
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  
  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=dim(W)[2],ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  #----------------------------------------------------
  wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
    trt = untrt = diff_vec = rep(NA,length(beta)) 
    names(trt) = names(untrt) = names(diff_vec) = varlist
    for(jj in 1:length(varlist)){ 
      this.var = paste("w",varlist[jj],sep="") 
      DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
      trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
      untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
      diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
    } 
    wdiff_vec = diff_vec * abs(beta) 
    wAMD = c( sum(wdiff_vec))
    ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
    return(ret) 
  }
  create_weights = function(fp,fA,fw){
    fw = (fp)^(-1)
    fw[fA==0] = (1 - fp[fA==0])^(-1)
    return(fw)
  }
  ######################################################################################
  #####  Run censored outcome adaptive lasso for each lambda value 
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  
  n <- length(A)
  
  for(lil in names(lambda_vec)){
  
    tryCatch({
      il = lambda_vec[lil]
      ig = gamma_vals[lil]
      
      ### create the censored outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      coal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
      ### run censored outcome-adaptive lasso model with appropriate penalty
      logit_coal = lqa(w.full.form, data=Data, penalty=coal_pen, family=binomial(logit))
      
      
      # generate propensity score
      Data[,paste("f.pA",lil,sep="")] = predict(logit_coal)$mu.new
      # save propensity score coefficients
      coeff_XA[var.list,lil] = coef(logit_coal)[var.list]
      # create inverse probability of treatment weights
      Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
      
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                    wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
      
    }, error=function(e){
      "Happy"   })
  } # close loop through lambda value
  tt=which.min(wAMD_vec)
  lil=names(tt)
  il = lambda_vec[lil]
  ig = gamma_vals[lil]
  coal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
  logit_coal = lqa(w.full.form, data=Data, penalty=coal_pen, family=binomial(logit))
  # save propensity score coefficients
  coeff_XA=coef(logit_coal)
  re=list(coef=coeff_XA,ps=predict(logit_coal)$mu.new)
  return(re)
}
