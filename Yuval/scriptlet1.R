####
##  Some additional code to be used with identification.Rmd
####

matchResults = function(preds, dataY, score,invCov = diag(length(usevox)), voxind=usevox, rankit = TRUE){
  npoints = nrow(preds)
  res <- list()
  for (i in 1:nrow(dataY)){
    if (score=='MSE'){
      Yi = matrix(rep(dataY[i,voxind],npoints),nr=npoints,byrow=T)
      ResYi = Yi - preds
      res[[i]] = diag(ResYi %*% invCov %*% t(ResYi))
    }
    else if (score=='COR'){
      res[[i]] = -apply(preds,1,cor,dataY[i,voxind])
    }
    else if (score=='COV'){
      res[[i]] = -apply(preds,1,cov,dataY[i,voxind])
    }
    if ((i%%20)==0){
      cat(',')
    }
  }
  res <- do.call(rbind, res)
  if (rankit) res <- t(apply(res, 1, rank))
  res
}

res <- matchResults(validpred,validY[1:120,],'MSE',invCov = ridgeinv)
res[1, ]

sum(diag(res) == 1)/120
