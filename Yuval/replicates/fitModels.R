
ncores = 13
B = 10
#Vs = 1:1250
Vs = 1:1250
first_seed = 100 
# Don't change - first_seed = 100

# Prepare features, voxel responses etc

split_train_test = function(seed, features, resp){
  # Remove features that have no variance
  standards = apply(features,2,sd)
  stopifnot(all(standards>0))

  set.seed(seed)
  train_samp = sample(1750,1500,replace = FALSE)
  test_samp = (1:1750)[-train_samp]
  trainF = features[train_samp,]
  trainY = v1$resp[train_samp,]
  testF = features[test_samp,]
  testY = v1$resp[test_samp,]
  data_st = list(trainF = trainF, trainY= trainY, 
                 testF = testF, testY = testY, 
                 seed = seed, train_samp = train_samp,
                 test_samp = test_samp)
}

library('glmnet')
run_glmnet = function(voxnum, trainF, trainY, alpha_val=0.8, nfolds = 6) {
  res = cv.glmnet(trainF, trainY[,voxnum],nfolds = 6,alpha = alpha_val)  
  best_lam = res$lambda.min
  best_ind = which(res$lambda==best_lam)
  pred_rule = res$glmnet.fit$beta[,best_ind]
  pred_score = res$cvm[best_ind]
  return(list(lambda = best_lam, alpha = alpha_val,beta_coef= pred_rule, score=pred_score))
}

getPreds = function(dataF,rules_struct,voxind){
  preds = matrix(0,nc = length(voxind),nr = nrow(dataF))
  for (i in 1:length(voxind)){
    preds[,i] = dataF %*% rules_struct[[voxind[i]]]$beta_coef
  }
  return (preds)
}


## Scale and center data
load('v1_data.RData')
standards = apply(fit_feat,2,sd)
constF = which(standards==0)
c_features = scale(fit_feat[,-constF],center = TRUE, scale = TRUE)
feat_mean_vec = attr(c_features,'scaled:center')
feat_scale_vec = attr(c_features,'scaled:center')

seeds = first_seed + 0:(B-1)

for(i in 1:B){
  require('parallel')
  data_st = split_train_test(seed = seeds[i], c_features, v1$resp)
  vox_prediction_rules = mclapply(Vs, run_glmnet, mc.cores = ncores,
                                  trainF=data_st$trainF, trainY = data_st$trainY) 
  attr(vox_prediction_rules,'seed') = data_st$seed
  attr(vox_prediction_rules,'train_samp') = data_st$train_samp
  attr(vox_prediction_rules,'test_samp') = data_st$test_samp
  save(vox_prediction_rules,file = sprintf('pred_rules_boots_%d.RData',data_st$seed))

  train_pred = getPreds(data_st$trainF,vox_prediction_rules,Vs)
  test_pred = getPreds(data_st$testF,vox_prediction_rules,Vs)
  vox_predictions = list(train_pred = train_pred, test_pred = test_pred,
                         seed = data_st$seed, train_samp = data_st$train_samp,
                         test_samp = data_st$test_samp)
  save(vox_predictions, file = sprintf('preds_boots_%d.RData',data_st$seed))
}

