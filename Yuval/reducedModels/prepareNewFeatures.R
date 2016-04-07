
# Prepare new feature set...

load("/Users/yuvalb/Studies/fMRI/V4Data/wavpyr.RData")
load("/Users/yuvalb/Teaching/StatMethNeuroW15/Assignments/Visual/stim.RData")


# generates the feature attributes vector.
# Orientation -> featAtt[1,] (1 vertical, 5 horizontal)  
# Pyramid level -> featAtt[2,]
# Vertical location -> featAtt[3,]  (vert/horz might be confused)
# Horizontal location -> featAtt[4,]

featAtt = function(){
  # featAttVec: row1 - orientation, row2 - pyr-level, row3 - locationx, row4 - lovationy
  pyr  = c(1,8,4*8,16*8,64*8,256*8,1024*8)
  lens = c(0,1,2,   4,   8,   16,    32,  64)
  featAttVec = matrix(0,nr=4,nc=sum(pyr))
  featAttVec[1,] = (0:(sum(pyr)-1)) %% 8
  cumsumpyr = cumsum(pyr)
  featAttVec[2,1] = 1;
  featAttVec[3:4,1] = 0;
  for (i in 1:(length(cumsumpyr)-1)) {
    k = 0
    for (j in seq((cumsumpyr[i]+1),cumsumpyr[i+1],8)) {
      featAttVec[2,j:(j+7)] = i+1
      featAttVec[3,j:(j+7)] = floor(k/(lens[i+1])) + 1
      featAttVec[4,j:(j+7)] = (k%%(lens[i+1])) + 1
      k = k + 1
    }
  }
  return(featAttVec)
}
featAll = featAtt()
stopifnot(ncol(wav.pyr)==ncol(featAll))

# Projection from large features to small features
library('Matrix')

n_oldfeat = ncol(wav.pyr)
n_newfeat = 1+(ncol(wav.pyr)-1)/8 
proj_matrix = Matrix(0,nr=n_oldfeat,nc=n_newfeat,sparse = TRUE)
proj_matrix[1,1]=1
for (i_oldfeat in 2:ncol(wav.pyr)){
  i_newfeat = 2 + floor((i_oldfeat - 1 - 1)/8)
  proj_matrix[i_oldfeat,i_newfeat] = 1
}

filter_res = abs(stim %*% wav.pyr)^2
filter_res_val = abs(val.stim %*% wav.pyr)^2

features_proj =  sqrt(filter_res %*% proj_matrix)
features_proj_val =  sqrt(filter_res_val %*% proj_matrix)

nor_features = sqrt(features_proj)  
nor_features_val = sqrt(features_proj_val)

nor_feat_ann = featAll[,featAll[1,]==0]
save(file = '~/Studies/fMRI/Charles/noOrientFeat.RData',nor_features, nor_features_val,nor_feat_ann)
