load("/home/snarles/stat312data/vidData1sec.RData")

dim(vidTrainF) ## 7128 sample size, 8556 features
dim(vidTrainY) ## 2088 dimensions Y
running
## running has 18 components, total number of components is equal to dim of Y (2088)
run_loc
## run_loc has 18 components
## all_samp has 2088

plot(vidTrainF[1:100, 1], type = "l")
plot(vidTrainF[, 6], type = "l")


## find places with large delta
diff <- function(v) v[-1] - v[-length(v)]
plot(sort(order(-abs(diff(vidTrainF[, 1])))[1:100]))
