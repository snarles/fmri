
library('R.matlab')

dat_1 = R.matlab::readMat('/Users/yuvalb/Dropbox/File requests/Matrices/Charles Zheng - All_Sub_Rest1.mat')
dat_2 = R.matlab::readMat('/Users/yuvalb/Dropbox/File requests/Matrices/Charles Zheng - All_Sub_Rest2.mat')

image(dat_1$All.Sub[1,,])

par(mfcol = c(2,1))
hist(dat_1$All.Sub[2,,]-dat_2$All.Sub[1,,],xlim = c(-1,1),
     breaks = seq(-1.2,1.2,0.05))
hist(dat_1$All.Sub[2,,]-dat_2$All.Sub[2,,],xlim = c(-1,1),
     breaks = seq(-1.2,1.2,0.05))

plot(dat_1$All.Sub[8,,],dat_2$All.Sub[100,,],xlim = c(-.5,1.6), ylim = c(-.5,1.6))

plot(dat_1$All.Sub[8,,],dat_2$All.Sub[100,,],xlim = c(-.5,1.6), ylim = c(-.5,1.6))

plot(dat_1$All.Sub[100,1,]+ dat_2$All.Sub[8,1,],
     dat_1$All.Sub[100,1,]-dat_2$All.Sub[8,1,],xlim = c(-.5,1.6), ylim = c(-.5,1.6))

dplot = function(a,b,c, newp = FALSE,...) {
    if (newp){
        plot(dat_1$All.Sub[a,c,],dat_2$All.Sub[b,c,],
             xlim = c(-.5,1.6), ylim = c(-.5,1.6),...)
    } else {
        points(dat_1$All.Sub[a,c,],dat_2$All.Sub[b,c,],
               xlim = c(-.5,1.6), ylim = c(-.5,1.6),...)
    }
}
dbplot = function(a,b,c, newp = FALSE,...) {
    if (newp){
        plot(dat_1$All.Sub[a,c,]+dat_2$All.Sub[b,c,],
             dat_1$All.Sub[a,c,]-dat_2$All.Sub[b,c,],
             xlim = c(-.5,1.6), ylim = c(-1.6,1.6),...)
    } else {
        points(dat_1$All.Sub[a,c,]+dat_2$All.Sub[b,c,],
               dat_1$All.Sub[a,c,]-dat_2$All.Sub[b,c,],
               xlim = c(-.5,1.6), ylim = c(-1.6,1.6),...)
    }
}

dplot(100,100,2,col = 0,newp = 1)
for (i in 1:10){
    dplot(100,100,i,col = i, pch = 20, cex = 0.3)
}

dbplot(100,100,2,col = 0,newp = 1)
for (i in 1:200){
    dbplot(8,8,i,col = 1, pch = 20, cex = 0.3)
}



cor(as.numeric(dat_1$All.Sub[7,,]),as.numeric(dat_2$All.Sub[8,,]))


