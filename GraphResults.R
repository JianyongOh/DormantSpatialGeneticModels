#why is rate estimate = E(difference)/T

latlong = env[,3:4]
dist = dist(latlong)
ponds = cmdscale(dist,k=1,eig=TRUE)
pov = ponds$eig[1]/sum(ponds$eig)
hist(dist)

#Read fasta file data and save it in a dataframe "newdata"
library("Biostrings")
Data <- readDNAStringSet("C:/Users/ohjia/Desktop/MORSE/URSS/INPonds.bac.final.0.03.fasta")
newdata <- as.data.frame(Data)
rownames(newdata) = c()

#finding out which location-pairs fit into which distance bracket
dist = as.matrix(dist)
indices = list()
for (i in seq(0,0.32,0.02)){
  temp = c()
  for (j in 1:49){
    for (k in j:49){
      if (dist[j,k]>i & dist[j,k]<i+0.02){
        temp = rbind(temp,c(j,k))
        indices[[i*50+1]] = temp
      }
    }
  }
}
indices

#calculate actual POI and log
POI_act = rep(0,17)
for (i in 1:17){
  for (j in 1:nrow(indices[[i]])){
    temp = indices[[i]][j,]
    POI_act[i] = POI_act[i] + sum(act.com3[temp[1],]*act.com3[temp[2],])
  }
  POI_act[i] = POI_act[i]/nrow(indices[[i]])
}
POI_act_log = log(POI_act)

#reading simulated data
for (i in 1:17){
  path = paste("C:/Users/ohjia/Desktop/MORSE/URSS/Data0.075/timedata",i,".csv",sep="")
  temp = read.csv(path)
  if (i == 1){
    df = temp
  }
  else{
    df = cbind(df,temp)
  }
}

#average time spent active at each separation 0-16
avg = c()
for (i in seq(2,51,3)){
  avg = append(avg,mean(df[,i]))
}

#calibrate rate = log(POI_act)/T
Pmut = -log(POI_act)/avg
Pmut = Pmut[5]

#calculate simulated POI and log
POI_sim_log = -Pmut*avg

#save results into a pdf
pdf("POI_sim0.075_log.pdf")
plot(0:16, POI_sim_log, type = 'l', ylim = c(-4.5,-3), xlab = "distance", ylab = "log(POI)")
par(new=TRUE)
plot(0:16, POI_act_log, type = 'l', ylim = c(-4.5,-3), xlab = "distance", ylab = "log(POI)", col = "red")
dev.off()

