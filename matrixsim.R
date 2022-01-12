#library('cultevo')

# 08/2021 Jonathan Fiorentino

# Change to your directory
wd <- '/Users/jonathan/Desktop/IIT/TIME_SERIES_PREDICTION/deep_rnn/FINAL_ANALYSIS/'
setwd(wd)

noise.mat <- as.matrix(read.csv(file = './matrices/noisematrix.csv',row.names = 1))
hub.mat <- as.matrix(read.csv(file = './matrices/matrix_HubScore.csv',row.names = 1))
bet.mat <- as.matrix(read.csv(file = './matrices/matrix_Betweenness.csv',row.names = 1))
clus.mat <- as.matrix(read.csv(file = './matrices/matrix_Clustering_Coefficient.csv',row.names = 1))

myLabels <- c('Noise_analysis',
                'HubScore',
                'Betweenness',
                'Clustering_Coefficient')

matList <- list(noise.mat,hub.mat,bet.mat,clus.mat)

pear.log <- data.frame(matrix(ncol = length(myLabels), nrow = length(myLabels)))
colnames(pear.log) <- myLabels
rownames(pear.log) <- myLabels

pear.pval.log <- data.frame(matrix(ncol = length(myLabels), nrow = length(myLabels)))
colnames(pear.pval.log) <- myLabels
rownames(pear.pval.log) <- myLabels

for(i in seq(1,length(myLabels),1)){
  for(j in seq(1,length(myLabels),1)){

    pear.log[i,j] <- cor(c(sign(matList[[i]])*log(1.0+abs(matList[[i]]))), c(sign(matList[[j]])*log(1.0+abs(matList[[j]]))))
    pear.pval.log[i,j] <- cor.test(c(sign(matList[[i]])*log(1.0+abs(matList[[i]]))), c(sign(matList[[j]])*log(1.0+abs(matList[[j]]))))$p.value
    
  
  }
}

# Save the dataframe to a csv file
write.csv(pear.log,"pearson_logmatrix.csv", row.names = TRUE)
write.csv(pear.pval.log,"pearson_logpval.csv", row.names = TRUE)