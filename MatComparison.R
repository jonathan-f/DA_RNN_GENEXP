# Edoardo Milanetti
# library
########################
require(corrplot)
require(igraph)
require(reshape2)
require(ggplot2)

# working directory
########################
dir_input <- "/Users/edo/Downloads/Attention_Matrices/"
dir_output <- "/Users/edo/Downloads/FigMat/"

# functions
########################
RowNormFun <- function(x){
  y <- x/max(x)
  return(y)
}

Local_Network_Paramters <- function(g){
  
  Betweenness_Centrality <- betweenness(g, v = V(g), directed = FALSE, nobigint = TRUE, normalized = TRUE,  weights = E(g)$weight)
  Closeness_Centrality <- closeness(g, vids = V(g), weights = E(g)$weight, normalized = TRUE)
  Strength_Value <- strength(g)
  Diversity_Index <- diversity(g, weights = E(g)$weight, vids = V(g))
  Mean_Shortest_Path <- apply(distances(g, algorithm = "dijkstra", weights = E(g)$weight), 1, mean)
  Hub_Score_aus <- hub_score(g, weights = E(g)$weight)
  Hub_Score <- Hub_Score_aus$vector
  Clustering_Coefficient <- transitivity(g, type="barrat", weights = E(g)$weight)
  Degree_Value <- degree(g, loops = FALSE, normalized = TRUE)
  
  Parameters_DF <- cbind(Betweenness_Centrality, Closeness_Centrality, Strength_Value, Diversity_Index, Mean_Shortest_Path, Hub_Score, Clustering_Coefficient, Degree_Value)
  #Parameters_DF <- cbind(Clustering_Coefficient, Strength_Value, Degree_Value, Betweenness_Centrality)
  
  Parameters_DF <- as.data.frame(Parameters_DF)
  
  return(Parameters_DF)
  
}

Network_Generation <- function(EdgeList){
  
  EdgeList[,1]=as.character(EdgeList[,1])
  EdgeList[,2]=as.character(EdgeList[,2])
  EdgeList <- as.matrix(EdgeList)
  g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
  E(g)$weight <- as.numeric(EdgeList[,3])
  
  return(g)
  
}

MatList <- list.files(dir_input)
IntMatList <- MatList[grep("InteractionMatrix",MatList)]
AttMatList <- MatList[grep("att",MatList)]

vet_max <- c()
for(m in 1:length(AttMatList)){
  name_mat <- AttMatList[m]
  mat_aus <- read.table(name_mat)
  vet_max <- c(vet_max, max(mat_aus))
}
max_ass <- max(vet_max)


##### ##### ##### ##### ##### 
      ##### MAIN ######
##### ##### ##### ##### ##### 


# matrices comparison
#####################################
setwd(dir_input)
for(m in 1:length(IntMatList)){
  name_int_mat <- MatList[m]
  mat_int_aus <- read.table(name_int_mat)
  colnames(mat_int_aus) <- 1:ncol(mat_int_aus)
  mat_int_aus_abs <- abs(mat_int_aus)
  
  if(!mean(unlist(mat_int_aus_abs))==1){
    
    sub_name <- unlist(strsplit(name_int_mat,"InteractionMatrix_Jij_"))[2]
    name_att_mat <- AttMatList[grep(sub_name,AttMatList)]
    mat_att_aus <- read.table(name_att_mat)
    colnames(mat_att_aus) <- 1:ncol(mat_att_aus)
    mat_att_aus_norm <- as.matrix(mat_att_aus/max_ass)
    #mat_att_aus_norm <- mat_att_aus/max(mat_att_aus)
    #mat_att_aus_norm <- t(apply(mat_att_aus, 1, RowNormFun))
    
    vet_lower <- c()
    vet_upper <- c()
    for(r in 1:nrow(mat_att_aus_norm)){
      for(c in r:ncol(mat_att_aus_norm)){
        if(!r==c){
          vet_upper <- c(vet_upper, mat_att_aus_norm[r,c])
          vet_lower <- c(vet_lower, mat_att_aus_norm[c,r])
          
        }
      }
    }
    
    cor_tri <- cor(vet_upper,vet_lower)
    
    for(r in 1:nrow(mat_att_aus_norm)){
      for(c in r:ncol(mat_att_aus_norm)){
        if(!r==c){
          mat_att_aus_norm[r,c] <- mat_att_aus_norm[c,r]
        }
      }
    }
    
    
    # best threshold binarization
    cutoff_min <- min(mat_att_aus_norm)
    cutoff <- cutoff_min
    deltaCutoff <- 0.05
    Nsteps <- trunc((1 - cutoff_min)/deltaCutoff)
    vet_cor <- c()
    vet_cutoff <- c()
    for(i in 1:Nsteps){
      mat_att_aus_norm_bin <- mat_att_aus_norm
      cutoff <- cutoff + deltaCutoff
      mat_att_aus_norm_bin[mat_att_aus_norm_bin < cutoff] <- 0
      mat_att_aus_norm_bin[mat_att_aus_norm_bin >= cutoff] <- 1
      #cor_aus <- cor(as.numeric(as.vector(mat_att_aus_norm_bin)),as.vector(as.matrix(mat_int_aus_abs)))
      cor_aus <- cor(as.vector(unlist(mat_att_aus_norm_bin)), as.vector(unlist(mat_int_aus_abs)))
      vet_cor <- c(vet_cor, cor_aus)
      vet_cutoff <- c(vet_cutoff, cutoff)
    }
    
    vet_cutoff <- vet_cutoff[!is.na(vet_cor)]
    vet_cor <- vet_cor[!is.na(vet_cor)]
    plot(vet_cutoff,vet_cor)
    
    optimal_cutoff <- vet_cutoff[which.max(vet_cor)]
    mat_att_aus_norm_bin_optimal <- mat_att_aus_norm
    mat_att_aus_norm_bin_optimal[mat_att_aus_norm_bin_optimal < optimal_cutoff] <- 0
    mat_att_aus_norm_bin_optimal[mat_att_aus_norm_bin_optimal >= optimal_cutoff] <- 1
    
    if(!sum(mat_int_aus_abs)== ncol(mat_int_aus_abs)*nrow(mat_int_aus_abs)){
      
      titolo <- paste("Comparison_Jij_Att_",unlist(strsplit(sub_name,".txt")),".pdf")
      pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
        par(mfrow=c(2,2))
        par(mgp = c(2.5, 1, 0))
        corrplot(as.matrix(mat_int_aus_abs), method = "square", main="Binary Jij Matrix", mar=c(0,0,10,0),tl.pos = 'n')
        corrplot(as.matrix(mat_att_aus_norm), method = "square", main="Normalized Attention Matrix", mar=c(0,0,10,0),tl.pos = 'n')
        par(mar=c(4,4,10,0))
        plot(vet_cutoff, vet_cor, type = "l", lwd=2, ylab="Pearson Correlation", xlab="cutoff", xlim=c(0,1))
        corrplot(as.matrix(mat_att_aus_norm_bin_optimal), method = "square", main="Optimal Binary Attention Matrix", mar=c(1,1,9,1.2),tl.pos = 'n')
      dev.off()
      
    }
  }
}

# matrices comparison (no normalization!)
###########################################
df_par_tot <- data.frame()
for(m in 1:length(AttMatList)){
  name_att_mat <- AttMatList[m]
  name_aus <- unlist(strsplit(unlist(strsplit(name_att_mat,"attMat_"))[2],".txt"))
  mat_att_aus <- read.table(name_att_mat)
  colnames(mat_att_aus) <- 1:ncol(mat_att_aus)
  #mat_att_aus_norm <- t(apply(mat_att_aus, 1, RowNormFun))
  mat_att_aus_norm <- as.matrix(mat_att_aus)/max(as.matrix(mat_att_aus))
  EdgeList_net <- melt(mat_att_aus_norm)
  net_aus <- Network_Generation(EdgeList_net)
  
  #plot(net_aus,  full = F, vertex.label.cex=0.7, vertex.size = 2,weights=0.7,
  #     edge.arrow.size=0.01,vertex.label=NA,edge.color=sample("green",180, replace = T))
  
  plot(net_aus,
       edge.color=rep(c("red","pink"),length(EdgeList_net$value)),           # Edge color
       edge.width=EdgeList_net$value,                        # Edge width, defaults to 1
       edge.arrow.size=1,                           # Arrow size, defaults to 1
       edge.arrow.width=1,                          # Arrow width, defaults to 1
       edge.lty=c("solid")                           # Line type, could be 0 or ???blank???, 1 or ???solid???, 2 or ???dashed???, 3 or ???dotted???, 4 or ???dotdash???, 5 or ???longdash???, 6 or ???twodash???
       #edge.curved=c(rep(0,5), rep(1,5))            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
  )
  
  df_local_par <- Local_Network_Paramters(net_aus)
  df_local_par$System <- name_aus
  df_local_par$Gene <- 1:ncol(mat_att_aus)
  
  df_par_tot <- as.matrix(rbind(df_par_tot,as.matrix(df_local_par)))
  
}

# dataframe with all parameters
###########################################
df_par_tot <- as.data.frame(df_par_tot)
rownames(df_par_tot) <- NULL
df_par_tot$Clustering_Coefficient <- as.numeric(as.vector(df_par_tot$Clustering_Coefficient))
df_par_tot$Strength_Value <- as.numeric(as.vector(df_par_tot$Strength_Value))
df_par_tot$Degree_Value <- as.numeric(as.vector(df_par_tot$Degree_Value))
df_par_tot$Betweenness_Centrality <- as.numeric(as.vector(df_par_tot$Betweenness_Centrality))
df_par_tot$Closeness_Centrality <- as.numeric(as.vector(df_par_tot$Closeness_Centrality))
df_par_tot$Diversity_Index <- as.numeric(as.vector(df_par_tot$Diversity_Index))
df_par_tot$Mean_Shortest_Path <- as.numeric(as.vector(df_par_tot$Mean_Shortest_Path))
df_par_tot$Hub_Score <- as.numeric(as.vector(df_par_tot$Hub_Score))

df_par_tot$System <- as.factor(df_par_tot$System)
df_par_tot$Gene <- as.numeric(as.vector(df_par_tot$Gene))

df_par_tot[,1:7] <- (df_par_tot[,1:7] - apply(df_par_tot[,1:7],2,mean))/apply(df_par_tot[,1:7],2,sd)

# plots
###########################################
cols <- c("darkred","red","orange","gray","cyan","blue", "steelblue","green","black","yellow","pink","darkgreen")

titolo <- paste("Clustering_Coefficient",".pdf", sep="")
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
p <- ggplot(df_par_tot, aes(x = Gene, y = Clustering_Coefficient)) + 
      geom_line(aes(color = System)) + 
      scale_color_manual(values = cols)
p + theme(panel.background = element_rect(fill='gray95', colour='black'))
dev.off()

titolo <- paste("Strength_Value",".pdf", sep="")
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
p <- ggplot(df_par_tot, aes(x = Gene, y = Strength_Value)) + 
  geom_line(aes(color = System)) + 
  scale_color_manual(values = cols)
p + theme(panel.background = element_rect(fill='gray95', colour='black'))
dev.off()

titolo <- paste("Betweenness_Centrality",".pdf", sep="")
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
p <- ggplot(df_par_tot, aes(x = Gene, y = Betweenness_Centrality)) + 
  geom_line(aes(color = System)) + 
  scale_color_manual(values = cols)
p + theme(panel.background = element_rect(fill='gray95', colour='black'))
dev.off()

titolo <- paste("Closeness_Centrality",".pdf", sep="")
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
p <- ggplot(df_par_tot, aes(x = Gene, y = Closeness_Centrality)) + 
  geom_line(aes(color = System)) + 
  scale_color_manual(values = cols)
p + theme(panel.background = element_rect(fill='gray95', colour='black'))
dev.off()

titolo <- paste("Diversity_Index",".pdf", sep="")
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
p <- ggplot(df_par_tot, aes(x = Gene, y = Diversity_Index)) + 
  geom_line(aes(color = System)) + 
  scale_color_manual(values = cols)
p + theme(panel.background = element_rect(fill='gray95', colour='black'))
dev.off()

titolo <- paste("Mean_Shortest_Path",".pdf", sep="")
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
p <- ggplot(df_par_tot, aes(x = Gene, y = Mean_Shortest_Path)) + 
  geom_line(aes(color = System)) + 
  scale_color_manual(values = cols)
p + theme(panel.background = element_rect(fill='gray95', colour='black'))
dev.off()

titolo <- paste("Hub_Score",".pdf", sep="")
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)
p <- ggplot(df_par_tot, aes(x = Gene, y = Hub_Score)) + 
  geom_line(aes(color = System)) + 
  scale_color_manual(values = cols)
p + theme(panel.background = element_rect(fill='gray95', colour='black'))
dev.off()

# barplot
################
df_par_tot$System <- as.character(df_par_tot$System)

titolo <- "ClusteringCoefficient_Barplot"
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)

par(mfrow=c(3,4))
df_Clustering_Coefficient <- data.frame()
for(i in 1:length(unique(df_par_tot$System))){
  title <-  unique(df_par_tot$System)[i]
  Ylab <- "Clustering Coefficient"
  
  barplot(sort(df_par_tot$Clustering_Coefficient[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T), ylim=c(-5,25), col="blue", main=title, ylab=Ylab)
  df_Clustering_Coefficient <- as.matrix(rbind(df_Clustering_Coefficient, sort(df_par_tot$Clustering_Coefficient[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T)))
}
dev.off()

df_Clustering_Coefficient <- as.data.frame(df_Clustering_Coefficient)
rownames(df_Clustering_Coefficient) <- unique(df_par_tot$System)
colnames(df_Clustering_Coefficient) <- NULL
write.csv(df_Clustering_Coefficient,paste(dir_output, "df_Clustering_Coefficient.csv", sep=""))

pdf(file = paste(dir_output,"Clustering_Clustering_Plot", ".pdf", sep=""), width = 20, height = 10)
plot(hclust(dist(df_Clustering_Coefficient)))
dev.off()

titolo <- "BetweennessCentrality_Barplot"
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)

par(mfrow=c(3,4))
df_Betweenness <- data.frame()
for(i in 1:length(unique(df_par_tot$System))){
  title <-  unique(df_par_tot$System)[i]
  Ylab <- "Betweenness Centrality"
  Ylim <- c(min(df_par_tot$Betweenness_Centrality)-0.10*min(df_par_tot$Betweenness_Centrality),max(df_par_tot$Betweenness_Centrality)+0.10*max(df_par_tot$Betweenness_Centrality))
  barplot(sort(df_par_tot$Betweenness_Centrality[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T), ylim=Ylim, col="blue", main=title, ylab=Ylab)
  df_Betweenness <- as.matrix(rbind(df_Betweenness, sort(df_par_tot$Betweenness_Centrality[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T)))
}
dev.off()

df_Betweenness <- as.data.frame(df_Betweenness)
rownames(df_Betweenness) <- unique(df_par_tot$System)
colnames(df_Betweenness) <- NULL

pdf(file = paste(dir_output,"Betweenness_Clustering_Plot", ".pdf", sep=""), width = 20, height = 10)
plot(hclust(dist(df_Betweenness)))
dev.off()


titolo <- "ClosenessCentrality_Barplot"
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)

par(mfrow=c(3,4))
df_Closeness <- data.frame()
for(i in 1:length(unique(df_par_tot$System))){
  title <-  unique(df_par_tot$System)[i]
  Ylab <- "Closeness Centrality"
  Ylim <- c(min(df_par_tot$Closeness_Centrality)-abs(0.10*min(df_par_tot$Closeness_Centrality)),max(df_par_tot$Closeness_Centrality)+abs(0.10*max(df_par_tot$Closeness_Centrality)))
  barplot(sort(df_par_tot$Closeness_Centrality[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T), ylim=Ylim, col="blue", main=title, ylab=Ylab)
  df_Closeness <- as.matrix(rbind(df_Closeness, sort(df_par_tot$Closeness_Centrality[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T)))
}
dev.off()

df_Closeness <- as.data.frame(df_Closeness)
rownames(df_Closeness) <- unique(df_par_tot$System)
colnames(df_Closeness) <- NULL

pdf(file = paste(dir_output,"Closeness_Clustering_Plot", ".pdf", sep=""), width = 20, height = 10)
plot(hclust(dist(df_Closeness)))
dev.off()


titolo <- "Strength_Barplot"
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)

par(mfrow=c(3,4))
df_Strength <- data.frame()
for(i in 1:length(unique(df_par_tot$System))){
  title <-  unique(df_par_tot$System)[i]
  Ylab <- "Strength"
  Ylim <- c(min(df_par_tot$Strength_Value)-abs(0.10*min(df_par_tot$Strength_Value)),max(df_par_tot$Strength_Value)+abs(0.10*max(df_par_tot$Strength_Value)))
  barplot(sort(df_par_tot$Strength_Value[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T), ylim=Ylim, col="blue", main=title, ylab=Ylab)
  df_Strength <- as.matrix(rbind(df_Strength, sort(df_par_tot$Strength_Value[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T)))
  
}
dev.off()

df_Strength <- as.data.frame(df_Strength)
rownames(df_Strength) <- unique(df_par_tot$System)
colnames(df_Strength) <- NULL

pdf(file = paste(dir_output,"Strength_Clustering_Plot", ".pdf", sep=""), width = 20, height = 10)
plot(hclust(dist(df_Strength)))
dev.off()


titolo <- "DiversityIndex_Barplot"
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)

par(mfrow=c(3,4))
df_Diversity <- data.frame()
for(i in 1:length(unique(df_par_tot$System))){
  title <-  unique(df_par_tot$System)[i]
  Ylab <- "Diversity Index"
  Ylim <- c(min(df_par_tot$Diversity_Index)-abs(0.10*min(df_par_tot$Diversity_Index)),max(df_par_tot$Diversity_Index)+abs(0.10*max(df_par_tot$Diversity_Index)))
  barplot(sort(df_par_tot$Diversity_Index[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T), ylim=Ylim, col="blue", main=title, ylab=Ylab)
  df_Diversity <- as.matrix(rbind(df_Diversity, sort(df_par_tot$Diversity_Index[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T)))
  
}
dev.off()

df_Diversity <- as.data.frame(df_Diversity)
rownames(df_Diversity) <- unique(df_par_tot$System)
colnames(df_Diversity) <- NULL

pdf(file = paste(dir_output,"Diversity_Clustering_Plot", ".pdf", sep=""), width = 20, height = 10)
plot(hclust(dist(df_Diversity)))
dev.off()

titolo <- "MeanShortestPath_Barplot"
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)

par(mfrow=c(3,4))
df_ShortestPath <- data.frame()
for(i in 1:length(unique(df_par_tot$System))){
  title <-  unique(df_par_tot$System)[i]
  Ylab <- "Mean Shortest Path"
  Ylim <- c(min(df_par_tot$Mean_Shortest_Path)-abs(0.10*min(df_par_tot$Mean_Shortest_Path)),max(df_par_tot$Mean_Shortest_Path)+abs(0.10*max(df_par_tot$Mean_Shortest_Path)))
  barplot(sort(df_par_tot$Mean_Shortest_Path[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T), ylim=Ylim, col="blue", main=title, ylab=Ylab)
  df_ShortestPath <- as.matrix(rbind(df_ShortestPath, sort(df_par_tot$Mean_Shortest_Path[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T)))
  
}
dev.off()

df_ShortestPath <- as.data.frame(df_ShortestPath)
rownames(df_ShortestPath) <- unique(df_par_tot$System)
colnames(df_ShortestPath) <- NULL

pdf(file = paste(dir_output,"ShortestPath_Clustering_Plot", ".pdf", sep=""), width = 20, height = 10)
plot(hclust(dist(df_ShortestPath)))
dev.off()


titolo <- "HubScore_Barplot"
pdf(file = paste(dir_output,titolo, ".pdf", sep=""), width = 20, height = 10)

par(mfrow=c(3,4))
df_HubScore <- data.frame()
for(i in 1:length(unique(df_par_tot$System))){
  title <-  unique(df_par_tot$System)[i]
  Ylab <- "HubScore"
  Ylim <- c(min(df_par_tot$Hub_Score)-abs(0.10*min(df_par_tot$Hub_Score)),max(df_par_tot$Hub_Score)+abs(0.10*max(df_par_tot$Hub_Score)))
  barplot(sort(df_par_tot$Hub_Score[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T), ylim=Ylim, col="blue", main=title, ylab=Ylab)
  df_HubScore <- as.matrix(rbind(df_HubScore, sort(df_par_tot$Hub_Score[df_par_tot$System %in% unique(df_par_tot$System)[i]], decreasing = T)))
  
}
dev.off

df_HubScore <- as.data.frame(df_HubScore)
rownames(df_HubScore) <- unique(df_par_tot$System)
colnames(df_HubScore) <- NULL

pdf(file = paste(dir_output,"HubScore_Clustering_Plot", ".pdf", sep=""), width = 20, height = 10)
plot(hclust(dist(df_HubScore)))
dev.off()


# matrix correlations
###########################################
setwd(dir_input)
df_Jij <- data.frame()
vet_name <- c()
for(m in 1:length(IntMatList)){
  name_int_mat <- MatList[m]
  mat_int_aus <- read.table(name_int_mat)
  name_int_mat <- unlist(strsplit(unlist(strsplit(name_int_mat,"InteractionMatrix_"))[2],".txt"))
  name_int_mat <- unlist(strsplit(name_int_mat,"Jij_"))[2]
  vet_name <- c(vet_name,name_int_mat)
  colnames(mat_int_aus) <- 1:ncol(mat_int_aus)
  df_Jij <- as.matrix(rbind(df_Jij, as.vector(unlist(mat_int_aus))))
}
control <- !(abs(apply(df_Jij, 1, mean))==0 | abs(apply(df_Jij, 1, mean))==1)
rownames(df_Jij) <- vet_name
df_Jij <- df_Jij[control,]
cor_Jij <- cor(t(df_Jij))

pdf(file = paste(dir_output,"density", ".pdf", sep=""), width = 20, height = 10)
plot(density(as.vector(cor_Jij[c(1:3,8:10),c(1:3,8:10)])[!as.vector(cor_Jij[c(1:3,8:10),c(1:3,8:10)])==1]), xlim=c(-0.15,0.15), ylim=c(0,15), main="", xlab="Pearson Correlation")
par(new=T)
plot(density(as.vector(cor_Jij[4:7,4:7])[!as.vector(cor_Jij[4:7,4:7])==1]), xlim=c(-0.15,0.15), ylim=c(0,15), col="red", main="", xlab="Pearson Correlation")
legend(0.07, 13, legend=c("Other Networks", "Oscillating Networks"),
       col=c("black", "red"), lty=1, cex=0.6)
dev.off()

att_corresponding <- rownames(cor_Jij)
df_Att <- data.frame()
vet_name <- c()
for(m in 1:length(AttMatList)){
  name_att_mat <- AttMatList[grep(att_corresponding[m],AttMatList)]
  if(length(name_att_mat)==1){
    name_aus <- unlist(strsplit(unlist(strsplit(name_att_mat,"attMat_"))[2],".txt"))
    vet_name <- c(vet_name,name_aus)
    
    mat_att_aus <- read.table(name_att_mat)
    df_Att <- as.matrix(rbind(df_Att, as.vector(unlist(mat_att_aus))))
  }
}

rownames(df_Att) <- vet_name
#colnames(df_Att) <- vet_name
cor_Att <- cor(t(df_Att))

cor_aus <- cor(as.vector(cor_Att)[!as.vector(cor_Att)==1], as.vector(cor_Jij)[!as.vector(cor_Jij)==1])
cor_test <- cor.test(as.vector(cor_Att)[!as.vector(cor_Att)==1], as.vector(cor_Jij)[!as.vector(cor_Jij)==1])
title <- paste("Pearson Correlation: ", round(cor_aus,3), " | p-value: ", round(cor_test$p.value,3), sep="")

pdf(file = paste(dir_output,"Correlations", ".pdf", sep=""), width = 20, height = 10)
plot(as.vector(cor_Att)[!as.vector(cor_Att)==1], as.vector(cor_Jij)[!as.vector(cor_Jij)==1], pch=20, cex.axis=1.5, cex=1.5, 
     ylab="Jij Matrix Correlation",  xlab="Attention Matrix Correlation", main=title, ylim=c(-0.15,0.15), xlim=c(-0.30,0.30))
abline(lm(as.vector(cor_Jij)[!as.vector(cor_Jij)==1] ~ as.vector(cor_Att)[!as.vector(cor_Att)==1], data = mtcars), col = "blue")
dev.off()

