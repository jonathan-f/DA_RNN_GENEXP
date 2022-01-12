# 08/2021 Jonathan Fiorentino

library('TreeDist')
library(ape)

# Load the trees from txt files
noisetree <- read.tree('./newick_trees/Noise_analysis_tree.txt')
hubtree <- read.tree('./newick_trees/HubScore_tree.txt')
betweennesstree<- read.tree('./newick_trees/Betweenness_tree.txt')
cluscoefftree<- read.tree('./newick_trees/Clustering_Coefficient_tree.txt')

# Put all the trees in a list
treeLabels <- c('Noise_analysis',
'HubScore',
'Betweenness',
'Clustering_Coefficient')

treeList <- c(noisetree,hubtree, betweennesstree,cluscoefftree)

# Define an empty dataframe with the labels of the trees on the rows and columns
distance.df <- data.frame(matrix(ncol = length(treeLabels), nrow = length(treeLabels)))
colnames(distance.df) <- treeLabels
rownames(distance.df) <- treeLabels

for(i in seq(1,length(treeLabels),1)){
  for(j in seq(1,length(treeLabels),1)){
    print(c(i,j))
    tree1 <- treeList[i]
    tree2 <- treeList[j]
    distance.df[i,j] <- TreeDistance(tree1, tree2)
  }
}

# Save the dataframe to a csv file
write.csv(distance.df,"dendrodist.csv", row.names = TRUE)