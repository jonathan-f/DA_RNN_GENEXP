# Dual Attention Recurrent Neural Network for the prediction of gene expression time series data

This repository contains the code and data associated to the manuscript [Monti*, Fiorentino*, et al, Prediction of gene expression time series and structural analysis of gene regulatory networks using recurrent neural networks, arXiv, 2021](https://arxiv.org/abs/2109.05849). We generate synthetic gene expression time series data from a set of archetypical gene regulatory networks and we use a parallel Dual Attention Recurrent Neural Network (DA-RNN) for the time series prediction. Moreover, using tools from graph theory we show that graph descriptors of the input attention layer allow to distinguish different gene regulatory network architectures. We also study the effect on the predcition of Gaussian noise addition on protein concentrations. More details can be found in the associated manuscript; please cite it if you use our code and/or data.

# 1. Generation of synthetic gene expression time series data

This is achieved using the Gillespie algorithm .... add details and instructions

The time series of protein concentrations and the ground truth gene interaction matrices are provided in the [DATA](/DATA/) folder for each gene regulatory network.

# 2. Setting up a python3 virtual environment

For all the following analyses a python3 virtual environment can be set up using the file [requirements.txt](/requirements.txt).

The virtual environment, called for instance 'darnn-venv', can be created through the following commands:

```
python -m venv darnn-venv
source darnn-venv/bin/activate
pip3 install ipykernel
ipython kernel install --user --name=darnn-venv
python -m pip install -r requirements.txt
```


# 3. Network training and analysis of the input attention layer

The python script [Parallel_DA_RNN_training.py](/Parallel_DA_RNN_training.py) allows to train a parallel DA-RNN on gene expression time series data for all the gene regulatory networks generated in section 1.. The network parameters obtained for the training are saved in dedicated folders; also the matrix of input attention is saved for downstream analyses (see below). They are already provided in the [DATA](/DATA/) folder. The python script [My_allFunctions.py](/My_allFunctions.py) contains all the functions needed for training the DA-RNN and for the subsequent time series prediction.

The R script XXX computes the graph theory descriptors (Clustering coefficient, Betweenness, Hub Score) of the input attention matrices for each gene regulatory network used in our study.

The analysis of Gaussian noise addition is provided in the python Jupyter notebook [DA_RNN_NOISE_ANALYSIS.ipynb](/DA_RNN_NOISE_ANALYSIS.ipynb); the results are in the [DATA](/DATA/) folder.

The Pearson correlation between the matrices obtained from the graph theory descriptors and from the noise analysis is computed using the R script [matrixsim.R](/matrixsim.R), with inout data provided in the [matrices](/matrices/) folder.

The PCA and clustering analyses of the matrices mentioned above is performed in the Jupyter notebook [DA_RNN_DOWNSTREAM_ANALYSIS.ipynb](/DA_RNN_DOWNSTREAM_ANALYSIS.ipynb). The dendrograms associated to the different clusterings are provided in the [newick_trees](/newick_trees/) folder in Newick tree format. and the distances between dendrograms are computed using the R script [treeanalysis.R](/treeanalysis.R).


