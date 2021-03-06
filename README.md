# Dual Attention Recurrent Neural Network for the prediction of gene expression time series data

This repository contains the code and data associated to the manuscript [Monti*, Fiorentino*, et al, Prediction of Time Series Gene Expression and Structural Analysis of Gene Regulatory Networks Using Recurrent Neural Networks, Entropy, 2022](https://www.mdpi.com/1454208). We generate synthetic gene expression time series data from a set of archetypical gene regulatory networks and we use a parallel Dual Attention Recurrent Neural Network (DA-RNN) for the time series prediction. Moreover, using tools from graph theory we show that graph descriptors of the input attention layer allow to distinguish different gene regulatory network architectures. We also study the effect on the predcition of Gaussian noise addition on protein concentrations. More details can be found in the associated manuscript; please cite it if you use our code and/or data.

# 1. Generation of synthetic gene expression time series data

Synthetic gene expression time series data are generated using the Gillespie algorithm. The code for performing the Gillespie simulations is provided in the folder [Gillespie_Simulation](/Gillespie_Simulation/), which also contains a [README](/Gillespie_Simulation/README.txt) file containing a description of the method and all the instructions needed to perform the simulation and visualize the results.

The time series of protein concentrations and the ground truth gene interaction matrices, for the gene regulatory netwroks used in our study, are provided in the [DATA](/DATA/) folder.

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

The python script [Parallel_DA_RNN_training.py](/Parallel_DA_RNN_training.py) allows to train a parallel DA-RNN on gene expression time series data for all the gene regulatory networks generated in section 1. The network parameters obtained for the training are saved in dedicated folders; also the matrix of the input attention is saved for downstream analyses (see below). The latter are also provided in the [Attention_Matrices](/DATA/Attention_Matrices) folder. The python script [My_allFunctions.py](/My_allFunctions.py) contains all the functions needed for training the DA-RNN and for the subsequent time series prediction.

The R script [MatComparison.R](/MatComparison.R) computes the graph theory descriptors (Clustering coefficient, Betweenness, Hub Score) of the input attention matrices for each gene regulatory network used in our study. The resulting matrices are provided in the folder [Graph_Descriptor_Matrices](/DATA/Graph_Descriptor_Matrices).

The analysis of Gaussian noise addition is provided in the python Jupyter notebook [DA_RNN_NOISE_ANALYSIS.ipynb](/DA_RNN_NOISE_ANALYSIS.ipynb); the results are in the [DATA](/DATA/) folder.

The Pearson correlation between the matrices obtained from the graph theory descriptors and from the noise analysis is computed using the R script [matrixsim.R](/matrixsim.R), with input data provided in the [matrices](/matrices/) folder.

The PCA and clustering analyses on the matrices mentioned above is performed in the Jupyter notebook [DA_RNN_DOWNSTREAM_ANALYSIS.ipynb](/DA_RNN_DOWNSTREAM_ANALYSIS.ipynb). The dendrograms associated to the different clusterings are provided in the [newick_trees](/newick_trees/) folder in Newick tree format. The distances between dendrograms are computed using the R script [treeanalysis.R](/treeanalysis.R).


