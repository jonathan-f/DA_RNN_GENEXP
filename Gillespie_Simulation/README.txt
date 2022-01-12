README:

Author Michele Monti 2021

This code simulates a multilayer Gene -  Protein Network. The formation of RNA's and its activity is not taken into consideration. We are working in the adiabatic approximation
where RNA translation to Protein is much faster then the Gene-RNA transcription, such that to skip the RNA messenger step.


The Biology behind:
- The expression of a gene - so a creation of an RNA - is the product of the Hill functions of the different TF or RF (transcription or repressor Factors) proteins
- Each gene is activated or repressed by a certain Hill function of the protein that act on it
- it is possible to directly choose the average number of activators and repressors in the System
- to each gene corresponds a protein that is diluted with a constant rate in solution

The Algorithm and how it is implemented:

I use a standard Gillespie algorithm with an external input dynamic (not present but implemented).

The idea is that there is a network of Hill functions that are used to simulate each process:

This code does not track intermediate products (like complexes of Protein RNA or protein proteins or even RNAs) that are not functional for the expression of the genes.
All these ingredients though are present but hidden in the Hill functions.


The Algorithm is a standard Gillespie implemented with matrix. This makes the code very fast.
The interaction matrix are generated randomly and saved out into files fr a given simulation.
The possibility of reading the interaction matrix from an external txt file is also implemented.


HOW TO USE THE CODE

The code generates a random matrix of interactions Jij (or it can be set with a specific structure by the user), in the file GeneProteinsDynamic.c 3 parameters can be set:

- N: The number of proteins (genes) in the simulation
- nr: percentage of repressors per protein
- na: percentage of activators per protein

Given N, nr and na the code generates a random Jij.

k and mu are the activation and dilution rates, respectively.

To run the simulation:

- Go to the folder Gillespie_Simulations
- Compile by typing make in the terminal
- Run the simulation typing sh running.sh

This opens a folder called run (the name can be changed in the file running.sh). In the folder the matrix of interactions Jij and the protein time traces are saved as txt files

The python scripts provide some plotting functions to visualise the results.