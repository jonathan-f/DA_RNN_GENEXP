# 07/2021 Michele Monti & Jonathan Fiorentino


#!/usr/bin/env python
# coding: utf-8

# # Discrimination of Gene Regulatory Nets architecture using a DA-RNN


import torch
import numpy as np
from matplotlib import pyplot as plt
import collections
import typing
from torch import nn
from torch.autograd import Variable
from torch.nn import functional as tf
import re

import logging
import os
import pickle

import typing
from typing import Tuple
import json

from torch import optim
import torch.multiprocessing as mp

from sklearn.preprocessing import StandardScaler
import joblib

import pandas as pd

from My_allFunctions import*

import time

# We consider time series data obtained through the Gillespie algorithm from gene regulatory networks with different architectures.
# 
# For each GRN, we train a DA-RNN for each target gene and save all the input attention vectors. We collect them in a matrix in which each row is the input attention vector for a gene. 
# 
# In the end we get a list of these matrices for all the GRNs used.

def findInteraction(attn,n_epochs, plotAll = False):
	att= []
	for i in attn:
		att.append(np.array(np.matrix((i.detach().numpy())).sum(axis=0)))
	totAtt=np.zeros(np.shape(att[0]))
	batch_x_epoch=int(len(attn)/n_epochs)
	for i in att[-batch_x_epoch:]:
		totAtt+=i[0]
	if plotAll:
		for i in att:
			plt.plot(i[0])
		plt.figure(2)
		plt.plot(totAtt[0]/batch_x_epoch)
	return totAtt/batch_x_epoch;

def SaveNetParams(folderName, model,prep,cols,scaler,att):
	da_rnn_kwargs = {"batch_size": 128, "T": 50}
	targ_cols = cols
	n_targs=len(targ_cols)
	encoder_hidden_size=64
	decoder_hidden_size=64
	T=50
	enc_kwargs = {"input_size": prep.feats.shape[1], "hidden_size": encoder_hidden_size, "T": T}
	with open(os.path.join(folderName, "enc_kwargs.json"), "w") as fi:
		json.dump(enc_kwargs, fi, indent=4)
	dec_kwargs = {"encoder_hidden_size": encoder_hidden_size,"decoder_hidden_size": decoder_hidden_size, "T": T, "out_feats": n_targs}
	decoder = Decoder(**dec_kwargs).to(device)
	with open(os.path.join(folderName, "dec_kwargs.json"), "w") as fi:
		json.dump(dec_kwargs, fi, indent=4)
	with open(os.path.join(folderName, "da_rnn_kwargs.json"), "w") as fi:
		json.dump(da_rnn_kwargs, fi, indent=4)
	joblib.dump(scaler, os.path.join(folderName, "scaler.pkl"))
	torch.save(model.encoder.state_dict(), os.path.join(folderName, "encoder.torch"))
	torch.save(model.decoder.state_dict(), os.path.join(folderName, "decoder.torch"))
	filehandler = open(folderName+"attn_list.att", 'wb') 
	pickle.dump(att, filehandler)		

def train_target(target,csvFile,cols,da_rnn_kwargs,sub,n_epochs):
	target=[target,]
	raw_data = pd.read_csv(csvFile, nrows=5000, usecols = cols)
	# Data preprocessing
	prep,scaler= my_preprocess_data(raw_data, target)
	
	#Train the DA-RNN net
	config, model = da_rnn(prep, n_targs=len(target), learning_rate=.001, **da_rnn_kwargs)
	iter_loss, epoch_loss, attn = train(model, prep, config, n_epochs=n_epochs, save_plots=False)
	
	mystring=re.sub(r'^.*?Protein', 'Protein', csvFile)
	mystring=mystring.replace('csv','')
	mystring=mystring.replace('.','')
	mystring=mystring.replace('/','_')
	if os.path.isdir(sub+'/'+mystring+'/')==False:
		os.mkdir(sub+'/'+mystring+'/')
	net_targ_dir=sub+'/'+mystring+'/'+target[0]+'/'
	
	if os.path.isdir(net_targ_dir)==False:
		os.mkdir(net_targ_dir)
	SaveNetParams(net_targ_dir, model,prep,cols,scaler,attn)
	attTot = findInteraction(attn,n_epochs, plotAll = False)
	
	return attTot;
    
FileListtxt=[]
FileListcsv=[]
subList=[]

# Get the names of the subfolders in DATA and the file with the protein time traces
#subfolders=[x[0] for x in os.walk('./DATA/')][1:]

#subfolders=['./DATA/FullyConnected',
#'./DATA/OscillatingNetworks_0123Clock',
#'./DATA/mediumConnection',
#'./DATA/MasterRegulator_gene1',
#'./DATA/FullyRepressed20Genes',
#'./DATA/ExternalSignal',
#'./DATA/SparseConnection']

subfolders=['./DATA/SparseConnection']

for sub in subfolders:
	print('----------- %s ------------' % (sub))
	for file in os.listdir(sub):
		if 'Protein' in file and 'txt' in file:
			print(file)
			subList.append(sub)
			FileListtxt.append(sub+'/'+file)
			FileListcsv.append(sub+'/'+file.replace('txt','csv'))

n_epochs=50
n_cores=1

print('Number of epochs: %d' % (n_epochs))
print('Number of cores used: %d' % (n_cores))

da_rnn_kwargs = {"batch_size": 128, "T": 50}

# List of attention matrices for all the networks used
attMatCollection=[]

print('Number of GRN',len(FileListtxt))

for i in range(len(FileListtxt)):
	print('Training GRN %d on each target' %(i+1))
	
	a=time.time()
	txtFile=FileListtxt[i]
	csvFile=FileListcsv[i]
	sub=subList[i]
	
	timelist = np.loadtxt(txtFile)[:,0]
	dt = timelist[3]-timelist[2]
	
	_, cols = txt_to_csv(txtFile, csvFile, ncols = [0])
	
	# Define the attention matrix for one net
	# The j-th row is the input attention vector for the j-th target gene
	attMat=np.zeros((len(cols),len(cols)))
	X = [(tg,csvFile,cols,da_rnn_kwargs,sub,n_epochs) for tg in cols]
	
	results = []
	if n_cores==None or n_cores > mp.cpu_count():
		print(mp.cpu_count())
		pool = mp.Pool(mp.cpu_count())
	else:
		print(n_cores)
		pool = mp.Pool(n_cores) # If the user wants to use a specified number of cores
		
	results = pool.starmap_async(train_target, X)
	results = results.get()
	results = [res[0] for res in results]
	b=time.time()
	print('Total training time for GRN %d: %f s' % (i+1,b-a))
	print(type(results))
	attMat=np.array(results)
	attMatCollection.append(attMat)
	
	att_matt_dir='./DATA/Attention_Matrices/'
	
	if os.path.isdir(att_matt_dir)==False:
		os.mkdir(att_matt_dir)
	np.savetxt(att_matt_dir+'attMat_'+sub.replace('/','_')+re.sub(r'^.*?Protein', 'Protein', txtFile),np.c_[attMat],fmt='%f',delimiter='\t')
