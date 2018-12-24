# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 14:51:32 2018

@author: katsuhisa
"""

# libraries
import numpy as np
import matplotlib.pyplot as plt
import uneye


########### Parameters ###########
sampfreq = 500 #Hz
weights_name = 'weights_synthetic'
min_sacc_dur = 6 # in ms
min_sacc_dist = 10 #in ms

########### YOUR DATA ############

datapath = 'Z:/Katsuhisa/headfree_project/dataset/eyes_dis_csv/' #example data
x_filename = '2015.12.16_x.csv'
y_filename = '2015.12.16_y.csv'

##################################

# load data
Xtest = np.loadtxt(datapath+x_filename,delimiter=',')
Ytest = np.loadtxt(datapath+y_filename,delimiter=',')


# Prediction
model = uneye.DNN(weights_name=weights_name,
                 sampfreq=sampfreq,
                 min_sacc_dur=min_sacc_dur,
                min_sacc_dist=min_sacc_dist)
Prediction,Probability = model.predict(Xtest,Ytest)