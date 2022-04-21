#!/usr/bin/env python
"""
Created on Fri Apr 23 17:04:58 2021
Plots the boundary layer at different locations

@author: ruloz
"""

# import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
description = """ Plot the boundary layer."""


mainRoute = os.getcwd()  # get current directory route
dir_list = os.listdir(mainRoute)  # get list of directories inside case path
# files = os.listdir(mainRoute + str(time[0])) #get file names

times = []
for i in range(len(dir_list)):
    if(dir_list[i].isdigit()):
        times.append(dir_list[i])  # get time directories
    else:
        try:
            float(dir_list[i])
            times.append(dir_list[i])
        except:
            print("Not a float or int")
   # if(dir_list[i].isdigit()):
   #      times.append(dir_list[i]) # get time directories
timename = max(times)  # get latest time

postProcessFiles = os.listdir(mainRoute + '/postProcessing/singleGraph/' + timename)

# dfs = list()
# def loadData():
    # for file in postProcessFiles:
        # df0 = pd.read_csv(filename, skiprows=0)
        # dfs.append(df0)
        # print(dfs.head)

 # initialize arrays
y = []
Ux = []
Uy = []
Uz = []
U = [] 
u_shyy = []
y_shyy = []

def loadData():

    for file in postProcessFiles:
    
        sourceRoute = mainRoute + '/postProcessing/singleGraph/' + timename +'/' + file #route for each file

        #laod data
        ydata, UxData, UyData, UzData = np.loadtxt(sourceRoute, delimiter=',', skiprows=1, unpack = True)
        # genfromtxt(sourceRoute, delimiter=',',skip_header = 1)

        #store data in arrays
        uData = np.sqrt(UxData**2 + UyData**2 + UzData **2)/5.0
        y.append(np.array(ydata))
        Ux.append(np.array(UxData))
        Uy.append(np.array(UyData))
        Uz.append(np.array(UzData))
        U.append(np.array(uData))

    return U, y
    
def loadValidationData():

    sourceRoute = mainRoute + '/shyyData.csv'
    
    for i in range(0,8,2):
    
        u, y = np.loadtxt(sourceRoute, delimiter=',', skiprows=1, usecols = (i,i+1), unpack = True)
        
        y_shyy.append(np.array(y)/1e3)
        u_shyy.append(np.array(u))
        
    return u_shyy, y_shyy    


def plotBl():

    U, y = loadData()
    u_shyy, y_shyy = loadValidationData()
    
    fig, ax = plt.subplots()
    labels = list(map(lambda sub:int(''.join(
    [ele for ele in sub if ele.isnumeric()])), postProcessFiles))
    
    for i in range(len(U)):
    
        ax.plot(U[i],y[i],label='ST%s' %labels[i])
        ax.scatter(u_shyy[i],y_shyy[i],label='ST%s Shyy' %labels[i])
        # plt.plot(U[0],y[0],label='x = 0 mm')
        # plt.plot(U[1],y[1],label='x = 3 mm')
        # plt.plot(U[2],y[2],label='x = 6 mm')
        plt.xlabel('u [m/s]')
        plt.ylabel("y [m]")
        plt.legend()
        
    plt.savefig('validation'+timename+'.png', dpi=90)
    plt.show()

plotBl()
