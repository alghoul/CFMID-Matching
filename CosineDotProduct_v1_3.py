# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 09:57:30 2018

@author: HALGhoul
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 10:11:13 2017

@author: Hussein Al Ghoul
"""

import os
import time
from matplotlib import pyplot as plt
from matplotlib import collections as matcoll
import pandas as pd
import numpy as np
from fractions import Fraction

m=0.6
n=3
def Commons(chunks=None,dfU=None,ppm_sl=0):
    print ("starting Commons")
    df_list=list()
    dfL_list=list()
    dfU['MASS_y'] = dfU['MASS'].round(6) 
    #dfU.rename(columns={'PMASS':'PMASS_y'},inplace=True) 
    dfU['MASS'] = dfU['MASS'].round(1)
    dfU['WEIGHTSM'] = (dfU['INTENSITY0M']**m)*(dfU['PMASS_y']**n)
    for chunk in chunks:
        df = None
        dfL = None
        dfL = chunk    
        dfL['MASS_x'] = dfL['MASS']
        #dfL.rename(columns={'PMASS':'PMASS_x'},inplace=True)   
        dfL['MASS'] = dfL['MASS'].round(1)
        dfL['WEIGHTSC'] = (dfL['INTENSITY0C']**m)*(dfL['PMASS_x']**n)
        dfL = dfL[(dfL['INTENSITY0C']<=100) & (dfL['INTENSITY0C']>0.01)]
        df = pd.merge(dfL,dfU,how='inner',on='MASS')  
        df['MATCHES'] = np.where((((abs(df.PMASS_x-df.PMASS_y)/df.PMASS_x)*1000000)<=ppm_sl),'1','0') 
        df.drop(df[df['MATCHES'] == '0'].index,inplace=True)
        df.sort_values(['DTXCID','ENERGY','PMASS_x','INTENSITY0C'],ascending=[True,True,True,False],inplace=True) 
        df_list.append(df)
        dfL_list.append(dfL)
    dft=pd.concat(df_list)
    dfLt=pd.concat(dfL_list)
    dfLt.to_csv("cfmid.csv",index=False)

    WLI = dfLt.groupby(['MASS_x','DTXCID','FORMULA','ENERGY'])['WEIGHTSC'].apply(list).to_dict()  
    #print WLI    
    WUI = dfU.groupby('MASS_y')['WEIGHTSM'].apply(list).to_dict() 
    print(WUI)
    #df.to_csv("Commons_Output.csv",index=False)
    #WL = ((df['INTENSITY0C']**m)*(df['PMASS_x']**n)).values.tolist()
    #WU = ((df['INTENSITY0M']**m)*(df['PMASS_y']**n)).values.tolist()
    WL = dft.groupby(['MASS_x','DTXCID','FORMULA','ENERGY'])['WEIGHTSC'].apply(list).to_dict()
    WU = dft.groupby(['MASS_x','DTXCID','FORMULA','ENERGY'])['WEIGHTSM'].apply(list).to_dict()
    print(len(WL))
    #print WUI
    W = list()
    W.append(WL)
    W.append(WU)
    W.append(WLI)
    W.append(WUI)
    return W

def FR(WL,WU):
    #print WL
    #print WU
    num =0.0
    den = 0.0
    SUM = 0.0
    for i in range(0,len(WL)):
        num = WL[i]*WU[i-1]
        den = WL[i-1]*WU[i]
        if (num/den) <= 1:
            l = 1
        else:
            l = -1
        SUM += (num/den)**l     
    F_R = (1.0/float(len(WL)))*SUM
    return F_R

def FD(WL,WU,WLI,WUI):
    #print WL
    #print WU
    SUMU = 0.0
    SUML = 0.0
    SUM = 0.0
    F_D = 0.0
    #print WUI
    for i in range(0,len(WUI)):
        #print WUI[i]
        SUMU += WUI[i]*WUI[i]
    #print SUMU
    for i in range(0,len(WLI)):
        SUML += WLI[i]*WLI[i]
    #print SUML
    for i in range(0,len(WL)):
        #print WU[i]
        SUM += WL[i]*WU[i]
    #print SUM
    F_D = (SUM*SUM)/(SUMU*SUML)
    #print F_D
    return F_D    
   
def Score(dfL=None,dfU=None,Mass=0.0,ppm_sl=0):
    DF=list()
    W = Commons(dfL,dfU,ppm_sl)
    WL=set(W[0])
    #print WL
    WLI=set(W[2])
    record = list()
    records = list()
    #print WLI
    for keys in WLI.intersection(WL):
        #print W[0][keys] 
        #print W[1][Mass]
        N_LU=0
        F_D=0.0
        F_R=0.0
        score=0.0
        N_LU = len(W[0][keys])
        N_U = len(W[3][Mass])
        F_D = FD(W[0][keys],W[1][keys],W[2][keys],W[3][Mass])
        F_R = FR(W[0][keys],W[1][keys])
        score = ((N_U*F_D) + (N_LU*F_R))/(N_U + N_LU)
        record = list(keys)
        record.append(F_D)
        record.append(score)
        records.append(record)
    #dfL_plot = dfL.loc[dfL['DTXCID'].isin([keys[1]])].reset_index()
    #plot(dfL_plot,dfU)
    df = pd.DataFrame.from_records(records,columns=['MASS','DTXCID','FORMULA','ENERGY','FD','SCORE'])
    df.sort_values(['ENERGY','SCORE'],ascending=[True,True],inplace=True)
    df['RANK'] = df.groupby(['FORMULA','ENERGY'])['SCORE'].rank(method='dense',ascending=True) # rank according to formula here by adding ['FORMULA','ENERGY']
    print (df)
    df.to_csv('Score_Alllevels.csv',index=False)
    dfs = df.groupby(['DTXCID','MASS','FORMULA'],as_index=False)[['SCORE','FD']].sum()
    dfs.reset_index()
    dfs['RANK'] = dfs.groupby(['MASS','FORMULA'])['SCORE'].rank(method='dense',ascending=True) # rank according to formula here by adding ['FORMULA','ENERGY']
    if not df.empty:
        dfs.sort_values('SCORE',ascending=False,inplace=True)    
    print (dfs)
    #dfs.to_csv("Score_Sum.csv",index=False)
    print ("Number of Matches: " + str(len(WL)))
    DF.append(df)
    DF.append(dfs)
    return DF


def plot(dfL=None,dfU=None):
    xc = dfL['PMASS_x']
    yc = dfL['INTENSITY0C']
    xn = dfU['PMASS_y']
    yn = dfU['INTENSITY0M']
    # intensity lines for every feature    
    linesc = []
    for i in range(len(xc)):
           pairc=[(xc[i],0), (xc[i], yc[i])]
           linesc.append(pairc)
    linecollc = matcoll.LineCollection(linesc,linewidths=(2,2,2,2))
    linesn = []
    for j in range(len(xn)):
           pairn=[(xn[j],0), (xn[j], yn[j])]
           #print pairn
           linesn.append(pairn)
    linecolln = matcoll.LineCollection(linesn,colors='red',linewidths=(2,2,2,2))

    # plot the cfmid spectrum and add lines for the acquired data
    ax = dfL.plot(x = 'PMASS_x',y='INTENSITY0C',kind = 'scatter',s=1,figsize=(20,6))
    ax.add_collection(linecollc)
    ax.add_collection(linecolln)    
    #plt.plot(x,y)
    ax.set_xlim(0.0,max(xc)*(1.1))
    ax.set_ylim(0,max(yc)*(1.1))
    plt.draw()    
    plt.show()















    
    
    
    
