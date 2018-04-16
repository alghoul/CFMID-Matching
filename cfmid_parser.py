# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 11:58:27 2018

@author: HALGhoul
"""
import os
import pandas as pd
import numpy as np

import CosineDotProduct_v1_3 as cpd
import pymysql as mysql
import time
#from sqlalchemy import create_engine
#cnx = create_engine('mysql://root:zahra_710@localhost/db')


dfc = None
dfn = None

Adduct_Mass = 1.007825



def parseCFMID(file=''):
    global dfc
    with open(file) as f:
        RESULT = list()
        for line_number, line in enumerate(f):
            if line_number > 450000:
                break
            #if line.startswith('# Name'):
            if line.startswith('# DTXCID'):
                result = {}
                titles = list()
                name = line.split(': ')[1]
                titles.append(name.rstrip())
                for line in f:
                    fields = line.split(': ')
                    #print fields
                    variable, value = fields
                    #if line.startswith('# InChI'):
                    if line.startswith('# RDMASS'):
                        titles.append(value.rstrip())
                        break
                    titles.append(value.rstrip())

            if line.startswith('PMASS'):
                e0_threads = list()
                e0_energy = list()
                e0_intensity = list()    
                e0 = line.split('\n')[0]
                for line in f:
                    if line.startswith('energy1'):
                        break
                    e0_threads.append(line.rstrip().split(' ')[0] + " | " + line.rstrip().split(' ')[1])
                    e0_energy.append(line.rstrip().split(' ')[0])
                    e0_intensity.append(line.rstrip().split(' ')[1])
                #result.append(e0_threads)
            if line.startswith('energy1'):
                e1_threads = list()
                e1 = line.split('\n')[0]
                for line in f:
                    if line.startswith('energy2'):
                        break
                    e1_threads.append(line.rstrip().split(' ')[0] + " | " + line.rstrip().split(' ')[1])
                #print e1_threads                
            if line.startswith('energy2'):
                e2_threads = list()
                e2 = line.split('\n')[0]
                for line in f:
                    if not line.strip():
                        break
                    e2_threads.append(line.rstrip().split(' ')[0] + " | " + line.rstrip().split(' ')[1])
                #result.append(e2_threads)
            if line.startswith('0'):
                msms_threads = list()
                msms_threads1 = line.rstrip()
                msms_threads.append(msms_threads1)
                for line in f:
                    if not line.strip():
                        break
                    msms_threads.append(line.rstrip())
                #result.append([('titles',titles), ('E0',e0_threads), ('E1',e1_threads), ('E2',e2_threads),('MSMS', msms_threads)])
                #result = { 'NAME':[titles[0]]*len(e0_threads), 'DTXCID':[titles[1]]*len(e0_threads), 'SMILES':[titles[2]]*len(e0_threads), 
                #    'CASRN':[titles[3]]*len(e0_threads), 'FORMULA':[titles[4]]*len(e0_threads), 'MASS':[titles[5]]*len(e0_threads), 'INCHI_KEY':[titles[6]],
                #    'PMASS':e0_energy, 'INTENSITY':e0_intensity, 'ENERGY1':e1_threads, 'ENERGY2':e2_threads}
                result = {'DTXCID':[titles[0]]*len(e0_threads), 'SMILES':[titles[1]]*len(e0_threads), 'MASS':[titles[2]]*len(e0_threads),
                    'PMASS':e0_energy, 'INTENSITY':e0_intensity, 'ENERGY1':e1_threads, 'ENERGY2':e2_threads}                
                RESULT.append(result)
    categories = [ "RUN %s" %i for i in range(0,len(RESULT))]
    dfc = pd.concat([pd.DataFrame.from_dict(d,orient='index').transpose() for d in RESULT],keys=categories)
    #dfc = dfc[['CASRN','MASS','PMASS','INTENSITY']].reset_index()
    dfc = dfc[['DTXCID','MASS','PMASS','INTENSITY']].reset_index()
    dfc['PMASS'] = dfc['PMASS'].astype('float64')
    dfc['MASS'] = dfc['MASS'].astype('float64')
    dfc['INTENSITY'] = dfc['INTENSITY'].astype('float64')
    dfc['INTENSITY0C'] = dfc.groupby(['MASS','DTXCID'])['INTENSITY'].apply(lambda x: (x/x.max())*100.0)
    dfc.loc[dfc['INTENSITY0C']>100,'INTENSITY0C']=0
    dfc.sort_values(['DTXCID','INTENSITY0C'],ascending=[True,False],inplace=True)   
    #print dfc
    #dfc.to_csv("cfmid.csv",index=False)
    #dfc.to_sql("CFMIDdb",con=cnx,if_exists='replace',index=False,chunksize=10000)    
    return dfc
