# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:11:13 2017

@author: Hussein Al Ghoul
"""

import os
import pandas as pd
import numpy as np

import CosineDotProduct_v1_2 as cpd
import pymysql as mysql
import time
#from sqlalchemy import create_engine
#cnx = create_engine('mysql://root:zahra_710@localhost/db')


dfc = None
dfn = None

Adduct_Mass = 1.007825


def read_NTA_data(file=''):  # read a csv file into a DataFrame
    global dfn
    dfn = pd.read_csv(file)
    #dfn.columns = dfn.columns.str.replace(' ','_')
    energy_column = 'Collision Energy A (Fragment Mass, Predicted Formula, Abundance, Abundance %)'
    regexI = "[\^,]*([0-9]+\.*[0-9]*)$" # a regex to find the score looking for a pattern of "=something_to_find ]" 
    regexE = "(.*?)\,{2}"
    dfn['INTENSITY0N'] = dfn[energy_column].str.extract(regexI,expand=True).astype('float64')
    dfn['PMASS'] = dfn[energy_column].str.extract(regexE,expand=True).astype('float64')
    dfn = dfn[['CASRN','PMASS','INTENSITY0N']]
    dfn.sort_values(['CASRN','INTENSITY0N'],ascending=[True,False],inplace=True)    
    dfn.to_csv("nta.csv",index=False)
    #print dfn
    return dfn





def parseMGF(file=''):
    NFile = file.rsplit('/',1)[-1]
    NewFile = NFile.rsplit('.',1)[0] + ".csv"
    with open(file) as f:
        RESULT = list()
        for line in f:
            if line.startswith('TITLE'):
                result = list()
                fields = line.split(' ')
                title, MS, of, pmass, charge, at, RT, mins, delimeter = fields
                #print (pmass, charge, RT)
                #result.append([pmass,RT])
            if line.startswith('RTINSECONDS'):
                RTS = line.split('=')[1]
                for line in f:
                    if line.split(' ')[0] == 'END':
                        break
                    a, b  = line.split('\t')
                    result.append([float(pmass), float(RT), charge, float(a),float(b)])
                RESULT.append(result)
        #print RESULT[0]
    categories = [ "RUN %s" %i for i in range(0,len(RESULT))]
    dfg = pd.concat([pd.DataFrame(d) for d in RESULT], keys=categories)
    dfg.columns = ["MASS", "RETENTION TIME", "CHARGE", "PMASS","INTENSITY"]
    dfg.sort_values(['MASS','RETENTION TIME'],ascending=[True,True],inplace=True) 
    df1 = dfg.reset_index()
    df1['TOTAL INTENSITY'] = df1.groupby(['MASS','RETENTION TIME'])['INTENSITY'].transform(sum)
    df1.sort_values(['MASS','TOTAL INTENSITY'],ascending=[True,True],inplace=True)
    df1 = df1.groupby('MASS').apply(lambda x: x[x['TOTAL INTENSITY'] == x['TOTAL INTENSITY'].max()])
    df1.loc[df1['PMASS']>df1['MASS'],'INTENSITY']=None
    df1.dropna(inplace=True)
    df1.sort_values(['MASS','INTENSITY'],ascending=[True,False],inplace=True)
    #df1['INTENSITY0M'] = df1.groupby('MASS')['INTENSITY'].apply(lambda x: (x/x.nlargest(2).min())*100.0)
    df1['INTENSITY0M'] = (df1['INTENSITY']/df1.groupby('MASS')['INTENSITY'].transform(np.max))*100.0
    df1.loc[df1['INTENSITY0M']>100,'INTENSITY0M']=None
    #df1.loc[df1['INTENSITY0M']<0.1,'INTENSITY0M']=None
    df1.reset_index(drop=True, inplace=True) # reset index
    df1.to_csv(NewFile,float_format='%.5f',index=False)
    '''
    x = dfmgf['ION MASS']
    y = dfmgf['INTENSITY']
    # intensity lines for every feature    
    lines = []
    for i in range(len(x)):
           pair=[(x[i],0), (x[i], y[i])]
        lines.append(pair)
    linecoll = matcoll.LineCollection(lines,linewidths=(2,2,2,2))
    ax = dfmgf.plot(x = 'ION MASS',y='INTENSITY',kind = 'scatter',s=1,figsize=(20,6))
    ax.add_collection(linecoll)
    plt.show()
    '''
    #print dfg
    dfg = df1
    #dfg.to_csv("CE10d_mgf.csv",index=False)
    return dfg


def spectrum_reader(file=''):
    dfg = pd.read_csv(file)
    dfg = dfg[(dfg['INTENSITY0M']<=100) & (dfg['INTENSITY0M']>0.01)]
    return dfg



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


def sqlCFMID(mass=0.0,ppm=0,mode=''):
    db = mysql.connect(host="host",
                   user="user",
                   passwd="pass",
                   db="db")
    cur = db.cursor()
    query= """select t1.dtxcid as DTXCID, t1.formula as FORMULA,t1.mass as MASS, t1.mz as PMASS, (t1.intensity/maxintensity)*100.0 as INTENSITY0C,
t1.energy as ENERGY 
from 
	(select c.*, p.* from peak p
		Inner Join job j on p.job_id=j.id
		Inner Join chemical c on j.dtxcid=c.dtxcid
		Inner Join spectra s on j.spectra_id = s.id
		Inner Join fragtool ft on j.fragtool_id=ft.id        
		where (abs(c.mass-"""+str(mass)+""")/"""+str(mass)+""")*1000000<"""+str(ppm)+ """
       and s.type='""" + mode + """') t1
left JOIN (select c.dtxcid, max(p.intensity) as maxintensity, p.energy from peak p
			Inner Join job j on p.job_id=j.id
			Inner Join chemical c on j.dtxcid=c.dtxcid
			Inner Join spectra s on j.spectra_id = s.id
			Inner Join fragtool ft on j.fragtool_id=ft.id
			where (abs(c.mass-"""+str(mass)+""")/"""+str(mass)+""")*1000000<"""+str(ppm)+ """
       and s.type='""" + mode + """'
			group by c.dtxcid, p.energy) t2
on t1.dtxcid=t2.dtxcid and t1.energy=t2.energy
order by DTXCID,ENERGY,INTENSITY0C desc;"""
    cur.execute(query)
    chunks=list()
    for chunk in pd.read_sql(query,db,chunksize=1000):
        chunks.append(chunk)
    #dfc.to_csv("cfmid.csv",index=False)
    return chunks
            



def compare_df(file='',ppm=10,ppm_sl=15,POSMODE=True):
    dfg = spectrum_reader(file)
    if POSMODE:
        mode='ESI-MSMS-pos'
        #CMass = Mass - Adduct_Mass
        dfg['MASS'] = dfg.groupby('RETENTION TIME')['MASS'].transform(lambda x: (x-1.007825))
        dfg['MASS'] = dfg['MASS'].round(6)
    else:
        mode='ESI-MSMS-neg'
        #CMass = Mass + Adduct_Mass
        dfg['MASS'] = dfg.groupby('RETENTION TIME')['MASS'].transform(lambda x: (x+1.007825)) 
        dfg['MASS'] = dfg['MASS'].round(6)
    #dfg.to_csv("dfg_aftercsv.csv",float_format='%.7f',index=False)  
    mass_list = dfg['MASS'].unique().tolist()
    mass_list = [196.077105]
    print mass_list
    dfAE_list=list()
    dfS_list=list()   
    for mass in mass_list:
        print("searching mass " + str(mass))
        dfcfmid = sqlCFMID(mass,ppm,mode)
        if not dfcfmid:
            print("No matches for this mass in CFMID library, consider changing the accuracy of the queried mass")
        else:    
            dfmgf = None
            df = None
            dfmgf = dfg[dfg['MASS'] == mass].reset_index()
            dfmgf.sort_values('MASS',ascending=True,inplace=True)
            df = cpd.Score(dfcfmid,dfmgf,mass,ppm_sl)
            dfAE_list.append(df[0]) #all energies scores
            dfS_list.append(df[1]) #sum of all energies score
            
    dfAE_total = pd.concat(dfAE_list) #all energies scores for all matches
    dfS_total = pd.concat(dfS_list) #Sum of scores for all matches
    th_dtxcid = dfS_list[0].at[0,'DTXCID'] #top hit dtxcid for plotting
    print(th_dtxcid) 
    dfcfm = pd.concat(dfcfmid)[(pd.concat(dfcfmid)['DTXCID'] == th_dtxcid) & (pd.concat(dfcfmid)['ENERGY'] == 'energy1')].reset_index()
    cpd.plot(dfcfm,dfmgf)
    print(dfAE_total)


    

#read_NTA_data('/home/hussein/Documents/NTA/Python_alt/ENTACT_DataReporting_EPA_MS2.csv')
#parseCFMID('/home/hussein/Documents/NTA/Python_alt/spectra_ESI-MSMS-neg_mass.dat')
#compare_df(183.057312)
#compare_df(183.058217)
#parseMGF(os.getcwd()+'/20180124_500_neg_MSMS_50ms.mgf') #<--Convert CSV to MGF
file= os.getcwd()+'/20180123_500_pos_MSMS_processed.csv'
#file= os.getcwd()+'/20180124_500_neg_MSMS_50ms.csv'   
t0=time.clock()
compare_df(file,10,15,True)
t1=time.clock()
print ("time to Process is " + str(t1-t0))
#sqlCFMID(144.081873,10)









