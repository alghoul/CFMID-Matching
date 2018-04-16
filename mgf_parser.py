# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:11:13 2017

@author: Hussein Al Ghoul
"""

import os
import pandas as pd
import numpy as np

import CosineDotProduct_v1_3 as cpd
import pymysql as mysql
import time



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
    dfg.columns = ["MASS", "RETENTION TIME", "CHARGE", "PMASS_y","INTENSITY"]
    dfg.sort_values(['MASS','RETENTION TIME'],ascending=[True,True],inplace=True) 
    df1 = dfg.reset_index()
    df1['TOTAL INTENSITY'] = df1.groupby(['MASS','RETENTION TIME'])['INTENSITY'].transform(sum)
    df1.sort_values(['MASS','TOTAL INTENSITY'],ascending=[True,True],inplace=True)
    df1 = df1.groupby('MASS').apply(lambda x: x[x['TOTAL INTENSITY'] == x['TOTAL INTENSITY'].max()])
    df1.loc[df1['PMASS_y']>df1['MASS'],'INTENSITY']=None
    df1.dropna(inplace=True)
    df1.sort_values(['MASS','INTENSITY'],ascending=[True,False],inplace=True)
    #df1['INTENSITY0M'] = df1.groupby('MASS')['INTENSITY'].apply(lambda x: (x/x.nlargest(2).min())*100.0)
    df1['INTENSITY0M'] = (df1['INTENSITY']/df1.groupby('MASS')['INTENSITY'].transform(np.max))*100.0
    df1.loc[df1['INTENSITY0M']>100,'INTENSITY0M']=None
    #df1.loc[df1['INTENSITY0M']<0.1,'INTENSITY0M']=None
    df1.reset_index(drop=True, inplace=True) # reset index
    df1.to_csv(NewFile,float_format='%.5f',index=False)
    dfg = df1
    #dfg.to_csv("CE10d_mgf.csv",index=False)
    return dfg



def spectrum_reader(file=''):
    dfg = pd.read_csv(file)
    dfg = dfg[(dfg['INTENSITY0M']<=100) & (dfg['INTENSITY0M']>0.01)]
    return dfg



def sqlCFMID(mass=0.0,ppm=0,mode=''):
    db = mysql.connect(host="host",
                   user="user",
                   passwd="pass",
                   db="db")
    cur = db.cursor()
    query= """select t1.dtxcid as DTXCID, t1.formula as FORMULA,t1.mass as MASS, t1.mz as PMASS_x, (t1.intensity/maxintensity)*100.0 as INTENSITY0C,
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
    return chunks
            

def list_maker(fpcdl='',dfg=None):
    # make a list of the MGF masses corresponding to the PCL Monoisotopic masses 
    dfpcdl = pd.read_csv(fpcdl)
    dfg['Mass_rounded'] = dfg['MASS'].round(1)
    dfpcdl['Mass_rounded'] = dfpcdl['Neutral Monoisotopic Mass'].round(1)
    df = pd.merge(dfpcdl,dfg,how='left',on='Mass_rounded') 
    df['MATCHES'] = np.where((((abs(df['Neutral Monoisotopic Mass']-df['MASS'])/df['MASS'])*1000000)<=15),'1','0') 
    df.drop(df[df['MATCHES'] == '0'].index,inplace=True)
    df=df[np.isfinite(df['MASS'])]
    mlist = df['MASS'].unique().tolist()
    print mlist
    return mlist


def compare_df(file='',fpcdl='',ppm=10,ppm_sl=15,POSMODE=True):
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
    #mass_list = dfg['MASS'].unique().tolist()
    #mass_list = [312.184525]
    mass_list = list_maker(fpcdl,dfg)
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
    #th_dtxcid = dfAE_list[0].at[0,'DTXCID'] #top hit dtxcid for plotting
    #print(th_dtxcid) 
    #dfcfm = pd.concat(dfcfmid)[(pd.concat(dfcfmid)['DTXCID'] == th_dtxcid) & (pd.concat(dfcfmid)['ENERGY'] == 'energy2')].reset_index()
    #cpd.plot(dfcfm,dfmgf)
    df_result = merge_pcdl(fpcdl,dfS_total)
    print(df_result)
    NFile = fpcdl.rsplit('/',1)[-1]
    NewFile = NFile.rsplit('.',1)[0] + "_compared_with_CFMID_onescore_wformula.csv"
    df_result.to_csv(NewFile,index=False)



def merge_pcdl(fpcdl='',df=None):
    
    dfpcdl = pd.read_csv(fpcdl)
    df = pd.merge(dfpcdl,df,how='left',on='DTXCID')
    return df

#compare_df(183.057312)
#compare_df(183.058217)
#parseMGF(os.getcwd()+'/20180124_500_neg_MSMS_50ms.mgf') #<--Convert CSV to MGF
file = os.getcwd()+'/20180123_500_pos_MSMS_processed.csv'
fpcdl = os.getcwd()+'/PCDL_Matches_500_v2_FINAL_CID.csv'
#file= os.getcwd()+'/20180124_500_neg_MSMS_50ms.csv'   
t0=time.clock()
compare_df(file,fpcdl,10,15,True)
t1=time.clock()
print ("time to Process is " + str(t1-t0))
#sqlCFMID(144.081873,10)









