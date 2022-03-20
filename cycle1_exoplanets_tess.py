#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 09:02:29 2021

@author: 
    
Plots updated to include which ones were discoverd by TESS that
will be observed in Cycle 1 JWST.
"""



import numpy as np
import requests
import matplotlib.pyplot as plt
import pandas as p
from IPython.display import display, HTML
import re



def get_planet_properties(planet_name):

    planeturl = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/"
    dvurl = "https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/"
    header={}
    
    if (planet_name[0:3] == 'TOI') & (re.match("b|c|d|e|f|g|h|j", planet_name[-1]) is not None):
        planet_name = planet_name.replace(' ', '-', 1)
     
    print(planet_name)
        
    url = planeturl + planet_name + "/properties/"

    r = requests.get(url = url, headers = header)

    planet_prop = r.json()
    
    #print(len(planet_prop))
    if len(planet_prop) == 2:
        return planet_prop[1]['Rp'], planet_prop[1]['Tp']
    elif len(planet_prop) == 1:
        return planet_prop[0]['Rp'], planet_prop[0]['Tp']
    else:
        return -1, -1
#%%
#Read in which objects will be 
ddir = "/Users/smullally/Python_Code/jwst_project/cycle1_science/"
jwst_files = ['cycle1_planets_transit.txt',
             'cycle1_planets_direct.txt',
             'gto_planets_transit.txt',
             'gto_planets_direct.txt']
jwst_types = ['GO','GO','GTO', 'GTO']
obs_types=['transit','direct','transit','direct']

previous_file = 'transitspec_2021-05-19_nexsci.csv'


all_dfs = []

for i,f in enumerate(jwst_files):
    
    df = p.read_csv(ddir+f, delimiter=",",names=['planet', 'inst'])
    df['calltype'] =jwst_types[i]
    df['obstype'] = obs_types[i]
    
    all_dfs.append(df)

jwstdf = p.concat(all_dfs,axis=0)

previousdf = p.read_csv(ddir+previous_file, comment="#")

unique_prev = previousdf[['plntname', 'instrument']].drop_duplicates()

#%%

Rp = []
Tp = []

for i,planet in jwstdf.iterrows():
    
    try:
        radii, temp = get_planet_properties(planet['planet'])
        #print(radii, temp)
    except:
        radii = -1
        temp = -1
        
    Rp.append(radii)
    Tp.append(temp)

jwstdf['Rp'] = Rp
jwstdf['Teq'] = Tp

jwstdf.to_csv(ddir + "jwst_cycle1_planet_prop.csv")

#%%
#Read in TESS List
tess_file = ddir + "tess_cycle1_planets.txt"

tess_list = p.read_csv(tess_file, comment="#", names=['planet'])

tessdf = p.merge(tess_list, jwstdf, how='inner', right_on="planet", left_on="planet")

#%%
#Print out cycle 1 transiting planet properites
Rp =[]
Tp = []
for i,planet in unique_prev.iterrows():
     
    try:
        radii, temp = get_planet_properties(planet['plntname'])
        print(radii, temp)
    except:
        radii = -1
        temp = -1
        
    Rp.append(radii)
    Tp.append(temp)

unique_prev['Rp'] = Rp
unique_prev['Teq'] = Tp

unique_prev.to_csv(ddir + "transitspec_planet_prop.csv")

#%%
inst = ['Wide Field Camera 3', 
        'Infared Array Camera (IRAC)', 
        'Space Telescope Imaging Spectrograph',
        'WFC3', 'IRAC', 'WFC',
        'Near Infrared Camera and Multi-Object Spectrometer'
        ]
aprev = unique_prev.copy()
cond = list(map(lambda x : x in inst, aprev['instrument']))
prev = aprev[cond]
go1_transit = jwstdf[(jwstdf['calltype']=="GO") & (jwstdf['obstype'] == "transit")]
#go1_direct = jwstdf[(jwstdf['calltype']=="GO") & (jwstdf['obstype'] == "direct")]
gto_transit = jwstdf[(jwstdf['calltype']=="GTO") & (jwstdf['obstype'] == "transit")]
#gto_direct = jwstdf[(jwstdf['calltype']=="GTO") & (jwstdf['obstype'] == "direct")]

go1 = jwstdf[(jwstdf['calltype']=="GO")]
gto = jwstdf[(jwstdf['calltype']=="GTO")]

Mj = 11.2089
#Create plots
plt.style.use('dark_background')

plt.figure(figsize=(8,8))
fnt=15



plt.plot(go1_transit['Teq'], go1_transit['Rp']*Mj, 's', color="darkorange",
             fillstyle='full', label="JWST GO1 Transit", ms=10, alpha=0.8)
#plt.plot(go1_direct['Teq'], go1_direct['Rp']*Mj, 's', color="darkorange",
#             fillstyle='none', label="JWST GO1 Direct", ms=10)

plt.plot(gto_transit['Teq'], gto_transit['Rp']*Mj, 'ro', fillstyle='full', 
         label="JWST GTO Transit", ms=8, alpha=0.7)
#plt.plot(gto_direct['Teq'], gto_direct['Rp']*Mj, 'ro', fillstyle='none', 
#         label="JWST GTO Direct", ms=8)
plt.plot(tessdf['Teq'], tessdf['Rp']*Mj, 's', color="white",
             fillstyle='none', label="TESS Discoveries", ms=10, alpha=1)

plt.plot(prev['Teq'], prev['Rp']*Mj,'w.', ms=5, label="HST/Spitzer < 2021")
plt.yscale('log', base=10)
plt.xscale('log', base=10)


plt.legend(fontsize=fnt-3)
xrange = [4000,51]
plt.xlim(xrange[0], xrange[1])
plt.ylim(.7, 31)
plt.xlabel('Equilibrium Temperature (K)', fontsize=fnt)
plt.ylabel('Planet Radius (Earth)', fontsize=fnt)
xpos = [1000, 300,100]
xlab = ["1000","300", "100"]
plt.xticks(xpos, xlab, fontsize=fnt)
ypos = [1,3,10,30]
ylab = ["1", "3", "10", "30"]
plt.yticks(ypos,ylab, fontsize=fnt)


plt.title("JWST Cycle 1 Transiting Exoplanets", fontsize=fnt)

#annotate

plt.plot(219, 1, 'o', color='cyan', alpha=1, ms=12)
plt.annotate('Earth', (205,.96), color='cyan', fontsize=fnt)

plt.plot(110, 11.2, 'o', color='khaki', alpha=1, ms=12)
plt.annotate('Jupiter', (102,11.1), color='khaki', fontsize=fnt)

plt.savefig(ddir+"cycle1_planets_tess.png")