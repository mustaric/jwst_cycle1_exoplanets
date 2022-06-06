#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 09:02:29 2021

@author: smullally
"""

import numpy as np
import requests
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display, HTML
import re

planeturl = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/"
dvurl = "https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/"
header={}

def get_planet_properties(planet_name):
    
    if (planet_name[0:3] == 'TOI') & (re.match("b|c|d|e|f|g|h|j", planet_name[-1]) is not None):
        planet_name = planet_name.replace(' ', '-', 1)
     
    print(planet_name)
        
    url = planeturl + planet_name + "/properties/"

    r = requests.get(url = url, headers = header)

    planet_prop = r.json()
    
    print(len(planet_prop))
    if len(planet_prop) == 2:
        return planet_prop[1]
    elif len(planet_prop) == 1:
        return planet_prop[0]
    else:
        return {'planet_name' : planet_name, 'catalog_name' : 'not found'}
#%%
#Create two dictionaries
# one for the jwst planets
# one for the planets previously observed by HST/Spitzer
ddir = "/Users/smullally/Python_Code/jwst_project/cycle1_science/jwst_cycle1_exoplanets/"
jwst_files = ['cycle1_planets_transit.txt',
             'cycle1_planets_direct.txt',
             'gto_planets_transit.txt',
             'gto_planets_direct.txt']
jwst_types = ['GO','GO','GTO', 'GTO']
obs_types=['transit','direct','transit','direct']

previous_file = 'transitspec_2021-05-19_nexsci.csv'


all_dfs = []

for i,f in enumerate(jwst_files):
    
    df = pd.read_csv(ddir+f, delimiter=",",names=['planet', 'inst'])
    df['calltype'] =jwst_types[i]
    df['obstype'] = obs_types[i]
    
    all_dfs.append(df)

jwstdf = pd.concat(all_dfs,axis=0)

previousdf = pd.read_csv(ddir+previous_file, comment="#")

unique_prev = previousdf[['plntname', 'instrument']].drop_duplicates()

#%%
#Gather information about the JWST targets
params = ['Rp', 'Tp', 'Mp','Teff', 'Kmag', 'orbital_period', 'Rs','a/Rs', 'catalog_name', 'planet_name']

planetdf  = pd.DataFrame(columns = params)

for i,planet in jwstdf.iterrows():
    
    try:
        planetdict = get_planet_properties(planet['planet'])
        shortdf = pd.DataFrame([planetdict], columns = params)
        print(shortdf)
    except:
        shortdf = pd.DataFrame([{'planet_name' : planet['planet'], 'catalog_name' : 'error'}])
    
    planetdf = pd.concat([planetdf, shortdf])    
    
final_jwstdf = pd.merge(jwstdf, planetdf, how='left', left_on='planet', 
                        right_on='planet_name').drop_duplicates()

final_jwstdf.to_csv(ddir + "jwst_cycle1_planet_prop_tess.csv")

#%%
params = ['Rp', 'Tp', 'Mp','Teff', 'Kmag', 'orbital_period', 'Rs','a/Rs', 'catalog_name', 'planet_name']

prevdf  = pd.DataFrame(columns = params)

for i,planet in unique_prev.iterrows():
     
    try:
        prevdict = get_planet_properties(planet['plntname'])
        shortdf = pd.DataFrame([prevdict], columns = params)
        print(shortdf)
        
    except:
        shortdf = pd.DataFrame([{'planet_name' : planet['plntname'], 'catalog_name' : 'error'}])
        
        
    prevdf = pd.concat([prevdf, shortdf])
    
final_uniqueprev = pd.merge(unique_prev, prevdf, how='left', 
                        left_on='plntname', right_on='planet_name').drop_duplicates()
final_uniqueprev.to_csv(ddir + "transitspec_planet_prop_tess.csv")


#%%
#Read in my tables, which were created above
#Adjust the parameters for one that is known to be in correct in exo.mast
final_uniqueprev = pd.read_csv(ddir+"transitspec_planet_prop_tess.csv")
final_jwstdf = pd.read_csv(ddir + "jwst_cycle1_planet_prop_tess.csv")

final_jwstdf['a/Rs'].iloc[47] = 11.7
final_jwstdf['Rs'].iloc[47] = 1.38
final_jwstdf['Teff'].iloc[47] = 5675

#%%


def flux_fromTp(Tp):
    boltzm = 5.670374419e-8  #Wm-2K-4
    try:
        I = 4.0 * boltzm * Tp**4
        return I
   
    except TypeError:
        return -1

#Add in instellation flux

Io = list( map( lambda x : flux_fromTp(x), final_jwstdf['Tp']))
final_jwstdf['Io'] = Io

Io = list( map( lambda x : flux_fromTp(x), final_uniqueprev['Tp']))
final_uniqueprev['Io'] = Io
#%%
def get_incidentFlux(aRstar,Teff):
    
    F = (aRstar*.0046505)**(-2) * Teff**4 / (5778)**4
    
    return F

inFlux = list(map (lambda x,y: get_incidentFlux(x,y), final_jwstdf['a/Rs'], final_jwstdf['Teff']))
final_jwstdf['inFlux'] = inFlux

inFlux = list( map( lambda x,y : get_incidentFlux(x,y), final_uniqueprev['a/Rs'], final_uniqueprev['Teff']))
final_uniqueprev['inFlux'] = inFlux


#%%
#Here is where we do the plotting.
#Rp vs Teq 

inst = ['Wide Field Camera 3', 
        'Infared Array Camera (IRAC)', 
        'Space Telescope Imaging Spectrograph',
        'WFC3', 'IRAC', 'WFC',
        'Near Infrared Camera and Multi-Object Spectrometer'
        ]

jwstdf = final_jwstdf
aprev = final_uniqueprev.copy()
cond = list(map(lambda x : x in inst, aprev['instrument']))
prev = aprev[cond]

#%%


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

plt.plot(prev['Tp'], prev['Rp']*Mj,'w.', ms=5, label="HST/Spitzer <2021")
plt.yscale('log', base=10)
plt.xscale('log', base=10)


plt.plot(gto_transit['Tp'], gto_transit['Rp']*Mj, 'r^', fillstyle='full', 
         label="JWST GTO Transit", ms=9, alpha=0.7)
#plt.plot(gto_direct['Teq'], gto_direct['Rp']*Mj, 'ro', fillstyle='none', 
#         label="JWST GTO Direct", ms=8)
plt.plot(go1_transit['Tp'], go1_transit['Rp']*Mj, 'h', color="darkorange",
             fillstyle='full', label="JWST GO1 Transit", ms=8, alpha=0.9)
#plt.plot(go1_direct['Teq'], go1_direct['Rp']*Mj, 's', color="darkorange",
#             fillstyle='none', label="JWST GO1 Direct", ms=10)

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


#%%
#Here is where we do the plotting.
# Mark TESS Planets

def is_tess(planet):
    tessplanets = ['HD 15337 c', 'L168-9 c', 'L98-59 c', 'LTT 1445 A b', 'TOI 260.01', 
                   'TOI 776 b', 'HD 15337 b',
                   'TOI 776 c', 'TOI 836.01', 'TOI 836.02', 'LTT 1445 A b',
                   'TOI 178 b', 'TOI 178 d','TOI 178 g','TOI 421 b','TOI 421 b',
                   'LHS 3844 b', 'WD 1856+534 b','HIP 67522 b', 'GJ 486 b',
                   'TOI 741.01', 'GJ 357 b', 'TOI 910.01',
                   'GJ 367 b', 'L 98-59 d', 'LP 791-18 c', 
                   'LTT 9779 b']
    
    if planet in tessplanets:
        return True
    else: 
        return False

tesspl = list( map( lambda x : is_tess(x), jwstdf['planet_name']))
jwstdf = final_jwstdf
jwstdf['tess'] = tesspl



#%%

go1_transit = jwstdf[(jwstdf['calltype']=="GO") & (jwstdf['obstype'] == "transit")]

gto_transit = jwstdf[(jwstdf['calltype']=="GTO") & (jwstdf['obstype'] == "transit")]

tess_transit = jwstdf[(jwstdf['tess']) & (jwstdf['obstype'] == 'transit')]

go1 = jwstdf[(jwstdf['calltype']=="GO")]
gto = jwstdf[(jwstdf['calltype']=="GTO")]

Mj = 11.2089
#Create plots
plt.style.use('dark_background')

plt.figure(figsize=(8,8))
fnt=15

plt.plot(gto_transit['Tp'], gto_transit['Rp']*Mj, 'r^', fillstyle='full', 
         label="JWST GTO Transit", ms=9, alpha=0.7)

plt.plot(go1_transit['Tp'], go1_transit['Rp']*Mj, 'h', color="darkorange",
             fillstyle='full', label="JWST GO1 Transit", ms=8, alpha=0.9)

plt.plot(tess_transit['Tp'], tess_transit['Rp']*Mj, 'D', color="white", 
             fillstyle='none', ms=10, label="TESS Planet")

plt.yscale('log', base=10)
plt.xscale('log', base=10)
plt.legend(fontsize=fnt-3)
xrange = [4000,51]
plt.xlim(xrange[0], xrange[1])
plt.ylim(.6, 31)
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

plt.savefig(ddir+"cycle1_planets_tesspl.png")

#%%
#Teq vs Radius and stellar temperature marking the TESS planets.
#if set ot True it will mark TESS planets
mark_tess = True

transit = jwstdf[ (jwstdf['obstype'] == "transit")]

tess_transit = jwstdf[(jwstdf['tess']) & (jwstdf['obstype'] == 'transit')]

#go1 = jwstdf[(jwstdf['calltype']=="GO")]
#gto = jwstdf[(jwstdf['calltype']=="GTO")]

Mj = 11.2089

#Create plots
plt.style.use('dark_background')

plt.figure(figsize=(9.3,8))
fnt=16


plt.scatter(transit['Tp'], transit['Rp']*Mj, marker='o', c=transit['Teff'], cmap='RdYlBu',
              label="JWST Transit", alpha=1, s=47)
cb = plt.colorbar(label='Stellar Teff (K)')
cb.set_label(label="Stellar Teff (K)", size=fnt)
cb.ax.tick_params(labelsize=fnt+1)

if mark_tess:
    plt.plot(tess_transit['Tp'], tess_transit['Rp']*Mj, 'D', color="white", 
             fillstyle='none', ms=11, label="TESS Planet",lw=2)

plt.yscale('log', base=10)
plt.xscale('log', base=10)
plt.legend(fontsize=fnt-3)
xrange = [4000,51]
plt.xlim(xrange[0], xrange[1])
plt.ylim(.6, 31)
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

if mark_tess:
    plt.savefig(ddir+"cycle1_planets_tess_star.png")
else:
    plt.savefig(ddir+"cycle1_planets_star.png")

#%%
#In Incident Flux and TESS planets marked.

mark_tess = False

transit = jwstdf[ (jwstdf['obstype'] == "transit")]

tess_transit = jwstdf[(jwstdf['tess']) & (jwstdf['obstype'] == 'transit')]

#go1 = jwstdf[(jwstdf['calltype']=="GO")]
#gto = jwstdf[(jwstdf['calltype']=="GTO")]

Mj = 11.2089

#Create plots
plt.style.use('dark_background')

plt.figure(figsize=(11,7))
fnt=15


plt.scatter(transit['inFlux'], transit['Rp']*Mj, marker='o', c=transit['Teff'], 
            cmap='RdYlBu',
            label="JWST Transit", alpha=1, s=48)
#plt.colorbar(label='Stellar Teff (K)')
cb = plt.colorbar(label='Stellar Teff (K)')
cb.set_label(label="Stellar Teff (K)", size=fnt)
cb.ax.tick_params(labelsize=fnt+1)

if mark_tess:
    plt.plot(tess_transit['inFlux'], tess_transit['Rp']*Mj, 'D', color="white", 
             fillstyle='none', ms=12, label="TESS Planet")

plt.yscale('log', base=10)
plt.xscale('log', base=10)
plt.legend(fontsize=fnt-3)
xrange = [10000,.01]
plt.xlim(xrange[0], xrange[1])
plt.ylim(.6, 31)
plt.xlabel('Incident Flux (Earth)', fontsize=fnt)
plt.ylabel('Planet Radius (Earth)', fontsize=fnt)
xpos = [1000, 300,100, 30, 1, .3, .1]
xlab = ["1000","300", "100","30", "1", "0.3", ".1"]
plt.xticks(xpos, xlab, fontsize=fnt)
ypos = [1,3,10,30]
ylab = ["1", "3", "10", "30"]
plt.yticks(ypos,ylab, fontsize=fnt)


plt.title("JWST Cycle 1 Transiting Exoplanets", fontsize=fnt)

#annotate

plt.plot(1, 1, 'o', color='cyan', alpha=1, ms=12)
plt.annotate('Earth', (1.1,0.83), color='cyan', fontsize=fnt)

plt.plot(.037, 11.2, 'o', color='khaki', alpha=1, ms=12)
plt.annotate('Jupiter', (.05,9.0), color='khaki', fontsize=fnt)

if mark_tess:
    plt.savefig(ddir+"cycle1_planets_tess_influx.png")
else:
    plt.savefig(ddir+"cycle1_planets_influx.png")
#%%
#In Incident Flux and wihtout TESS planets

transit = jwstdf[ (jwstdf['obstype'] == "transit")]

tess_transit = jwstdf[(jwstdf['tess']) & (jwstdf['obstype'] == 'transit')]

#go1 = jwstdf[(jwstdf['calltype']=="GO")]
#gto = jwstdf[(jwstdf['calltype']=="GTO")]

Mj = 11.2089

#Create plots
plt.style.use('dark_background')

plt.figure(figsize=(11,7))
fnt=16


plt.scatter(transit['inFlux'], transit['Rp']*Mj, marker='o', c=transit['Teff'], cmap='RdYlBu',
              label="JWST Transit", alpha=1, s=45)
plt.colorbar(label='Stellar Teff (K)')

#plt.plot(tess_transit['inFlux'], tess_transit['Rp']*Mj, 'D', color="white", 
#             fillstyle='none', ms=10, label="TESS Planet")

plt.yscale('log', base=10)
plt.xscale('log', base=10)
plt.legend(fontsize=fnt-3)
xrange = [10000,.01]
plt.xlim(xrange[0], xrange[1])
plt.ylim(.6, 31)
plt.xlabel('Incident Flux (Earth)', fontsize=fnt)
plt.ylabel('Planet Radius (Earth)', fontsize=fnt)
xpos = [1000, 300,100, 30, 1, .3, .1]
xlab = ["1000","300", "100","30", "1", "0.3", ".1"]
plt.xticks(xpos, xlab, fontsize=fnt)
ypos = [1,3,10,30]
ylab = ["1", "3", "10", "30"]
plt.yticks(ypos,ylab, fontsize=fnt)


plt.title("JWST Cycle 1 Transiting Exoplanets", fontsize=fnt)

#annotate

plt.plot(1, 1, 'o', color='cyan', alpha=1, ms=12)
plt.annotate('Earth', (1.1,0.83), color='cyan', fontsize=fnt)

plt.plot(.037, 11.2, 'o', color='khaki', alpha=1, ms=12)
plt.annotate('Jupiter', (.05,9.0), color='khaki', fontsize=fnt)

plt.savefig(ddir+"cycle1_planets_influx.png")
