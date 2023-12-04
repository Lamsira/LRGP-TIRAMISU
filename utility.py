# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 17:01:31 2023

@author: fages1
"""
import time
import datetime
import pickle 
import os 
import read_model_PAH as rm
import make_translation as mt
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell
import re

# output_gen='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/dump/'
output_gen='E:/Chemkin/Exgas_ALLNL_Sensi/'

# path_mech='C:/Users/fages1/OneDrive - Universite de Lorraine/Mecanisme/LLNL_PAH23/LLNL_PAH23.mech'
# path_mech='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/dump/cyc5h71--_main.mech'
# path_mech='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/LLNL_PAH23/LLNL_PAH23_PAH_mech.txt'
# path_therm='C:/Users/fages1/OneDrive - Universite de Lorraine/Mecanisme/LLNL_PAH23/LLNL_PAH23.therm'
# path_trans='C:/Users/fages1/OneDrive - Universite de Lorraine/Mecanisme/LLNL_PAH23/LLNL_PAH23.trans'

path_mech='E:/Chemkin/Exgas_ALLNL_Sensi/Exgas/CHEMKIN_exgas_modified.inp'
path_therm='E:/Chemkin/Exgas_ALLNL_Sensi/Exgas/CHEMKIN_exgas_modified.inp'
mech1=rm.Model(path_mech, path_therm, output_gen, '')
# path_mech='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/dump/temp.mech'
# path_mech='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/+Tsilla/ALLNLV2_B_N+Tsilla.mech'
# path_therm='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/+Tsilla/ALLNLV2_B_N+Tsilla.therm'
# path_trans='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/+Tsilla/ALLNLV2_B_N+Tsilla.trans'
# path_mech='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/ALLNLV2_B_N.mech'
# path_therm='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/ALLNLV2_B_N.therm'
# path_trans='C:/Users/fages1/OneDrive - Universite de Lorraine/Code_these/PAH/ALLNLV2_B_N/Iterations/ALLNLV2_B_N.trans'
# path_mech='C:/Users/fages1/OneDrive - Universite de Lorraine/Mecanisme/Tsilla/Tsilla.mech'
# path_therm='C:/Users/fages1/OneDrive - Universite de Lorraine/Mecanisme/Tsilla/Tsilla.therm'
# path_trans='C:/Users/fages1/OneDrive - Universite de Lorraine/Mecanisme/Tsilla/Tsilla.trans'
path_mech='E:/Chemkin/Exgas_ALLNL_Sensi/ALLNL/AlkLLNL.inp'
path_therm='E:/Chemkin/Exgas_ALLNL_Sensi/ALLNL/AlkLLNL_THERMO.dat'
path_trans='E:/Chemkin/Exgas_ALLNL_Sensi/ALLNL/AlkLLNL_TRANSPORT.dat'


mech2=rm.Model(path_mech, path_therm, output_gen, path_trans)

#%%

# spe_list=['CYC5H71-3','CYC5H71-4']
spe_list=['C3H6O1-2']
spe_list=[x.lower() for x in spe_list]

aze=mech2

write=aze.micro_mech(spe_list)

with open(output_gen+'temp'+'.mech', 'w') as out :
    """ ELEMENTS """

    out.write('elements\n'.upper())
    elements=set(aze.reactions[0])
    for element in elements :
        out.write(f'{element}\n'.upper())

    """ SPECIES """

    out.write('end\nspecies\n'.upper())
    for species in write[0] :
        out.write(f'{species:<40}!\n'.upper())
        
    out.write('end\nreactions\n'.upper())

count=0
for reaction in write[1] :
    count+=1
    reaction.print_reaction(output_gen+'temp.mech')
    
with open(output_gen+'temp'+'.mech', 'a') as out :
    out.write('end'.upper())
    
print(count)
#%%
import matplotlib.pyplot as plt
from cycler import cycler
line=['-','-' ,'-' ,'--', '--', '--']
color=['b','r']#, 'b', 'y', 'm']
cycle=cycler(linestyle=line)
therm1=['R3OOH']
therm2=['HO2']
name_list=[['EXGAS'], ['ALLNL']]
count=0
fig1, axcp = plt.subplots()
fig2,axh = plt.subplots()
fig3,axs = plt.subplots()
fig4,axg = plt.subplots()

fig5,ax5 = plt.subplots()
# ax5.set_prop_cycle(cycle)
ax5.set_xlabel('T (K)', fontsize=14)
ax5.set_ylabel('H (kcal/mol), S (cal/mol/K)')
ax6=ax5.twinx()
ax6.set_ylabel('Cp (cal/mol/K)')
axcp.set_prop_cycle(cycle)
axcp.set_xlabel('T (K)', fontsize=14)
# axcp.set_ylabel('cp/R', fontsize=14)
axcp.set_ylabel('Cp', fontsize=14)
axh.set_prop_cycle(cycle)
axh.set_xlabel('T (K)', fontsize=14)
# axh.set_ylabel('h/RT', fontsize=14)
axh.set_ylabel('H', fontsize=14)
axs.set_prop_cycle(cycle)
axs.set_xlabel('T (K)', fontsize=14)
# axs.set_ylabel('s/R', fontsize=14)
axs.set_ylabel('S', fontsize=14)
axg.set_prop_cycle(cycle)
axg.set_xlabel('T (K)', fontsize=14)
# axg.set_ylabel('h/RT-s/R', fontsize=14)
axg.set_ylabel('G', fontsize=14)

h1=[]
h2=[]

for ind,therm in enumerate(therm1):
    name=name_list[0][ind]
    data=mech1.therm[therm.lower()].plot_data()
    axcp.plot(data['T'], data['cp'], label = name+' Cp', color=color[count])
    axh.plot(data['T'], data['h'], label = name+' H', color=color[count])
    axs.plot(data['T'], data['s'], label = name+' S', color=color[count])
    # axg.plot(*mech1.therm[therm.lower()].G(), label = therm, color=color[count])
    axg.plot(data['T'], data['g'], label = name+' G', color=color[count])
    
    
    ax5.plot(data['T'], data['hk'], label = name+' H', color=color[count], linestyle='--')
    ax5.plot(data['T'], data['s'], label = name+' S', color=color[count], linestyle=':')
    
    ax6.plot(data['T'], data['cp'], label = name+' Cp', color=color[count], linestyle='-')
    # ax5.plot(data['T'], data['mg10k'], label = name+' -G/10000', color=color[count], linestyle= 'dashdot')
    h1.append(data['h_value'])
    count+=1
    if count == 5 :
        count=0
for ind,therm in enumerate(therm2):
    name=name_list[1][ind]
    data=mech2.therm[therm.lower()].plot_data()
    axcp.plot(data['T'], data['cp'], label = name+' Cp', color=color[count])
    axh.plot(data['T'], data['h'], label = name+' H', color=color[count])
    axs.plot(data['T'], data['s'], label = name+' S', color=color[count])
    # axg.plot(*mech2.therm[therm.lower()].G(), label = therm, color=color[count])
    axg.plot(data['T'], data['g'], label = name+' G', color=color[count])
    
    ax6.plot(data['T'], data['cp'], label = name+' Cp', color=color[count], linestyle='-')
    ax5.plot(data['T'], data['hk'], label = name+' H', color=color[count], linestyle='--')
    ax5.plot(data['T'], data['s'], label = name+' S', color=color[count], linestyle=':')
    # ax5.plot(data['T'], data['mg10k'], label = name+' -G/10000', color=color[count], linestyle= 'dashdot')
    h2.append(data['h_value'])
    count+=1
    if count == 5 :
        count=0
        
handles1, labels1 = ax5.get_legend_handles_labels()
handles2, labels2 = ax6.get_legend_handles_labels()

# handles = [handles2[0], handles1[0], handles1[1], handles1[2], handles1[1], handles1[3], handles1[4], handles1[5]]
# labels = [labels2[0], labels1[0], labels1[1], labels1[2], labels2[1], labels1[3], labels1[4], labels1[5]]
handles = [handles2[0], handles1[0], handles1[1], handles2[1], handles1[2], handles1[3]]
labels = [labels2[0], labels1[0], labels1[1], labels2[1], labels1[2], labels1[3]]

axcp.legend()
axh.legend()
axs.legend()
axg.legend()
# ax5.legend()
# axh.set_xlim(300,305)
# axh.set_ylim(-5,50)
ax5.set_title('Thermodynamic Data of HO2')
# ax5.set_xlim(280, 1200)
fig5.legend(handles, labels ,loc='lower right', bbox_to_anchor= (0.9,0.12))
val_tem=[300,500,700]
for idj,j in enumerate(h1) :
    for idk,k in enumerate(h2) :
        prt=f'{therm1[idj]}-{therm2[idk]} :\n'
        for ind,i in enumerate(val_tem) :
            prt+=f'{i}K : {j[ind]-k[ind]:<.3f}  {(j[ind]-k[ind])/j[ind]*100: <.1f}%  ({j[ind]:<.1f})\n'
        print(prt)
fig5.savefig(('C:/Users/fages1/OneDrive - Universite de Lorraine/Presentations/thermo_HO2.png'))
#%%
stock=mt.models_translation(mech1, mech2, output_gen, 350, 150000)
#%%
print_list=[]
for species in mech1.reactions[1] :
    write=mech1.micro_mech([species])
    if len(write[1]) <= 1 :
        print_list.append([species, len(write[1])])
    
for spe in print_list :
    print(f'species = {spe[0]} nb_rea = {spe[1]}')
    
print(len(print_list))