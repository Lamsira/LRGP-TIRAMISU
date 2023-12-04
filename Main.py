import time
import datetime
import pickle 
    
import read_model_PAH as rm
import make_translation as mt
import add_and_write_mechanism as awm

#%%

"""-----------------------    INPUT 1    -----------------------"""

# pick=True
pick=False

output='E:/Chemkin/Auto_simu/ALLNLV2_4B/1/'

path_mech2='E:/Chemkin/Auto_simu/ALLNLV2-1/ALLNLV2-1.mech'
path_therm2='E:/Chemkin/Auto_simu/ALLNLV2-1/ALLNLV2-1.therm'
path_trans2='E:/Chemkin/Auto_simu/ALLNLV2-1/ALLNLV2-1.trans'
# path_mech2='E:/Chemkin/Auto_simu/ALLNLV2-4/ALLNLV2-4.mech'
# path_therm2='E:/Chemkin/Auto_simu/ALLNLV2-4/ALLNLV2-4.therm'
# path_trans2='E:/Chemkin/Auto_simu/ALLNLV2-4/ALLNLV2-4.trans'


path_mech1='E:/Chemkin/Auto_simu/4-B/H_But_Mech_(FS).txt'
path_therm1='E:/Chemkin/Auto_simu/4-B/But_Therm.txt'
path_trans1='E:/Chemkin/Auto_simu/4-B/But_Trans.txt'

n=100000

print(f'NB PERMUTATIONS : {n}')

"""-----------------------      END      -----------------------"""
if pick == True :
    time_start = time.time()
    
    model1=rm.Model(path_mech1, path_therm1, output, path_trans1)
    model2=rm.Model(path_mech2, path_therm2, output, path_trans2)
    
    time_end=time.time()
    
    print(f'\nTotal time = {time_end-time_start} seconds or {datetime.timedelta(seconds=(time_end-time_start))} ')
    
    print('\nREGARDE ICI, LA : ',len(model1.reactions[1]), len(model1.therm))
    print('REGARDE ICI, LA : ',len(model2.reactions[1]), len(model2.therm))
    with open(f'{output}model1.pick', 'wb') as m1 :
        pickle.dump(model1, m1)
    with open(f'{output}model2.pick', 'wb') as m2 :
        pickle.dump(model2, m2)


else :
    with open(f'{output}model1.pick', 'rb') as m1 :
        model1 = pickle.load(m1)
    with open(f'{output}model2.pick', 'rb') as m2 :
        model2 = pickle.load(m2)

    

time_start = time.time()

# stock=mt.models_translation(model1, model2, output, 350, n)

time_end=time.time()



print(f'\nTotal time = {time_end-time_start} seconds or {datetime.timedelta(seconds=(time_end-time_start))} ')



#%%


print(len(model2.reactions[1]), len(model2.reactions[2]))
"""-----------------------    INPUT 2    -----------------------"""

time_start = time.time()

""" TO CHANGE START """

out= f'{output}ALLNLV2-1_4B'
# out= f'{output}ALLNLV2-4_4B'

md1='4-BUTANOLS'

md2='ALLNLV2'

depth=1e10

introduction=f'! New mechanism (V2) made from LLNL and Alkylaromatic with 4-Butanols IButanol sub-mechanism.\n! n={n} depth={depth}\n! {datetime.date.today().strftime("%d/%m/%Y")}\n'

list_to_add=['ic4h9oh']

""" END """

inp=f'{output}tmp_species_traduction.inp'
# inp='E:/Chemkin/Auto_simu/ALLNLV2_4B/1/tmp_species_traduction.inp'

translation, warning=mt.parse_file(inp)

############################## VERIF
count=0
for element in translation :
    if translation[element] not in model2.reactions[1] :
        print(f'WTF {element} {translation[element]} not in species !')
    else :
        count+=1
    # if translation[element] != element :
    #     print(element, translation[element])
print('\n\n1/2 = nb of translation (different if bad translation !), 3 = nb of species of model1, 4 = nb of species of model2')
print(count, len(translation), len(model1.reactions[1]), len(model2.reactions[1]))
##############################

model1_translated = mt.translate(model1, translation)

to_add = awm.recursive_addition(list_to_add, model1_translated, model2, warning, depth)

if to_add == '' :
    print('STOP') 

############################## VERIF
print('\n\nLength of to_add : (reactions/species)')
print(len(to_add[0]), len(to_add[1]))
# print('\nList of species that are not added to the new mechanism')
# for species in model1_translated.reactions[1] :
#     if species not in model2.reactions[1] and species not in to_add[1] :
#         print(species)
allspe=set(model2.reactions[1])
print('\nShow all reactions that are in model1 but not in model2 if all of their species are already in model2')
for reaction in model1_translated.reactions[2] :
    if reaction.species[0] == reaction.species[1] : 
        print(f'UN-LUMPING : {reaction.species[0]} = {reaction.species[1]}')
    # tmp_spe=set([*reaction.species[0], *reaction.species[1]])
    # if tmp_spe.issubset(allspe) :
    #     if reaction not in model2.reactions[2] :
    #         print(reaction.species)
############################## VERIF


mechanism=awm.mechanism_writting(out, md1, md2, model1_translated, model2, to_add, introduction)

time_end=time.time()

print(f'\nTotal time = {time_end-time_start} seconds or {datetime.timedelta(seconds=(time_end-time_start))} ')