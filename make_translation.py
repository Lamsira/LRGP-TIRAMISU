import numpy as np
import matplotlib.pyplot as plt
import itertools
import copy

class IsomerGroup:
    def __init__(self, output, n, max_permut, atomic_composition=dict(), model1_species=list(),
                 model2_species=list()):
        self.output = output
        self.atomic_composition = atomic_composition
        self.model1_species = model1_species
        self.model2_species = model2_species
        self.n = n
        self.max_permut = max_permut
        self.translations = dict()
        self.r2_permutations = list()
        self.model_save=[therm.species for therm in self.model2_species]
        self.flat_errors=[]

    def nb_permutation(self,x,y) :
        return(np.math.factorial(x)/(np.math.factorial(x-y)))

    def trad(self):
        # print("REGARDE   ",self.model2_species)
        model1_species = [therm.species for therm in self.model1_species]
        model2_species = [therm.species for therm in self.model2_species]
        
        if self.atomic_composition == {'c':4, 'h':9, 'o' : 1}:
            print(model1_species, model2_species)
        
        trad_errors = np.zeros((len(model1_species), len(model2_species)))
        n_rows, n_columns = trad_errors.shape
        flat_errors = np.zeros((len(model1_species), len(model2_species)))
        n_rows, n_columns = flat_errors.shape
        
        """-----------------------------------------------------------------"""
        
        model1_tmp_therm=[]
        for model1_therm in self.model1_species :
            model1_tmp_therm.append(calc_identical_therm(model1_therm))
        model2_tmp_therm=[]
        for model2_therm in self.model2_species :
            model2_tmp_therm.append(calc_identical_therm(model2_therm))
            
        for n_mod, mod in enumerate([model1_tmp_therm, model2_tmp_therm]) :
            
            count = 0
            indices_dict = {}
            
            for subl in range(len(mod)):
                sublist = tuple(mod[subl])
                if sublist not in indices_dict:
                    indices_dict[sublist] = [subl]
                else:
                    indices_dict[sublist].append(subl)
            
            
            for sublist, indices in indices_dict.items():
                if len(indices) > 1:
                    count += 1
                    # print(f"model{n_mod+1} : Identical sublist {sublist} found at indices: {indices}")
                    if n_mod == 0 :
                        with open(self.output+'similar_therm_model1.txt', 'a') as rt :
                            # print(f'{[model1_species[x] for x in indices]}')
                            rt.write(f'{[model1_species[x] for x in indices]}\n')
                    else :
                        with open(self.output+'similar_therm_model2.txt', 'a') as rt :
                            # print(f'{[model2_species[x] for x in indices]}')
                            rt.write(f'{[model2_species[x] for x in indices]}\n')
                    
            
            # if count != 0:
            #     print("There are", count, "groups of identical sublists.")
        
        """-----------------------------------------------------------------"""
        
        for i, model1_therm in enumerate(self.model1_species):
            for j, model2_therm in enumerate(self.model2_species):
                trad_errors[i][j], flat_errors[i][j] = calc_error(model1_therm, model2_therm, self.n)
                
                
        # if self.atomic_composition == {'c':2, 'o': 2, 'h': 2} :
        #     print(f'TRAD ERRORS :\n{trad_errors}')
            
        # if trad_errors != [] :
        #     print(trad_errors)
        """------------------------------------------------------------""" 
        is_sup=0
        is_shape=(len(model1_species),len(model2_species))
        if len(model1_species)*len(model2_species) > 1 : 
            for spe in self.model1_species :
                self.translations[str(spe)]=[[None, None, None, [42,0,0]], False, 0]
            is_sup=2
            maxp=max(is_shape)
            minp=min(is_shape)
            n_permut=self.nb_permutation(maxp,minp)
            if n_permut < self.max_permut :
                print(f'n_permut : {n_permut}')
                is_sup=1
                unique_combinations = []

                # Getting all permutations of list_1
                # with length of list_2
                if len(self.model1_species) != len(self.model2_species) :
                    max_spe=max(self.model1_species, self.model2_species, key=len)
                    min_spe=min(self.model1_species, self.model2_species, key=len)
                else :
                    max_spe=self.model1_species
                    min_spe=self.model2_species
                permut = itertools.permutations(max_spe, len(min_spe))

                # zip() is called to pair each permutation
                # and shorter list element into combination
                for comb in permut:
                    zipped = zip(comb, min_spe)
                    unique_combinations.append(list(zipped))
                if min_spe == self.model1_species :
                    for i,to_reverse in enumerate(unique_combinations.copy()) :
                        for j,trad in enumerate(to_reverse) :
                            unique_combinations[i][j]=tuple(reversed(unique_combinations[i][j]))

                unique_errors = []
                
                for combination in unique_combinations :
                    tmp_errors=[]
                    for tmp_translation in combination :
                        tmp_errors.append(calc_error(tmp_translation[0], tmp_translation[1], self.n)[0])
                    unique_errors.append(tmp_errors)

                r2_errors = []
                for errors in unique_errors :
                    r2_errors.append(np.sqrt(np.sum(np.array(errors)**2)))
                r2_min=min([abs(r2) for r2 in r2_errors])
                permutation_translation=unique_combinations[r2_errors.index(r2_min)]
                # if self.atomic_composition == {'C':2, 'H':3, 'O':2} :
                #     print(f'PERMUTATION :\n{permutation_translation}')
                # else : 
                #     print(self.atomic_composition)
                tmp_argsort=np.argsort(r2_errors)
                permutation_translation2=unique_combinations[tmp_argsort[1]]
                
                for spe in self.model1_species :
                    self.translations[str(spe)]=[[None, None, None, [0]+sorted(r2_errors)[0:2]], False, 0]
                for trad in permutation_translation :
                    # print(f'Trad : {trad}')
                    # self.translations[str(trad[0])]=[[None, str(trad[1])],[False]]
                    self.translations[str(trad[0])][0][1]=str(trad[1])
                    # print(f'LA TRAD : {str(trad[1])}')
                    
                for trad in permutation_translation2 :
                    try :
                        # print(str(trad[1]), sorted(r2_errors)[0:2])
                        self.translations[str(trad[0])][0][2]=str(trad[1])
                        # self.translations[str(trad[0])][0][3]=sorted(r2_errors)[0:2]
                    except KeyError :
                        print("!!!!!!!!!!!!!!!! SHOULDN'T HAPPEN !!!!!!!!!!!!!!")
                        self.translations[str(trad[0])]=[[None, None, str(trad[1]), sorted(r2_errors)[0:2]], False, 'NaN']
                        # self.translations[str(trad[0])][0][3]=sorted(r2_errors)[0:2]
                self.r2_permutations=sorted(r2_errors)[:10]
        """------------------------------------------------------------"""


        while 0 not in trad_errors.shape:
            ambiguous_key=False
            i_min, j_min = np.unravel_index(trad_errors.argmin(), trad_errors.shape)
            if trad_errors[i_min][j_min] > 1:
                ambiguous_key=True
                print("{} --> {} : {} %".format(model1_species[i_min],
                                              model2_species[j_min],
                                              trad_errors[i_min][j_min]))
                
            if is_sup==0 :
                self.translations[model1_species.pop(i_min)] = [model2_species.pop(j_min), ambiguous_key, [trad_errors[i_min][j_min]], flat_errors[i_min][j_min]]
                
            elif is_sup==2 :
                ambiguous_key=True
                self.translations[model1_species.pop(i_min)] = [model2_species.pop(j_min), ambiguous_key, [trad_errors[i_min][j_min]], flat_errors[i_min][j_min]]
                
            elif model1_species[i_min] in self.translations and self.translations[model1_species[i_min]][0][3][0] != 42 :
                
                self.translations[model1_species[i_min]][0][3][0] = trad_errors[i_min][j_min]
                self.translations[model1_species[i_min]][1] = ambiguous_key
                self.translations[model1_species[i_min]][2] = flat_errors[i_min][j_min] 
                self.translations[model1_species.pop(i_min)][0][0] = model2_species.pop(j_min)
                
                # if self.atomic_composition == {'c':2, 'o': 2, 'h': 2} :
                #     i_min, j_min = np.unravel_index(trad_errors.argmin(), trad_errors.shape)
                #     print(trad_errors[i_min][j_min], trad_errors[''])
                
            else :
                # self.translations[model1_species.pop(i_min)] = [[model2_species.pop(j_min), None],[ambiguous_key]]
                print(f'\n Species never trad before method 1 (which is the second applied ... i know.). Watchout !\n{model1_species[i_min]} ')
                self.translations[model1_species.pop(i_min)] = [[model2_species.pop(j_min), None, None, [1000,1000,1000]],[ambiguous_key]]
            trad_errors = np.delete(trad_errors, i_min, 0)
            trad_errors = np.delete(trad_errors, j_min, 1)
            flat_errors = np.delete(flat_errors, i_min, 0)
            flat_errors = np.delete(flat_errors, j_min, 1)

        # if is_sup == 1 :  
        #     print(self.translations,'\n')
        # for species in model1_species :
        #     'xar'
        #     # print(f"{{'{species}':''}}")
        # for species in model2_species :
        #     'xar'
        #     # print(f"{{'':'{species}'}}")

    def plot(self):
        for therm in self.model1_species:
            plt.plot(*therm.G(), label = therm.species + " 1")
        for therm in self.model2_species:
            plt.plot(*therm.G(), label = therm.species + " 2")
        plt.legend()
        
class translation:
    def __init__(self, erreur=float(), species=dict()):
        self.erreur = erreur
        self.species =  species
        

"""-------------------------------------------------------------------------
                               NECESSARY FUNCTIONS
-------------------------------------------------------------------------"""
    
def calc_error(model1, model2, n):
    low_temp = max(model1.temp_list[0], model2.temp_list[0])
    high_temp = min(model1.temp_list[2], model1.temp_list[2])
    temp = np.linspace(low_temp, high_temp, n)
    G_model1 = model1.G_new(temp) ###
    G_model2 = model2.G_new(temp) ###
    norm = min(np.max(G_model1**2) - np.min(G_model1**2), np.max(G_model2**2) - np.min(G_model2**2))
    # norm = 1
    maxn = np.max(np.sqrt((G_model1-G_model2)**2))
    error = np.mean((G_model1-G_model2)**2)*100/norm
    return error, maxn

def calc_identical_therm(model) :
    low_temp = 1000
    high_temp = 2000
    temp = np.linspace(low_temp, high_temp, 2)
    G_model = model.G_new(temp)
    return(G_model)

def get_isomer_group(isomer_groups, therm):
    for isomer_group in isomer_groups:
        if isomer_group.atomic_composition == therm.atomic_composition:
            return isomer_group

def get_atomic_number(atom, atomic_composition):
    try:
        return atomic_composition[atom]
    except KeyError:
        return 0

def comp_sort_key(atomic_composition):
    atoms_priority_list = ["c", "o", "n", "h", "he", "ar"]
    return [get_atomic_number(atom, atomic_composition) for atom in atoms_priority_list]

def dict_sort_key(atom):
    atoms_priority_list = ["c", "o", "n", "h", "he", "ar"]
    return atoms_priority_list.index(atom[0])

"""-------------------------------------------------------------------------
                                   JOB FUNCTIONS
-------------------------------------------------------------------------"""
def coherence_check(model1, model2, spe_to_trad, translat) :
    translat=[str(x) for x in translat]
    translat=set(translat)
    translat-={None}
    toprint='!'
    if translat == set() :
        return('')
    else :
        watch_reaction=[]
        to_set_w=[]
        for reaction in model1.reactions[2] :
            if spe_to_trad in [*reaction.species[0], *reaction.species[1]] :
                tmp_rea=copy.deepcopy(reaction.species)
                tmp_set=[set(x) for x in tmp_rea]
                if tmp_set not in to_set_w and tmp_set[::-1] not in to_set_w :
                    for ind_i,side in enumerate(tmp_rea) :
                        for ind_j,spe in enumerate(side) :
                            if spe != spe_to_trad :
                                tmp_rea[ind_i][ind_j]=model1.therm[spe].atomic_composition
                    watch_reaction.append(tmp_rea)
        # for i in watch_reaction :
            # print(i)
        
        for trad in translat :
            
            test_reaction=[]
            to_set=[]
            for reaction in model2.reactions[2] :
                if trad in [*reaction.species[0], *reaction.species[1]] :
                    tmp_rea=copy.deepcopy(reaction.species)
                    tmp_set=[set(x) for x in tmp_rea]
                    if tmp_set not in to_set and tmp_set[::-1] not in to_set :
                        for ind_i,side in enumerate(tmp_rea) :
                            for ind_j,spe in enumerate(side) :
                                if spe != trad :
                                    tmp_rea[ind_i][ind_j]=model2.therm[spe].atomic_composition
                        to_set.append(tmp_set)
                        test_reaction.append(tmp_rea)
            
            tmp_wr=copy.deepcopy(watch_reaction)
            for ind_i,reaction in enumerate(tmp_wr.copy()) :
                for ind_j,side in enumerate(reaction) :
                    for ind_k, spe in enumerate(side) :
                        if spe == spe_to_trad :
                            tmp_wr[ind_i][ind_j][ind_k] = trad
            
            # for i in tmp_wr :
                # print(i)
            count=0
            for rea1 in tmp_wr :
                for rea2 in test_reaction :
                    if rea1 == rea2 or rea1[::-1] == rea2 :
                        count+=1
                        break
            
            toprint+=f'     {trad} : {count} {len(watch_reaction)}/{len(test_reaction)}'
    toprint+='\n'
    return(toprint)
                        
                
            
    
def models_translation(model1, model2, output, n, max_permut):
    
    with open(output+'similar_therm_model1.txt', 'w') as c1, open(output+'similar_therm_model2.txt', 'w') as c2 :
        ''

    
    atomic_compositions = []
    for therm in model1.therm :
        if model1.therm[therm].atomic_composition not in atomic_compositions:
            atomic_compositions.append(model1.therm[therm].atomic_composition)
    for therm in model2.therm :
        if model2.therm[therm].atomic_composition not in atomic_compositions:
            atomic_compositions.append(model2.therm[therm].atomic_composition)
    
    
    atomic_compositions = [dict(sorted(atomic_composition.items(), key=dict_sort_key))
                       for atomic_composition in atomic_compositions]
    atomic_compositions = sorted(atomic_compositions, key=comp_sort_key)
    
    
    isomer_groups = [IsomerGroup(output, n, max_permut,atomic_composition=atomic_composition,
                                 model1_species=list(), model2_species=list()) for atomic_composition in atomic_compositions]
    
    
    for therm1 in model1.therm:
        get_isomer_group(isomer_groups, model1.therm[therm1]).model1_species.append(model1.therm[therm1])
    for therm2 in model2.therm:
        get_isomer_group(isomer_groups, model2.therm[therm2]).model2_species.append(model2.therm[therm2])
        
        
    with open(output+'all_comparisons.log','w') as comparisons :
        for i in isomer_groups :
            name1=model1.path_mech[model1.path_mech.rfind("/")+1:model1.path_mech.rfind(".")]
            name2=model2.path_mech[model2.path_mech.rfind("/")+1:model2.path_mech.rfind(".")]
            # print(i.atomic_composition,f'\n{name1:<20} : {" ".join(str(x) for x in i.model1_species)}'\
            #       f'\n{name2 :<20} : {" ".join(str(x) for x in i.model2_species)}\n')
            comparisons.write((f'{i.atomic_composition}\n{name1:<20} : {" ".join(str(x) for x in i.model1_species)}'\
                  f'\n{name2 :<20} : {" ".join(str(x) for x in i.model2_species)}\n\n'))
    
    
    count=0
#     def parallelize(isomer_group):
#         # code to be executed in parallel goes here
#         isomer_group.trad()
#         if isomer_group.translations != {} :
#             count+=1
#             for iso in isomer_group.translations :
#                 if type(isomer_group.translations[iso][0]) is list :
#                     if isomer_group.translations[iso][0][0] != isomer_group.translations[iso][0][1] :
#                         print('DIFFERENT TRAD SPECIES BETWEEN THE TWO TECHNICS')
#                         print(isomer_group.r2_permutations)
#                         break
    
    
#     with ThreadPoolExecutor(max_workers=8) as executor:
#         for isomer_group in isomer_groups[:]:
#             # print(f'\n{isomer_group.atomic_composition}')
#             executor.submit(parallelize, isomer_group)
    
    for isomer_group in isomer_groups :
        # print(f'\n{isomer_group.atomic_composition}')
        isomer_group.trad()
        if isomer_group.translations != {} :
            count+=1
            # print(f'\n{isomer_group.atomic_composition}')
            # 'xar'
            # print(isomer_group.translations)
            for iso in isomer_group.translations :
                
                if type(isomer_group.translations[iso][0]) is list :
                    if isomer_group.translations[iso][0][0] != isomer_group.translations[iso][0][1] :
                        print('DIFFERENT TRAD SPECIES BETWEEN THE TWO TECHNICS')
                        print(isomer_group.r2_permutations)
                        break
        
    print(f'Nb of isomer_groups to trad : {count}')
    with open(output+'tmp_species_traduction.tmp','w') as result :
        result.write('! TEMPORARY TRADUCTION FILE\n! Choose an isomer for every line commented by a "!CHECK" and rename the file into'+\
                    f' tmp_speces_traduction.inp\n! Made with {n} permutations allowed.\n! [X, Y, Z] : \n! X = error between species to trad and the traduction'+\
                    '\n! Y = error of the best fit for all isomers\n! Z = error of the second best fit for all isomers\n')
        for isomer_group in isomer_groups :
            if isomer_group.atomic_composition == {'c':4, 'h':9, 'o':1} :
                print(isomer_group.translations)
            if isomer_group.translations != {} : 
                
                cool = [f'{x}{isomer_group.atomic_composition[x]} ' for x in isomer_group.atomic_composition]
                result.write(f'!\n! {"".join(cool)}\n!{isomer_group.model2_species}\n!\n')
                for iso in isomer_group.translations :
                    # if iso in [x.species for x in isomer_group.model2_species] :
                    #     print(iso, [x.species for x in isomer_group.model2_species])
                    translation=isomer_group.translations[iso]
                    if type(translation[0]) is list :
                        # print(translation)
                        check_error=translation[0][3]
                        if check_error[1] == 0 :
                            result.write(f'{iso} = {translation[0][1]}\n')
                        elif check_error[0] == 0 and translation[0][0] != None :
                            first_part=f'{iso} = {translation[0][0]}'
                            result.write(f'{first_part:<50}{" ":10}!WATCHOUT {translation[0][3]}\n')
                        else :
                            diff=check_error[1]/check_error[2]
                            if diff >= 0.8 :
                                
                                sec_part=coherence_check(model1, model2, iso, isomer_group.model2_species)
                                first_part=f'!CHECK {iso} = {translation[0][0:3]}'
                                result.write(f'{first_part:<50}{" ":10}!{translation[0][3]} {translation[2]:.1f}\n{sec_part}')
                            elif diff >= 0.1 :
                                if translation[0][0] != translation[0][1] :
                                    sec_part=coherence_check(model1, model2, iso, isomer_group.model2_species)
                                    first_part=f'!CHECK {iso} = {translation[0][0:3]}'
                                    result.write(f'{first_part:<50}{" ":10}!{translation[0][3]} {translation[2]:.1f}\n{sec_part}')
                                else :
                                    if translation[0][1] != translation[0][2] :
                                        first_part=f'{iso} = {translation[0][1]}'
                                        result.write(f'{first_part:<50}{" ":10}!Should be ok but look : {translation[0][0:3]} {translation[0][3]} {translation[2]:.1f}\n')
                                    else :
                                        first_part=f'{iso} = {translation[0][1]}'
                                        result.write(f'{first_part:<50}{" ":10}!Should be ok {translation[2]:.1f}\n')
                            else :
                                first_part=f'{iso} = {translation[0][1]}'
                                result.write(f'{first_part:<50}{" ":10}!diff<0.1\n')
                    else :
                        if translation[1] is False :
                            result.write(f'{iso} = {translation[0]}\n')
                        elif translation[1] is True :
                            if translation[2][0] == 0 :
                                result.write(f'{iso} = {translation[0]}\n')
                            else :
                                sec_part=coherence_check(model1, model2, iso, isomer_group.model2_species)
                                first_part=f'!CHECK {iso} = {translation[0]}'
                                result.write(f'{first_part:<50}{" ":10}!Maybe too much isomers. Be careful. Error = {translation[2]} {translation[3]:.1f}\n{sec_part}')
                        else :
                            print(f'ERROR : {translation} not True nor False')
                
        # if isomer_group.atomic_composition == {'C': 3, 'H': 2}:
        #     isomer_group.plot()
    return(isomer_groups)

def parse_file(filename):
    result = {}
    warning = []
    with open(filename, 'r') as f:
        # do_break=False
        # skip_all=False
        line_count=0
        for line in f:
            # Ignore lines that start with '!' or are empty
            # if line.startswith('!CHECK') and skip_all == False:
            if line.startswith('!CHECK') :
                line_count+=1
                # if skip_all == True :
                    
                #     line_is=line[6:]
                #     key, trash = line_is.split('=', 1)
                #     key = key.strip()
                #     # print(f'{key} should have been added to warning')
                #     warning.append(key)
                #     continue
                # else :
                #     line_is=line[6:]
                #     do_break=False
                    
                #     while True : 
                #         check=input(f"You didn't choose an option for {line_is}Are you sure you want to pass it (y/y_all/n)?")
                #         if check == 'n' or check == "N" :
                #             do_break=True
                #             break
                #         elif check == 'y' or check == "Y" :
                #             line_is=line[6:]
                #             key, trash = line_is.split('=', 1)
                #             key = key.strip()
                #             warning.append(key)
                #             break
                #         elif check == 'y_all' :
                #             skip_all=True
                #             line_is=line[6:]
                #             key, trash = line_is.split('=', 1)
                #             key = key.strip()
                #             warning.append(key)
                #             break
                #         else :
                #             print("You potato ! You didn't put an acceptable answer !")
                            
                # if do_break==True :
                #     print("FINAL WORD")
                #     result={}
                #     break
                    
            if line.startswith('!') or line.strip() == '':
                continue

            # Split the line on the first '=' character
            key, value = line.split('=', 1)

            # Strip any leading or trailing white space from the key and value
            key = key.strip()
            value = value.strip()

            # Check if the value has a comment following it
            if '!' in value:
                # Split the value on the '!' character to separate the value from the comment
                value, comment = value.split('!', 1)
                # Strip any leading or trailing white space from the value and comment
                value = value.strip()
                comment = comment.strip()

            # Add the key-value pair to the dictionary
            if value == "None" :
                # print(f'{key} is None')
                continue
            result[key] = value
    # print(warning)
    print(f'Number of !CHECK : {line_count}')
    return(result, warning)


def translate(model1, translation):
    
    translated = copy.deepcopy(model1)
    
    for species in translation :
        if species not in model1.reactions[1] :
            raise Exception(f'You probably erased a letter. {species} (trad : {translation[species]}) not in the model to translate.')
    
    for ind,species in enumerate(translated.reactions[1].copy()) :
        if species in translation :
            translated.reactions[1][ind]=translation[species]

    for ind,species in enumerate(translated.therm.copy()) :
        if species in translation :
            translated.therm[translation[species]]=translated.therm.pop(species)
            translated.therm[translation[species]].species=translation[species]

    for ind,species in enumerate(translated.trans.copy()) :
        if species in translation :
            translated.trans[translation[species]]=translated.trans.pop(species)
            translated.trans[translation[species]].species=translation[species]
    for ind,reaction in enumerate(translated.reactions[2].copy()) :
        for ind2,side in enumerate(reaction.species) :
            for ind3,spe in enumerate(side) :
                if spe in translation :
                    translated.reactions[2][ind].species[ind2][ind3] = translation[spe]
                    
    return(translated)