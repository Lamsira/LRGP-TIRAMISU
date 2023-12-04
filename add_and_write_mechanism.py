import copy
import re

def recursive_addition(to_add_mech, model1_translated, model2, warning, depth=1e10) :
    warn_to_check=set()
    new_species=[]
    skip=0
    for species in to_add_mech :
        if species not in model1_translated.reactions[1] :
            print(f'ERROR : {species} not in model1')
            print(species, to_add_mech)
            print(model1_translated.reactions[1])
            return('')
        if species in model2.reactions[1] :
            
            while True and skip !=1 :
                conti=input(f'I think the mech {species} is already in model2. Do you want to continue (y/y_all/n)?')
                if conti == 'y' or conti == 'Y':
                    break
                if conti == 'y_all' :
                    skip = 1
                    break
                elif conti == 'n' or conti == 'N' :
                    return('ERROR : already in model2')
        new_species.append(species)
    new_reactions=[]
    
    n_species=copy.deepcopy(new_species)
    count_depth=0
    last_reactions=[]
    while (n_species != [] and count_depth < depth) :
        last2_reactions = copy.deepcopy(last_reactions)
        last_reactions=[]
        print(n_species)
        count_depth+=1
        tmp_n_species=copy.deepcopy(n_species)
        n_species = []
        if tmp_n_species != [] :
            for n_spe in tmp_n_species :
                count=0
                for reaction in model1_translated.reactions[2] :

                    if n_spe in reaction.species[0] or n_spe in reaction.species[1] :

                        if reaction in model2.reactions[2] :
                            print(reaction.species, "Error ?")

                        else :

                            if reaction not in new_reactions :
                                new_reactions.append(reaction)
                            else :
                                
                                len1=len([0 for x in new_reactions[new_reactions.index(reaction)].supinfo if x.strip().startswith('dup')])
                                len2=len([0 for x in reaction.supinfo if x.strip().startswith('dup')])   
                                if len1 == 0 and len2 != 0 :
                                    print(f'\nERROR LUMPING : {new_reactions[new_reactions.index(reaction)].__dict__} {reaction.__dict__}\n')
                                    new_reactions.remove(reaction)
                                    new_reactions.append(reaction)
                                    last_reactions.append(reaction)
                                elif len1 != 0 and len2 != 0 :
                                    xor=[rea.speed for rea in new_reactions if rea==reaction]
                                    # print(f'Faux xor : {xor}')
                                    if reaction.speed not in xor :
                                        new_reactions.append(reaction)
                                        last_reactions.append(reaction)
                                    # print(f'Xor : {xor}')
                            # print(f'{" + ".join(reaction.species[0])} = {" + ".join(reaction.species[1]):30}   {reaction.third_body.upper() if reaction.third_body!=0 else ""}')
                            for side in reaction.species :
                                for species in side :

                                    if species == "hv" :
                                        print(f'REGARDE ! UN PHOTON ! {reaction.species}')
                                    if species != "hv" and species not in new_species \
                                    and species not in model2.reactions[1] :
                                        if species in warning :
                                            print(f'WARNING : {species} ENCOUNTERED IN REACTION AND NOT TRADUCTED (!CHECK)')
                                            warn_to_check.add(species)

                                        n_species.append(species)
                                        new_species.append(species)

                        count+=1
        else :
            break
            # list_to_check=[]
            # for reaction in last_reactions :
            #     for species in [*reaction.species[0],reaction.species[1]]
            
            
            # print(count)
    for species in new_species.copy() :
        if species in model2.reactions[1] :
            print('Already there ....', species)
            new_species.pop(species)
    for reaction in new_reactions :
        len1=len([0 for x in reaction.supinfo if x.strip().startswith('dup')])
        if len1 != 0 :
            len2=len([0 for x in new_reactions if x==reaction])
            if len2 == 1 :
                print(f'\nA reaction have been duplicated :\n{reaction.__dict__}\n')
                new_reactions.append(reaction)
            elif len2 < 1 :
                print('FATAL ERROR')
            
        
    print(f'DEPTH = {count_depth}')
    print(f'\nNumber of species that can cause issue : {len(warn_to_check)}\n {warn_to_check}')
    return(new_reactions, new_species)

def mechanism_writting(outp, md1, md2, model1_translated, model2, to_add, custom='') :
    
    """ MECH """
    
    with open(outp+'.mech', 'w') as out :

        """ INTRODUCTION """

        out.write(custom)

        """ ELEMENTS """

        out.write('elements\n'.upper())
        elements=set(model1_translated.reactions[0]+model2.reactions[0])
        for element in elements :
            out.write(f'{element}\n'.upper())

        """ SPECIES 2 """

        out.write(f'end\n{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" SPECIES":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\nspecies\n'.upper())
        for species in model2.reactions[1] :
            out.write(f'{species}\n'.upper())

        """ SPECIES 1 """

        out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" SPECIES":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
        if to_add != '' :
            for species in to_add[1] :
                out.write(f'{species}\n'.upper())

        """ REACTIONS 2 """

        out.write(f'end\n{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" REACTIONS":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\nreactions\n'.upper())
        for reaction in model2.reactions[2] :
            if reaction.is_rev == True :
                equal= '='
            elif reaction.is_rev == False :
                equal= '=>'
            else :
                print(f'ERROR, NO REV {reaction.__dict__}')

            if reaction.third_body != 0 :
                tb=f"+{reaction.third_body}"
                if reaction.supinfo != [] :
                    low_check=["low/" for x in reaction.supinfo if 'low/' in x or 'low /' in x]
                    tcheb_check=["tcheb/" for x in reaction.supinfo if 'tcheb/' in x or 'tcheb /' in x]
                    if len(low_check) >= 1 or len(tcheb_check) >= 1 :
                        tb=f"(+{reaction.third_body})"
            else :
                tb=''

            reactant=[]
            product=[]
            # print(reaction.speed[1])
            for s, c in zip(reaction.species[0], reaction.coef[0]) :
                if c == 1 :
                    c=''
                elif int(c) == c :
                    c=int(c)
                reactant.append(str(c) + str (s))

            for s, c in zip(reaction.species[1], reaction.coef[1]) :
                if c == 1 :
                    c=''
                elif int(c) == c :
                    c=int(c)
                product.append(str(c) + str (s))

            final_reaction=f'{"+".join(reactant)}{tb if tb!= "" else ""}{equal}{"+".join(product)}{tb if tb!= "" else ""}'
            # to_write=f'{final_reaction:48} {reaction.speed[0]}{reaction.speed[1]:<+7.3f}{reaction.speed[2]:+11.3e}'
            to_write='{:<48}{:>+12.3e}{:>+9.3f}{:>+11.1f}'.format\
            (final_reaction, reaction.speed[0], reaction.speed[1], reaction.speed[2])
            if len(to_write)>80 :
                print(f'Reaction from base model too big ! :\n{to_write}')
                to_write=re.sub("\s+", " ", to_write)
                print(to_write)
            out.write((to_write+'\n').upper())
            # print(to_write,len(to_write))
            for sup in reaction.supinfo :
                out.write(('\t'+sup+'\n').upper())
    
            """ REACTIONS 1 """

        out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" REACTIONS":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
        if to_add != '' :
            
            for reaction in to_add[0] :

                if reaction.is_rev == True :
                    equal= '='
                elif reaction.is_rev == False :
                    equal= '=>'
                else :
                    print(f'ERROR, NO REV {reaction.__dict__}')

                if reaction.third_body != 0 :
                    tb=f"+{reaction.third_body}"
                    if reaction.supinfo != [] :
                        low_check=["low/" for x in reaction.supinfo if 'low/' in x or 'low /' in x]
                        tcheb_check=["tcheb/" for x in reaction.supinfo if 'tcheb/' in x or 'tcheb /' in x]
                        if len(low_check) >= 1 or len(tcheb_check) >= 1 :
                            tb=f"(+{reaction.third_body})"
                else :
                    tb=''

                reactant=[]
                product=[]
                # print(reaction.speed[1])
                for s, c in zip(reaction.species[0], reaction.coef[0]) :
                    if c == 1 :
                        c=''
                    elif int(c) == c :
                        c=int(c)
                    reactant.append(str(c) + str (s))

                for s, c in zip(reaction.species[1], reaction.coef[1]) :
                    if c == 1 :
                        c=''
                    elif int(c) == c :
                        c=int(c)
                    product.append(str(c) + str (s))

                final_reaction=f'{"+".join(reactant)}{tb if tb!= "" else ""}{equal}{"+".join(product)}{tb if tb!= "" else ""}'
                to_write='{:<48}{:>+12.3e}{:>+9.3f}{:>+11.1f}'.format\
                (final_reaction, reaction.speed[0], reaction.speed[1], reaction.speed[2])
                # to_write=f'{final_reaction:48} {reaction.speed[0]:<+12.3e}{reaction.speed[1]:<+7.3f}{reaction.speed[2]:+11.3e}'
                if len(to_write)>80 :
                    print(f'Reaction from added model too big !\n{to_write}')
                    to_write=re.sub("\s+", " ", to_write)
                    print(to_write)
                                     
                out.write((to_write+'\n').upper())
                # print(to_write,len(to_write))
                for sup in reaction.supinfo :
                    out.write(('\t'+sup+'\n').upper())
        out.write('end'.upper())
    
    """ THERMO """
    
    with open(outp+'.therm', 'w') as out :
        
        """ INTRODUCTION """

        out.write(custom)
        
        """ START """
        
        out.write(f'thermo\n{300:>10.3f}{1000:>10.3f}{5000:>10.3f}\n'.upper())
        
        out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" THERMO":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
        
        for species in model2.reactions[1] :
            name=model2.therm[species]
            ato_comp=name.atomic_composition
            string=''
            for  ato in ato_comp :
                string+=ato+' '*(5-len(str(ato_comp[ato]))-len(ato))+str(ato_comp[ato])
#           
            temp=name.temp_list
            
            out.write(f'{name.species:<16}{" "*8}{string:<20}G{temp[0]:>10.2f}{temp[2]:>10.2f}{temp[1]:>8.2f}{" "*6}1\n'.upper())
            for idx,item in enumerate(name.nasa_list) :
                if idx in [4,9,13] :
                    
                    if idx == 4 :
                        strsup='    2'
                    if idx == 9 :
                        strsup='    3'
                    if idx == 13 :
                        strsup=' '*15+'    4'
                        
                    out.write(f'{item:>+15.8e}{strsup}\n'.upper())
                else :
                    out.write(f'{item:>+15.8e}'.upper())
        
        out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" THERMO":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
        if to_add != '' :
            for species in to_add[1] :
                name=model1_translated.therm[species]
                ato_comp=name.atomic_composition
                string=''
                for  ato in ato_comp :
                    string+=ato+' '*(5-len(str(ato_comp[ato]))-len(ato))+str(ato_comp[ato])
    #           
                temp=name.temp_list

                out.write(f'{name.species:<16}{" "*8}{string:<20}G{temp[0]:>10.2f}{temp[2]:>10.2f}{temp[1]:>8.2f}{" "*6}1\n'.upper())
                for idx,item in enumerate(name.nasa_list) :
                    if idx in [4,9,13] :

                        if idx == 4 :
                            strsup='    2'
                        if idx == 9 :
                            strsup='    3'
                        if idx == 13 :
                            strsup=' '*15+'    4'

                        out.write(f'{item:>+15.8e}{strsup}\n'.upper())
                    else :
                        out.write(f'{item:>+15.8e}'.upper())
                    
        out.write('end'.upper())
        
    """ TRANSPORT """
    
    with open(outp+'.trans', 'w') as out :
        
        """ INTRODUCTION """
        
        out.write(custom)
        
        """ START """
        
        out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" TRANSPORT":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
        
        for species in  model2.reactions[1] :
            
            trans=model2.trans[species]
            out.write(f'{trans.species:<20}{trans.geo:>10.0f}{trans.lj_p:>10.3f}{trans.lj_cd:>10.3f}{trans.dp:>10.3f}{trans.p:>10.3f}{trans.rr:>10.3f}\n'.upper())
            
        out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" THERMO":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
        if to_add !='' and model1_translated.path_trans != '':
            for species in  to_add[1] :

                trans=model1_translated.trans[species]
                out.write(f'{trans.species:<20}{trans.geo:>10.0f}{trans.lj_p:>10.3f}{trans.lj_cd:>10.3f}{trans.dp:>10.3f}{trans.p:>10.3f}{trans.rr:>10.3f}\n'.upper())
        else :
            if model1_translated.path_trans == '' :
                out.write('NOT COMPLETE\nNOT COMPLETE\n')
        out.write('end'.upper())

def mechanism_writting_LLNL(outp, md1, md2, model1_translated, model2, order, custom='') :

    def write_species(out, model, mlist, nb) :
        if nb == 1 :
            for species in model.reactions[1] :
                out.write(f'{species}\n'.upper())

        elif nb == 2 :
            for species in model.reactions[1] :
                if species not in mlist :
                    out.write(f'{species}\n'.upper())
        else :
            raise Exception(f'Error in spe nb given : {nb}')

    def write_reactions(out, model, mlist, nb) :
        if nb == 1 :
            for reaction in model.reactions[2] :
                if reaction.is_rev == True :
                    equal= '='
                elif reaction.is_rev == False :
                    equal= '=>'
                else :
                    print(f'ERROR, NO REV {reaction.__dict__}')
    
                if reaction.third_body != 0 :
                    tb=f"+{reaction.third_body}"
                    if reaction.supinfo != [] :
                        low_check=["low/" for x in reaction.supinfo if 'low/' in x or 'low /' in x]
                        tcheb_check=["tcheb/" for x in reaction.supinfo if 'tcheb/' in x or 'tcheb /' in x]
                        if len(low_check) >= 1 or len(tcheb_check) >= 1 :
                            tb=f"(+{reaction.third_body})"
                else :
                    tb=''
    
                reactant=[]
                product=[]
                for s, c in zip(reaction.species[0], reaction.coef[0]) :
                    if c == 1 :
                        c=''
                    elif int(c) == c :
                        c=int(c)
                    reactant.append(str(c) + str (s))
    
                for s, c in zip(reaction.species[1], reaction.coef[1]) :
                    if c == 1 :
                        c=''
                    elif int(c) == c :
                        c=int(c)
                    product.append(str(c) + str (s))
    
                final_reaction=f'{"+".join(reactant)}{tb if tb!= "" else ""}{equal}{"+".join(product)}{tb if tb!= "" else ""}'
                # to_write=f'{final_reaction:48} {reaction.speed[0]}{reaction.speed[1]:<+7.3f}{reaction.speed[2]:+11.3e}'
                to_write='{:<48}{:>+12.3e}{:>+9.3f}{:>+11.1f}'.format\
                (final_reaction, reaction.speed[0], reaction.speed[1], reaction.speed[2])
                if len(to_write)>80 :
                    print(f'Reaction from base model too big ! :\n{to_write}')
                    to_write=re.sub("\s+", " ", to_write)
                    print(to_write)
                out.write((to_write+'\n').upper())
                # print(to_write,len(to_write))
                for sup in reaction.supinfo :
                    out.write(('\t'+sup+'\n').upper())
        elif nb == 2 :
            for reaction in model.reactions[2] :
                if reaction not in mlist :
                    if reaction.is_rev == True :
                        equal= '='
                    elif reaction.is_rev == False :
                        equal= '=>'
                    else :
                        print(f'ERROR, NO REV {reaction.__dict__}')
        
                    if reaction.third_body != 0 :
                        tb=f"+{reaction.third_body}"
                        if reaction.supinfo != [] :
                            low_check=["low/" for x in reaction.supinfo if 'low/' in x or 'low /' in x]
                            tcheb_check=["tcheb/" for x in reaction.supinfo if 'tcheb/' in x or 'tcheb /' in x]
                            if len(low_check) >= 1 or len(tcheb_check) >= 1 :
                                tb=f"(+{reaction.third_body})"
                    else :
                        tb=''
        
                    reactant=[]
                    product=[]
                    for s, c in zip(reaction.species[0], reaction.coef[0]) :
                        if c == 1 :
                            c=''
                        elif int(c) == c :
                            c=int(c)
                        reactant.append(str(c) + str (s))
        
                    for s, c in zip(reaction.species[1], reaction.coef[1]) :
                        if c == 1 :
                            c=''
                        elif int(c) == c :
                            c=int(c)
                        product.append(str(c) + str (s))
        
                    final_reaction=f'{"+".join(reactant)}{tb if tb!= "" else ""}{equal}{"+".join(product)}{tb if tb!= "" else ""}'
                    # to_write=f'{final_reaction:48} {reaction.speed[0]}{reaction.speed[1]:<+7.3f}{reaction.speed[2]:+11.3e}'
                    to_write='{:<48}{:>+12.3e}{:>+9.3f}{:>+11.1f}'.format\
                    (final_reaction, reaction.speed[0], reaction.speed[1], reaction.speed[2])
                    if len(to_write)>80 :
                        print(f'Reaction from added model too big ! :\n{to_write}')
                        to_write=re.sub("\s+", " ", to_write)
                        print(to_write)
                    out.write((to_write+'\n').upper())
                    # print(to_write,len(to_write))
                    for sup in reaction.supinfo :
                        out.write(('\t'+sup+'\n').upper())
        else :
            raise Exception(f'Error in rea nb given : {nb}')


    def write_therm(out, model, mlist, nb) :
        if nb == 1 :
            for species in model.reactions[1] :
                name=model.therm[species]
                ato_comp=name.atomic_composition
                string=''
                for  ato in ato_comp :
                    string+=ato+' '*(5-len(str(ato_comp[ato]))-len(ato))+str(ato_comp[ato])
    
                temp=name.temp_list
                
                out.write(f'{name.species:<16}{" "*8}{string:<20}G{temp[0]:>10.2f}{temp[2]:>10.2f}{temp[1]:>8.2f}{" "*6}1\n'.upper())
    
                for idx,item in enumerate(name.nasa_list) :
                    if idx in [4,9,13] :
                        
                        if idx == 4 :
                            strsup='    2'
                        if idx == 9 :
                            strsup='    3'
                        if idx == 13 :
                            strsup=' '*15+'    4'
                            
                        out.write(f'{item:>+15.8e}{strsup}\n'.upper())
                    else :
                        out.write(f'{item:>+15.8e}'.upper())
        elif nb == 2 :
            for species in model.reactions[1] :
                if species not in mlist :
                    name=model.therm[species]
                    ato_comp=name.atomic_composition
                    string=''
                    for  ato in ato_comp :
                        string+=ato+' '*(5-len(str(ato_comp[ato]))-len(ato))+str(ato_comp[ato])
        
                    temp=name.temp_list
                    
                    out.write(f'{name.species:<16}{" "*8}{string:<20}G{temp[0]:>10.2f}{temp[2]:>10.2f}{temp[1]:>8.2f}{" "*6}1\n'.upper())
        
                    for idx,item in enumerate(name.nasa_list) :
                        if idx in [4,9,13] :
                            
                            if idx == 4 :
                                strsup='    2'
                            if idx == 9 :
                                strsup='    3'
                            if idx == 13 :
                                strsup=' '*15+'    4'
                                
                            out.write(f'{item:>+15.8e}{strsup}\n'.upper())
                        else :
                            out.write(f'{item:>+15.8e}'.upper())

        else :
            raise Exception(f'Error in therm nb given : {nb}')


    def write_trans(out, model, mlist, nb) :
        if nb == 1 :
            for species in model.reactions[1] :
                trans=model.trans[species]
                out.write(f'{trans.species:<20}{trans.geo:>10.0f}{trans.lj_p:>10.3f}{trans.lj_cd:>10.3f}{trans.dp:>10.3f}{trans.p:>10.3f}{trans.rr:>10.3f}\n'.upper())
        elif nb == 2 :
            for species in model.reactions[1] :
                if species not in mlist :
                    trans=model.trans[species]
                    out.write(f'{trans.species:<20}{trans.geo:>10.0f}{trans.lj_p:>10.3f}{trans.lj_cd:>10.3f}{trans.dp:>10.3f}{trans.p:>10.3f}{trans.rr:>10.3f}\n'.upper())

        else :
            raise Exception(f'Error in trans nb given : {nb}')


    """ MECH """
    
    with open(outp+'.mech', 'w') as out :

        """ INTRODUCTION """

        out.write(custom)
        out.write('!'+order+'\n')

        """ ELEMENTS """

        out.write('elements\n'.upper())
        elements=set(model1_translated.reactions[0]+model2.reactions[0])
        for element in elements :
            out.write(f'{element}\n'.upper())

        """ SPECIES """
        if order[0] == '1' :
            """ SPECIES """
            out.write(f'end\n{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" SPECIES":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\nspecies\n'.upper())
            write_species(out, model1_translated, model2.reactions[1], 1)
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" SPECIES":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_species(out, model2, model1_translated.reactions[1], 2)

            """ REACTIONS """
            out.write(f'end\n{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" REACTIONS":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\nreactions\n'.upper())
            write_reactions(out, model1_translated, model2.reactions[2], 1)
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" REACTIONS":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_reactions(out, model2, model1_translated.reactions[2], 2)

        elif order[0] == '2' :
            """ SPECIES """
            out.write(f'end\n{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" SPECIES":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\nspecies\n'.upper())
            write_species(out, model2, model1_translated.reactions[1], 1)
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" SPECIES":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_species(out, model1_translated, model2.reactions[1], 2)

            """ REACTIONS """
            out.write(f'end\n{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" REACTIONS":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\nreactions\n'.upper())
            write_reactions(out, model2, model1_translated.reactions[2], 1)
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" REACTIONS":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_reactions(out, model1_translated, model2.reactions[2], 2)

        else :
            raise Exception(f'Error in mech order given : {order}')
        out.write('end'.upper())



    """ THERMO """

    with open(outp+'.therm', 'w') as out :

        """ INTRODUCTION """

        out.write(custom)
        out.write('!'+order+'\n')

        """ START """
        
        out.write(f'thermo\n{300:>10.3f}{1000:>10.3f}{5000:>10.3f}\n'.upper())

        if order[1] == '1' :
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" THERMO":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_therm(out, model1_translated, model2.reactions[1], 1)
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" THERMO":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_therm(out, model2, model1_translated.reactions[1], 2)

        elif order[1] == '2' :
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" THERMO":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_therm(out, model2, model1_translated.reactions[1], 1)
            out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" THERMO":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
            write_therm(out, model1_translated, model2.reactions[1], 2)

        else :
            raise Exception(f'Error in therm order given : {order}')

        out.write('end'.upper())
        
    """ TRANSPORT """
    if model1_translated.path_trans != '' and model2.path_trans != '' :
        with open(outp+'.trans', 'w') as out :
            
            """ INTRODUCTION """
            
            out.write(custom)
            out.write('!'+order+'\n')
    
            """ START """
    
            if order[1] == '1' :
                out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" TRANSPORT":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
                write_trans(out, model1_translated, model2.reactions[1], 1)
                out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" TRANSPORT":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
                write_trans(out, model2, model1_translated.reactions[1], 2)
    
            elif order[1] == '2' :
                out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md2+" TRANSPORT":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
                write_trans(out, model2, model1_translated.reactions[1], 1)
                out.write(f'{"!"*80}\n!!!!!{" "*70}!!!!!\n!!!!!{md1+" TRANSPORT":^70}!!!!!\n!!!!!{" "*70}!!!!!\n{"!"*80}\n'.upper())
                write_trans(out, model1_translated, model2.reactions[1], 2)
    
            else :
                raise Exception(f'Error in trans order given : {order}')


            out.write('end'.upper())