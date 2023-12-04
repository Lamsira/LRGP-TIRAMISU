import numpy as np
import os
import re
from collections import Counter

class Model :
    def __init__(self, path_mech, path_therm, out, path_trans=''):
        self.path_mech = path_mech
        self.path_therm = path_therm
        self.out = out
        self.reactions = self.read_model()
        self.therm = self.read_therm()
        if path_trans != '' :
            self.path_trans = path_trans
            self.trans = self.read_trans()
        else :
            self.path_trans = ''
        self.check = self.atomic_conservation()
            
    def select_reactions(self, which=[]) :
        if which == [] :
            return([x for x in self.reactions[2]])
        else :
            to_return=[]
            for typ in which :
                to_return+=[x for x in self.reactions[2] if typ in x.reaction_types]
                # to_return+=[x.reaction_types for x in self.reactions[2]]
                # print(to_return)
            return(to_return)
    
    def read_model(self) :

        reactions=[]

        # print('\n',f'\n {self.path_mech}')
        
        """---------------------------------------------------------------------
                                   READING MECHANISM   
        ---------------------------------------------------------------------"""
        
        with open (self.path_mech,'r') as open_mech :
            read_mech=open_mech.readlines()
        index_line=0
        
        """---------------------------------------------------------------------
                                   PURGING MECHANISM  
        ---------------------------------------------------------------------"""

        #### Replacing float number like .x by 0.x (.3 into 0.3)
        purged_mech=[re.sub('\+\.','+0.',re.sub('=\.','=0.',
        re.sub('>\.','>0.' ,re.sub('\s\.', ' 0.', line)))) 
                     for line in read_mech]
        
        #### Removing every comments (!***) and blank lines
        purged_mech=[line.strip() for line in purged_mech if line.strip() != ''
                     and line.strip()[0] != '!']
        
        #### Removing every useless white space (in reactions and multiples 
        #### white spaces)
        # purged_mech=[re.sub("\s+", " ", line.split('!')[0]) 
        #              for line in purged_mech]     
        purged_mech=[re.sub("\s+", " ", line) 
                      for line in purged_mech] 
        
        purged_mech=[re.sub('> ','>',re.sub(
            ' <', '<',re.sub('[= ][= ]+','=',re.sub('\+ ', '+', line)))) 
                     for line in purged_mech]
        
        #### Split Comment
        purged_mech=[line.split('!') for line in purged_mech]
        
        for entireline in purged_mech.copy() :
            line=entireline[0]
            listplus=(re.findall(' \+', line))
            tmpline=line
            for iterplus in listplus :
                tmpnumb=tmpline[tmpline.index(' +')+2:].split(' ')[0]
                try : 
                    float(tmpnumb) 
                    break
                except ValueError :
                    tmpline=tmpline.replace(' +', '+', 1)
            purged_mech[purged_mech.index(entireline)][0]=tmpline
            
        for entireline in purged_mech.copy() :
            
            line=entireline[0]
            listplus=(re.findall(' \(\+', line))
            tmpline=line
            for iterplus in listplus :
                tmpnumb=tmpline[tmpline.index(' (+')+3:].split(' ')[0]
                try : 
                    float(tmpnumb)
                    break
                except ValueError :
                    tmpline=tmpline.replace(' (+', '(+', 1)
            purged_mech[purged_mech.index(entireline)][0]=tmpline
        
        
        for ind,line in enumerate(purged_mech.copy()) :
            purged_mech[ind][0]=purged_mech[ind][0].lower()
            
        for x in purged_mech :
            if '!' in x[0] : 
                print(f'warning {x[0]}')
      
        """------------------------------------------------------------------"""

        len_mech = len(purged_mech)
        elements = []
        species = []

        index_line=0

        """---------------------------------------------------------------------
                                   READING ELEMENTS  
        ---------------------------------------------------------------------"""
        print(f'\n------------------------------{os.path.basename(self.path_mech)}------------------------------\n')
        while purged_mech[index_line][0].strip() != 'elements' :
            index_line+=1
        print(f'LINE FOR START OF ELEMENTS : {purged_mech[index_line][0]}')
        index_line+=1
        while purged_mech[index_line][0].strip() != 'end' :
            tmplist=purged_mech[index_line][0].split()
            for ele in tmplist :
                elements.append(ele)
            index_line+=1
        print(f'LINE FOR END OF ELEMENTS : {purged_mech[index_line][0]}')
        index_line+=1

        #### check if blank
        [print('BLANK ELEMENT !!!!!!!') for ele in elements if ele=='']
        
        """---------------------------------------------------------------------
                                   READING SPECIES  
        ---------------------------------------------------------------------"""

        while purged_mech[index_line][0].strip() != 'species' :
            index_line+=1
        print(f'LINE FOR START OF SPECIES : {purged_mech[index_line][0]}')
        index_line+=1
        while 'end' not in purged_mech[index_line][0].strip() != 'end' :
            tmplist=purged_mech[index_line][0].split()
            for spe in tmplist :
                if spe in species :
                    print('Mech duplicate :',spe, index_line, purged_mech[index_line])
                else :
                    species.append(spe)
            index_line+=1
        print(f'LINE FOR END OF SPECIES : {purged_mech[index_line][0]}')
        index_line+=1


        #### check if blank
        if self.path_mech!=self.path_therm : 
            [print('BLANK SPECIES !!!!!!!') for spe in species if spe=='']
        else :
            print('Mech == Therm')
        """---------------------------------------------------------------------
                                   READING REACTIONS  
        ---------------------------------------------------------------------"""
        
        while purged_mech[index_line][0].strip() != 'reactions' and 'maxsp' not in purged_mech[index_line][0]:
            index_line+=1
        
        if 'maxsp' in purged_mech[index_line][0] :
            print('WARNING : MAXSP HERE')
        print(f'LINE FOR START OF REACTIONS : {purged_mech[index_line][0]}')
        
        count_of_not_int = 0
        index_line+=1

        from_model=''

        initiation=0
        
        identifier=0
        
        while 'end' not in purged_mech[index_line][0] :
            line=purged_mech[index_line][0]

            #### if '=' in line, it's a reaction line
            if '=' in line :
               
                if initiation!=0 :
                    reaction_types=[]
                    
                    for info in sup_i :
                        
                        if 'plog' in info and 'plog' not in reaction_types :
                            reaction_types.append('plog')
                        elif ('low' in info or 'troe' in info) and 'low/troe' not in reaction_types :
                            if 'low' in info and 'low' not in reaction_types : 
                                reaction_types.append('low')
                            else :
                                if 'low' in reaction_types :
                                    reaction_types.remove('low')
                                reaction_types.append('low/troe')
                        elif 'dup' in info and 'duplicate' not in reaction_types :
                            reaction_types.append('duplicate')
                        elif 'rev' in info and 'rev_keyword' not in reaction_types :
                            reaction_types.append('rev_keyword')
                        elif 'cheb' in info and 'cheb' not in reaction_types :
                            reaction_types.append('cheb')
                        
                        else :
                            nb_occurence=0
                            temp_inf=info.replace('/', ' ')
                            temp_inf=temp_inf.split()
                            for word in temp_inf :
                                if word in ['plog', 'low', 'dup', 'rev', 'cheb'] :
                                    nb_occurence+=1
                                    
                            if nb_occurence>= 2 :
                                print(f'WARNING IN SUPINFO : {info}')
                            if nb_occurence == 0 :
                                for word in temp_inf :
                                    if word in ['colleff','exci','fit1','ford','high','jan','lt','mome','rlt','rord', 'sri', 'tdep', 'usrprog', 'units'] :
                                        raise Exception(f'KeyWordError : {info}')
                        
                    
                        
                    reactions.append(Reactions([reactant,product],
                                stoechio_coef, rev, speed, sup_i, third_body, comment, reaction_types, identifier, from_model, 
                                os.path.splitext(os.path.basename(self.path_mech))[0]))
                    identifier+=1
                    #appel classe
               
                if len(purged_mech[index_line])>1 :
                    try :
                        comment='!'.join(purged_mech[index_line][1:])
                        if '*:_:*' in comment :
                            from_model = comment.split('*:_:*')[1]
                        else :
                            from_model=''
                    except TypeError :
                        print(f'INT OBJECT IS NOT SUBSCRIPTABLE : {purged_mech[index_line][1:]}')
                        comment='!'.join([str(x) for x in purged_mech[index_line][1:]])
                        print(comment)
                else :
                    comment=''
                    from_model=''
                
                
    
                sup_i=[]

            #### splitting reaction line into [reaction, A, n, Ea]
                third_body=0
                line=line.split()
                if len(line)!=4 :
                    line=[''.join(line[0:-3]),line[-3],line[-2],line[-1]]
                    print(f'New line = {line}')
                    # while len(line) != 4 :
                    #     print('\n##### ERROR :',line,'\n')
                    #     line[0]=line[0]+line[1]
                    #     del line[1]
                    #     while True :
                    #         conti=input(f'Is the reaction correct ? (y/n)?\n{line}')
                    #         print([conti])
                    #         print(len(line))
                    #         if conti == 'n' :
                    #             raise Exception(f'Reaction error : {line}')
                    #         elif conti == 'y' :
                    #             print('should break')
                    #             break
                    #         else :
                    #             print(f'Error input : {conti}')

                reaction=line[0]
                
                if '.' in reaction : count_of_not_int+=1 #print('Not int species : ', reaction)
                speed=[]
                for speed_coef in line[1:] :
            #### if one of speed_coef is badly written (not a number), append a
            #### 0 instead. (created due to some numbers written 0.0000+03 
            #### instead of 0.0000E+03 in LLNL) 
                    try :
                        speed.append(float(speed_coef.replace('d', 'e')))
                        if 'd' in speed_coef :
                            print(f'WARNING IN THE MODEL {speed_coef} (should work)')
                    except ValueError :
                        speed.append(0)
                        print(f'One of them is transformed into 0 : {line[1:]}')
            #### check if the reaction is reversible or not and split into
            #### [reactant, product] (and check if a reaction as been written 
            #### "in reverse" with a <=)
                if '<=>' in reaction : 
                    rev=True
                    reaction=reaction.split('<=>')
                elif '<=' in reaction :
                    rev=False
                    reaction=reaction.split('<=')

                    print('Reaction written in a wrong way ....', reaction)
                    reaction=[reaction[1],reaction[0]]
                    print('Reaction written in a good way ....', reaction)

                elif '=>' in reaction :
                    rev=False
                    reaction=reaction.split('=>')

                else :
                    rev=True
                    reaction=reaction.split('=')

            #### protect third body (like (+M) or (+H2O)) 
                if '(+' in reaction.copy()[0] :
                    third_body = reaction[0][reaction[0].index('(+')+2 : -1]
                    reaction[0] = reaction[0][:reaction[0].index('(+')]
                    reaction[1] = reaction[1][:reaction[1].index('(+')]

            #### split species

                reactant=reaction[0].split('+')

                product=reaction[1].split('+')
            #### removing third body disguised in species ....
                if 'm' in reactant :
                    if 'm' not in product or third_body != 0 :
                        print('Third_body ERROR !! Either third_body'\
                              'protection as failed or a third_body is only in'\
                              'one side of the reaction (reactant) ',
                              third_body, reaction)
                    else :
                        third_body = 'm'
                        reactant.remove('m')
                        product.remove('m')

                if 'm' in product :
                    print('Third_body ERROR !!A third_body is only in one side'\
                          'of the reaction (product)', third_body, reaction)


            #### reading the stoechiometric coefficient before each species 
                stoechio_coef=[[],[]]

            #### for reactant
                for rea in reactant.copy() :
                    stoe_rea=''
                    tmp_rea_index=0
                    while rea[tmp_rea_index].isnumeric() == True \
                    or rea[tmp_rea_index]=='.' :

                        stoe_rea+=rea[tmp_rea_index]
                        tmp_rea_index+=1
                    
                    if stoe_rea=='' :
                        stoe_rea=1
                    else :
#                         print(reactant)
                        reactant[reactant.index(rea)]=rea[tmp_rea_index:]
#                         print(reactant)
                    stoechio_coef[0].append(float(stoe_rea))
            #### for product
                for prod in product.copy() :
                    stoe_prod=''
                    tmp_prod_index=0
                    try :
                        while prod[tmp_prod_index].isnumeric() == True \
                        or prod[tmp_prod_index]=='.' :
                            
                            stoe_prod+=prod[tmp_prod_index]
                            tmp_prod_index+=1
                    except IndexError :
                        raise Exception(f'IndexError in : {prod} {product} {reaction} {tmp_prod_index}')
                        
                    if stoe_prod=='' :
                        stoe_prod=1
                    else :
#                         print(product)
                        product[product.index(prod)]=prod[tmp_prod_index:]
#                         print(product)
                    stoechio_coef[1].append(float(stoe_prod))

            #### removing multiple species occurence and summing their 
            #### stoechiometric coeffiscient for reactant
                stop_loop1=0
                while [1 for spe in Counter(reactant) if \
                       Counter(reactant)[spe]>1] != [] :

                    stop_loop1+=1
                    tmp_counter= [spe for spe in Counter(reactant) \
                                  if Counter(reactant)[spe]>1]
                    
                    species_idx = [i for i, item in enumerate(reactant) \
                                   if tmp_counter[0] == item]
                    
                    stoechio_coef[0][species_idx[0]]=sum(\
                        [x for num,x in enumerate(stoechio_coef[0])\
                         if num in species_idx])
                    
                    stoechio_coef[0]=[x for num,x in \
                    enumerate(stoechio_coef[0]) if num not in species_idx[1:]]
                    reactant=[x for num,x in enumerate(reactant) \
                              if num not in species_idx[1:]]

                    if stop_loop1>=4 :
                        print('FATAL ERROR IN WHILE LOOP FOR stoechio_coef '+
                              'MERGE')
                        break

            #### for product
                stop_loop2=0        
                while [1 for spe in Counter(product) 
                       if Counter(product)[spe]>1] != [] :

                    stop_loop2+=1
                    tmp_counter= [spe for spe in Counter(product) 
                                  if Counter(product)[spe]>1]
                    
                    species_idx = [i for i, item in enumerate(product) 
                                   if str(tmp_counter[0]) == str(item)]
#                     print('\n')
#                     for item in product :
#                         print([tmp_counter[0]], [item], re.fullmatch(tmp_counter[0], item))
                    
#                     print(f'{reactant} = {product},\nSt_coef {stoechio_coef[1]},\nSp_idx {species_idx},\ntmp_count {tmp_counter}')
                    
                    stoechio_coef[1][species_idx[0]]=sum([x for num,x in enumerate(stoechio_coef[1])\
                                                          if num in species_idx])
                    
                    stoechio_coef[1]=[x for num,x in enumerate(stoechio_coef[1]) \
                                      if num not in species_idx[1:]]
                    
                    product=[x for num,x in enumerate(product) if num not in species_idx[1:]]

                    if stop_loop2>=4 :
                        print('FATAL ERROR IN WHILE LOOP FOR stoechio_coef MERGE')
                        break


                initiation=1

            else :
                sup_i.append(line)
            index_line+=1
        
        for info in sup_i :
            if 'plog' in info and 'plog' not in reaction_types :
                reaction_types.append('plog')
            if ('low' in info or 'troe' in info) and 'low/troe' not in reaction_types :
                reaction_types.append('low_troe')
            if 'dup' in info and 'duplicate' not in reaction_types :
                reaction_types.append('duplicate')
            if 'rev' in info and 'rev_keyword' not in reaction_types :
                reaction_types.append('rev_keyword')
            if 'cheb' in info and 'cheb' not in reaction_types :
                reaction_types.append('cheb')
                
        reactions.append(Reactions([reactant,product], stoechio_coef, rev, speed, sup_i, third_body, comment, reaction_types, identifier,
                                   from_model, os.path.splitext(os.path.basename(self.path_mech))[0]))
    
        #### if one species start with "(", the code doesn't work anymore
        for spe in species :
            if spe[0] == '(' :
                print('NEED TO CHANGE THE CODE BECAUSE "(" AT THE START OF A SPECIES : ',spe)
        
        ############## WRITING LOGS #################   
        root, extension = os.path.splitext(self.path_mech)
        with open(self.out+os.path.basename(root)+"_mech.log", 'w') as log:
        # with open(self.path_mech+'.log', 'w') as log:
            for rea in reactions :
                log.write(f'Species : {rea.species} Coef : {rea.coef} is_rev : {rea.is_rev}'
            +f' third_body : {rea.third_body} len speed&sup_i {len(rea.speed)} {len(rea.supinfo)}\n')
        
        
        is_in_rea=set()
        for reaction in reactions :
            for side in reaction.species :
                for spe in side :
                    is_in_rea.add(spe)
        
        for spe in species.copy() :
            if spe not in is_in_rea and spe not in elements :
                print(f'{spe} NOT IN REACTIONS')
                species.remove(spe)
        ############## RETURNING RESULTS ################# 
        print(f'\nNot-int species in {count_of_not_int} reactions')
        return(elements, species, reactions)
    
    """-------------------------------------------------------------------------
    ----------------------------------------------------------------------------
    
                                       THERMO  
                                   
    ----------------------------------------------------------------------------
    -------------------------------------------------------------------------"""

    
    def read_therm(self) :
        
        not_in_mech=[]
        therm_duplicate=[]
        therm_data={}
        
        """---------------------------------------------------------------------
                                   READING THERM FILE  
        ---------------------------------------------------------------------"""
    
        with open(self.path_therm, 'r') as therm :
            read_therm=therm.readlines()
            
        """---------------------------------------------------------------------
                                   PURGING THERM FILE  
        ---------------------------------------------------------------------"""

        purged_therm=[line for line in read_therm if line.strip() != '' and line.strip()[0] != '!']
        purged_therm=[line.split('!')[0] for line in purged_therm]
        index=0
        while 'thermo' not in purged_therm[index].lower() :
            index+=1
            
        if index != 0 :
            print(f'Number of lines before THERMO or THERMO ALL : {index}')
            purged_therm=purged_therm[index:]
        if 'thermo' not in purged_therm[0].lower() :
            print('ERROR IN THERMO')
        temperature_save=purged_therm[1]
        purged_therm=purged_therm[2:]
        purged_therm=[x.lower() for x in purged_therm]
        index_line=0

        """---------------------------------------------------------------------
                                   SAVING BASE TEMP  
        ---------------------------------------------------------------------"""

        base_temp=temperature_save.strip().split()
        base_temp=list(map(float, base_temp))
        base_temp=sorted(base_temp)

        multiple_list=[]
        
        """---------------------------------------------------------------------
                                READING SPECIES AND NASA 
        ---------------------------------------------------------------------"""
        # print(self.reactions[1])
        while 'end' not in purged_therm[index_line][:10] :
            line=purged_therm[index_line]
            line=line[:80]  

        #### reading species
            species=line.strip()[:line.strip().index(' ')]
        #### counting atomes of each species

            if species not in multiple_list and line[79]=='1' and species in self.reactions[1] :
                multiple_list.append(species)

                
                hidden_formula=line[24:44]
                hidden_formula=re.sub('\.',' ',hidden_formula)
                hidden_formula=[hidden_formula[x:x+5] for x in range(0,len(hidden_formula),5)]

                for num,atomedata in enumerate(hidden_formula) :
                    hidden_formula[num]=atomedata.strip()

                hidden_formula=set(hidden_formula)
                hidden_formula=list(hidden_formula)

                if '0' in hidden_formula or '00' in hidden_formula :
                    hidden_formula=list(filter(lambda a: a != '0' and a != '00', hidden_formula))

                if '' in hidden_formula :
                    hidden_formula=list(filter(lambda a: a != '', hidden_formula))

                tmp_elements_list={}
                
                # print('HD :',hidden_formula)
                for num,atome in enumerate(hidden_formula) :
                    element=atome.split()[0]
                    # print('atome :',atome)
                    if atome.split()[1] == '0':
                        continue
                    else :
                        nb=atome.split()[1]
                        tmp_elements_list[element]=int(nb)

     

            #### actual temp integration
                test_temp=0
                actual_temp=line[45:76].split()
                actual_temp=list(map(float, actual_temp))
                actual_temp=sorted(actual_temp)

                if len(actual_temp) == 2 :
                    test_temp=actual_temp.copy()
                    actual_temp.append(base_temp[1])
                    
                elif len(actual_temp) == 0 :
                    actual_temp=base_temp.copy()
                    
                elif len(actual_temp) == 1 : 
                    print(f'ERROR : TOO FEW THERMO TEMPERATURE {len(actual_temp)} FOR {species}')
                actual_temp=sorted(actual_temp)
                
                if actual_temp[1]-actual_temp[0]<300 or actual_temp[2]-actual_temp[1]<500 :
                    print(f'Thermo temperature diff too low !!! : {species}'\
                          f' {actual_temp[1]-actual_temp[0]}, {actual_temp[2]-actual_temp[1]}'\
                          , actual_temp, test_temp)
                    
                split_nasa=[]

            #### NASA reading
                for data_line in range(1, 4) :

                    if data_line==3 :
                        line_nasa=purged_therm[index_line+data_line][:60]
                        line_nasa=re.sub('e 0', 'e+0', line_nasa)
                        split_nasa=split_nasa+[float(line_nasa[x:x+15]) for x in range(0, len(line_nasa), 15)]
                        index_line+=4

                        if len(split_nasa) != 14 : print('SPLIT NASA LENGTH = ', len(split_nasa))
                        therm_data[species]=(Thermo(species, tmp_elements_list, actual_temp, split_nasa))
                        

                    else :
                        line_nasa=purged_therm[index_line+data_line][:75]
                        line_nasa=re.sub('e 0', 'e+0', line_nasa)
                        split_nasa=split_nasa+[float(line_nasa[x:x+15]) for x in range(0, len(line_nasa), 15)]
            else :
                if species not in self.reactions[1] :
                    # print(f'THERMO DATA NOT IN MECH : {species}')
                    not_in_mech.append(f'THERMO DATA NOT IN MECH : {species}\n')
                else :
                    # print('THERM DUPLICATE : ', species)
                    therm_duplicate.append(f'THERM DUPLICATE : {species}\n')

                index_line+=1

                for k_iter in range(index_line, len(purged_therm)-1) :

                    if len(purged_therm[index_line][0:80]) < 80 :
                        index_line+=1
                    elif purged_therm[index_line][0:80][-1] != '1' :

                        index_line+=1
                    else :
                        break

        
        root, extension = os.path.splitext(self.path_therm)
        with open(self.out+os.path.basename(root)+"_therm.log", 'w') as log:
        # with open(self.path_therm+'.log', 'w') as log :
            for therm in therm_data :
                log.write(f'Temp : {therm_data[therm].temp_list} {therm_data[therm].species} {therm_data[therm].atomic_composition} \n')
            log.write("".join(not_in_mech))
            log.write("".join(therm_duplicate))
                
        return(therm_data)

    
    """-------------------------------------------------------------------------
    ----------------------------------------------------------------------------
    
                                    TRANSPORT  
                                   
    ----------------------------------------------------------------------------
    -------------------------------------------------------------------------"""
                
    def read_trans(self) :
        
        trans_data={}
        
        """---------------------------------------------------------------------
                                   READING TRANSPORT   
        ---------------------------------------------------------------------"""
        
        with open (self.path_trans, 'r', encoding='latin-1') as open_trans :
            read_trans=open_trans.readlines()
        index_line=0
        
        """---------------------------------------------------------------------
                                   PURGING TRANSPORT  
        ---------------------------------------------------------------------"""

        
        #### Removing every comments (!***) and blank lines
        purged_trans=[line.strip() for line in read_trans if line.strip() != ''
                     and line.strip()[0] != '!']

        #### Removing every useless white space (in reactions and multiples 
        #### white spaces)
        purged_trans=[re.sub("\s+", " ", line.split('!')[0]) 
                     for line in purged_trans]     
            
        purged_trans=[x.lower() for x in purged_trans]
        if 'end' in purged_trans[-1] :
            purged_trans=purged_trans[:-1]
            
        not_in_mech=[]
        
        for data_line in purged_trans :
            data = data_line.split()
            if len(data) != 7 :
                print(data_line,'\n', data)
                
            else :
                if data[0] in self.reactions[1] :
                    trans_data[data[0]]=(Transport(data[0], float(data[1]), float(data[2]), float(data[3]), float(data[4]), float(data[5]), float(data[6])))
                else :
                    # print(f'TRANS DATA NOT IN MECH : {data[0]}')
                    not_in_mech.append(f'TRANS DATA NOT IN MECH : {data[0]}\n')
                    
        root, extension = os.path.splitext(self.path_trans)
        with open(self.out+os.path.basename(root)+"_trans.log", 'w') as log:
            for trans in trans_data :
                t_d=trans_data[trans]
                log.write(f'{t_d.species} : {t_d.geo} {t_d.lj_p} {t_d.lj_cd} {t_d.dp} {t_d.p} {t_d.rr}\n')
            log.write("".join(not_in_mech))
        return(trans_data)
                
    
    """-------------------------------------------------------------------------
    ----------------------------------------------------------------------------
                            READING SPECIES AND NASA 
    ----------------------------------------------------------------------------
    -------------------------------------------------------------------------"""

    def atomic_conservation(self) :
        iso_reactions={}
        for reaction in self.reactions[2] :
#             print(reaction.__dict__)
            reactant=[reaction.species[0],reaction.coef[0]]
            product=[reaction.species[1],reaction.coef[1]]
            rea_dict={}
            prod_dict={}
            for num,rea in enumerate(reactant[0]) :
                if rea != 'hv' :
                    rea_tmp=self.therm[rea].atomic_composition
                    for atome in rea_tmp :
                        if atome not in rea_dict :
                            rea_dict[atome]=rea_tmp[atome]*reactant[1][num]
                        else :
                            rea_dict[atome]=rea_dict[atome]+rea_tmp[atome]*reactant[1][num]

            for num,prod in enumerate(product[0]) :
#                 print(prod, num, product)
                if prod != 'hv' :
                    prod_tmp=self.therm[prod].atomic_composition
                
                    for atome in prod_tmp :
                        if atome not in prod_dict :
                            prod_dict[atome]=prod_tmp[atome]*product[1][num]
                        else :
                            prod_dict[atome]=prod_dict[atome]+prod_tmp[atome]*product[1][num]
            
            
            balance_error=False
            if rea_dict != prod_dict :
                is_different=False
                for atome in rea_dict :
                
                    if rea_dict[atome]*0.99 >= prod_dict[atome] or rea_dict[atome] <= 0.99*prod_dict[atome] :
                        is_different=True
                        balance_error=True
                if is_different==True :
                    print('ERREUR : ', rea_dict, prod_dict)
                    print(reaction.__dict__)
                # else : 
                #     print('Rounded : ', rea_dict, prod_dict, reaction.__dict__)
            if balance_error==True :
                print(f'\n{"#"*100}\n{"#"*35} FATAL ERROR IN ATOMS BALANCE {"#"*35}\n{"#"*100}')
                reaction.print_reaction()
                raise Exception('FATAL ERROR IN ATOMS BALANCE')
            
            order=sorted([x+str(int(rea_dict[x])) for x in rea_dict])
            iso=''.join(order)
            if iso in iso_reactions : 
                iso_reactions[iso].append(reaction)
            else :
                iso_reactions[iso] = [reaction]
            reaction.iso=iso
        return(iso_reactions)
                
    def micro_mech(self, list_species) :
        reaction_list=[]
        ind_rea_list=[]
        inreaction_species=[]
        for species in list_species :
            for reaction in self.reactions[2] :
                sperea=[*reaction.species[0], *reaction.species[1]]
                if species in sperea and reaction.identifier not in ind_rea_list :
                    reaction_list.append(reaction)
                    for spe in sperea :
                        if spe not in inreaction_species :
                            inreaction_species.append(spe)
        return(inreaction_species,reaction_list)
        
        
class Thermo :
    def __init__(self, species, atomic_composition, temp_list, nasa_list) :
        self.species=species
        self.atomic_composition=atomic_composition
        self.temp_list=temp_list
        self.nasa_list=nasa_list
    
    ######FUCKING DELETE THIS

    ###########################
    def __repr__(self):
        return self.species

    """-------------------------------------------------------------------------
                       NASA COMPUTED THERMODYNAMIC FUNCTIONS
    -------------------------------------------------------------------------"""
    
    def cp(nasa_coeffs, temp):
        a = nasa_coeffs
        t = temp
        return a[0] + a[1]*t + a[2]*t**2 + a[3]*t**3 + a[4]*t**4

    def h(nasa_coeffs, temp):
        a = nasa_coeffs
        t = temp
        return a[0] + a[1]*t/2 + a[2]*t**2/3 + a[3]*t**3/4 + a[4]*t**4/5 + a[5]/t

    def s(nasa_coeffs, temp):
        a = nasa_coeffs
        t = temp
        return a[0]*np.log(t) + a[1]*t + a[2]*t**2/2 + a[3]*t**3/3 + a[4]*t**4/4 + a[6]
    
    def calc_therm(self, f):
        n=350
        low_temps = np.linspace(self.temp_list[0], self.temp_list[1], n)
        high_temps = np.linspace(self.temp_list[1], self.temp_list[2], n)
        high_temp_f = f(self.nasa_list[:7], high_temps)
        low_temp_f = f(self.nasa_list[7:], low_temps)
        return np.hstack((low_temps, high_temps)), np.hstack((low_temp_f, high_temp_f))
    
    def G(self):
        temp, a = self.calc_therm(Thermo.s)
        temp, b = self.calc_therm(Thermo.h)
        return temp, a - b
    
    def new_calc_therm(self, temp, f):
        if self.temp_list[1]<temp[0] :
            if len(temp) != 2 :    
                print(f'Could happend with other, len(temp) should be 2 : {len(temp)}')
            mid_temp_index = 0
        else :
            mid_temp_index = np.argmin(np.abs(temp - self.temp_list[1]))
        high_temps = temp[mid_temp_index:]
        low_temps = temp[:mid_temp_index]
        high_temp_cp = f(self.nasa_list[:7], high_temps)
        if mid_temp_index != 0 :
            low_temp_cp = f(self.nasa_list[7:], low_temps)
            return np.hstack((low_temp_cp, high_temp_cp))
        else :
            return np.hstack((high_temp_cp))
    
    def G_new(self, temp):
        return self.new_calc_therm(temp, Thermo.s) - self.new_calc_therm(temp, Thermo.h)
    
    def plot_data(self) :
        # n=350
       
        T,data_cp= self.calc_therm(Thermo.cp)
        T,data_h= self.calc_therm(Thermo.h)
        T,data_s= self.calc_therm(Thermo.s)
        h_dif=[Thermo.h(self.nasa_list[7:], 300)*1.987*300, Thermo.h(self.nasa_list[7:], 500)*1.987*500, Thermo.h(self.nasa_list[7:], 700)*1.987*700]
        
        true_cp=[x*1.987 for x in data_cp]
        true_h=[x*1.987*t for x,t in zip(data_h,T)]
        true_s=[x*1.987 for x in data_s]
        true_g=[h-t*s for h,t,s in zip(true_h,T,true_s)]
        true_hk=[x/1000 for x in true_h]
        true_mg10k=[-x/10000  for x in true_g]
        return {'T':T, 'cp/R':data_cp,
                'h/RT':data_h, 's/R':data_s, 'h_value':h_dif, 'cp' : true_cp,
                'h' : true_h, 's' : true_s, 'g' : true_g, 'hk' : true_hk, 'mg10k' : true_mg10k}
        


class Reactions :
    def __init__(self, species, coef, is_rev, speed, supinfo, third_body, comment, reaction_types, identifier, from_model, model, iso=0) :
        self.species = species
        self.coef = coef
        self.is_rev = is_rev
        self.speed = speed
        self.supinfo = supinfo
        self.third_body = third_body
        self.comment = comment
        self.reaction_types = reaction_types
        self.identifier=identifier
        self.from_model=from_model
        self.model=model
        self.iso=iso
    
    def __eq__(self, other) :
        tmp_set1=[set(x) for x in self.species]
        tmp_set2=[set(x) for x in other.species]
        tmp_coef1=[set(x) for x in self.coef]
        tmp_coef2=[set(x) for x in other.coef]
        
        if self.third_body == other.third_body :
            if (tmp_set1 == tmp_set2 or tmp_set1 == tmp_set2[::-1]) and (tmp_coef1 == tmp_coef2 or tmp_coef1==tmp_coef2[::-1])  :
                if self.is_rev != other.is_rev :
                    if self.model == other.model :
                        print('IMPORTANT WARNING : REV DIFFERENT IN SAME MECH !')
                        self.print_reaction()
                        other.print_reaction()
                        return(False)
                    else :
                        print('REV DIFFERENT IN DIFFERENT MECH ! (so that you know, but it will be considered as the same reactions)')
                        self.print_reaction()
                        other.print_reaction()
                        return(True)
                    
                else :
                    if self.is_rev == False:
                        if tmp_set1 == tmp_set2[::-1] :
                            print(f'THAT IS THE EXCEPTION')
                            return(False)
                        elif tmp_set1 == tmp_set2 :
                            return(True)
                        else :
                            print('WTF ????????')
                    else :        
                        return(True)
            
            else :
                if tmp_set1 == tmp_set2 or tmp_set1 == tmp_set2[::-1] :
                    print('MAYBE FATALA ERROR, LOOK : ')
                    self.print_reaction()
                    other.print_reaction()
                return(False)
        else :
            return(False)
        
        
    def print_reaction(self, path='') :
        if self.is_rev == True :
            equal= '='
        elif self.is_rev == False :
            equal= '=>'
        else :
            print(f'ERROR, NO REV {self.__dict__}')

        if self.third_body != 0 :
            tb=f"+{self.third_body}"
            if self.supinfo != [] :
                low_check=["low/" for x in self.supinfo if 'low/' in x or 'low /' in x]
                tcheb_check=["tcheb/" for x in self.supinfo if 'tcheb/' in x or 'tcheb /' in x]
                if len(low_check) >= 1 or len(tcheb_check) >= 1 :
                    tb=f"(+{self.third_body})"
        else :
            tb=''

        reactant=[]
        product=[]
        # print(reaction.speed[1])
        for s, c in zip(self.species[0], self.coef[0]) :
            if c == 1 :
                c=''
            elif int(c) == c :
                c=int(c)
            reactant.append(str(c) + str (s))

        for s, c in zip(self.species[1], self.coef[1]) :
            if c == 1 :
                c=''
            elif int(c) == c :
                c=int(c)
            product.append(str(c) + str (s))

        final_reaction=f'{"+".join(reactant)}{tb if tb!= "" else ""}{equal}{"+".join(product)}{tb if tb!= "" else ""}'
        # to_write=f'{final_reaction:48} {reaction.speed[0]}{reaction.speed[1]:<+7.3f}{reaction.speed[2]:+11.3e}'
        to_write='{:<48}{:>+12.3e}{:>+9.3f}{:>+11.1f}'.format\
        (final_reaction, self.speed[0], self.speed[1], self.speed[2])
        if len(to_write)>80 :
            print(f'Reaction from base model too big ! :\n{to_write}')
            to_write=re.sub("\s+", " ", to_write)
            # print(to_write)
        print((to_write).upper())
            
        # print(to_write,len(to_write))
        for sup in self.supinfo :
            print(('\t'+sup).upper())
            
        if path != '' :
            with open(path, 'a') as extract :
                extract.write((to_write+f'\t{self.comment}'+'\n').upper())
                for sup in self.supinfo :
                    extract.write(('\t'+sup+'\n').upper())
            
class Transport :
    def __init__(self, species, geo, lj_p, lj_cd, dp, p, rr) :
        self.species = species
        self.geo = geo
        self.lj_p = lj_p
        self.lj_cd = lj_cd
        self.dp = dp
        self.p = p
        self.rr = rr
        
