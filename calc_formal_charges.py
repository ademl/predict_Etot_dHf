# Calculates formal charges of species for given chemical formula
# input=vals dictionary with key=elements, value=[stoich, space_holder]
# output=vals dictionary with key=elements, value=[stoich, formal_chg]
# Limited ability to calculate charges for compounds with multiple TMs and/or
#   col 14/15 elements which can be cations or anions. If charges undetermined,
#   returns dict with value=[stoich, nan]
# 2015-10-05 A.Deml

#######################
# assign charges for elements with a single charge state
# assign NaN for elements with multiple possible charge states
fml_charges = {
    'Li': 1,
    'Be': 2,
    'N': -3,
    'O': -2,
    'F': -1,
    'Na': 1,
    'Mg': 2,
    'Al': 3,
    'Si': float('NaN'),
    'P': -3,
    'S': -2,
    'Cl': -1,
    'K': 1,
    'Ca': 2,
    'Sc': 3,
    'Ti': float('NaN'),
    'V': float('NaN'),
    'Cr': float('NaN'),
    'Mn': float('NaN'),
    'Fe': float('NaN'),
    'Co': float('NaN'),
    'Ni': float('NaN'),
    'Cu': float('NaN'),
    'Zn': 2,
    'Ga': 3,
    'Ge': float('NaN'),
    'As': -3,
    'Se': -2,
    'Br': -1,
    'Rb': 1,
    'Sr': 2,
    'Y': 3,
    'Zr': 4,
    'Nb': float('NaN'),
    'Mo': float('NaN'),
    'Rh': float('NaN'),
    'Pd': float('NaN'),
    'Ag': 1,
    'Cd': 2,
    'In': 3,
    'Sn': float('NaN'),
    'Sb': float('NaN'),
    'Te': -2,
    'I': -1,
    'Cs': 1,
    'Ba': 2,
    'Hf': 4,
    'Ta': float('NaN'),
    'W' : float('NaN'),
    'Ir': float('NaN'),
    'Pt': float('NaN'),
    'Au': 1,
    'Hg': 2,
    'Pb': float('NaN'),
    'Bi': float('NaN'),
    'La': 3}

# Col 14-15 elements with multiple charge options including cation or anion
cation_anion=['Si', 'Ge', 'Sn', 'Sb', 'Bi']

# List of charge options for elements with multiple charge options (TMs, Col 14-15)
# For TMs allow all integers in reasonable range even if uncommon
# Do not change order of col 14-15 charges!
fml_multi_charges = {'Si': [-4,2,4],
                     'Ge': [-4,2,4],
                     'Sn': [-4,2,4],
                     'Sb': [-3,3,5],
                     'Bi': [-3,3,5],

                     'Ti': [2,3,4],
                     'V' : [2,3,4,5],
                     'Cr': [2,3,4,5,6],
                     'Mn': [2,3,4,5,6,7],
                     'Fe': [2,3,4],
                     'Co': [2,3,4],
                     'Ni': [2,3,4],
                     'Cu': [1,2],
                     'Nb': [2,3,4,5],
                     'Mo': [2,3,4,5,6],
                     'Rh': [2,3,4],
                     'Pd': [2,3,4],
                     'Ta': [2,3,4,5],
                     'W' : [2,3,4,5,6],
                     'Ir': [2,3,4,5,6],
                     'Pt': [2,3,4],
                     'Pb': [2,3,4]}

#########################################
# Calculates formal charges with simple algorithm
# Electronegativities used to assign charges for cmpds w/2 cations with multiple charge options
# input=vals dictionary with key=elements, value=[stoich, space_holder]
# output=vals dictionary with key=elements, value=[stoich, formal_chg]

def calc_chg(vals=None):
    from elemental_properties import pauling
    import numpy as np

    # create list of elements
    elements=[]
    for elem in vals:
        elements.append(elem)

    # assign charges for elements with single charge state
    for elem in elements:
        if elem in fml_charges.keys():
            vals[elem][1]=fml_charges[elem]
        else:
            vals[elem][1]=float('NaN')

    # count num of elements with multi charges
    num_multichg=np.count_nonzero([np.isnan(vals[elem][1]) for elem in elements])

    # if more than 2 unknown, assign all to NaN -- cannot determine charges
    if num_multichg>2:
        for elem in elements:
            vals[elem][1]=float('NaN')

    # calc formal charge if two elements unknown
    if num_multichg==2:
        # if one element from Col 14-15 (cation_anion)
        if np.count_nonzero([elem in cation_anion for elem in elements])==1:
            # if an anion has already been assigned, assume col 14-15 elem is cation
            if any([vals[elem][1]<0 for elem in elements]):
                # get elements with unknown charges
                for elem in elements:
                    if elem in cation_anion:
                        ctn_an=elem
                    if elem not in cation_anion and elem in fml_multi_charges.keys():
                        TM=elem
                # compare EN of TM and cation_anion to determine charge of cation_anion
                EN={}
                EN[ctn_an]=float(pauling[ctn_an])
                EN[TM]=float(pauling[TM])
                # if cation_anion EN is smaller, will give more e- for more pos chg
                if EN[ctn_an]<=EN[TM]:
                    vals[ctn_an][1]=max(fml_multi_charges[ctn_an])
                else:
                    vals[ctn_an][1]=fml_multi_charges[ctn_an][1]
                # calc chg of TM below since only 1 unknown remains

            else: # no anion assigned
                for elem in elements:
                    if elem in cation_anion:
                        ctn_an=elem
                vals[ctn_an][1]=fml_multi_charges[ctn_an][0]
                # calc chg of other elem below since only 1 unknown remains

        # if two elements from Col 14-15
        elif np.count_nonzero([elem in cation_anion for elem in elements])==2:
            # if another anion has already been assigned, assume col 14-15 elem are cations
            if any([vals[elem][1]<0 for elem in elements]):
                # get elements with unknown charges
                unknown_elem=[]
                for elem in elements:
                    if elem in cation_anion:
                        unknown_elem.append(elem)
                # arbitarily assign one to act as 'TM'
                ctn_an=unknown_elem[1]
                TM=unknown_elem[0]
                # compare EN to determine charge of cation_anion
                EN={}
                EN[ctn_an]=float(pauling[ctn_an])
                EN[TM]=float(pauling[TM])
                # if cation_anion EN is smaller, will give more e- for more pos chg
                if EN[ctn_an]<=EN[TM]:
                    vals[ctn_an][1]=max(fml_multi_charges[ctn_an])
                else:
                    vals[ctn_an][1]=fml_multi_charges[ctn_an][1]
                # calc chg of other col 14-15 elem below since only 1 unknown remains

            else: # no anion assigned -- cannot determine charges
                for elem in elements:
                    vals[elem][1]=float('NaN')

        # if two TMs, then >2 cation charge options for each -- cannot determine charges
        elif np.count_nonzero([elem in cation_anion for elem in elements])==0:
            for elem in elements:
                vals[elem][1]=float('NaN')

    # count remaining num of elements with multi charges
    num_multichg=np.count_nonzero([np.isnan(vals[elem][1]) for elem in elements])

    # calc formal charge if one element unknown
    if num_multichg==1:
        # sum all other charges * stoichiometry
        sumchg=0.
        for elem in elements:
            if np.isnan(vals[elem][1])==False:
                sumchg=sumchg+vals[elem][0]*vals[elem][1]
            else:
                unknown_elem=elem
        # calc charge of unknown elem that gives zero net charge
        calc_chg=-sumchg/vals[unknown_elem][0]
        
        # assign charge if calculated charge matches options given in fml_multi_charges
        for chg in fml_multi_charges[unknown_elem]:
            if calc_chg==chg:
                vals[unknown_elem][1]=calc_chg
        # if no match found, set all to NaN
        if np.count_nonzero([np.isnan(vals[elem][1]) for elem in elements])==1:
            for elem in elements:
                vals[elem][1]=float('NaN')

    # count remaining num of elements with multi charges
    num_multichg=np.count_nonzero([np.isnan(vals[elem][1]) for elem in elements])

    # if zero unknown charges, check that charges sum to 0
    if num_multichg==0:
        sumchg=0.
        for elem in elements:
            sumchg=sumchg+vals[elem][0]*vals[elem][1]
            
        # if charges do not sum to zero, assign all as NaN
        if sumchg != 0:
            for elem in elements:
                vals[elem][1]=float('NaN')

    # change all charges to float for later calculations
    for elem in elements:
        vals[elem][1]=float(vals[elem][1])

    return (vals)
