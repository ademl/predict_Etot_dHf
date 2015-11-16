# python definition to predict Etot and dHf from input composition
#   input: composition, e.g. 'Al2O3'
#   returns: Etot(eV/atom), dHf(eV/atom), error messages
#   2015-10-11 A.Deml

######################################
from elemental_properties import *
from calc_formal_charges import *
import numpy as np
import re
from copy import deepcopy

def predict_Etot_dHf(chem_form=None):
    # define a class to calc statists for elemental properties
    class elem_prop_stats:
        def __init__(self):
            stoich=np.array(vals.values())[:,0]
            prop=np.array(vals.values())[:,1]
            stoich_prop=np.array(vals.values())[:,0]* \
                        np.array(vals.values())[:,1]
            self.max=max(prop)
            self.min=min(prop)
            self.range=self.max-self.min
            self.std=np.std(prop) #population stdev
            self.mean=np.sum(stoich_prop) / np.sum(stoich)

    # start calculating...
    # get elements and stoich
    err_msg=None
    elements=filter(None, re.split('[0-9]', chem_form))
    stoich=map(int, filter(None, re.split('\D', chem_form)))

    # exit if compound is an element or if not all elements in FERE list
    if len(elements)==1:
        return('','','Error: %s is an element.' %chem_form)
    elif not all([elem in mus_fere for elem in elements]):
        return('','','Error: %s contains an element without a FERE value.' %chem_form)
    else: 
        # make dict with key=elements, value=[stoich, elemental property]
        # elemental property initially set to 'NaN'
        vals={}
        for elem in elements:
            vals[elem]=[float(stoich[elements.index(elem)]), float('NaN')]

        # calculate formal charges
        # send vals dict to function that returns new vals dict with calculated charges
        vals=calc_chg(vals)

        # exit if charges not known
        if any(np.isnan(vals[elem][1])==True for elem in vals):
            return('','','Error: Formal charges could not be determined for %s.' %chem_form)
        else:
            # save vals dict with formal charges for later
            charges=deepcopy(vals)
            #print 'charges=', charges

            # create separate lists for cations and anions
            cation_elements=[]
            anion_elements=[]
            cation_vals={}
            anion_vals={}
            for elem in elements:
                # identify cations
                if charges[elem][1]>0:
                    cation_elements.append(elem)
                    cation_vals[elem]=[vals[elem][0],charges[elem][1]]
                # identify anion(s)
                else:
                    anion_elements.append(elem)
                    anion_vals[elem]=[vals[elem][0],charges[elem][1]]

            ##############
            # calculate elemental properties and related stats
            # atoms per formula unit
            atoms_per_fu=sum(v[0] for v in vals.values())

            # fraction of TMs
            num_magn=0
            for elem in elements:
                if elem in magn_elem:
                    num_magn=num_magn+vals[elem][0]
            fraction_TMs=float(num_magn)/atoms_per_fu
            #print('fraction TM='),fraction_TMs

            # atomic number
            for elem in elements:
                vals[elem][1]=atom_num[elem]
            atomic_number=elem_prop_stats()

            #print 'atomic_number.min=', atomic_number.min
            #print 'atomic_number.max=', atomic_number.max
            #print 'atomic_number.range=', atomic_number.range
            #print 'atomic_number.std=', atomic_number.std
            #print 'atomic_number.mean=', atomic_number.mean

            # atomic mass
            for elem in elements:
                vals[elem][1]=atom_mass[elem]
            atomic_mass=elem_prop_stats()
            #print 'atomic_mass.mean=', atomic_mass.mean

            # row number
            for elem in elements:
                vals[elem][1]=row_num[elem]
            row_number=elem_prop_stats()
            #print 'row_number.mean=',row_number.mean

            # column number
            for elem in elements:
                vals[elem][1]=col_num[elem]
            columns=deepcopy(vals)
            column_number=elem_prop_stats()
            #print 'column_number.mean=', column_number.mean

            # ave number of s,p,d valence electrons
            sum_s=0; sum_p=0; sum_d=0
            for elem in elements: 
                col=columns[elem][1]
                sum_s=sum_s+valence_e[col]['s']*vals[elem][0]
                sum_p=sum_p+valence_e[col]['p']*vals[elem][0]
                sum_d=sum_d+valence_e[col]['d']*vals[elem][0]
            num_s_valence_e=float(sum_s/atoms_per_fu)
            num_p_valence_e=float(sum_p/atoms_per_fu)
            num_d_valence_e=float(sum_d/atoms_per_fu)
            #print 'num_s_valence_e=', num_s_valence_e
            #print 'num_p_valence_e=', num_p_valence_e
            #print 'num_d_valence_e=', num_d_valence_e

            # fraction of s,p,d valence electrons
            tot_val_e=sum_s+sum_p+sum_d
            fraction_s_valence_e=float(sum_s/tot_val_e)
            fraction_p_valence_e=float(sum_p/tot_val_e)
            fraction_d_valence_e=float(sum_d/tot_val_e)
            #print 'fraction_s_valence_e=', fraction_s_valence_e
            #print 'fraction_p_valence_e=', fraction_p_valence_e
            #print 'fraction_d_valence_e=', fraction_d_valence_e

            # atomic radius
            for elem in elements:
                vals[elem][1]=atom_radius[elem]
            atomic_radius=elem_prop_stats()
            #print 'atomic_radius.mean=', atomic_radius.mean

            # molar volume
            for elem in elements:
                vals[elem][1]=molar_vol[elem]
            molar_volume=elem_prop_stats()
            #print 'molar_volume.mean=', molar_volume.mean

            # latent heat of fusion
            for elem in elements:
                vals[elem][1]=heat_fusion[elem]
            latent_heat_fusion=elem_prop_stats()
            #print 'latent_heat_fusion.mean=', latent_heat_fusion.mean

            # melting point
            for elem in elements:
                vals[elem][1]=mp[elem]
            melting_point=elem_prop_stats()
            #print 'melting_point.mean=', melting_point.mean

            # boiling point
            for elem in elements:
                vals[elem][1]=bp[elem]
            boiling_point=elem_prop_stats()
            #print 'boiling_point.mean=', boiling_point.mean

            # heat capacity
            for elem in elements:
                vals[elem][1]=heat_cap[elem]
            heat_capacity=elem_prop_stats()
            #print 'heat_capacity.mean=', heat_capacity.mean

            # 1st ionization energy (includes cations and anions)
            for elem in elements:
                vals[elem][1]=ionization_en[elem][0]/1000.
            first_ionization_energy=elem_prop_stats()
            #print 'first_ionization_energy.mean=', first_ionization_energy.mean

            # cumulative ionization energy of cations
            for elem in cation_elements:
                # check if data table has all ionization energies up to formal charge
                if len(ionization_en[elem])>=charges[elem][1]:
                    sum_ioniz=sum(ionization_en[elem][i] \
                              for i in range(0,int(charges[elem][1])))
                    cation_vals[elem][1]=sum_ioniz
                else:
                    cation_vals[elem][1]=float('NaN')
            # save copy of original vals to get stats for cation_vals (clumsy!) 
            vals_orig=deepcopy(vals)
            vals=cation_vals
            cumulative_ionization_energy=elem_prop_stats()
            # change weighted mean to per atom, not per cation
            cumulative_ionization_energy.mean=cumulative_ionization_energy.mean* \
                sum(cation_vals[elem][0] for elem in cation_vals)/atoms_per_fu
            vals=vals_orig
            #print 'cumulative_ionization_energy.mean=', cumulative_ionization_energy.mean

            # average electron affinity of anions
            # multiply by formal charge to flip sign and scale for 2-,3-,...
            # no stats since most have only 1 anion
            sum_elect_affin=0
            for elem in anion_elements:
                sum_elect_affin=sum_elect_affin+electron_affin[elem]*charges[elem][1]*vals[elem][0]
            electron_affinity=sum_elect_affin/atoms_per_fu
            #print 'electron_affinity=', electron_affinity

            # Pauling electronegativity
            for elem in elements:
                vals[elem][1]=pauling[elem]
            Pauling_electronegativity=elem_prop_stats()
            #print 'Pauling_electronegativity.mean=', Pauling_electronegativity.mean

            # Pauling electronegativity differences between cations and anions
            # assume each cation has 1 anion neighbor 
            #  (or fractional sums of 2+ anions to make 1 anion neighbor)
            for cation in cation_elements:
                val=0
                cation_EN=pauling[cation]
                for anion in anion_elements:
                    fraction_anion=anion_vals[anion][0]/sum([anion_vals[anion][0] for anion in anion_elements])
                    # cummulative sum of EN diff = ave weighted by fraction of each anion
                    anion_EN=pauling[anion]
                    val=val+abs(cation_EN-anion_EN)*fraction_anion
                cation_vals[cation][1]=val
            # again send cation_vals for stats
            vals=cation_vals
            Pauling_electronegativity_diff=elem_prop_stats()
            # change weighted mean to per atom, not per cation
            Pauling_electronegativity_diff.mean=Pauling_electronegativity_diff.mean* \
                sum(cation_vals[elem][0] for elem in cation_vals)/atoms_per_fu
            vals=vals_orig
            #print 'Pauling_electronegativity_diff.mean=', Pauling_electronegativity_diff.mean

            # formal charges
            vals=deepcopy(charges)
            formal_charges=elem_prop_stats()
            #print 'formal_charges.mean=', formal_charges.mean

            # crystal field splitting
            for elem in elements:
                if elem in xtal_field_split:
                    if charges[elem][1] in xtal_field_split[elem]:
                        vals[elem][1]=xtal_field_split[elem][charges[elem][1]]
                    else: vals[elem][1]=0
                else: vals[elem][1]=0
            crystal_field_splitting=elem_prop_stats()
            #print 'crystal_field_splitting.mean=', crystal_field_splitting.mean

            # magnetic moment
            for elem in elements:
                if elem in magn_moment:
                    if charges[elem][1] in magn_moment[elem]:
                        vals[elem][1]=magn_moment[elem][charges[elem][1]]
                    else: vals[elem][1]=0
                else: vals[elem][1]=0
            magnetic_moment=elem_prop_stats()
            #print 'magnetic_moment.mean=', magnetic_moment.mean

            # spin orbit coupling
            for elem in elements:
                if elem in so_coupling:
                    if charges[elem][1] in so_coupling[elem]:
                        vals[elem][1]=so_coupling[elem][charges[elem][1]]
                    else: vals[elem][1]=0
                else: vals[elem][1]=0
            spin_orbit_coupling=elem_prop_stats()
            #print 'spin_orbit_coupling.mean=', spin_orbit_coupling.mean

            # max total electron spin
            for elem in elements:
                if elem in sat_magn:
                    if charges[elem][1] in sat_magn[elem]:
                        vals[elem][1]=sat_magn[elem][charges[elem][1]]
                    else: vals[elem][1]=0
            max_up_spin=elem_prop_stats()
            #print 'max_up_spin.mean=', max_up_spin.mean

            # electric polarizability
            for elem in elements:
                vals[elem][1]=electric_pol[elem]
            electric_polarizability=elem_prop_stats()
            #print 'electric_polarizability.mean=', electric_polarizability.mean

            # GGA+U total energies of elements
            for elem in elements:
                vals[elem][1]=GGAU_Etot[elem]
            GGAU_elemental_energy=elem_prop_stats()
            #print 'GGAU_elemental_energy.mean=', GGAU_elemental_energy.mean

            # FERE values, fitted elemental reference phase energies
            for elem in elements:
                vals[elem][1]=mus_fere[elem]
            FERE_elemental_energy=elem_prop_stats()
            #print 'FERE_elemental_energy.mean=', FERE_elemental_energy.mean

            # Difference between GGA+U and FERE elemental energies
            for elem in elements:
                vals[elem][1]=mus_fere[elem]-GGAU_Etot[elem]
            FERE_correction=elem_prop_stats()
            #print 'FERE_correction.mean=', FERE_correction.mean

            ##########################################
            # Predict Etot/dHf from Model 1
            Etot=8.12606055439434+ \
                -0.00516599516895853*np.sqrt(abs(electron_affinity))+ \
                -2.16552982624354*num_s_valence_e+ \
                0.0000003472765154249*(cumulative_ionization_energy.mean)+ \
                -4.60856942619469*1/(row_number.mean+0.0000000001)+ \
                7.03070901279658*1/(electric_polarizability.mean+0.0000000001)+ \
                ((FERE_elemental_energy.mean)-(-3.49292008244051))*((np.sqrt(abs(atomic_radius.mean))-10.98011879974)*-0.173834402047588)+ \
                (1/(electric_polarizability.mean+0.0000000001)-0.125512187547373)*((np.sqrt(abs(electric_polarizability.mean))-3.13425119907364)*3.25849616319473)+ \
                -1.51163499103685*np.sqrt(abs(Pauling_electronegativity_diff.mean))+ \
                ((GGAU_elemental_energy.mean)-(-3.65976180832511))*((np.sqrt(abs(spin_orbit_coupling.mean))-1.40633111650202)*-0.0224093137364463)+ \
                -0.175503248790532*(formal_charges.max)+ \
                (fraction_p_valence_e-0.443292841388864)*((np.sqrt(abs(Pauling_electronegativity_diff.mean))-0.839558464165695)*-3.68833219690241)+ \
                ((heat_capacity.mean)-27.0029100427661)*((np.sqrt(abs(electron_affinity))-389.26142930946)*-0.000602792922611928)+ \
                (1/(fraction_s_valence_e+0.0000000001)-2.72948120067356)*((np.sqrt(abs(fraction_p_valence_e))-0.656497458174226)*1.18841896557305)+ \
                (1/(fraction_TMs+0.0000000001)-6097889539.5562)*((np.sqrt(abs(fraction_TMs))-0.188297973750786)*-1.57502208979947e-10)+ \
                -0.608230891317656*(Pauling_electronegativity_diff.min)+ \
                (1/(magnetic_moment.mean+0.0000000001)-8616973507.11494)*((np.sqrt(abs(boiling_point.mean))-35.3083340280017)*3.12113639050643e-12)+ \
                (fraction_TMs-0.0970381498118542)*((np.sqrt(abs(cumulative_ionization_energy.mean))-1247.16387554414)*-0.00121909314326555)+ \
                0.748910543304164*(GGAU_elemental_energy.mean)+ \
                -0.0631895270977302*(GGAU_elemental_energy.min)+ \
                ((GGAU_elemental_energy.mean)-(-3.65976180832511))*(((molar_volume.mean)-18.920755815203)*0.0173534584853804)+ \
                (1/(atomic_radius.mean+0.0000000001)-0.0085551764503817)*((1/(cumulative_ionization_energy.mean+0.0000000001)-0.0000008393214301751)*104301741.674977)+ \
                ((boiling_point.mean)-1279.45531437311)*((1/(first_ionization_energy.mean+0.0000000001)-0.00111357431746744)*1.86246897134558)+ \
                -0.00476168064672303*(atomic_radius.range)+ \
                0.656773058959826*(FERE_elemental_energy.mean)+ \
                ((cumulative_ionization_energy.mean)-1672776.17500611)*((1/(atomic_radius.mean+0.0000000001)-0.0085551764503817)*0.0000465316751773129)+ \
                (1/(fraction_s_valence_e+0.0000000001)-2.72948120067356)*((np.sqrt(abs(fraction_d_valence_e))-0.319036888509652)*1.48388053501385)+ \
                ((Pauling_electronegativity_diff.range)-0.643551863493488)*(((Pauling_electronegativity_diff.range)-0.643551863493488)*0.263863294072938)+ \
                (np.sqrt(abs(boiling_point.mean))-35.3083340280017)*((np.sqrt(abs(first_ionization_energy.mean))-30.3462301670948)*0.0130432483965199)+ \
                ((FERE_correction.mean)-0.16684172587427)*((1/(atomic_radius.mean+0.0000000001)-0.0085551764503817)*-331.385575641025)+ \
                ((melting_point.mean)-617.337892133765)*((1/(atomic_mass.mean+0.0000000001)-0.0234330603282443)*0.0248804586680415)+ \
                (electron_affinity-(-159236.476609832))*((1/(num_s_valence_e+0.0000000001)-0.536630873602602)*-0.000017185622956708)+ \
                (np.sqrt(abs(crystal_field_splitting.mean))-0.263067373872025)*((np.sqrt(abs(magnetic_moment.mean))-0.145892764929502)*0.0964552351900768)+ \
                (1/(melting_point.min+0.0000000001)-0.0110718424382578)*((1/(melting_point.min+0.0000000001)-0.0110718424382578)*4314.84724939908)+ \
                (1/(cumulative_ionization_energy.mean+0.0000000001)-0.0000008393214301751)*((np.sqrt(abs(num_d_valence_e))-0.773934796161201)*-136752.773664298)+ \
                ((GGAU_elemental_energy.mean)-(-3.65976180832511))*((fraction_d_valence_e-0.172015354157611)*0.282357610480358)+ \
                (1/(atomic_mass.min+0.0000000001)-0.0476143524921413)*((1/(atomic_mass.min+0.0000000001)-0.0476143524921413)*-123.676053865014)+ \
                (1/(molar_volume.min+0.0000000001)-0.098814862277504)*((1/(molar_volume.min+0.0000000001)-0.098814862277504)*41.8038054696031)+ \
                (1/(column_number.range+0.0000000001)-0.0987004027543787)*((1/(column_number.range+0.0000000001)-0.0987004027543787)*-0.339874727945705)+ \
                -8.41317801871503*1/(molar_volume.mean+0.0000000001)+ \
                0.010288012742297*(electric_polarizability.range)+ \
                (1/(Pauling_electronegativity.mean+0.0000000001)-0.462906601985631)*((np.sqrt(abs(first_ionization_energy.mean))-30.3462301670948)*0.292390863275612)+ \
                ((electric_polarizability.std)-11.3655347533718)*(((electric_polarizability.std)-11.3655347533718)*0.00210660002843539)+ \
                ((molar_volume.std)-8.80623161110011)*(((molar_volume.std)-8.80623161110011)*-0.00192287569262687)+ \
                ((melting_point.mean)-617.337892133765)*((1/(melting_point.mean+0.0000000001)-0.00191958560574764)*-0.133763013677746)+ \
                ((row_number.mean)-3.4314253306547)*((1/(boiling_point.mean+0.0000000001)-0.000870652143242029)*-256.475129692838)+ \
                ((first_ionization_energy.mean)-928.495400698205)*((np.sqrt(abs(electron_affinity))-389.26142930946)*-0.0000415791303319574)+ \
                (np.sqrt(abs(atomic_radius.mean))-10.98011879974)*((np.sqrt(abs(Pauling_electronegativity_diff.mean))-0.839558464165695)*0.477159656486591)+ \
                0.0000548417431654954*np.sqrt(abs(cumulative_ionization_energy.range))+ \
                ((FERE_correction.min)-(-0.0740497188325089))*(((FERE_correction.min)-(-0.0740497188325089))*-0.182359700916887)+ \
                (1/(heat_capacity.min+0.0000000001)-0.0423072605922765)*((1/(heat_capacity.min+0.0000000001)-0.0423072605922765)*-1275.35134527302)+ \
                0.516502813094666*np.sqrt(abs(Pauling_electronegativity_diff.min))+ \
                ((FERE_correction.mean)-0.16684172587427)*((1/(first_ionization_energy.mean+0.0000000001)-0.00111357431746744)*-1688.85230525828)+ \
                -0.136607885652447*(FERE_correction.max)+ \
                ((Pauling_electronegativity_diff.mean)-0.733951924899417)*((np.sqrt(abs(Pauling_electronegativity.mean))-1.49519614790166)*-1.83619494544683)+ \
                (fraction_s_valence_e-0.384691804511449)*((fraction_d_valence_e-0.172015354157611)*8.52478035604471)+ \
                0.000175216652418157*(melting_point.std)+ \
                (np.sqrt(abs(first_ionization_energy.mean))-30.3462301670948)*((np.sqrt(abs(electron_affinity))-389.26142930946)*0.00199965533444257)+ \
                (num_s_valence_e-1.88320741991693)*((1/(atomic_number.mean+0.0000000001)-0.050790478194881)*-91.109417316705)+ \
                ((cumulative_ionization_energy.mean)-1672776.17500611)*((np.sqrt(abs(molar_volume.mean))-4.31439645320296)*0.0000000757855817766)+ \
                ((FERE_correction.mean)-0.16684172587427)*((1/(latent_heat_fusion.mean+0.0000000001)-0.00019552666519982)*1515.82030927662)+ \
                ((Pauling_electronegativity.mean)-2.25936651545308)*(((electric_polarizability.mean)-10.5095049891352)*-0.0218619175889465)+ \
                0.0427776971965048*(molar_volume.std)+ \
                (1/(fraction_TMs+0.0000000001)-6097889539.5562)*((np.sqrt(abs(fraction_p_valence_e))-0.656497458174226)*-5.40219191610727e-11)+ \
                -0.000697095589975313*np.sqrt(abs(latent_heat_fusion.max))+ \
                (fraction_TMs-0.0970381498118542)*((1/(Pauling_electronegativity.mean+0.0000000001)-0.462906601985631)*-2.15714242697945)+ \
                ((Pauling_electronegativity.mean)-2.25936651545308)*((1/(atoms_per_fu+0.0000000001)-0.161837673109115)*0.35420427868522)+ \
                (1/(electric_polarizability.max+0.0000000001)-0.0497358596376293)*((1/(electric_polarizability.max+0.0000000001)-0.0497358596376293)*-8.61072322708728)+ \
                ((atomic_number.range)-30.7467444993264)*(((atomic_number.range)-30.7467444993264)*-0.0000451112348726259)+ \
                -0.188680461393955*np.sqrt(abs(molar_volume.std))+ \
                (1/(atomic_mass.mean+0.0000000001)-0.0234330603282443)*((1/(magnetic_moment.mean+0.0000000001)-8616973507.11494)*0.0000000005632523931)+ \
                ((Pauling_electronegativity.mean)-2.25936651545308)*((np.sqrt(abs(boiling_point.mean))-35.3083340280017)*-0.023680130402303)+ \
                (1/(electron_affinity+0.0000000001)-(-0.0000112317548585542))*((np.sqrt(abs(max_up_spin.mean))-0.110945824556354)*2051.51721939395)+ \
                ((molar_volume.mean)-18.920755815203)*((1/(cumulative_ionization_energy.mean+0.0000000001)-0.0000008393214301751)*6760.29318832862)+ \
                (num_s_valence_e-1.88320741991693)*((1/(atomic_mass.mean+0.0000000001)-0.0234330603282443)*141.830343576902)+ \
                (1/(cumulative_ionization_energy.min+0.0000000001)-0.0000012419221686125)*((1/(cumulative_ionization_energy.min+0.0000000001)-0.0000012419221686125)*-55146914540.3623)+ \
                (1/(column_number.max+0.0000000001)-0.0630675507965874)*((1/(column_number.max+0.0000000001)-0.0630675507965874)*7143.82899247505)+ \
                ((latent_heat_fusion.mean)-7247.33713798697)*((fraction_p_valence_e-0.443292841388864)*-0.0000295531194358355)+ \
                ((FERE_correction.max)-0.455707500459811)*(((FERE_correction.max)-0.455707500459811)*-0.273633753565192)+ \
                (1/(num_p_valence_e+0.0000000001)-0.517465002920072)*((np.sqrt(abs(Pauling_electronegativity.mean))-1.49519614790166)*-0.799249155571807)+ \
                ((electric_polarizability.mean)-10.5095049891352)*((np.sqrt(abs(electron_affinity))-389.26142930946)*-0.0000537922783063893)+ \
                ((atomic_radius.std)-44.8774760812231)*(((atomic_radius.std)-44.8774760812231)*-0.0000589970677042493)+ \
                ((boiling_point.mean)-1279.45531437311)*(((cumulative_ionization_energy.mean)-1672776.17500611)*-3.99823606897642e-11) 

            # add known corrections for elements not included in model build
            if 'Mo' in elements:
                Etot=Etot-0.3
            if 'Be' in elements:
                Etot=Etot+0.45

            # add warnings for specific element charge states
            # Col 14 elements that can be either cations or anions, but FERE was fit only for cations
            # Anionic Si compounds were included in fitting the Etot model (not true for Ge,Sn)
            concerns=['Si','Ge','Sn']
            if any(elem in elements for elem in concerns):
                check_chg=list(set(elements).intersection(concerns))
                for elem in check_chg:
                    if charges[elem][1]<0:
                        err_msg=('Warning: Etot/dHf may not be accurate for %s; '
                                 'FERE and the Etot model were not fit for anionic %s.'
                                 %(chem_form,elem))
            # Col 15 elements that can be either cations or anions, but FERE was fit only for anions
            # No cationic Col 15 compounds were included in fitting the Etot model
            concerns=['Sb','Bi']
            if any(elem in elements for elem in concerns):
                check_chg=list(set(elements).intersection(concerns))
                for elem in check_chg:
                    if charges[elem][1]>0:
                        err_msg=('Warning: Etot/dHf may not be accurate for %s; '
                                 'FERE and the Etot model were not fit for cationic %s.'
                                 %(chem_form,elem))

            # calculate dHf from Etot 
            mus_sum=0.
            for elem in elements:
                mus_sum=mus_sum+mus_fere[elem]*vals[elem][0]
            dHf=Etot-mus_sum/atoms_per_fu

            # return error message if Etot could not be calculated due to missing values
            # occurs mostly for charge dependent TM properties
            if np.isnan(Etot)==True:
                return('','','Error: Not all elemental properties are available for %s.' %chem_form)
            # check for error message
            if err_msg is None:
                err_msg=''

            # return Etot/dHf and any error messages
            return (Etot, dHf, err_msg)



################# RUN CODE #################

# read list of chemical formulas from text file
f=open('chemical_formulas.txt')
chemical_formulas=[]
for line in f:
    chemical_formulas.extend(line.split())
f.close()

# output results to new text file
f=open('Etot_dHf.txt','w')
f.write('Predicted Etot and dHf from Model 1 [submitted, 2015]\n'
        'chem_form\t Etot\t dHf\t Errors/Warnings\n'
        '\t eV/atom\t eV/atom\t \n'
        '\n')

for chem_form in chemical_formulas:
    print chem_form

    # get predicted Etot, dHf, error messages
    predicted_vals=predict_Etot_dHf(chem_form)
    Etot=predicted_vals[0]
    if isinstance(Etot,float): Etot=('%.4f' %Etot)
    dHf=predicted_vals[1]
    if isinstance(dHf,float): dHf=('%.4f' %dHf)
    err_msg=predicted_vals[2]

    f.write('%s\t %s\t %s\t %s\n' %(chem_form, Etot, dHf, err_msg))
    f.flush()
f.close()
