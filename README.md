# predict_Etot_dHf
Predicts the total energies and enthalpies of formation of metal-nonmetal compounds from their composition. 

Predicting Etot/dHf from Model 1,
Submitted to Phys Rev B (2015)

This python code contains a linear model to predict the Etot, and 
consequently dHf, of metal-nonmetal compounds from properties 
describing the compound compositions and their elemental constituents.

Additional information regarding development of the model and
evaluation of its predictive accuracy is available in the published 
article above. 


Required files:

1. chemical_formulas.txt (user input)
   User provided list of chemical formulas of the form AnBm... 
   	where A,B are elements and n,m are stoichiometric coefficients.

2. predict_Etot_dHf.py (executable)
   Provides a definition for predicting Etot and dHf: def predict_Etot_dHf(chem_form=None)
   Requires only an input chemical formulas.
   
   The definition returns Etot in eV/atom, dHf in eV/atom, and
   any error messages.
   
   This code calls 3 and 4 below.
   
   This code reads the user specified list of chemical compositions from chemical_formulas.txt
   	and writes the predicted Etot and dHf to a file: Etot_dHf.txt

3. calc_formal_charges.py
   Provides a definition to calculate formal charges on the ions: def calc_chg(vals=None)
   	Requires an input dictionary containing the elements and their stoichiometries. 
   	Format specified in predict_Etot_dHf.py
   	
   	Returns a dictionary with the formal elements and their formal charges. 
   	If charges could not be determined, all values are returned as nan.

4. elemental_properties.py
   Provides data tables of the elemental properties used in predicting Etot and dHf.

Procedure:
1. Specify list of chemical formulas in chemical_formulas.txt.
2. Run python predict_Etot_dHf.py to generate Etot_dHf.txt with predicted
   Etot and dHf.
