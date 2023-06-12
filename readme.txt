1. FBAwEB
【Step 1; Python】Get Ensemble biomass 

This module generates an ensemble of biomass equations, with the number specified by the user.


    1. Input
        organism=''         # ecoli / yeast / cho
                            # This 'organism' option is introduced to accommodate the variations in the metabolic networks of different models, particularly 
							regarding lipid synthesis.
                            For example, in the E. coli model, lipid synthesis is limited to fatty acids.
                            On the other hand, the CHO cells model involves a broader range of lipid types, such as 'chsterol' and 'clpn_cho'.
                            These lipid metabolites are produced through a synthesis network that incorporates metabolites denoted as "Rtotal" and "Rtotalcoa",
                            which represent the relative composition of fatty acids.
                            Each model has its unique metabolic structure and intermediate metabolites involved in biomass synthesis.
                            Hence, we offer the organism option, enabling users to select the most relevant model according to their specific requirements
                            and obtain the biomass equation accordingly.

        file2read=''        # Directory of biomass excel file (Ex: 'Ecoli test1.xlsx')
        sampling_n=         # The number of biomass equations to generate (Ex: 5000)
        macro_cols="A:W"    # This is the column range to use in the file2read (Ex: "Ecoli test1.xlsx').
                            Use the range "A:W" unless you comprehensively understand the codes involved
                            and intend to modify both the biomass excel file and the associated code accordingly.

    2. Output
        Output Directory = "{0}_Ensemble biomass_macro&FA +-2STDEV ... .xlsx' "  #{0} is organism name

        In the Output result,
            "Random coefficient" Sheet:
                                        "Data Average" and "Data Stdev" refer to the average and standard deviation, respectively, of the biomass composition data 					   obtained from the biomass file (e.g., "Ecoli test1.xlsx").
                                        "Norm_Data Average" represents the normalized biomass composition data used as a reference for the ensemble of biomass equations.
                                        "Random Average" and "Random Stdec" represent the average and standard deviation of the biomass composition data for the resulting 5000 ensemble biomass equations.
                                        Columns from "0" to "4999" represent the different biomass compositions used for the 5000 biomass equations. These compositions were randomized within specific ranges based on the coefficient of variations (CVs).

            "ref" Sheet:
                                        The reference biomass equations formulated based on the biomass composition of "Norm_Data Average" column in the "Random coefiiceint" sheet.

            "PROTsyn", "DNAsyn", ... "Biomass" Sheets:
                                        "PROTsyn" sheet provides 5000 distinct protein synthesis equations, with each row index corresponding to the column header 					   in the "Random coefficient" sheet.
                                        Similarly, the same concept applies to other columns such as "DNAsyn" and "RNAsyn."
                                        These columns also contain 5000 different equations for DNA synthesis and RNA synthesis, respectively.
                                        In both cases, the row indices correspond to the column headers in the "Random coefficient" sheet.
                                        The stoichiometric coefficients in all sheets from "PROTsyn" to "Biomass" have a unit of "mmol/gDCW/h"

    3. Usage Example
        testfile1='Ecoli test1.xlsx'
        # Generate Ensemble biomass equations
        a=randomBiomass(organism='ecoli',file2read=testfile1,sampling_n=5000,macro_cols="A:W")
        # Save the result
        b=a.exportBiomassEqns()


【Step 2; Python/Matlab】 run (p)FBAwEB

-Python-

This function generates the results of the Flux Balance Analysis (FBA), where each column corresponds to one of the 5000 biomass equations.


    1. Input
        organism=''     # ecoli / yeast / cho
                        # This 'organism' option is introduced to accommodate the variations in the metabolic networks of different models, particularly regarding lipid synthesis.
                        For example, in the E. coli model, lipid synthesis is limited to fatty acids.
                        On the other hand, the CHO cells model involves a broader range of lipid types, such as 'chsterol' and 'clpn_cho'.
                        These lipid metabolites are produced through a synthesis network that incorporates metabolites denoted as "Rtotal" and "Rtotalcoa",
                        which represent the relative composition of fatty acids.
                        Each model has its unique metabolic structure and intermediate metabolites involved in biomass synthesis.
                        Hence, we offer the organism option, enabling users to select the most relevant model according to their specific requirements
                        and obtain the biomass equation accordingly.
                        
        modelspath=''   # GEMs_FBAwEB.mat directory
        model_variable=''
                        # This option is for model conditions
                          ecoli_model_FBAwEB                # E.coli aerobic condition
                          ecoli_model_anaerobic_FBAwEB      # E.coli anaerobic condition
                          yeast_model_800_glucose_FBAwEB    # S.cerevisiae glucose condition
                          yeast_model_800_ethanol_FBAwEB    # S.cerevisiae ethanol condition
                          yeast_model_800_xylose_FBAwEB     # S.cerevisiae xylose condition
                          cho_model_control_FBAwEB          # CHO cells control condition
                          cho_model_he_FBAwEB               # CHO cells high expression condition
                          cho_model_le_FBAwEB               # CHO cells low expression condition
        ensemblepath='' # '$$$_Ensemble biomass_macro&FA +-2STDEV ... .xlsx' directory (Output file of formulate_biomass.py)
        solver='gurobi'

    2. Output
        Output Directory = "FBAwEB result\{0}_pfBAwEB result_{1}.xlsx"  #{0} is organism name and {1} is date.


    3. Usage Example
        FBAwEB(ecoli,'GEMs_FBAwEB.mat','ecoli_model_FBAwEB','ECOLI_Ensemble biomass_macro&FA +-2STDEV_Feb04 15;24','gurobi')

   
-MATLAB-
(Example) 
---------------------------------------
initCobraToolbox(false)
load('Model/GEMs_FBAwEB_c13.mat') $ or load('Model/GEMs_FBAwEB.mat')
biomass_excel_name='Script/main/RandomBiomassExcel/ECOLI_Ensemble biomass_macro&FA +-2STDEV_date.xlsx'
[model_modified]=FBAwEB("ECOLI",iML1515_aerobic,biomass_excel_name,5000,0,'output file name,xlsx',false)

$function [model_modi1] = FBAwEB(species,model_name,biomass_excel_name,number_of_biomass,starting_point,filename2save,FVAoption);
---------------------------------------
(Note)

# Input :
	species = applicable for "ECOLI" or "YEAST" or "CHO"
	model_name = GEM model file (.mat format) to read 
	biomass_excel_name = sampled biomass equations excel file to import
	number_of_biomass = the number of sampled biomass equations that you want to implement FBAwEB 
	starting_point = if implementation failed during MATLAB run, you can restart from the failure point
	filename2save = name of exporting excel file ("FBAwEB_datetime+filenae2save(.xlsx)") 
	[ex. "FBAwEB_2020-01-01-00:00filename2save.xlsx"]
# Directory :
	Script/main/FBAwEB.mat
# Output :
	"FBAwEB_datetime+filenae2save.xlsx"


----------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------

2. Sensitivity analysis
【Step 1; Python】 Get biomass equations

(Example)
---------------------------------------
import Formulate_Eco_Macro_minmax as ema
import Formulate_Eco_mono_minmax as emo

testfile='Script/Example/Ecoli test1.xlsx'
a1=ema.makeBiomass(file2read=testfile)
b1=a1.Formulate_min_max()

a2=emo.makeBiomass(file2read=testfile)
b2=a2.Formulate_min_max()
-----------------------------------------
(Note)

# Input :
	file2read=biomass composition data file (excel) 
	[ex. 'test/Ecoli test1.xlsx']

# Output :
	'Script/Analysis/Sensitivity/RandomBiomassExcel/Ecoli_macro_+-25% biomass_date.xlsx'
	'Script/Analysis/Sensitivity/RandomBiomassExcel/Ecoli_mono_+-25% biomass_date.xlsx'
	
	
【Step 2; MATLAB/Python】 Sensitivity analysis monomer/macromolecular

(Example)
---------------------------------------
initCobraToolbox(false)
load('Model/GEMs_FBAwEB_c13.mat') $ or load('Model/GEMs_FBAwEB.mat')

[model_modified, flux_table, summary_table] = monoSensitivity("CHO",iCHO2291_control,'Script\Analysis\Sensitivitiy\RandomBiomassExcel\Cho_mono_+-25% biomass_date.xlsx','cho.xlsx',false)
[model_modified1, flux_table1, summary_table1] = macroSensitivity("CHO",iCHO2291_control,'Script\Analysis\Sensitivitiy\RandomBiomassExcel\Cho_macro_+-25% biomass_date.xlsx','cho.xlsx',false)
-----------------------------------------
(Note)

# Input :
	species = applicable for "ECOLI" or "YEAST" or "CHO"
	model_name = GEM model file (.mat format) to read 
	biomass_excel_name = macromolecular or monomeric varied biomass equations excel file to import
	filename2save = name of exporting excel file ("MACRO_datetime+filenae2save(.xlsx)" or "MONO_datetime+filenae2save(.xlsx)") 
	[ex. "MACRO/MONO_2020-01-01-00:00filename2save.xlsx"]

# Output :
	'MACRO 2022-01-01_00_00_00cho.xlsx'
	'MONO 2022-01-01_00_00_00cho.xlsx'

