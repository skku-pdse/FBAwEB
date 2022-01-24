1. FBAwEB
【Step 1; Python】Get Ensemble biomass 

(Example)
---------------------------------------
import formulate_biomass as fb

testfile='Script/Example/Ecoli test1.xlsx'
a=fb.randomBiomass(organism='ecoli',file2read=testfile,sampling_n=5000,macro_cols='A:W')
b=a.exportBiomassEqns()
---------------------------------------
(Note)
# Output :
	Script/main/RandomBiomassExcel/ECOLI_Ensemble biomass_macro&FA +-2STDEV_date.xlsx
# Option :
	organism = 'ecoli' or 'yeast' or 'cho'
	file2read = biomass composition data file (excel)
	sampling_n = number of biomass equations to make
	macro_cols = a column range to read in a sheet 'Overall' 
# Directory :
	Script/main/formulate_biomass.py

【Step 2; MATLAB】 run (p)FBAwEB

(Example)
---------------------------------------
initCobraToolbox(false)
load('Model/GEMs_FBAwEB_c13.mat') $ or load('Model/GEMs_FBAwEB.mat')
biomass_excel_name='Script/main/RandomBiomassExcel/ECOLI_Ensemble biomass_macro&FA +-2STDEV_date.xlsx'
[model_modified]=FBAwEB("ECOLI",iML1515_aerobic,biomass_excel_name,5000,0,'output file name,xlsx',false)

$function [model_modi1] = FBAwEB(species,model_name,biomass_excel_name,number_of_biomass,starting_point,filename2save,FVAoption);
---------------------------------------
(Note)
# Output :
	"FBAwEB_datetime+filenae2save.xlsx"
# Option :
	species = applicable for "ECOLI" or "YEAST" or "CHO"
	model_name = GEM model file (.mat format) to read 
	biomass_excel_name = sampled biomass equations excel file to import
	number_of_biomass = the number of sampled biomass equations that you want to implement FBAwEB 
	starting_point = if implementation failed during MATLAB run, you can restart from the failure point
	filename2save = name of exporting excel file ("FBAwEB_datetime+filenae2save(.xlsx)") 
	[ex. "FBAwEB_2020-01-01-00:00filename2save.xlsx"]
# Directory :
	Script/main/FBAwEB.mat







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
# Output :
	'Script/Analysis/Sensitivity/RandomBiomassExcel/Ecoli_macro_+-25% biomass_date.xlsx'
	'Script/Analysis/Sensitivity/RandomBiomassExcel/Ecoli_mono_+-25% biomass_date.xlsx'
# Option :
	file2read=biomass composition data file (excel) 
	[ex. 'test/Ecoli test1.xlsx']


【Step 2; MATLAB】 Sensitivity analysis monomer/macromolecular

(Example)
---------------------------------------
initCobraToolbox(false)
load('Model/GEMs_FBAwEB_c13.mat') $ or load('Model/GEMs_FBAwEB.mat')

[model_modified, flux_table, summary_table] = monoSensitivity("CHO",iCHO2291_control,'Script\Analysis\Sensitivitiy\RandomBiomassExcel\Cho_mono_+-25% biomass_date.xlsx','cho.xlsx',false)
[model_modified1, flux_table1, summary_table1] = macroSensitivity("CHO",iCHO2291_control,'Script\Analysis\Sensitivitiy\RandomBiomassExcel\Cho_macro_+-25% biomass_date.xlsx','cho.xlsx',false)
-----------------------------------------
(Note)
# Output :
	'MACRO 2022-01-01_00_00_00cho.xlsx'
	'MONO 2022-01-01_00_00_00cho.xlsx'
# Option :
	species = applicable for "ECOLI" or "YEAST" or "CHO"
	model_name = GEM model file (.mat format) to read 
	biomass_excel_name = macromolecular or monomeric varied biomass equations excel file to import
	filename2save = name of exporting excel file ("MACRO_datetime+filenae2save(.xlsx)" or "MONO_datetime+filenae2save(.xlsx)") 
	[ex. "MACRO/MONO_2020-01-01-00:00filename2save.xlsx"]


