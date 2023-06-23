import cobra
import pandas as pd
from pandas import ExcelWriter
from os.path import join
from pathlib import Path
from datetime import datetime
from tqdm import tqdm

def FBAwEB(modelspath,model_variable,ensemblepath,solver):
    """
    This function is a python version of 'Script/main/FBAwEB.m'
    https://github.com/skku-pdse/FBAwEB/blob/3e7799f19256c6a174b1a62389620b4a2d807eef/Script/main/FBAwEB.m

    This function generates the results of the Flux Balance Analysis (FBA), where each column corresponds to one of the 5000 biomass equations.


    1. Input

        modelspath=''   # GEMs_FBAwEB.mat directory
        model_variable=''
                        # This option is for model conditions
                          ecoli_model_FBAwEB                # E.coli aerobic condition
                          ecoli_model_anaerobic_FBAwEB      # E.coli anaerobic condition

        ensemblepath='' # '$$$_Ensemble biomass_macro&FA +-2STDEV ... .xlsx' directory (Output file of formulate_biomass.py)
        solver='gurobi'

    2. Output
        Output Directory = "FBAwEB result\{0}_pfBAwEB result_{1}.xlsx"  #{0} is organism name and {1} is date.


    3. Usage Example
        FBAwEB('GEMs_FBAwEB.mat','ecoli_model_FBAwEB','ECOLI_Ensemble biomass_macro&FA +-2STDEV_Feb04 15;24','gurobi')

    """



    # result file name
    saveDir = "FBAwEB result"
    Path(saveDir).mkdir(parents=True, exist_ok=True)
    file2save='{0}_pfBAwEB result_{1}.xlsx'.format(model_variable,datetime.now().strftime("%b%d %H;%M"))
    file2save = join(saveDir, file2save)

    #load matlab model
    m_model=cobra.io.load_matlab_model(modelspath,variable_name=model_variable)
    m_model.solver=solver

    #read Ensemble biomass
    Ensemblexls=pd.ExcelFile(ensemblepath)
    refdf=pd.read_excel(Ensemblexls,'ref',index_col=0, header=0)
    proteindf=pd.read_excel(Ensemblexls, 'PROTsyn',index_col=0, header=0)
    dnadf=pd.read_excel(Ensemblexls,'DNAsyn',index_col=0, header=0)
    rnadf=pd.read_excel(Ensemblexls,'RNAsyn',index_col=0, header=0)
    carbdf = pd.read_excel(Ensemblexls, 'CARBsyn', index_col=0, header=0)
    lipiddf=pd.read_excel(Ensemblexls,'LIPIDsyn',index_col=0, header=0)
    biomassdf=pd.read_excel(Ensemblexls,'biomass',index_col=0, header=0)

    # ADD metabolite

    m_prot=cobra.Metabolite('prot[c]',compartment='c')
    m_dna=cobra.Metabolite('dna[c]',compartment='c')
    m_rna = cobra.Metabolite('rna[c]', compartment='c')
    m_carb = cobra.Metabolite('carb[c]', compartment='c')
    m_lipid = cobra.Metabolite('lipid[c]', compartment='c')
    m_model.add_metabolites([m_prot,m_dna,m_rna,m_carb,m_lipid])



    # Add cobra reaction
    PROTsyn = cobra.Reaction(id="PROTsyn", name='PROTsyn')
    DNAsyn = cobra.Reaction(id="DNAsyn", name='DNAsyn')
    RNAsyn = cobra.Reaction(id="RNAsyn", name='RNAsyn')
    CARBsyn = cobra.Reaction(id="CARBsyn", name='CARBsyn')
    LIPIDsyn = cobra.Reaction(id="LIPIDsyn", name='LIPIDsyn')
    Biomass = cobra.Reaction(id="Biomass", name='Biomass')
    m_model.add_reactions([PROTsyn, DNAsyn, RNAsyn, CARBsyn,LIPIDsyn, Biomass])

    ###################
    # Reference equation pFBA

    ref_model = m_model

    ref_model.reactions.PROTsyn.build_reaction_from_string(refdf.iloc[0, 0])
    ref_model.reactions.DNAsyn.build_reaction_from_string(refdf.iloc[1, 0])
    ref_model.reactions.RNAsyn.build_reaction_from_string(refdf.iloc[2, 0])
    ref_model.reactions.CARBsyn.build_reaction_from_string(refdf.iloc[3, 0])
    ref_model.reactions.LIPIDsyn.build_reaction_from_string(refdf.iloc[4, 0])
    ref_model.reactions.Biomass.build_reaction_from_string(refdf.iloc[5, 0])

    # change Obective
    ref_model.objective = "Biomass"

    # pFBA
    ref_pfba_solution = cobra.flux_analysis.pfba(ref_model)

    # save solution
    ref_sol_df = pd.DataFrame()
    ref_sol_df[0] = ref_pfba_solution.fluxes


    ###################
    # FBAwEB
    eqnDict={}
    for i in tqdm(range(len(biomassdf))):

        i_model=m_model

        # Add string reaction eqn
        i_model.reactions.PROTsyn.build_reaction_from_string(proteindf.iloc[i, 0])
        i_model.reactions.DNAsyn.build_reaction_from_string(dnadf.iloc[i,0])
        i_model.reactions.RNAsyn.build_reaction_from_string(rnadf.iloc[i, 0])
        i_model.reactions.CARBsyn.build_reaction_from_string(carbdf.iloc[i, 0])
        i_model.reactions.LIPIDsyn.build_reaction_from_string(lipiddf.iloc[i, 0])
        i_model.reactions.Biomass.build_reaction_from_string(biomassdf.iloc[i, 0])


        # change Obective
        i_model.objective="Biomass"

        # pFBA
        pfba_solution = cobra.flux_analysis.pfba(i_model)

        # save solution
        if i==0:
            sol_df=pd.DataFrame()
            for j in range(len(i_model.reactions)):
                rxn_j = m_model.reactions[j]
                rxn_id = rxn_j.id
                rxn_eqn = rxn_j.reaction
                rxn_sys= rxn_j.subsystem
                eqnDict[rxn_id]=[rxn_eqn,rxn_sys]
        sol_df[i] = pfba_solution.fluxes

    eqndf=pd.DataFrame.from_dict(eqnDict,orient='index',columns=['RxnFormula','Subsystem'])
    eqndf.index.name = "Rxn"

    # flux table sheet
    f_df=pd.concat([eqndf, sol_df], axis=1)
    ref_sol_df=pd.concat([eqndf,ref_sol_df],axis=1)

    # summary table sheet
    flux_median=f_df.iloc[:,3:].median(axis='columns')
    flux_std=f_df.iloc[:,3:].std(axis='columns')
    flux_max=f_df.iloc[:,3:].max(axis='columns')
    flux_min = f_df.iloc[:, 3:].min(axis='columns')

    s_df=f_df.iloc[:,:3]
    s_df['flux_median']=flux_median
    s_df['flux_std'] = flux_std
    s_df['flux_max'] = flux_max
    s_df['flux_min'] = flux_min
    print (s_df)

    with ExcelWriter(file2save) as writer :
        ref_sol_df.to_excel(writer,sheet_name='ref flux')
        f_df.to_excel(writer,sheet_name='flux table')
        s_df.to_excel(writer,sheet_name='summary table')
        writer.save()



# if __name__ == '__main__':
#     modelspath=""
#     model_variable=""
#     ensemblepath=""
#     solver='gurobi'
#     FBAwEB(modelspath,model_variable,ensemblepath,solver)

