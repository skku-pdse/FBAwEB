import cobra
import pandas as pd
from pandas import ExcelWriter
from os.path import join
from pathlib import Path
from datetime import datetime
import re
def FBAwEB(organism,modelspath,model_variable,ensemblepath,solver):
    # This function is a python version of 'Script/main/FBAwEB.m'
    # https://github.com/skku-pdse/FBAwEB/blob/3e7799f19256c6a174b1a62389620b4a2d807eef/Script/main/FBAwEB.m

    # -----organism-----
    # ecoli / yeast / cho

    # -----models-----
    # GEMs_FBAwEB.mat directory

    # -----Model variable-----
    # ecoli_model_FBAwEB
    # ecoli_model_anaerobic_FBAwEB
    # yeast_model_800_glucose_FBAwEB
    # yeast_model_800_ethanol_FBAwEB
    # yeast_model_800_xylose_FBAwEB
    # cho_model_control_FBAwEB
    # cho_model_he_FBAwEB
    # cho_model_le_FBAwEB

    # -----ensemblepath-----
    # '$$$_Ensemble biomass_macro&FA +-2STDEV ... .xlsx' directory

    # -----solver-----
    # we used solver gurobi

    # result file name
    saveDir = "FBAwEB result"
    Path(saveDir).mkdir(parents=True, exist_ok=True)
    file2save='{0}_pfBAwEB result_{1}.xlsx'.format(organism,datetime.now().strftime("%b%d %H;%M"))
    file2save = join(saveDir, file2save)

    #load matlab model
    m_model=cobra.io.load_matlab_model(modelspath,variable_name=model_variable)
    m_model.solver=solver

    #read Ensemble biomass
    Ensemblexls=pd.ExcelFile(ensemblepath)
    proteindf=pd.read_excel(Ensemblexls, 'PROTsyn',index_col=0, header=0)
    dnadf=pd.read_excel(Ensemblexls,'DNAsyn',index_col=0, header=0)
    rnadf=pd.read_excel(Ensemblexls,'RNAsyn',index_col=0, header=0)
    lipiddf=pd.read_excel(Ensemblexls,'LIPIDsyn',index_col=0, header=0)
    biomassdf=pd.read_excel(Ensemblexls,'biomass',index_col=0, header=0)
    if organism == "yeast":
        carbdf = pd.read_excel(Ensemblexls, 'CARBsyn', index_col=0, header=0)
        fadf = pd.read_excel(Ensemblexls, 'FAsyn', index_col=0, header=0)
    elif organism == 'cho':
        carbdf = pd.read_excel(Ensemblexls, 'CARBsyn', index_col=0, header=0)
        fadf = pd.read_excel(Ensemblexls, 'FAsyn', index_col=0, header=0)
        facoadf=pd.read_excel(Ensemblexls, 'FACOAsyn', index_col=0, header=0)

    # ADD metabolite
    if (organism== o for o in ['ecoli','cho']) :
        m_prot=cobra.Metabolite('prot[c]',compartment='c')
        m_dna=cobra.Metabolite('dna[c]',compartment='c')
        m_rna = cobra.Metabolite('rna[c]', compartment='c')
        m_lipid = cobra.Metabolite('lipid[c]', compartment='c')
        m_model.add_metabolites([m_prot,m_dna,m_rna,m_lipid])
    elif organism == "yeast":
        m_prot=cobra.Metabolite('protein[c]',compartment='c')
        m_dna=cobra.Metabolite('DNA[c]',compartment='c')
        m_rna = cobra.Metabolite('RNA[c]', compartment='c')
        m_lipid = cobra.Metabolite('lipid[c]', compartment='c')
        m_carb = cobra.Metabolite('carbohydrate[c]',compartment='c')
        m_fa = cobra.Metabolite('fatty_acid[c]',compartment='c')
        m_model.add_metabolites([m_prot, m_dna, m_rna, m_lipid,m_carb,m_fa])

    if organism == "cho":
        m_carb = cobra.Metabolite('carbohydrate[c]', compartment='c' )
        m_rtotal_c = cobra.Metabolite('Rtotal[c]', compartment='c')
        m_rtotalcoa_c = cobra.Metabolite('Rtotalcoa[c]',compartment='c')
        m_rtotal_e = cobra.Metabolite('Rtotal[e]', compartment='e')
        m_rtotalcoa_e = cobra.Metabolite('Rtotalcoa[e]',compartment='e')
        m_rtotal_m = cobra.Metabolite('Rtotal[m]', compartment='m')
        m_rtotalcoa_m = cobra.Metabolite('Rtotalcoa[m]',compartment='m')
        m_rtotal_x = cobra.Metabolite('Rtotal[x]', compartment='x')
        m_rtotalcoa_x = cobra.Metabolite('Rtotalcoa[x]',compartment='x')
        m_model.add_metabolites([m_carb,m_rtotal_c,m_rtotalcoa_c,m_rtotal_e,m_rtotalcoa_e,m_rtotal_m,m_rtotalcoa_m,m_rtotal_x,m_rtotalcoa_x])


    # Add cobra reaction
    PROTsyn = cobra.Reaction(id="PROTsyn", name='PROTsyn')
    DNAsyn = cobra.Reaction(id="DNAsyn", name='DNAsyn')
    RNAsyn = cobra.Reaction(id="RNAsyn", name='RNAsyn')
    LIPIDsyn = cobra.Reaction(id="LIPIDsyn", name='LIPIDsyn')
    Biomass = cobra.Reaction(id="Biomass", name='Biomass')
    m_model.add_reactions([PROTsyn, DNAsyn, RNAsyn, LIPIDsyn, Biomass])
    if organism == "yeast" :
        CARBsyn = cobra.Reaction(id="CARBsyn", name='CARBsyn')
        FAsyn = cobra.Reaction(id='FAsyn', name='FAsyn')
        m_model.add_reactions([CARBsyn,FAsyn])
    elif organism == 'cho':
        CARBsyn = cobra.Reaction(id="CARBsyn", name='CARBsyn')

        FATTYACIDsyn_c = cobra.Reaction(id='FATTYACIDsyn_c', name='FATTYACIDsyn_c')
        FATTYACIDCOAsyn_c = cobra.Reaction(id='FATTYACIDCOAsyn_c', name='FATTYACIDCOAsyn_c')
        FATTYACIDsyn_e = cobra.Reaction(id='FATTYACIDsyn_e', name='FATTYACIDsyn_e')
        # FATTYACIDCOAsyn_e = cobra.Reaction(id='FACOAsyn_e', name='FACOAsyn_e')
        FATTYACIDsyn_m = cobra.Reaction(id='FATTYACIDsyn_m', name='FATTYACIDsyn_m')
        FATTYACIDCOAsyn_m = cobra.Reaction(id='FATTYACIDCOAsyn_m', name='FATTYACIDCOAsyn_m')
        FATTYACIDsyn_x = cobra.Reaction(id='FATTYACIDsyn_x', name='FATTYACIDsyn_x')
        FATTYACIDCOAsyn_x = cobra.Reaction(id='FATTYACIDCOAsyn_x', name='FATTYACIDCOAsyn_x')
        m_model.add_reactions([CARBsyn,FATTYACIDsyn_c,FATTYACIDCOAsyn_c,FATTYACIDsyn_e,FATTYACIDsyn_m,
                               FATTYACIDCOAsyn_m,FATTYACIDsyn_x,FATTYACIDCOAsyn_x])
    #rxn id, equation
    eqnDict={}
    for i in range(len(biomassdf)):
        i_model=m_model

        # Add string reaction eqn
        i_model.reactions.PROTsyn.build_reaction_from_string(proteindf.iloc[i, 0])
        i_model.reactions.DNAsyn.build_reaction_from_string(dnadf.iloc[i,0])
        i_model.reactions.RNAsyn.build_reaction_from_string(rnadf.iloc[i, 0])
        i_model.reactions.LIPIDsyn.build_reaction_from_string(lipiddf.iloc[i, 0])
        i_model.reactions.Biomass.build_reaction_from_string(biomassdf.iloc[i, 0])
        if organism == "yeast":
            i_model.remove_reactions(['r_FATTYACIDsyn_modi','r_2108','r_4041','r_4047','r_4048','r_4049','r_4050'])
            i_model.reactions.CARBsyn.build_reaction_from_string(carbdf.iloc[i, 0],fwd_arrow=" -> ")
            i_model.reactions.FAsyn.build_reaction_from_string(fadf.iloc[i,0])
        elif organism == 'cho':
            # delete existing FA rxns
            i_model.remove_reactions(['FATTYACIDSyn_c','FATTYACIDSyn_e','FATTYACIDSyn_m','FATTYACIDSyn_x',
                                      'FATTYACIDCOASyn_c','FATTYACIDCOASyn_m','FATTYACIDCOASyn_x'])

            # add new FA rxns
            i_model.reactions.CARBsyn.build_reaction_from_string(carbdf.iloc[i, 0])
            i_model.reactions.FATTYACIDsyn_c.build_reaction_from_string(fadf.iloc[i, 0])
            i_model.reactions.FATTYACIDCOAsyn_c.build_reaction_from_string(facoadf.iloc[i, 0])

            i_model.reactions.FATTYACIDsyn_e.build_reaction_from_string(fadf.iloc[i,0].replace('[c]','[e]'))
            # i_model.reactions.FATTYACIDCOAsyn_e.build_reaction_from_string(facoadf.iloc[i, 0].replace('[c]','[e]'))

            i_model.reactions.FATTYACIDsyn_m.build_reaction_from_string(fadf.iloc[i,0].replace('[c]','[m]'))
            i_model.reactions.FATTYACIDCOAsyn_m.build_reaction_from_string(facoadf.iloc[i, 0].replace('[c]','[m]'))

            i_model.reactions.FATTYACIDsyn_x.build_reaction_from_string(fadf.iloc[i,0].replace('[c]','[x]'))
            i_model.reactions.FATTYACIDCOAsyn_x.build_reaction_from_string(facoadf.iloc[i, 0].replace('[c]','[x]'))


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
        f_df.to_excel(writer,sheet_name='flux table')
        s_df.to_excel(writer,sheet_name='summary table')
        writer.save()




if __name__ == '__main__':
    organism=''
    modelspath=''
    model_variable=""
    ensemblepath=''
    solver='gurobi'

    FBAwEB(organism,modelspath,model_variable,ensemblepath,solver)

