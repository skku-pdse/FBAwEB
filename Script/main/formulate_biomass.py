
import pandas as pd
from pandas import ExcelWriter
from datetime import date, datetime
import random
import numpy as np
import os.path as path
from pathlib import Path
import cobra
import os
from tqdm import tqdm

class randomBiomass():
    """

    This module generates an ensemble of biomass equations, with the number specified by the user.


    1. Input

        file2read=''        # Directory of biomass excel file (Ex: 'Ecoli test1.xlsx')
                            # The "Reactant" component of lipids must conform to the prescribed format provided in the LIPIDsyn sheet.
                            # To determine the coefficient of each 'Reactant', such as pe160, the fatty acids and lipid kinds undergo individual normalization processes.
        sampling_n=         # The number of biomass equations to generate (Ex: 5000)
        macro_cols="A:W"    # This is the column range to use in the file2read (Ex: "Ecoli test1.xlsx').
                            Use the range "A:W" unless you comprehensively understand the codes involved
                            and intend to modify both the biomass excel file and the associated code accordingly.


    2. Output
        Output Directory = "{0}_Ensemble biomass_macro&FA +-2STDEV ... .xlsx' "  #{0} is organism name

        In the Output result,
            "Random coefficient" Sheet:
                                        "Data Average" and "Data Stdev" refer to the average and standard deviation, respectively, of the biomass composition data obtained from the biomass file (e.g., "Ecoli test1.xlsx").
                                        "Norm_Data Average" represents the normalized biomass composition data used as a reference for the ensemble of biomass equations.
                                        "Random Average" and "Random Stdec" represent the average and standard deviation of the biomass composition data for the resulting 5000 ensemble biomass equations.
                                        Columns from "0" to "4999" represent the different biomass compositions used for the 5000 biomass equations. These compositions were randomized within specific ranges based on the coefficient of variations (CVs).

            "ref" Sheet:
                                        The reference biomass equations formulated based on the biomass composition of "Norm_Data Average" column in the "Random coefiiceint" sheet.

            "PROTsyn", "DNAsyn", ... "Biomass" Sheets:
                                        "PROTsyn" sheet provides 5000 distinct protein synthesis equations, with each row index corresponding to the column header in the "Random coefficient" sheet.
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

    """

    def __init__(self,file2read,sampling_n,macro_cols):

        self.date = datetime.now().strftime("%b%d %H;%M")
        self._commonReadCols="A:H"
        self.file2read =file2read
        self.sampling_n =sampling_n
        self.macro_cols=macro_cols
        self.file2save="Ensemble biomass_macro&FA +-2STDEV_{0}.xlsx".format(self.date)
        self.BiomassSets=self.RandomBiomass()


    def _processBiomassFrame(self,_df,_in_or_out):
        if _in_or_out == 0 :
            if len(_df["Reactant"]) >=  len(_df["Product"]) :
                df_processed=_df[:len(_df["Reactant"])]
            else :
                df_processed= _df[:len(_df["Product"])]

        elif _in_or_out == -1 :
            df_processed =_df.iloc[:,:list(_df.columns).index("Product")]
            for i in range(len(_df)) :
                if pd.isna(_df["Reactant"][i]) == True :
                    df_processed=df_processed[:i]
                    break
        elif _in_or_out == 1:
            df_processed = _df.iloc[:, list(_df.columns).index("Product"):]
            for i in range(len(_df)):
                if pd.isna(_df["Product"][i]) == True:
                    df_processed = df_processed[:i]
                    break
        return df_processed

    def _returnSynthesisEquation(self,df_in,df_out):

        for i in range(len(df_in)) :
            coeff = df_in.iloc[i, 4]
            species = df_in.iloc[i, 2]
            compart = df_in.iloc[i, 3]
            if i == 0:
                sub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                sub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

        for j in range(len(df_out)):
            coeff = df_out.iloc[j, 2]
            species = df_out.iloc[j, 0]
            compart = df_out.iloc[j, 1]
            if j == 0:
                pro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                pro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

            # Synthesis equation
            syn = sub + " -> " + pro
        return syn

    def _makeRef(self,randomdf):
        refDict=dict()
        refCol="Norm_Data Average"
        refDict["PROTsyn_ref"]=self.PROTEIN(randomdf.loc['Protein',refCol])
        refDict["DNAsyn_ref"]=self.DNA_or_RNA(randomdf.loc['Dna',refCol],DNA_or_RNA="DNA")
        refDict["RNAsyn_ref"]=self.DNA_or_RNA(randomdf.loc['Rna',refCol],DNA_or_RNA="RNA")
        refDict["CARBsyn_ref"]=self.CARB(randomdf.loc['Carb',refCol])
        refDict["LIPIDsyn_ref"]=self.LIPID(randomdf.loc['Lipid',refCol],FA_CV=0.7145)
        refDict["BIOMASS_ref"]=self.BIOMASS(randomdf.loc['Lps',refCol], randomdf.loc['Murein',refCol])
        print (randomdf)

        return refDict


    def exportBiomassEqns(self):

        func=self.BiomassSets
        saveDir="RandomBiomassExcel"
        Path(saveDir).mkdir(parents=True,exist_ok=True)
        self.file2save=path.join(saveDir,self.file2save)

        ref=self._makeRef(randomdf=func[0])
        with ExcelWriter(self.file2save) as writer:
            func[0].to_excel(writer,sheet_name="Random coefficient")
            pd.DataFrame(ref,index=[0]).T.to_excel(writer, sheet_name='ref')
            pd.DataFrame(func[2]).to_excel(writer, sheet_name='PROTsyn')
            pd.DataFrame(func[3]).to_excel(writer, sheet_name='DNAsyn')
            pd.DataFrame(func[4]).to_excel(writer, sheet_name='RNAsyn')
            pd.DataFrame(func[5]).to_excel(writer, sheet_name='CARBsyn')
            pd.DataFrame(func[6]).to_excel(writer, sheet_name='LIPIDsyn')
            pd.DataFrame(func[7]).to_excel(writer, sheet_name='FAsyn')
            pd.DataFrame(func[8]).to_excel(writer, sheet_name='biomass')
            writer.save()
        print ("Result has been saved in this file: ",self.file2save)

    def RandomBiomass(self):
        pd.set_option("display.max_columns",999)

        _PROTlist=list()
        _DNAlist=list()
        _RNAlist=list()
        _CARBlist=list()
        _LIPIDlist = list()
        _BIOMASSlist=list()
        #Overall composition
        _all = pd.read_excel(self.file2read, sheet_name='Overall',usecols=self.macro_cols,convert_float=True)
        _all=_all.set_index(_all.columns[0])

        #USE average and stdev
        _all["mean"]=_all.mean(axis=1)
        _all["stdev"]=_all.iloc[:,:-1].std(axis=1)

        component_list_default=_all.index.values
        component_list = [ x for x in component_list_default if x not in ("Sum", "sum", "Total", "total", "Summation", "summation")]
        component_list.append("Sum")
        randomcoeffpd_default=pd.DataFrame({"Component":component_list,"Data Average":[np.nan]*len(component_list),"Data Stdev":[np.nan]*len(component_list),
                                    "Norm_Data Average": [np.nan] * len(component_list),
                                    "Random Average": [np.nan] * len(component_list), "Random Stdev": [np.nan] * len(component_list)})

        randomcoeffpd_default=randomcoeffpd_default.set_index(randomcoeffpd_default.columns[0])
        randomcoeffpd = randomcoeffpd_default.copy()

        for m in component_list[:-1]:
            randomcoeffpd.loc[m, "Data Average"] = _all.loc[m, "mean"]
            randomcoeffpd.loc[m, "Data Stdev"] = _all.loc[m, "stdev"]

        randomcoeffpd.index = randomcoeffpd.index.str.capitalize()
        new_index = ["Protein", "Dna", "Rna", "Carb","Lipid", "Lps", "Murein", "Inorganic", "Soluble pool", "Sum"]
        for i, ind in enumerate(randomcoeffpd.index.values):
            if ind.capitalize() not in new_index :
                if any(name.lower() in ind.lower() for name in ("ion", "metal", "debris")):
                    randomcoeffpd.index.values[i] = "Ion"
                elif "dna" in ind.lower():
                    randomcoeffpd.index.values[i] = "Dna"
                elif "rna" in ind.lower():
                    randomcoeffpd.index.values[i] = "Rna"
                elif "carbohydrate" in ind.lower():
                    randomcoeffpd.index.values[i] = "Carb"
        randomcoeffpd.reindex(new_index)

        # Amount of DNA and RNA should be expressed as ratio to Protein. DNA/Protein, RNA/Protein
        for m_i,m in enumerate(new_index[:-1]):
            if "dna" in m.lower():
                randomcoeffpd.loc[m,"Data Average"] = randomcoeffpd["Data Average"][m_i] * randomcoeffpd.loc["Protein","Data Average"]
            # for DNA/Protein , RNA/Protein
            elif "rna" in m.lower():
                randomcoeffpd.loc[m, "Data Average"] = randomcoeffpd["Data Average"][m_i] * randomcoeffpd.loc["Protein", "Data Average"]

        randomcoeffpd.loc["Sum", "Data Average"]= randomcoeffpd["Data Average"][:len(component_list)-1].sum(axis=0)


        fixed_coeff = randomcoeffpd.loc[["Inorganic", "Soluble pool"],"Data Average"].sum(axis=0)
        varing_coeff = randomcoeffpd.loc["Sum", "Data Average"]-fixed_coeff

        # Normalization of average value
        for i in randomcoeffpd.index.values:
            if i in ["Inorganic", "Soluble pool"] :
                randomcoeffpd.loc[i, "Norm_Data Average"] = randomcoeffpd.loc[i, "Data Average"]
            else :
                randomcoeffpd.loc[i, "Norm_Data Average"] = randomcoeffpd.loc[i, "Data Average"] / varing_coeff * (1 - fixed_coeff)

        randomcoeffpd.loc["Sum","Norm_Data Average"] =randomcoeffpd["Norm_Data Average"][:len(component_list)-1].sum(axis=0)

        col_k=0
        while len(randomcoeffpd.columns) < self.sampling_n + len(randomcoeffpd_default.columns) :
            coe_min_list=list()
            coe_max_list=list()
            for r in range(5) : #Prot, DNA, RNA, Carb, Lipid
                coe_min=randomcoeffpd["Norm_Data Average"][r] - 2 * randomcoeffpd["Data Stdev"][r]
                coe_max = randomcoeffpd["Norm_Data Average"][r] + 2 * randomcoeffpd["Data Stdev"][r]
                if coe_min < 0 :
                    coe_min =0
                coe_min_list.append(coe_min)
                coe_max_list.append(coe_max)

            PROTcoe = random.uniform(coe_min_list[0], coe_max_list[0])
            DNAcoe = random.uniform(coe_min_list[1], coe_max_list[1])
            RNAcoe = random.uniform(coe_min_list[2], coe_max_list[2])
            CARBcoe =random.uniform(coe_min_list[3], coe_max_list[3])
            LIPIDcoe = random.uniform(coe_min_list[4], coe_max_list[4])
            LPScoe= randomcoeffpd["Norm_Data Average"][5]
            Mureincoe=randomcoeffpd["Norm_Data Average"][6]
            Inorganiccoe=randomcoeffpd["Norm_Data Average"][7]
            Solublecoe=randomcoeffpd["Norm_Data Average"][8]

            varing_coe=PROTcoe+DNAcoe+RNAcoe+CARBcoe+LIPIDcoe
            fixed_coe=LPScoe+Mureincoe+Inorganiccoe+ Solublecoe

            PROTcoe_nor= PROTcoe/varing_coe*(1-fixed_coe)
            DNAcoe_nor=DNAcoe/varing_coe*(1-fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe*(1-fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe*(1-fixed_coe)

            totalwt_nor=PROTcoe_nor+DNAcoe_nor+RNAcoe_nor+CARBcoe_nor+LIPIDcoe_nor+fixed_coe

            randomcoeffpd[str(col_k)] = [PROTcoe_nor, DNAcoe_nor, RNAcoe_nor, CARBcoe_nor, LIPIDcoe_nor, LPScoe, Mureincoe,
                                         Inorganiccoe, Solublecoe, totalwt_nor]
            col_k += 1


        # Random result summary
        for i, r in enumerate(randomcoeffpd.index.values):
            randomcoeffpd.loc[r, "Random Average"] = randomcoeffpd.iloc[i, len(randomcoeffpd_default.columns):].mean(
                skipna=True)
            randomcoeffpd.loc[r, "Random Stdev"] = randomcoeffpd.iloc[i, len(randomcoeffpd_default.columns):].std(
                skipna=True)

        # multiple overall composition loop
        for col_i in tqdm(range(len(randomcoeffpd.columns))):

            if col_i >= len(randomcoeffpd_default.columns) :
                randomcoeffpd.iloc[len(component_list)-1,col_i]= randomcoeffpd.iloc[:len(component_list)-1, col_i].sum(axis=0)

                _PROTsyn=self.PROTEIN(randomcoeffpd.iloc[0,col_i])
                _DNAsyn=self.DNA_or_RNA(randomcoeffpd.iloc[1,col_i],DNA_or_RNA="DNA")
                _RNAsyn=self.DNA_or_RNA(randomcoeffpd.iloc[2,col_i],DNA_or_RNA="RNA")
                _CARBsyn = self.CARB(randomcoeffpd.iloc[3, col_i])
                _LIPIDsyn=self.LIPID(randomcoeffpd.iloc[4,col_i],FA_CV=0.7145)
                _BIOMASSsyn= self.BIOMASS(randomcoeffpd.iloc[5,col_i],randomcoeffpd.iloc[6,col_i])
                _PROTlist.append(_PROTsyn)
                _DNAlist.append(_DNAsyn)
                _RNAlist.append(_RNAsyn)
                _CARBlist.append(_CARBsyn)
                _LIPIDlist.append(_LIPIDsyn)
                _BIOMASSlist.append(_BIOMASSsyn)

        return randomcoeffpd,["NA"],_PROTlist,_DNAlist,_RNAlist,_CARBlist,_LIPIDlist,["NA"],_BIOMASSlist


    def PROTEIN(self,Prot_wt_per_DCW):

        PROTpd= pd.read_excel(self.file2read, sheet_name='PROTsyn', usecols=self._commonReadCols ,convert_float=True)
        PROTin= self._processBiomassFrame(_df=PROTpd,_in_or_out=-1)
        ##################

        PROTin["mmol/gDCW"]=Prot_wt_per_DCW*PROTin.iloc[:,0]/PROTin.iloc[:,1]*1000

        #PROT out
        PROTout = self._processBiomassFrame(_df=PROTpd, _in_or_out=1)
        # PROTout coeff cal
        for k_i,k in enumerate(PROTout["Product"]):
            if k.lower() == 'h2o':
                PROTout.at[k_i, "mmol/gDCW.1"] = PROTin.iloc[:20, 4].sum(axis=0)
            elif k.lower() in ("protein", "prot") :
                PROTout.at[k_i, "mmol/gDCW.1"] = 1
            else : # amino acids
                for ami_i, ami in enumerate(PROTin["Reactant"]):
                    if k.lower().replace('trna','') in ami.lower():
                        PROTout.at[k_i, "mmol/gDCW.1"]= PROTin.iloc[ami_i ,4]
                        break

        PROTsyn= self._returnSynthesisEquation(PROTin,PROTout)
        return PROTsyn

    def DNA_or_RNA(self,wt_per_DCW, DNA_or_RNA):

        DNApd = pd.read_excel(self.file2read, sheet_name=DNA_or_RNA+"syn", usecols=self._commonReadCols, convert_float=True)
        DNAin = self._processBiomassFrame(DNApd, -1)
        # Nucleotides in
        for i in range(len(DNAin)):
            gDNApergDCW = DNAin.iloc[i, 0]
            MW = DNAin.iloc[i, 1]
            DNAin.at[i, "mmol/gDCW"] = wt_per_DCW * gDNApergDCW / MW * 1000

        # Nucleotides out
        DNAout = self._processBiomassFrame(DNApd, 1)

        for p_ind,p in enumerate(DNAout["Product"]):
            if any( c == p.lower() for c in ("ppi","pi","pp","p") ):
                DNAout.at[p_ind, "mmol/gDCW.1"] = DNAin.iloc[:, 4].sum(axis=0)
            elif any (d == p.lower() for d in ("dna","rna") ) :
                DNAout.iloc[p_ind,2] = 1
        DNAsyn = self._returnSynthesisEquation(DNAin, DNAout)
        return DNAsyn

    def CARB(self,Carb_per_DCW):
        CARBpd = pd.read_excel(self.file2read, sheet_name='CARBsyn', usecols=self._commonReadCols, convert_float=True)
        CARBin = self._processBiomassFrame(CARBpd, -1)
        # Note
        #     gpergCARB = CARBin.iloc[:, 0]
        #     MW = CARBin.iloc[:, 1]
        #     CARBin["mmol/gDCW"] = Carb_per_DCW * gpergCARB / MW * 1000
        CARBin["mmol/gDCW"] = Carb_per_DCW * CARBin.iloc[:, 0] / CARBin.iloc[:, 1] * 1000
        CARBout = self._processBiomassFrame(CARBpd, 1)
        for p_ind,p in enumerate(CARBout["Product"]):
            if any( c == p.lower() for c in ("carbohydrate","carb","carbo")) :
                CARBout.at[p_ind,"mmol/gDCW.1"] = 1

        CARBsyn=self._returnSynthesisEquation(CARBin,CARBout)
        return CARBsyn

    def LIPID(self,LIPID_wt_per_DCW,FA_CV):
        LIPIDpd = pd.read_excel(self.file2read, sheet_name='LIPIDsyn', usecols=self._commonReadCols, convert_float=True)

        LIPIDin = self._processBiomassFrame(_df=LIPIDpd, _in_or_out=-1)

        _Lipids = LIPIDin["Reactant"].str.extract(r'(\D+)').values.tolist()
        _FAs = LIPIDin["Reactant"].str.extract(r'(\d+)').values.tolist()
        Lipids = []
        FAs = []
        for i in _Lipids:
            if i[0] not in Lipids:
                Lipids.append(i[0])
        for j in _FAs:
            if j[0] not in FAs:
                FAs.append(j[0])

        # FATTY ACID
        Lipid_mono_ratio_df = pd.DataFrame(columns=Lipids, index=FAs)

        for l_i, l in enumerate(Lipids):
            for l_j, j in enumerate(FAs):
                lipid_fa = l + j
                gpergLipid = LIPIDin["g/g Lipid"][LIPIDin['Reactant'] == lipid_fa].values[0]
                Lipid_mono_ratio_df.iloc[l_j, l_i] = gpergLipid / LIPIDin["g/g Lipid"].sum()

        Lipid_mono_ratio_df["FA_Sum"] = Lipid_mono_ratio_df.sum(axis=1)
        Lipid_mono_ratio_df.loc["Sum"] = Lipid_mono_ratio_df.sum(axis=0)
        Fatty_acid_ratio = pd.DataFrame({"ref": list(Lipid_mono_ratio_df["FA_Sum"][:-1])}, index=FAs)

        new_fa = list()
        # Fatty acid random sampling

        fa_i = 0
        while len(new_fa) < len(Fatty_acid_ratio):  # 3
            ref_ratio = Fatty_acid_ratio['ref'][fa_i]
            coe_min = ref_ratio * (1 - 2 * FA_CV)
            coe_max = ref_ratio * (1 + 2 * FA_CV)
            if coe_min < 0:
                coe_min = 0
            rand_coe = random.uniform(coe_min, coe_max)
            new_fa.append(rand_coe)
            fa_i += 1

        # Fatty acid ratio normalization
        norm_new_fa = np.divide(new_fa, sum(new_fa))
        Fatty_acid_ratio["normalized_fa_coe"] = norm_new_fa

        # End
        Lipid_kind_ratio = pd.DataFrame(columns=Lipids, index=["ref"])
        for i, i_c in enumerate(Lipid_kind_ratio.columns):
            Lipid_kind_ratio[i_c] = list(Lipid_mono_ratio_df.loc["Sum"])[:-1][i]

        for l_i, l_kind in enumerate(Lipids):
            for f_j, f_kind in enumerate(FAs):
                Lipid_mono_ratio_df.iloc[f_j, l_i] = Fatty_acid_ratio["normalized_fa_coe"].loc[f_kind] * \
                                                     Lipid_kind_ratio[l_kind]

        lipid_coeff_dict = dict()
        for c in Lipids:
            for f in FAs:
                lipid_mono = c + f
                monomer_g = Lipid_mono_ratio_df.loc[f, c]
                lipid_coeff_dict[lipid_mono] = monomer_g

        for i in range(len(LIPIDin)):
            lipid_comp = LIPIDin["Reactant"][i]
            gpergP = lipid_coeff_dict[lipid_comp]
            MW = LIPIDin.iloc[i, 1]
            LIPIDin.loc[i, "mmol/gDCW"] = gpergP / MW * 1000 * LIPID_wt_per_DCW

        # LIPID out
        LIPIDout = self._processBiomassFrame(_df=LIPIDpd, _in_or_out=1)
        LIPIDsyn = self._returnSynthesisEquation(LIPIDin, LIPIDout)

        return LIPIDsyn

    # def LIPID(self,LIPID_wt_per_DCW,FA_CV):
    #     LIPIDpd = pd.read_excel(self.file2read, sheet_name='LIPIDsyn', usecols=self._commonReadCols, convert_float=True)
    #
    #     LIPIDin = self._processBiomassFrame(_df=LIPIDpd, _in_or_out=-1)
    #
    #     # FATTY ACID#####################################
    #     Lipid_mono_ratio_df = pd.DataFrame({"pe": np.nan * 3, "pg": np.nan * 3, "clpn": np.nan * 3},
    #                                        index=["160", "161", "181"])
    #
    #     for l_i, l in enumerate(LIPIDin["Reactant"]):
    #         if "pe" in l:
    #             col_loc = 0
    #         elif "pg" in l:
    #             col_loc = 1
    #         elif "clpn" in l:
    #             col_loc = 2
    #
    #         if "160" in l:
    #             row_loc = 0
    #         elif "161" in l:
    #             row_loc = 1
    #         elif "181" in l:
    #             row_loc = 2
    #         Lipid_mono_ratio_df.iloc[row_loc, col_loc] = LIPIDin["g/g Lipid"][l_i] / LIPIDin["g/g Lipid"].sum()
    #
    #     Lipid_mono_ratio_df["FA_Sum"] = Lipid_mono_ratio_df.sum(axis=1)
    #     Lipid_mono_ratio_df.loc["Sum"] = Lipid_mono_ratio_df.sum(axis=0)
    #     Fatty_acid_ratio = pd.DataFrame({"ref": list(Lipid_mono_ratio_df["FA_Sum"][:-1])}, index=["160", "161", "181"])
    #     new_fa=list()
    #     #Fatty acid random sampling
    #
    #     fa_i=0
    #     while len(new_fa) < len(Fatty_acid_ratio): #3
    #         ref_ratio=Fatty_acid_ratio['ref'][fa_i]
    #         coe_min = ref_ratio*(1 - 2 * FA_CV)
    #         coe_max = ref_ratio*(1 + 2 * FA_CV)
    #         if coe_min < 0 :
    #             coe_min=0
    #         rand_coe=random.uniform(coe_min,coe_max)
    #         new_fa.append(rand_coe)
    #         fa_i+=1
    #
    #     #Fatty acid ratio normalization
    #     norm_new_fa=np.divide(new_fa, sum(new_fa))
    #     Fatty_acid_ratio["normalized_fa_coe"] = norm_new_fa
    #
    #     #End
    #     Lipid_kind_ratio = pd.DataFrame(
    #         {"pe": [Lipid_mono_ratio_df.loc["Sum", "pe"]], "pg": [Lipid_mono_ratio_df.loc["Sum", "pg"]],
    #          "clpn": [Lipid_mono_ratio_df.loc["Sum", "clpn"]]}, index=["ref"])
    #
    #     for l_i, l_kind in enumerate(["pe", "pg", "clpn"]):
    #         Lipid_mono_ratio_df.iloc[0, l_i] = Fatty_acid_ratio["normalized_fa_coe"][0] * Lipid_kind_ratio[l_kind][0]
    #         Lipid_mono_ratio_df.iloc[1, l_i] = Fatty_acid_ratio["normalized_fa_coe"][1] * Lipid_kind_ratio[l_kind][0]
    #         Lipid_mono_ratio_df.iloc[2, l_i] = Fatty_acid_ratio["normalized_fa_coe"][2] * Lipid_kind_ratio[l_kind][0]
    #
    #     lipid_coeff_dict = dict()
    #     for c in ["pe", "pg", "clpn"]:
    #         for f in ["160", "161", "181"]:
    #             lipid_mono = c + f
    #             monomer_g = Lipid_mono_ratio_df.loc[f, c]
    #             lipid_coeff_dict[lipid_mono] = monomer_g
    #
    #     for i in range(len(LIPIDin)):
    #         lipid_comp = LIPIDin["Reactant"][i]
    #         gpergP = lipid_coeff_dict[lipid_comp]
    #         MW = LIPIDin.iloc[i, 1]
    #         LIPIDin.loc[i, "mmol/gDCW"] = gpergP / MW * 1000 * LIPID_wt_per_DCW
    #
    #     # LIPID out
    #     LIPIDout = self._processBiomassFrame(_df=LIPIDpd, _in_or_out=1)
    #     LIPIDsyn = self._returnSynthesisEquation(LIPIDin, LIPIDout)
    #
    #     return LIPIDsyn

    def BIOMASS(self, LPS_per_DCW, Murein_per_DCW):
        BIOMASSpd = pd.read_excel(self.file2read, sheet_name='Biomass', usecols=self._commonReadCols, convert_float=True)
        BIOMASSin = self._processBiomassFrame(BIOMASSpd, -1)
        for s_i,s in enumerate(BIOMASSin["Reactant"]):
            if pd.isna(BIOMASSin["mmol/gDCW"][s_i]) == True:
                if any (r in s.lower() for r in ("lipid4","colipa")) :
                    Molecularweight = BIOMASSin["MW"][s_i]
                    BIOMASSin.loc[s_i,"mmol/gDCW"] = LPS_per_DCW / Molecularweight * 1000
                elif "murein" in s.lower():
                    Molecularweight = BIOMASSin["MW"][s_i]
                    BIOMASSin.loc[s_i,"mmol/gDCW"] = Murein_per_DCW / Molecularweight * 1000
            else:
                pass

        BIOMASSout = self._processBiomassFrame(BIOMASSpd, 1)
        BIOMASSsyn=self._returnSynthesisEquation(BIOMASSin,BIOMASSout)

        return BIOMASSsyn


class randomFBA():

    def __init__(self,model,RandomBiomassExcelFile,n_Biomass):
        if type(model)==cobra.core.model.Model :
            self.model=model
        elif '.mat' in model :
            self.model= cobra.io.load_matlab_model(model)
        elif '.xml' in model :
            self.model=cobra.io.read_sbml_model(model)

        self.model_modi = self.model.copy()
        self.RandomBiomassExcelFile=RandomBiomassExcelFile
        self.n_Biomass = n_Biomass

        biomass_xls=pd.ExcelFile(self.RandomBiomassExcelFile)
        self.protSheet=pd.read_excel(biomass_xls,"PROTsyn").drop(columns="Unnamed: 0")
        self.dnaSheet=pd.read_excel(biomass_xls,"DNAsyn").drop(columns="Unnamed: 0")
        self.rnaSheet=pd.read_excel(biomass_xls,"RNAsyn").drop(columns="Unnamed: 0")
        self.carbSheet = pd.read_excel(biomass_xls, "CARBsyn").drop(columns="Unnamed: 0")
        self.lipidSheet=pd.read_excel(biomass_xls,"LIPIDsyn").drop(columns="Unnamed: 0")
        self.biomassSheet=pd.read_excel(biomass_xls,"biomass").drop(columns="Unnamed: 0")


    def saveResult(self):

        for i in self.model_modi.reactions:
            if "biomass" in i.id:
                self.model_modi = self.model_modi.remove_reactions([self.model_modi.reactions.get_by_id(i.id)])

        for r in range(self.n_Biomass):

            protEqn=self.protSheet.iloc[r,0]
            dnaEqn=self.dnaSheet.iloc[r,0]
            rnaEqn=self.rnaSheet.iloc[r,0]
            carbEqn=self.carbSheet.iloc[r,0]
            lipidEqn=self.lipidSheet.iloc[r,0]
            biomassEqn=self.biomassSheet.iloc[r,0]
            eqn_list=[protEqn,dnaEqn,rnaEqn,carbEqn,lipidEqn,biomassEqn]

            rxn2add=dict()
            for m_i, m in enumerate(["PROTsyn_modi",'DNAsyn_modi','RNAsyn_modi','CARBsyn_modi','LIPIDsyn_modi','biomass_modi']):
                rxn2add[m]=eqn_list[m_i]
                rxn=cobra.Reaction(m)
                rxn.name=m
                rxn.subsystem="Biomass"
                rxn.lower_bound=0
                rxn.upper_bound=1000
                rxn.reaction= rxn2add[m]
                self.model_modi=self.model_modi.add_reactions([rxn])


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    testfile1= "Ecoli test1.xlsx"

    a=randomBiomass(file2read=testfile1,sampling_n=5000,macro_cols="A:W")
    b=a.exportBiomassEqns()



