#-*-coding:utf-8-*-
import pandas as pd
from pandas import ExcelWriter
from datetime import datetime
import numpy as np
from pathlib import Path
import os


class makeBiomass():

    def __init__(self, file2read):
        date = datetime.now().strftime("%b%d %H;%M")
        self.file2read=file2read        
        saveDir="RandomBiomassExcel"
        Path(saveDir).mkdir(parents=True,exist_ok=True)
        savefile = 'Cho_macro_+-25% biomass_{}.xlsx'.format(date)
        self.file2save=os.path.join(saveDir,savefile)

    def _processBiomassFrame(self, df, in_or_out):
        if in_or_out == 0:
            if len(df["Reactant"]) >= len(df["Product"]):
                df_processed = df[:len(df["Reactant"])]
            else:
                df_processed = df[:len(df["Product"])]

        elif in_or_out == -1:
            df_processed = df.iloc[:, :list(df.columns).index("Product")]
            for i in range(len(df)):
                if pd.isna(df["Reactant"][i]) == True:
                    df_processed = df_processed[:i]
                    break

        elif in_or_out == 1:
            df_processed = df.iloc[:, list(df.columns).index("Product"):]
            for i in range(len(df)):
                if pd.isna(df["Product"][i]) == True:
                    df_processed = df_processed[:i]
                    break
        return df_processed

    def _returnSynthesisEquation(self, df_in, df_out):

        for i in range(len(df_in)):
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

    def RandomBiomass(self, Protein, DNA, RNA, Carbohydrate, Lipid):

        all = pd.read_excel(self.file2read, sheet_name='Overall', usecols="A:Q", convert_float=True,index_col=0)

        #USE average and stdev
        all["mean"]=all.mean(axis=1)
        all["stdev"]=all.iloc[:,:-1].std(axis=1)

        componentList=["Protein", "DNA", "RNA", "Carbohydrate", "Lipid", "Sum"]
        randomcoeffpd_default=pd.DataFrame({"Component":componentList,"Data Average":[np.nan]*len(componentList),"Data Stdev":[np.nan]*len(componentList),
                                    "Norm_Data Average": [np.nan] * len(componentList),
                                    "Random Average": [np.nan] * len(componentList), "Random Stdev": [np.nan] * len(componentList)})

        randomcoeffpd=randomcoeffpd_default.copy()
        randomcoeffpd_default.set_index(randomcoeffpd_default.columns[0])
        for n,m in enumerate(componentList[:-1]):
            macro_mean_mass=all.loc[m,"mean"]
            macro_stdev=all.loc[m,"stdev"]
            randomcoeffpd.at[n,"Data Average"]=macro_mean_mass
            randomcoeffpd.at[n,"Data Stdev"] = macro_stdev

        randomcoeffpd.at[len(componentList)-1, "Data Average"] = randomcoeffpd["Data Average"][:len(componentList)-1].sum(axis=0)

        # Normalization of average value
        for i in range(len(componentList)):
            randomcoeffpd.at[i, "Norm_Data Average"] = randomcoeffpd["Data Average"][i] / randomcoeffpd["Data Average"][len(componentList)-1]

        randomcoeffpd.at[len(componentList)-1, "Norm_Data Average"] = randomcoeffpd["Norm_Data Average"][:len(componentList)-1].sum(axis=0)

        if Protein == 'max':
            PROTcoe = randomcoeffpd["Norm_Data Average"][0] * 1.25
        elif Protein == 'min':
            PROTcoe = randomcoeffpd["Norm_Data Average"][0] * 0.75
        elif Protein == "mean":
            PROTcoe = randomcoeffpd["Norm_Data Average"][0]

        if DNA == 'max':
            DNAcoe = randomcoeffpd["Norm_Data Average"][1] * 1.25
        elif DNA == 'min':
            DNAcoe = randomcoeffpd["Norm_Data Average"][1]  * 0.75
        elif DNA == 'mean':
            DNAcoe = randomcoeffpd["Norm_Data Average"][1]

        if RNA == 'max':
            RNAcoe = randomcoeffpd["Norm_Data Average"][2] * 1.25
        elif RNA == 'min':
            RNAcoe =randomcoeffpd["Norm_Data Average"][2] * 0.75
        elif RNA == 'mean':
            RNAcoe = randomcoeffpd["Norm_Data Average"][2]


        if Carbohydrate == 'max':
            CARBcoe = randomcoeffpd["Norm_Data Average"][3] * 1.25
        elif Carbohydrate == 'min':
            CARBcoe = randomcoeffpd["Norm_Data Average"][3] * 0.75
        elif Carbohydrate == 'mean':
            CARBcoe = randomcoeffpd["Norm_Data Average"][3]

        if Lipid == 'max':
            LIPIDcoe = randomcoeffpd["Norm_Data Average"][4] * 1.25
        elif Lipid == 'min':
            LIPIDcoe = randomcoeffpd["Norm_Data Average"][4] * 0.75
        elif Lipid == 'mean':
            LIPIDcoe = randomcoeffpd["Norm_Data Average"][4]

        if (Protein == 'max') or (Protein == 'min'):
            varing_coe = DNAcoe + RNAcoe + CARBcoe + LIPIDcoe
            fixed_coe = PROTcoe
            PROTcoe_nor = PROTcoe
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)

            varing_coe_nor = DNAcoe_nor + RNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor
        elif (DNA == 'max') or (DNA == 'min'):
            varing_coe = PROTcoe + RNAcoe + CARBcoe+ LIPIDcoe
            fixed_coe = DNAcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)

            varing_coe_nor = PROTcoe_nor + RNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor
        elif (RNA == 'max') or (RNA == 'min'):
            varing_coe = PROTcoe + DNAcoe + CARBcoe+ LIPIDcoe
            fixed_coe = RNAcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)

            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor

        elif (Carbohydrate == 'max') or (Carbohydrate == 'min'):
            varing_coe = PROTcoe + DNAcoe + RNAcoe + LIPIDcoe
            fixed_coe = CARBcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor + LIPIDcoe_nor


        elif (Lipid == 'max') or (Lipid == 'min'):
            varing_coe = PROTcoe + DNAcoe + RNAcoe + CARBcoe
            fixed_coe = LIPIDcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor + CARBcoe_nor

        else : # all mean
            varing_coe = PROTcoe + DNAcoe + RNAcoe + CARBcoe + LIPIDcoe
            fixed_coe = 0
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor


        # totalwt_nor=PROTcoe_nor+DNAcoe_nor+RNAcoe_nor+GLYCOcoe_nor+LIPIDcoe_nor + fixed_coe
        totalwt_nor = varing_coe_nor + fixed_coe

        # randomcoeffpd[str(col_k)]=[PROTcoe_nor,DNAcoe_nor,RNAcoe_nor,GLYCOcoe_nor,LIPIDcoe_nor, LPScoe, MUREINcoe, INORGANICcoe,SOLUBLEPcoe,totalwt_nor]
        Final_coe_list = [PROTcoe_nor, DNAcoe_nor, RNAcoe_nor, CARBcoe_nor, LIPIDcoe_nor, totalwt_nor]

        return Final_coe_list

    def Formulate_min_max(self):

        minmax_list = ['Prot_min', "Prot_max", "DNA_min", "DNA_max", "RNA_min", "RNA_max", "Carb_min", "Carb_max", "Lipid_min", "Lipid,max","Ref"]
        coeff_index = ["PROTcoe_nor", "DNAcoe_nor", "RNAcoe_nor", "CARBcoe_nor", "LIPIDcoe_nor","totalwt_nor"]
        coeff_pd = pd.DataFrame(index=coeff_index, columns=minmax_list)
        rxn_list=["ProtSyn","DNASyn","RNASyn","CarbSyn", "LipidSyn","FATTYACIDSyn_c","FATTYACIDSyn_e","FATTYACIDSyn_m","FATTYACIDSyn_x",
                  "FATTYACUDCOASyn_c","FATTYACUDCOASyn_m","FATTYACUDCOASyn_x","BiomassSyn"]
        minmax_pd=pd.DataFrame(index=rxn_list, columns=minmax_list)

        for col_i, head in enumerate(minmax_list):
            if col_i == 0:
                coeff_list = self.RandomBiomass("min", "mean", "mean", "mean", "mean")
            elif col_i == 1:
                coeff_list = self.RandomBiomass("max", "mean", "mean", "mean", "mean")
            elif col_i == 2:
                coeff_list = self.RandomBiomass("mean", "min", "mean", "mean", "mean")
            elif col_i == 3:
                coeff_list = self.RandomBiomass("mean", "max", "mean", "mean", "mean")
            elif col_i == 4:
                coeff_list = self.RandomBiomass("mean", "mean", "min", "mean", "mean")
            elif col_i == 5:
                coeff_list = self.RandomBiomass("mean", "mean", "max", "mean", "mean")
            elif col_i == 6:
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "min", "mean")
            elif col_i == 7:
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "max", "mean")
            elif col_i == 8:
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "mean", "min")
            elif col_i == 9:
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "mean", "max")
            elif col_i == 10:
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "mean", "mean")
            minmax_pd[head][0] = self.PROTEIN(coeff_list[0])
            minmax_pd[head][1] = self.DNA_or_RNA(coeff_list[1], DNA_or_RNA="DNA")
            minmax_pd[head][2] = self.DNA_or_RNA(coeff_list[2], DNA_or_RNA="RNA")
            minmax_pd[head][3] = self.CARB(coeff_list[3])
            minmax_pd[head][4] = self.LIPID(coeff_list[4],FA_variation=False)[0]
            minmax_pd[head][5] = self.LIPID(coeff_list[4], FA_variation=False)[1]["ref"][0]
            minmax_pd[head][6] = self.LIPID(coeff_list[4], FA_variation=False)[1]["ref"][1]
            minmax_pd[head][7] = self.LIPID(coeff_list[4], FA_variation=False)[1]["ref"][2]
            minmax_pd[head][8] = self.LIPID(coeff_list[4], FA_variation=False)[1]["ref"][3]
            minmax_pd[head][9] = self.LIPID(coeff_list[4], FA_variation=False)[2]["ref"][0]
            minmax_pd[head][10] = self.LIPID(coeff_list[4], FA_variation=False)[2]["ref"][1]
            minmax_pd[head][11] = self.LIPID(coeff_list[4], FA_variation=False)[2]["ref"][2]
            minmax_pd[head][12] = self.BIOMASS()

            for k in range(len(coeff_index)):
                coeff_pd[head][k] = coeff_list[k]

        with ExcelWriter(self.file2save) as writer:
            coeff_pd.to_excel(writer,sheet_name="Random coefficient")
            minmax_pd.to_excel(writer, sheet_name="25% minmax biomass")
            writer.save()

    def PROTEIN(self,Prot_wt_per_DCW):

        PROTpd= pd.read_excel(self.file2read, sheet_name='PROTsyn', usecols="A:H" ,convert_float=True)
        PROTin=self._processBiomassFrame(df=PROTpd,in_or_out=-1)

        PROTin["mmol/gDCW"] = Prot_wt_per_DCW * PROTin.iloc[:, 0] / PROTin.iloc[:, 1] * 1000

        # PROT out
        PROTout = self._processBiomassFrame(df=PROTpd, in_or_out=1)
        # PROTout coeff cal
        for k_i, k in enumerate(PROTout["Product"]):
            if k.lower() == 'h2o':
                PROTout.at[k_i, "mmol/gDCW.1"] = PROTin.iloc[:20, 4].sum(axis=0)
            elif k.lower() in ("protein", "prot"):
                PROTout.at[k_i, "mmol/gDCW.1"] = 1
            else:  # amino acids
                for ami_i, ami in enumerate(PROTin["Reactant"]):
                    if k.lower().replace('trna', '') in ami:
                        PROTout.at[k_i, "mmol/gDCW.1"] = PROTin.iloc[ami_i, 4]
                        break
        PROTsyn = self._returnSynthesisEquation(PROTin, PROTout)
        return PROTsyn

    def DNA_or_RNA(self,wt_per_DCW, DNA_or_RNA):

        DNApd = pd.read_excel(self.file2read, sheet_name=DNA_or_RNA+"syn", usecols="A:H", convert_float=True)
        DNAin = self._processBiomassFrame(df=DNApd,in_or_out=-1)

        #
        for i in range(len(DNAin)):
            gDNApergDCW = DNAin.iloc[i, 0]
            MW = DNAin.iloc[i, 1]
            DNAin.at[i, "mmol/gDCW"] = wt_per_DCW * gDNApergDCW / MW * 1000

        DNAout = self._processBiomassFrame(df=DNApd,in_or_out=1)

        # ppi
        for p_ind, p in enumerate(DNAout["Product"]):
            if DNAout["Product"][p_ind].lower() == "ppi":
                DNAout.at[p_ind, "mmol/gDCW.1"] = DNAin.iloc[:, 4].sum(axis=0)


        # PROTEIN BIOMASS
        DNAsyn = self._returnSynthesisEquation(df_in=DNAin,df_out=DNAout)
        return DNAsyn

    def CARB(self, Carb_per_DCW):
        CARBpd = pd.read_excel(self.file2read, sheet_name='CARBsyn', usecols="A:H", convert_float=True)
        CARBin = self._processBiomassFrame(df=CARBpd,in_or_out=-1)
        for i in range(len(CARBin)):
            gpergCARB = CARBin.iloc[i, 0]
            MW = CARBin.iloc[i, 1]
            CARBin.at[i, "mmol/gDCW"] = Carb_per_DCW * gpergCARB / MW * 1000

        CARBout=self._processBiomassFrame(df=CARBpd,in_or_out=1)

        for p_ind, p in enumerate(CARBout["Product"]):
            if any(c == p.lower() for c in ("carbohydrate", "carb", "carbo")):
                CARBout.at[p_ind, "mmol/gDCW.1"] = 1

        CARBsyn=self._returnSynthesisEquation(df_in=CARBin,df_out=CARBout)
        return CARBsyn

    def LIPID(self,LIPID_wt_per_DCW,FA_variation):
        LIPIDpd= pd.read_excel(self.file2read, sheet_name='LIPIDsyn', usecols="A:H" ,convert_float=True)

        #Lipid in
        LIPIDin=self._processBiomassFrame(LIPIDpd,-1)

        for i in range(len(LIPIDin)):
            gpergLipid=LIPIDin.iloc[i,0]
            MW=LIPIDin.iloc[i,1]
            LIPIDin.at[i,"mmol/gDCW"]= LIPID_wt_per_DCW*gpergLipid/MW*1000

        #Lipid out
        LIPIDout = self._processBiomassFrame(LIPIDpd, 1)
        LIPIDsyn = self._returnSynthesisEquation(df_in=LIPIDin,df_out=LIPIDout)

        FA_dict=dict()
        FAcoa_dict=dict()

        FApd = pd.read_excel(self.file2read, sheet_name='FATTYACIDsyn', usecols="A:H", convert_float=True)
        FAcoapd = pd.read_excel(self.file2read, sheet_name='FATTYACIDCOAsyn', usecols="A:H", convert_float=True)
        FAin = self._processBiomassFrame(FApd, -1)
        FAcoain = self._processBiomassFrame(FAcoapd, -1)
        FA_list_default = FAin["Reactant"]
        FA_list = [x for x in FA_list_default if
                   x not in ("Sum", "sum", "Total", "total", "Summation", "summation")]
        FA_list.append("Sum")
        mean_list = FAin.iloc[:, 0].tolist()
        mean_list.append(0)
        FA_minmax_pd = pd.DataFrame({"Component": FA_list, "mean": mean_list})
        FA_minmax_pd = FA_minmax_pd.set_index(FA_minmax_pd.columns[0])
        FA_minmax_pd.loc["Sum", "mean"] = FA_minmax_pd["mean"][:-1].sum(axis=0)
        FA_minmax_pd["norm_mean"] = FA_minmax_pd["mean"] / FA_minmax_pd.loc["Sum", "mean"]
        FA_list.remove("Sum")
        FA_component_list = FA_list
        if FA_variation == True:
            max_cv = 0.25
            FA_minmax_pd["min"] = FA_minmax_pd["norm_mean"] * (1 - max_cv)  ##추가
            FA_minmax_pd["max"] = FA_minmax_pd["norm_mean"] * (1 + max_cv)
            for c_i, c in enumerate(FA_component_list):
                for min_or_max in ["min", "max"]:
                    new_col_name = c + "_{}".format(min_or_max)
                    FA_minmax_pd[new_col_name] = np.nan * len(FA_minmax_pd)
                    for row_i, fa_name in enumerate(FA_minmax_pd.index.values):
                        if fa_name == c:
                            FA_minmax_pd.loc[fa_name, new_col_name] = FA_minmax_pd.loc[c, min_or_max]
                        else:
                            FA_minmax_pd.loc[fa_name, new_col_name] = FA_minmax_pd.loc[fa_name, "norm_mean"] / \
                                                                      (FA_minmax_pd.loc["Sum", "norm_mean"] -
                                                                       FA_minmax_pd.loc[c, "norm_mean"]) * \
                                                                      (1 - FA_minmax_pd.loc[c, min_or_max])
                    FA_minmax_pd.loc["Sum", new_col_name] = FA_minmax_pd[new_col_name][:-1].sum(axis=0)

            min_max_col_name = [x for x in FA_minmax_pd.columns.values if any(m in x for m in ("_min", "_max"))]
            min_max_col_name.append("norm_mean")

            for c_iter in min_max_col_name:
                gpergFA_Series = FA_minmax_pd[c_iter][:-1]
                molpergFA_Series = gpergFA_Series / FAin["MW"].values
                molpergFA_sum = molpergFA_Series.sum()
                molpermolFA_Series = molpergFA_Series / molpergFA_sum
                FA_minmax_pd[c_iter + "mol%"] = molpermolFA_Series
                for fa in range(len(FA_component_list)):
                    coeff = molpermolFA_Series[fa]  # mol/molFA
                    species = FAin.iloc[fa, 2]
                    compart = FAin.iloc[fa, 3]
                    if fa == 0:
                        FAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                    else:
                        FAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                for facoa in range(len(FAcoain["Reactant"])) :
                    coeff = molpermolFA_Series[facoa]
                    species = FAcoain.iloc[facoa, 2]
                    compart = FAcoain.iloc[facoa, 3]
                    if facoa == 0:
                        FAcoasub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                    else:
                        FAcoasub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                FAsyn_c = FAsub + " -> 1.0 Rtotal[c]"
                FAsyn_e = FAsyn_c.replace("[c]","[e]")
                FAsyn_m = FAsyn_c.replace("[c]","[m]")
                FAsyn_x= FAsyn_c.replace("[c]","[x]")

                FAcoasyn_c = FAcoasub + " -> 1.0 Rtotalcoa[c]"
                FAcoasyn_m = FAcoasyn_c.replace("[c]","[m]")
                FAcoasyn_x = FAcoasyn_c.replace("[c]", "[x]")

                FA_dict[c_iter] = [FAsyn_c, FAsyn_e, FAsyn_m, FAsyn_x]
                FAcoa_dict[c_iter] = [FAcoasyn_c, FAcoasyn_m, FAcoasyn_x]

        else:
            gpergFA_Series = FA_minmax_pd["norm_mean"][:-1]
            molpergFA_Series = gpergFA_Series / FAin["MW"].values
            molpergFA_sum = molpergFA_Series.sum()
            molpermolFA_Series = molpergFA_Series / molpergFA_sum
            FA_minmax_pd["norm_mean_mol%"] = molpermolFA_Series
            for fa in range(len(FA_component_list)):
                coeff = molpermolFA_Series[fa]  # mol/molFA
                species = FAin.iloc[fa, 2]
                compart = FAin.iloc[fa, 3]
                if fa == 0:
                    FAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                else:
                    FAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            for facoa in range(len(FAcoain["Reactant"])) :
                coeff = molpermolFA_Series[facoa] 
                species = FAcoain.iloc[facoa, 2]
                compart = FAcoain.iloc[facoa, 3]
                if facoa == 0:
                    FAcoasub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                else:
                    FAcoasub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

            FAsyn_c = FAsub + " -> 1.0 Rtotal[c]"
            FAsyn_e = FAsyn_c.replace("[c]", "[e]")
            FAsyn_m = FAsyn_c.replace("[c]", "[m]")
            FAsyn_x = FAsyn_c.replace("[c]", "[x]")

            FAcoasyn_c = FAcoasub + " -> 1.0 Rtotalcoa[c]"
            FAcoasyn_m = FAcoasyn_c.replace("[c]", "[m]")
            FAcoasyn_x = FAcoasyn_c.replace("[c]", "[x]")

            FA_dict["ref"] = [FAsyn_c, FAsyn_e, FAsyn_m, FAsyn_x]
            FAcoa_dict["ref"]= [FAcoasyn_c, FAcoasyn_m, FAcoasyn_x]

        return LIPIDsyn, FA_dict, FAcoa_dict

    def BIOMASS(self):
        BIOMASSpd = pd.read_excel(self.file2read, sheet_name='Biomass', usecols="A:H", convert_float=True)
        BIOMASSin = self._processBiomassFrame(df=BIOMASSpd,in_or_out=-1)

        # out
        BIOMASSout = self._processBiomassFrame(df=BIOMASSpd,in_or_out=1)

        # PROTEIN BIOMASS
        BIOMASSsyn = self._returnSynthesisEquation(df_in=BIOMASSin,df_out=BIOMASSout)
        return BIOMASSsyn


#
# if __name__ == '__main__':
#     # file2read=testfile1
#     fileloc = ''
#     filename = fileloc + '\Cho test1.xlsx'
#
#     a = makeBiomass(file2read=filename)
#     b = a.Formulate_min_max()



