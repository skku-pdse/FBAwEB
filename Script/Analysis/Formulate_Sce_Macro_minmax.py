#-*-coding:utf-8-*-
import pandas as pd
from pandas import ExcelWriter
from datetime import date, datetime
import random
import numpy as np
import os.path as path


class makeBiomass():

    def __init__(self, file2read):
        date = datetime.now().strftime("%b%d %H;%M")
        self.file2read=file2read
        self.file2save='Yeast_macro_+-25% biomass_{}.xlsx'.format(date)

        # if path.isfile(file2save) == True :
        #     self.file2save=file2save + "(1)"

    def processBiomassFrame(self, df, in_or_out):
        self.df = df
        self.in_or_out = in_or_out

        if in_or_out == 0:
            if len(self.df["Reactant"]) >= len(self.df["Product"]):
                df_processed = self.df[:len(df["Reactant"])]
            else:
                df_processed = self.df[:len(df["Product"])]

        elif in_or_out == -1:
            df_processed = self.df.iloc[:, :list(self.df.columns).index("Product")]
            for i in range(len(self.df)):
                if pd.isna(self.df["Reactant"][i]) == True:
                    df_processed = df_processed[:i]
                    break

        elif in_or_out == 1:
            df_processed = self.df.iloc[:, list(self.df.columns).index("Product"):]
            for i in range(len(self.df)):
                if pd.isna(self.df["Product"][i]) == True:
                    df_processed = df_processed[:i]
                    break
        return df_processed


    def RandomBiomass(self, Protein, DNA, RNA, Carbohydrate, Lipid):
        pd.set_option("display.max_columns",999)

        #Overall composition
        all = pd.read_excel(self.file2read, sheet_name='Overall', usecols="A:T", convert_float=True,index_col=0)
        print(all)

        #USE average and stdev
        all["mean"]=all.mean(axis=1)
        all["stdev"]=all.iloc[:,:-1].std(axis=1)


        componentList=["Protein","DNA","RNA","Carbohydrate","Lipid","Metals + sulfate","Sum"]
        randomcoeffpd_default=pd.DataFrame({"Component":componentList,"Data Average":[np.nan]*len(componentList),"Data Stdev":[np.nan]*len(componentList),
                                    "Norm_Data Average": [np.nan]*len(componentList),
                                    "Random Average":[np.nan]*len(componentList), "Random Stdev": [np.nan]*len(componentList)})

        randomcoeffpd = randomcoeffpd_default.copy()
        randomcoeffpd_default.set_index(randomcoeffpd_default.columns[0])
        for n,m in enumerate(componentList[:-1]):
            macro_mean_mass=all.loc[m,"mean"]
            macro_stdev=all.loc[m,"stdev"]
            randomcoeffpd.at[n,"Data Average"]=macro_mean_mass
            randomcoeffpd.at[n,"Data Stdev"] = macro_stdev

        randomcoeffpd.at[len(componentList)-1, "Data Average"]= randomcoeffpd["Data Average"][:len(componentList)-1].sum(axis=0)

        # Normalization of average value
        for i in range(len(componentList)):
            randomcoeffpd.at[i, "Norm_Data Average"] = randomcoeffpd["Data Average"][i] / randomcoeffpd["Data Average"][len(componentList)-1]

        randomcoeffpd.at[len(componentList)-1, "Norm_Data Average"] =randomcoeffpd["Norm_Data Average"][:len(componentList)-1].sum(axis=0)



        #biomass normalization factor
        norm_fac =  randomcoeffpd["Norm_Data Average"][5] / 0.029

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

            # GLYCOcoe = random.uniform(randomcoeffpd["Norm_Data Average"][3] - randomcoeffpd["Data Stdev"][3],randomcoeffpd["Norm_Data Average"][3]+randomcoeffpd["Data Stdev"][3])
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

        Metalcoe= randomcoeffpd["Norm_Data Average"][5]

        if (Protein == 'max') or (Protein == 'min'):
            varing_coe = DNAcoe + RNAcoe + CARBcoe + LIPIDcoe
            fixed_coe = PROTcoe + Metalcoe
            PROTcoe_nor = PROTcoe
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)

            varing_coe_nor = DNAcoe_nor + RNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor
        elif (DNA == 'max') or (DNA == 'min'):
            varing_coe = PROTcoe + RNAcoe + CARBcoe+ LIPIDcoe
            fixed_coe = DNAcoe+ Metalcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)

            varing_coe_nor = PROTcoe_nor + RNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor
        elif (RNA == 'max') or (RNA == 'min'):
            varing_coe = PROTcoe + DNAcoe + CARBcoe+ LIPIDcoe
            fixed_coe = RNAcoe+ Metalcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)

            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor

        elif (Carbohydrate == 'max') or (Carbohydrate == 'min'):
            varing_coe = PROTcoe + DNAcoe + RNAcoe + LIPIDcoe
            fixed_coe = CARBcoe+ Metalcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor + LIPIDcoe_nor


        elif (Lipid == 'max') or (Lipid == 'min'):
            varing_coe = PROTcoe + DNAcoe + RNAcoe + CARBcoe
            fixed_coe = LIPIDcoe+ Metalcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            CARBcoe_nor = CARBcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor + CARBcoe_nor

        else : # all mean
            fixed_coe =  Metalcoe
            PROTcoe_nor = PROTcoe
            DNAcoe_nor = DNAcoe
            RNAcoe_nor = RNAcoe
            CARBcoe_nor = CARBcoe
            LIPIDcoe_nor = LIPIDcoe
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor + CARBcoe_nor + LIPIDcoe_nor


        # totalwt_nor=PROTcoe_nor+DNAcoe_nor+RNAcoe_nor+GLYCOcoe_nor+LIPIDcoe_nor + fixed_coe
        totalwt_nor = varing_coe_nor + fixed_coe

        # randomcoeffpd[str(col_k)]=[PROTcoe_nor,DNAcoe_nor,RNAcoe_nor,GLYCOcoe_nor,LIPIDcoe_nor, LPScoe, MUREINcoe, INORGANICcoe,SOLUBLEPcoe,totalwt_nor]
        Final_coe_list = [PROTcoe_nor, DNAcoe_nor, RNAcoe_nor, CARBcoe_nor, LIPIDcoe_nor, Metalcoe, totalwt_nor, norm_fac]

        return Final_coe_list

    def Formulate_min_max(self):

        minmax_list = ['Prot_min', "Prot_max", "DNA_min", "DNA_max", "RNA_min", "RNA_max", "Carb_min", "Carb_max",
                       "Lipid_min", "Lipid,max", "Ref"]
        coeff_index = ["PROTcoe_nor", "DNAcoe_nor", "RNAcoe_nor", "CARBcoe_nor", "LIPIDcoe_nor", "Metal_coe", "totalwt_nor"]
        coeff_pd = pd.DataFrame(index=coeff_index, columns=minmax_list)
        rxn_list = ["ProtSyn", "DNASyn", "RNASyn", "CarbSyn", "LipidSyn", "FATTYACIDSyn","BiomassSyn"]
        minmax_pd = pd.DataFrame(index=rxn_list, columns=minmax_list)

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

            minmax_pd[head][0] = self.PROTEIN(coeff_list[0])[0]
            minmax_pd[head][1] = self.DNA_or_RNA(coeff_list[1], DNA_or_RNA="DNA")
            minmax_pd[head][2] = self.DNA_or_RNA(coeff_list[2], DNA_or_RNA="RNA")
            minmax_pd[head][3] = self.CARB(coeff_list[3])
            minmax_pd[head][4] = self.LIPID(coeff_list[4],FA_variation=False)[0]
            minmax_pd[head][5] = self.LIPID(coeff_list[4], FA_variation=False)[1]["ref"]
            minmax_pd[head][6] = self.BIOMASS(coeff_list[7])

            for k in range(len(coeff_index)):
                coeff_pd[head][k] = coeff_list[k]

        with ExcelWriter(self.file2save) as writer:
            coeff_pd.to_excel(writer, sheet_name="Random coefficient")
            self.PROTEIN(coeff_list[0])[1].to_excel(writer, sheet_name="AA_composition")
            minmax_pd.to_excel(writer, sheet_name="25% minmax biomass")
            writer.save()


    def PROTEIN(self,Prot_wt_per_DCW):

        PROTpd= pd.read_excel(self.file2read, sheet_name='PROTsyn', usecols="A:H" ,convert_float=True)
        PROTin=self.processBiomassFrame(df=PROTpd,in_or_out=-1)

        AApd = pd.read_excel(self.file2read, sheet_name="PROTsyn", usecols="J:U",convert_float=True)
        AApd["mean"] = AApd.mean(axis=1)
        AApd.at[20,"mean"]=AApd["mean"][:20].sum(axis=0)

        AApd["stdev"] = AApd.iloc[:,:-1].std(axis=1)
        ##################

        #####################

        AAmean_sum=AApd["mean"][20]
        AApd["mean_normalized"]=AApd["mean"]/AAmean_sum
        AApd.at[20,"mean_normalized"]=AApd["mean_normalized"][:20].sum(axis=0)
        for i in range(len(PROTin)):
            gpergP=AApd["mean_normalized"][i]
            MW=PROTin.iloc[i,1]
            PROTin.at[i,"mmol/gDCW"]= Prot_wt_per_DCW*gpergP/MW*1000

            coeff= PROTin.iloc[i,4]
            species= PROTin.iloc[i,2]
            compart= PROTin.iloc[i,3]

            if i==0 :
                PROTsub= str(round(coeff,6)) + " "+ species + "["+compart+"]"
            else :
                PROTsub += " + " + str(round(coeff,6)) + " " + species +"["+compart+"]"

            #PROT out
            PROTout=self.processBiomassFrame(df=PROTpd, in_or_out=1)

            # PROTout coeff cal
            for k in range(len(PROTin)):
                PROTout.at[k, "mmol/gDCW.1"] = PROTin.iloc[k, 4]
            # H2O
            PROTout.at[20, "mmol/gDCW.1"] = PROTout.iloc[:20, 2].sum(axis=0)

            for j in range(len(PROTout)):
                coeff = PROTout.iloc[j, 2]
                species = PROTout.iloc[j, 0]
                compart = PROTout.iloc[j, 1]
                if j == 0:
                    PROTpro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                else:
                    PROTpro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

            #PROTEIN BIOMASS
            PROTsyn=PROTsub + " -> " + PROTpro

        return PROTsyn, AApd

    def DNA_or_RNA(self,wt_per_DCW, DNA_or_RNA):

        DNApd = pd.read_excel(self.file2read, sheet_name=DNA_or_RNA+"syn", usecols="A:H", convert_float=True)
        DNAin = self.processBiomassFrame(df=DNApd,in_or_out=-1)

        # Nucleotide in
        for i in range(len(DNAin)):
            gDNApergDCW = DNAin.iloc[i, 0]
            MW = DNAin.iloc[i, 1]
            DNAin.at[i, "mmol/gDCW"] = wt_per_DCW * gDNApergDCW / MW * 1000

            coeff = DNAin.iloc[i, 4]
            species = DNAin.iloc[i, 2]
            compart = DNAin.iloc[i, 3]

            if i==0 :
                DNAsub = str(round(coeff,6)) + " " + species + "[" + compart + "]"
            else:
                DNAsub += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"


        # Nucleotide out
        DNAout = self.processBiomassFrame(df=DNApd,in_or_out=1)
        # ppi
        for p_ind,p in enumerate(DNAout["Product"]):
            if DNAout["Product"][p_ind].lower()=="ppi" :
                DNAout.at[p_ind, "mmol/gDCW.1"] = DNAin.iloc[:,4].sum(axis=0)
        for j in range(len(DNAout)):
            coeff = DNAout.iloc[j, 2]
            species = DNAout.iloc[j, 0]
            compart = DNAout.iloc[j, 1]
            if j==0 :
                DNApro = str(round(coeff,6)) + " " + species + "[" + compart + "]"
            else:
                DNApro += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"
        #  Nucleotide  BIOMASS
        DNAsyn = DNAsub + " -> " + DNApro
        return DNAsyn

    def CARB(self,Carb_per_DCW):
        CARBpd = pd.read_excel(self.file2read, sheet_name='CARBsyn', usecols="A:H", convert_float=True)
        CARBin = self.processBiomassFrame(df=CARBpd,in_or_out=-1)

        for i in range(len(CARBin)):
            gpergCARB = CARBin.iloc[i, 0]
            MW = CARBin.iloc[i, 1]
            CARBin.at[i, "mmol/gDCW"] = Carb_per_DCW * gpergCARB / MW * 1000

            coeff = CARBin.iloc[i, 4]
            species = CARBin.iloc[i, 2]
            compart = CARBin.iloc[i, 3]
            if i == 0:
                CARBsub = str(round(coeff,6)) + " " + species + "[" + compart + "]"
            else:
                CARBsub += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"

        CARBout=self.processBiomassFrame(df=CARBpd,in_or_out=1)
        for p_ind, p in enumerate(CARBout["Product"]):
            if any(c == p.lower() for c in ("carbohydrate", "carb", "carbo")):
                CARBout.loc[p, "mmol/gDCW.1"] = 1

            coeff = CARBout.iloc[p_ind, 2]
            species = CARBout.iloc[p_ind, 0]
            compart = CARBout.iloc[p_ind, 1]
            if p_ind == 0:
                CARBpro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                CARBpro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"
        # PROTEIN BIOMASS
        CARBsyn = CARBsub + " -> " + CARBpro
        return CARBsyn

    def LIPID(self,LIPID_wt_per_DCW,FA_variation):
        LIPIDpd= pd.read_excel(self.file2read, sheet_name='LIPIDsyn', usecols="A:H" ,convert_float=True)

        #Lipid in
        LIPIDin=self.processBiomassFrame(df=LIPIDpd,in_or_out=-1)
        AverageLipidMW=LIPIDin.iloc[0,1]
        for i in range(len(LIPIDin)):
            molpermolLipid=LIPIDin.iloc[i,0] # when g lipid/g DCW=0.042788756

            LIPIDin.at[i,"mmol/gDCW"]= molpermolLipid/AverageLipidMW*LIPID_wt_per_DCW*1000

            coeff= LIPIDin.iloc[i,4]
            species= LIPIDin.iloc[i,2]
            compart= LIPIDin.iloc[i,3]

            if i==0 :
                LIPIDsub= str(round(coeff,6)) + " "+ species + "["+compart+"]"
            else :
                LIPIDsub += " + " + str(round(coeff,6)) + " " + species +"["+compart+"]"

        #Lipid out
        LIPIDout=self.processBiomassFrame(LIPIDpd, 1)
        for p_ind,p in enumerate(LIPIDout["Product"]):
            coeff = LIPIDout.iloc[p_ind, 2]
            species = LIPIDout.iloc[p_ind, 0]
            compart = LIPIDout.iloc[p_ind, 1]
            if p_ind==0 :
                LIPIDpro = str(round(coeff,6)) + " " + species + "[" + compart + "]"
            else:
                LIPIDpro += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"
        # LIPID BIOMASS
        LIPIDsyn = LIPIDsub + " -> " + LIPIDpro

        FA_dict = dict()
        FApd = pd.read_excel(self.file2read, sheet_name='FattyAcidSyn', usecols="A:H", convert_float=True)
        FAin = self.processBiomassFrame(FApd, -1)
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
                print(molpergFA_sum)
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

                FAsyn = FAsub + " -> 1.0 fatty_acid[c]"

                FA_dict[c_iter] = FAsyn

        else:
            gpergFA_Series = FA_minmax_pd["norm_mean"][:-1]
            molpergFA_Series = gpergFA_Series / FAin["MW"].values
            molpergFA_sum = molpergFA_Series.sum()
            print(molpergFA_sum)
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

            FAsyn = FAsub + " -> 1.0 fatty_acid[c]"
            FA_dict["ref"] = FAsyn

        return LIPIDsyn, FA_dict


    def BIOMASS(self,norm_fac):
        BIOMASSpd = pd.read_excel(self.file2read, sheet_name='Biomass', usecols="A:H", convert_float=True)
        BIOMASSin = self.processBiomassFrame(df=BIOMASSpd,in_or_out=-1)

        for i in range(len(BIOMASSin)):
            if (BIOMASSin["Reactant"][i].lower() in ("atp","h2o","protein","dna","rna","carbohydrate","lipid")):
                coeff = BIOMASSin.iloc[i, 4]
            else:
                coeff = BIOMASSin.iloc[i, 4] * norm_fac
            species = BIOMASSin.iloc[i, 2]
            compart = BIOMASSin.iloc[i, 3]
            if i == 0:
                BIOMASSsub = str(round(coeff,6)) + " " + species + "[" + compart + "]"
            else:
                BIOMASSsub += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"
        # out
        BIOMASSout = self.processBiomassFrame(df=BIOMASSpd,in_or_out=1)
        for j in range(len(BIOMASSout)):
            coeff = BIOMASSout.iloc[j, 2]
            species = BIOMASSout.iloc[j, 0]
            compart = BIOMASSout.iloc[j, 1]
            if j == 0:
                BIOMASSpro = str(round(coeff,6)) + " " + species + "[" + compart + "]"

            else:
                BIOMASSpro += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"

        # PROTEIN BIOMASS
        BIOMASSsyn = BIOMASSsub + " -> " + BIOMASSpro
        return BIOMASSsyn


if __name__ == '__main__':
    
    #Directory = "&&& for python.xlsx"
    #file2read= 
    a=makeBiomass(file2read=file2read)
    b=a.Formulate_min_max()


