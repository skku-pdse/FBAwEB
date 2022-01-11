#-*-coding:utf-8-*-
import pandas as pd
from pandas import ExcelWriter
from datetime import date, datetime
import random
import numpy as np
import os.path as path



class makeBiomass():

    def __init__(self, file2read ):
        date = datetime.now().strftime("%b%d %H;%M")
        self.file2read=file2read
        saveDir="RandomBiomassExcel"
        Path(saveDir).mkdir(parents=True,exist_ok=True)
        self.file2save=os.path.join(saveDir,'Ecoli_macro_+-25% biomass_{}.xlsx'.format(date))

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


    def RandomBiomass(self,Protein, DNA, RNA, Lipid):
        pd.set_option("display.max_columns",999)
        all = pd.read_excel(self.file2read, sheet_name='Overall', usecols="A:W", convert_float=True,index_col=0)


        #USE average and stdev
        all["mean"]=all.mean(axis=1,skipna=True)
        all["stdev"]=all.iloc[:,:-1].std(axis=1,skipna=True)

        componentList = ["Protein","DNA(/Protein)", "RNA(/Protein)","Lipid","LPS", "Murein", "Inorganic", "Soluble pool","Sum"]

        # DNA/Protein or RNA/Protein -> For "Data Average","Norm_Data Average", "Data Stdev" columns : DNA g / Protein g, or RNA g/ Protein g.
        #                            -> For Random coefficient : DNA g/gDCW , RNA g/gDCW
        randomcoeffpd_default=pd.DataFrame({"Component":componentList,
                                    "Data Average":[np.nan]*len(componentList),"Data Stdev":[np.nan]*len(componentList),
                                    "Norm_Data Average" : [np.nan]*len(componentList),
                                    "Random Average": [np.nan]*len(componentList),"Random Stdev" : [np.nan]*len(componentList) })

        randomcoeffpd = randomcoeffpd_default.copy()
        randomcoeffpd_default.set_index(randomcoeffpd_default.columns[0])

        for n, m in enumerate(componentList[:-1]):
            if n==1 :
                m= "DNA/Protein"
            elif n==2 :
                m= "RNA/Protein"
            macro_mean_mass = all.loc[m, "mean"]
            macro_stdev = all.loc[m, "stdev"]
            randomcoeffpd.at[n,"Data Average"]=macro_mean_mass
            randomcoeffpd.at[n,"Data Stdev"] = macro_stdev

        #DNA RNA
        randomcoeffpd["Data Average"][1] = randomcoeffpd["Data Average"][0] * randomcoeffpd["Data Average"][1]
        randomcoeffpd["Data Average"][2] = randomcoeffpd["Data Average"][0] * randomcoeffpd["Data Average"][2]

        randomcoeffpd.at[len(componentList)-1,"Data Average"] = randomcoeffpd["Data Average"][len(componentList)-1].sum()

        varing_coeff = randomcoeffpd["Data Average"][:-3].sum(axis=0)
        fixed_coeff = randomcoeffpd[["Data Average"]][6:-1].sum(axis=0)
        #Normalization of average value
        for i in range(len(componentList)):
            if i < 6:
                randomcoeffpd["Norm_Data Average"][i] = randomcoeffpd["Data Average"][i] / varing_coeff * (
                            1 - fixed_coeff)
            else:
                randomcoeffpd["Norm_Data Average"][i] = randomcoeffpd["Data Average"][i]
        randomcoeffpd["Norm_Data Average"][-1] = randomcoeffpd["Norm_Data Average"][:-1].sum(axis=0)

        #biomass normalization factor
        norm_fac=randomcoeffpd["Norm_Data Average"][len(componentList)-1] / randomcoeffpd["Data Average"][len(componentList)-1]


        if Protein == 'max' :
            PROTcoe=randomcoeffpd["Norm_Data Average"][0] * 1.25
        elif Protein == 'min' :
            PROTcoe = randomcoeffpd["Norm_Data Average"][0] * 0.75
        elif Protein == "mean" :
            PROTcoe = randomcoeffpd["Norm_Data Average"][0]

        if DNA == 'max' :
            DNAcoe = randomcoeffpd["Norm_Data Average"][1] * 1.25
        elif DNA == 'min' :
            DNAcoe = randomcoeffpd["Norm_Data Average"][1] * 0.75
        elif DNA == 'mean' :
            DNAcoe = randomcoeffpd["Norm_Data Average"][1]

        if RNA == 'max':
            RNAcoe = randomcoeffpd["Norm_Data Average"][2] * 1.25
        elif RNA == 'min':
            RNAcoe = randomcoeffpd["Norm_Data Average"][2] * 0.75
        elif RNA == 'mean':
            RNAcoe = randomcoeffpd["Norm_Data Average"][2]

        # GLYCOcoe = random.uniform(randomcoeffpd["Norm_Data Average"][3] - randomcoeffpd["Data Stdev"][3],randomcoeffpd["Norm_Data Average"][3]+randomcoeffpd["Data Stdev"][3])

        if Lipid  == 'max':
            LIPIDcoe = randomcoeffpd["Norm_Data Average"][3] *1.25
        elif Lipid == 'min':
            LIPIDcoe = randomcoeffpd["Norm_Data Average"][3] * 0.75
        elif Lipid == 'mean':
            LIPIDcoe = randomcoeffpd["Norm_Data Average"][3]

        LPScoe= randomcoeffpd["Norm_Data Average"][4]
        MUREINcoe= randomcoeffpd["Norm_Data Average"][5]
        INORGANICcoe= randomcoeffpd["Norm_Data Average"][6]
        SOLUBLEPcoe= randomcoeffpd["Norm_Data Average"][7]

        if (Protein == 'max') or (Protein == 'min') :
            varing_coe = DNAcoe + RNAcoe + LIPIDcoe
            fixed_coe= PROTcoe + LPScoe+MUREINcoe+INORGANICcoe+SOLUBLEPcoe
            PROTcoe_nor = PROTcoe
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)
            varing_coe_nor= DNAcoe_nor + RNAcoe_nor + LIPIDcoe_nor
        elif (DNA == 'max') or (DNA == 'min') :
            varing_coe = PROTcoe + RNAcoe + LIPIDcoe
            fixed_coe= DNAcoe + LPScoe+MUREINcoe+INORGANICcoe+SOLUBLEPcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)
            varing_coe_nor = PROTcoe_nor + RNAcoe_nor + LIPIDcoe_nor
        elif (RNA == 'max') or (RNA == 'min') :
            varing_coe = PROTcoe + DNAcoe +  LIPIDcoe
            fixed_coe=RNAcoe +LPScoe+MUREINcoe+INORGANICcoe+SOLUBLEPcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + LIPIDcoe_nor
        elif (Lipid == 'max') or (Lipid == 'min') :
            varing_coe = PROTcoe + DNAcoe + RNAcoe
            fixed_coe =  LIPIDcoe + LPScoe + MUREINcoe + INORGANICcoe + SOLUBLEPcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor

        else : # all mean
            varing_coe = PROTcoe + DNAcoe + RNAcoe + LIPIDcoe
            fixed_coe =   LPScoe + MUREINcoe + INORGANICcoe + SOLUBLEPcoe
            PROTcoe_nor = PROTcoe / varing_coe * (1 - fixed_coe)
            DNAcoe_nor = DNAcoe / varing_coe * (1 - fixed_coe)
            RNAcoe_nor = RNAcoe / varing_coe * (1 - fixed_coe)
            LIPIDcoe_nor = LIPIDcoe / varing_coe * (1 - fixed_coe)
            varing_coe_nor = PROTcoe_nor + DNAcoe_nor + RNAcoe_nor + LIPIDcoe_nor

        # totalwt_nor=PROTcoe_nor+DNAcoe_nor+RNAcoe_nor+GLYCOcoe_nor+LIPIDcoe_nor + fixed_coe
        totalwt_nor = varing_coe_nor + fixed_coe

        # randomcoeffpd[str(col_k)]=[PROTcoe_nor,DNAcoe_nor,RNAcoe_nor,GLYCOcoe_nor,LIPIDcoe_nor, LPScoe, MUREINcoe, INORGANICcoe,SOLUBLEPcoe,totalwt_nor]
        Final_coe_list = [PROTcoe_nor, DNAcoe_nor, RNAcoe_nor, LIPIDcoe_nor, LPScoe, MUREINcoe, INORGANICcoe, SOLUBLEPcoe, totalwt_nor, norm_fac]

        return Final_coe_list
        # # Random result summary
        # print(randomcoeffpd)
        # for r in range(len(randomcoeffpd)):
        #     randomcoeffpd.at[r, "Random Average"] = randomcoeffpd.iloc[r, len(randomcoeffpd_default.columns):].mean( skipna=True)
        #     randomcoeffpd.at[r, "Random Stdev"] = randomcoeffpd.iloc[r, len(randomcoeffpd_default.columns):].std(skipna=True)
        #
        # randomcoeffpd.at[len(componentList)-1, "Random Average"] = randomcoeffpd["Random Average"][:8].sum(axis=0)

            # multiple overall composition loop
    def Formulate_min_max(self):

        minmax_list = ['Prot_min', "Prot_max", "DNA_min", "DNA_max", "RNA_min", "RNA_max", "Lipid_min", "Lipid,max", "Ref"]
        coeff_index = ["PROTcoe_nor", "DNAcoe_nor", "RNAcoe_nor", "LIPIDcoe_nor", "LPScoe", "MUREINcoe", "INORGANICcoe",
                       "SOLUBLEPcoe","totalwt_nor", "norm_fac"]
        coeff_pd = pd.DataFrame(index=coeff_index, columns=minmax_list)
        rxn_list=["ProtSyn","DNASyn","RNASyn","LipidSyn","BiomassSyn"]
        minmax_pd=pd.DataFrame(index=rxn_list, columns=minmax_list)


        for col_i,head in enumerate(minmax_list):
            if col_i == 0 :
                coeff_list = self.RandomBiomass("min", "mean", "mean", "mean")
            elif col_i == 1 :
                coeff_list = self.RandomBiomass("max", "mean", "mean", "mean")
            elif col_i == 2 :
                coeff_list = self.RandomBiomass("mean", "min", "mean", "mean")
            elif col_i == 3 :
                coeff_list = self.RandomBiomass("mean", "max", "mean", "mean")
            elif col_i == 4 :
                coeff_list = self.RandomBiomass("mean", "mean", "min", "mean")
            elif col_i == 5 :
                coeff_list = self.RandomBiomass("mean", "mean", "max", "mean")
            elif col_i == 6 :
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "min")
            elif col_i == 7 :
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "max")
            elif col_i == 8:
                coeff_list = self.RandomBiomass("mean", "mean", "mean", "mean")

            minmax_pd[head][0]= self.PROTEIN(coeff_list[0])[0]
            minmax_pd[head][1] = self.DNA_or_RNA(coeff_list[1],DNA_or_RNA="DNA")
            minmax_pd[head][2] = self.DNA_or_RNA(coeff_list[2],DNA_or_RNA="RNA")
            minmax_pd[head][3] = self. LIPID(coeff_list[3])
            minmax_pd[head][4] = self. BIOMASS(coeff_list[4], coeff_list[5])
            for k in range(len(coeff_index)):
                coeff_pd[head][k] = coeff_list[k]

        with ExcelWriter(self.file2save) as writer:
            coeff_pd.to_excel(writer,sheet_name="Random coefficient")
            self.PROTEIN(coeff_list[0])[1].to_excel(writer, sheet_name="AA_composition")
            minmax_pd.to_excel(writer, sheet_name="25% minmax biomass")
            writer.save()

    def PROTEIN(self,Prot_wt_per_DCW):

        PROTpd= pd.read_excel(self.file2read, sheet_name='PROTsyn', usecols="A:H" ,convert_float=True)
        PROTin=self.processBiomassFrame(df=PROTpd, in_or_out=-1)

        AApd = pd.read_excel(self.file2read, sheet_name="PROTsyn", usecols="J:W",convert_float=True)
        AApd["mean"] = AApd.mean(axis=1)
        AApd.at[20,"mean"]=AApd["mean"][:20].sum(axis=0)
        AApd["stdev"] = AApd.iloc[:,:-1].std(axis=1)
        #################

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
            PROTout=self.processBiomassFrame(PROTpd,1)

            # PROTout coeff cal
            # H2O
            PROTout.at[0,"mmol/gDCW.1"]=PROTin.iloc[:20,4].sum(axis=0)

            for j in range(len(PROTout)):
                coeff= PROTout.iloc[j,2]
                species= PROTout.iloc[j,0]
                compart= PROTout.iloc[j,1]
                if j==0 :
                    PROTpro= str(round(coeff,6)) + " "+ species + "["+ compart +"]"

                else:
                    PROTpro += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"

            #PROTEIN BIOMASS
            PROTsyn=PROTsub + " -> " + PROTpro

        return PROTsyn, AApd

    def DNA_or_RNA(self,wt_per_DCW, DNA_or_RNA):

        DNApd = pd.read_excel(self.file2read, sheet_name=DNA_or_RNA+"syn", usecols="A:H", convert_float=True)
        DNAin = self.processBiomassFrame(df=DNApd, in_or_out=-1)

        # PROT in
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


        # DNA out
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

        # PROTEIN BIOMASS
        DNAsyn = DNAsub + " -> " + DNApro
        return DNAsyn

    # def CARB(self,Carb_per_DCW):
    #     CARBpd = pd.read_excel(self.filename, sheet_name='CARBsyn', usecols="A:H", convert_float=True)
    #     CARBin = CARBpd.iloc[0, :5]
    #     gpergCARB = CARBin.iloc[0]
    #
    #     MW = CARBin.iloc[1]
    #     CARBin.at["mmol/gDCW"] = Carb_per_DCW * gpergCARB / MW * 1000
    #
    #     coeff = CARBin.iloc[4]
    #     species = CARBin.iloc[2]
    #     compart = CARBin.iloc[3]
    #
    #     CARBsub = str(round(coeff,6)) + " " + species + "[" + compart + "]"
    #
    #     # PROTEIN BIOMASS
    #     CARBsyn = CARBsub + " -> 1.0 carbohydrate[c]"
    #     return CARBsyn

    def LIPID(self,LIPID_wt_per_DCW):
        LIPIDpd= pd.read_excel(self.file2read, sheet_name='LIPIDsyn', usecols="A:H" ,convert_float=True)

        #Lipid in
        LIPIDin=self.processBiomassFrame(df=LIPIDpd,in_or_out=-1)

        # for i in range(len(LIPIDin)):
        #     gLIPIDpergDCW=LIPIDin.iloc[i,0]
        #     MW=LIPIDin.iloc[i,1]
        #     LIPIDin.at[i,"mmol/gDCW"]= LIPID_wt_per_DCW*gLIPIDpergDCW/MW*1000
        #
        #     coeff= LIPIDin.iloc[i,4]
        #     species= LIPIDin.iloc[i,2]
        #     compart= LIPIDin.iloc[i,3]
        #
        #     if i==0 :
        #         LIPIDsub= str(round(coeff,6)) + " "+ species + "["+compart+"]"
        #     else :
        #         LIPIDsub += " + " + str(round(coeff,6)) + " " + species +"["+compart+"]"


# FATTY ACID#####################################
        mono_norm_fac = LIPIDin["g/g Lipid"].sum()
        Lipid_mono_ratio_df = pd.DataFrame({"pe": np.nan * 3, "pg": np.nan * 3, "clpn": np.nan * 3},
                                           index=["160", "161", "181"])

        for l_i, l in enumerate(LIPIDin["Reactant"]):
            if "pe" in l:
                col_loc = 0
            elif "pg" in l:
                col_loc = 1
            elif "clpn" in l:
                col_loc = 2

            if "160" in l:
                row_loc = 0
            elif "161" in l:
                row_loc = 1
            elif "181" in l:
                row_loc = 2
            Lipid_mono_ratio_df.iloc[row_loc, col_loc] = LIPIDin["g/g Lipid"][l_i] / mono_norm_fac
        print(Lipid_mono_ratio_df)

        Lipid_mono_ratio_df["FA_Sum"] = Lipid_mono_ratio_df.sum(axis=1)
        Lipid_mono_ratio_df.loc["Sum"] = Lipid_mono_ratio_df.sum(axis=0)
        print(Lipid_mono_ratio_df)
        Fatty_acid_ratio = pd.DataFrame({"ref": list(Lipid_mono_ratio_df["FA_Sum"][:-1])}, index=["160", "161", "181"])
        print(Fatty_acid_ratio)
        Lipid_kind_ratio = pd.DataFrame(
            {"pe": [Lipid_mono_ratio_df.loc["Sum", "pe"]], "pg": [Lipid_mono_ratio_df.loc["Sum", "pg"]],
             "clpn": [Lipid_mono_ratio_df.loc["Sum", "clpn"]]},index=["ref"])
        print(Lipid_kind_ratio)

# Ref relative mass recal
#         LIPIDDict = dict()

        for l_i, l_kind in enumerate(["pe", "pg", "clpn"]):
            Lipid_mono_ratio_df.iloc[0, l_i] = Fatty_acid_ratio.iloc[0,0] * Lipid_kind_ratio[l_kind][0]
            Lipid_mono_ratio_df.iloc[1, l_i] = Fatty_acid_ratio.iloc[1,0] * Lipid_kind_ratio[l_kind][0]
            Lipid_mono_ratio_df.iloc[2, l_i] = Fatty_acid_ratio.iloc[2,0] * Lipid_kind_ratio[l_kind][0]

# Option 2 (No Use): No fatty acid variation, other lipid monomer variation
                # Option 2 : No fatty acid variation, other lipid monomer variation
        # for lip in ["pe", "pg", "clpn"]:
        #     for min_or_max in ["min", "max"]:
        #         if min_or_max == "min":
        #             fa_min_max_val = Lipid_kind_ratio.loc["ref", lip] * 0.75
        #         elif min_or_max == "max":
        #             fa_min_max_val = Lipid_kind_ratio.loc["ref", lip] * 1.25
        #
        #         Lipid_kind_ratio.loc[lip + "_" + min_or_max] = np.nan * len(Lipid_kind_ratio.columns)
        #         for r in range(len(Lipid_kind_ratio.columns)):
        #             Lipid_kind_ratio.loc[Lipid_kind_ratio, r] = Lipid_kind_ratio["ref"][r] / (
        #                     1 - Lipid_kind_ratio.loc["ref", lip]) * (1 - fa_min_max_val)
        #         Lipid_kind_ratio.loc[lip + "_" + min_or_max, lip] = fa_min_max_val
        #
        # LIPIDDict = dict()
        # for lip_mm_i in range(len(Lipid_kind_ratio.index)):
        #     for l_i, l_kind in enumerate(["pe", "pg", "clpn"]):
        #         Lipid_mono_ratio_df.iloc[0, l_i] = Fatty_acid_ratio.iloc[0,0] * Lipid_kind_ratio[l_kind][lip_mm_i]
        #         Lipid_mono_ratio_df.iloc[1, l_i] = Fatty_acid_ratio.iloc[1,0] * Lipid_kind_ratio[l_kind][lip_mm_i]
        #         Lipid_mono_ratio_df.iloc[2, l_i] = Fatty_acid_ratio.iloc[2,0] * Lipid_kind_ratio[l_kind][lip_mm_i]
##################################
# Should be run from here when Option 1 or 2 are activated
        lipid_coeff_dict = dict()
        for c in ["pe", "pg", "clpn"]:
            for f in ["160", "161", "181"]:
                lipid_mono = c + f
                monomer_g = Lipid_mono_ratio_df.loc[f, c]
                lipid_coeff_dict[lipid_mono] = monomer_g
        print("this should be 1", sum(list(lipid_coeff_dict.values())))

        for i in range(len(LIPIDin)):
            lipid_comp = LIPIDin["Reactant"][i]
            gpergP = lipid_coeff_dict[lipid_comp]
            MW = LIPIDin.iloc[i, 1]
            LIPIDin.at[i, "mmol/gDCW"] = gpergP / MW * 1000 * LIPID_wt_per_DCW
            coeff = LIPIDin.iloc[i, 4]
            species = LIPIDin.iloc[i, 2]
            compart = LIPIDin.iloc[i, 3]
            if i == 0:
                LIPIDsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                LIPIDsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"
        # LIPID out
        LIPIDout=self.processBiomassFrame(df=LIPIDpd,in_or_out=1)
        for p_ind, p in enumerate(LIPIDout["Product"]):
            coeff = LIPIDout.iloc[p_ind, 2]
            species = LIPIDout.iloc[p_ind, 0]
            compart = LIPIDout.iloc[p_ind, 1]
            if p_ind == 0:
                LIPIDpro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                LIPIDpro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"
        # LIPID BIOMASS
        LIPIDsyn = LIPIDsub + " -> " + LIPIDpro

        return LIPIDsyn


    def BIOMASS(self,LPS_per_DCW,Murein_per_DCW):
        BIOMASSpd = pd.read_excel(self.file2read, sheet_name='Biomass', usecols="A:H", convert_float=True)
        BIOMASSin = self.processBiomassFrame(df=BIOMASSpd,in_or_out=-1)


        #LPS
        LPSMW=BIOMASSin.iloc[0,1]
        BIOMASSin.at[0,"mmol/gDCW"]= LPS_per_DCW / LPSMW * 1000

        #Murein
        MureinMW=BIOMASSin.iloc[1,1]
        BIOMASSin.at[1,"mmol/gDCW"]= Murein_per_DCW / MureinMW *1000

        for i in range(len(BIOMASSin)):
            #No norm_fac in case of E.coli
            coeff = BIOMASSin.iloc[i, 4]
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

    # file2read= "&&& for python.xlsx"
    a=makeBiomass(file2read=file2read)
    b=a.Formulate_min_max()



