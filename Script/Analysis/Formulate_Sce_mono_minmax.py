#-*-coding:utf-8-*-
import pandas as pd
from pandas import ExcelWriter
from datetime import date, datetime
import random
import numpy as np
import os.path as path
import operator


class makeBiomass():

    def __init__(self, filename):
        date = datetime.now().strftime("%b%d %H;%M")
        self.filename = filename
        self.file2save='Yeast_mono_+-25% biomass_{}.xlsx'.format(date)
        if path.isfile(self.file2save) == True :
            self.file2save=file2save.format(date + "(1)")


    def RandomBiomass(self):
        pd.set_option("display.max_columns",999)
        CARBlist=list()
        BIOMASSlist=list()

        #Overall composition
        all = pd.read_excel(self.filename, sheet_name='Overall',usecols="A:T",convert_float=True)
        print(all)

        #USE average and stdev
        all["mean"]=all.mean(axis=1)
        all["stdev"]=all.iloc[:,:-1].std(axis=1)

        componentList = ["Protein", "DNA", "RNA", "Carbohydrate", "Lipid", "Metals + sulfate", "Sum"]

        randomcoeffpd=pd.DataFrame({"Component":componentList,"Data Average":list(all["mean"]),"Data Stdev":list(all["stdev"]),
                                    "Norm_Data Average": [np.nan] * 7,
                                    "Random Average": [np.nan] * 7, "Random Stdev": [np.nan] * 7})

        randomcoeffpd.at[6, "Data Average"]= randomcoeffpd["Data Average"][:6].sum(axis=0)

        # Normalization of average value
        for i in range(len(randomcoeffpd)):
            randomcoeffpd.at[i, "Norm_Data Average"] = randomcoeffpd["Data Average"][i] / randomcoeffpd.at[6, "Data Average"]

        randomcoeffpd.at[6, "Norm_Data Average"] =randomcoeffpd["Norm_Data Average"][:6].sum(axis=0)

        norm_fac = randomcoeffpd["Norm_Data Average"][5] / 0.029

        PROT_sensitivity = self.PROTEIN(randomcoeffpd["Norm_Data Average"][0])
        LIPID_sensitivity = self.LIPID(randomcoeffpd["Norm_Data Average"][4])
        DNA_sensitivity=self.DNA_or_RNA(randomcoeffpd["Norm_Data Average"][1],DNA_or_RNA="DNA")
        RNA_sensitivity=self.DNA_or_RNA(randomcoeffpd["Norm_Data Average"][2],DNA_or_RNA="RNA")
        CARBsyn=self.CARB(randomcoeffpd["Norm_Data Average"][3])

        BIOMASSsyn= self.BIOMASS(norm_fac)
        CARBlist.append(CARBsyn)
        BIOMASSlist.append(BIOMASSsyn)


        with ExcelWriter(self.file2save, mode='w', engine='openpyxl') as writer:

            randomcoeffpd.to_excel(writer,sheet_name="Random coefficient")

            pd.DataFrame(PROT_sensitivity[0].items()).to_excel(writer, sheet_name='PROTsyn')
            pd.DataFrame(LIPID_sensitivity[0].items()).to_excel(writer, sheet_name='LIPIDsyn')
            pd.DataFrame(LIPID_sensitivity[1].items()).to_excel(writer, sheet_name='FATTYACIDsyn')
            pd.DataFrame(LIPID_sensitivity[2]).to_excel(writer, sheet_name="Fatty acid data")
            # pd.DataFrame(PROTlist).to_excel(writer, index=False, sheet_name='PROTsyn_new')
            pd.DataFrame(DNA_sensitivity[0].items()).to_excel(writer, sheet_name='DNAsyn')
            pd.DataFrame(RNA_sensitivity[0].items()).to_excel(writer, sheet_name='RNAsyn')
            pd.DataFrame(CARBlist).to_excel(writer, sheet_name='CARBsyn')

            pd.DataFrame(BIOMASSlist).to_excel(writer, sheet_name='biomass_yeast')
            writer.save()

    def PROTEIN(self,wt_per_DCW):
        PROTDict={}
        PROTpd= pd.read_excel(self.filename, sheet_name='PROTsyn', usecols="A:H" ,convert_float=True)
        PROTin=PROTpd.iloc[:20,0:5]

        AApd = pd.read_excel(self.filename, sheet_name="PROTsyn", usecols="J:U",convert_float=True)
        AApd["mean"] = AApd.mean(axis=1)
        AApd.at[20,"mean"]=AApd["mean"][:20].sum(axis=0)
        max_cv = 0.25
        AApd["min"] = AApd["mean"][:20] - (max_cv * AApd["mean"][:20])
        AApd["max"] = AApd["mean"][:20] + (max_cv * AApd["mean"][:20])



        for aa_i in range(20):
            aa_num_range = range(20)
            aa_num_list = [x for i, x in enumerate(aa_num_range) if i != aa_i]
            TOTAL_1 = AApd["mean"][aa_num_list].sum()

            for min_or_max in ["min","max"]:
                AA_norm = (AApd["mean"])/TOTAL_1*(1 - AApd[min_or_max][aa_i])
                AA_norm[aa_i]=AApd[min_or_max][aa_i]
                minmax_col="aacoeff_" + str(AApd["g/g"][aa_i]) + "_" + str(min_or_max)
                AApd[minmax_col]=AA_norm
                AApd.at[20,minmax_col] = AApd[minmax_col][:20].sum()
        recal_AA=AApd.filter(like="aacoeff",axis=1)


        for i in range(len(recal_AA.columns)):
            for coeff_i in range(len(PROTin)):
                gpergP = recal_AA.iloc[coeff_i,i]
                MW = PROTin.iloc[coeff_i, 1]
                PROTin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                coeff = PROTin.iloc[coeff_i, 4]
                species = PROTin.iloc[coeff_i, 2]
                compart = PROTin.iloc[coeff_i, 3]
                if coeff_i == 0:
                    PROTsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                else:
                    PROTsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

        # PROT out
            PROTout = PROTpd.iloc[:22, 5:8]
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
            # PROTEIN BIOMASS
            PROTsyn = PROTsub + " -> " + PROTpro
            PROTDict[recal_AA.columns[i]]=PROTsyn


        ##mean_ref
        AA_norm_ref = (AApd["mean"])/AApd["mean"][0:20].sum()
        mean_col="aacoeff_" + "mean"
        AApd[mean_col]=AA_norm_ref
        AApd.at[20,mean_col] = AApd[mean_col][:20].sum()
        recal_AA_ref=AApd.filter(like="aacoeff_mean",axis=1)


        for coeff_i in range(len(PROTin)):
            gpergP = recal_AA_ref.iloc[coeff_i,0]
            MW = PROTin.iloc[coeff_i, 1]
            PROTin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
            coeff = PROTin.iloc[coeff_i, 4]
            species = PROTin.iloc[coeff_i, 2]
            compart = PROTin.iloc[coeff_i, 3]
            if coeff_i == 0:
                PROTsub_ref = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                PROTsub_ref += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"


        # PROT_REF out
        PROTout = PROTpd.iloc[:22, 5:8]
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
                PROTpro_ref = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                PROTpro_ref += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"
        #REF_BIOMASS
        PROTsyn_ref = PROTsub_ref + " -> " + PROTpro_ref
        for ref_i in range(44):
            PROTDict["ref"+str(ref_i)] = PROTsyn_ref

        return PROTDict, recal_AA


    def DNA_or_RNA(self,wt_per_DCW, DNA_or_RNA):
        RNAorDNA_Dict = {}

        # Nucleotide in
        #RNA
        if DNA_or_RNA == "RNA":
            RNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            RNAin = RNApd.iloc[:4, 0:5]

            RNA_comp_pd = pd.read_excel(self.filename, sheet_name='RNAsyn', usecols="A:G", skiprows=22 ,convert_float=True)
            RNA_comp_pd = RNA_comp_pd.iloc[:5]
            RNA_comp_pd["mean"] = RNA_comp_pd.loc[:, "A":].mean(axis=1)
            RNA_comp_pd.at[4, "mean"] = RNA_comp_pd["mean"][:4].sum(axis=0)
            RNA_comp_pd["mean"] = RNA_comp_pd["mean"] / RNA_comp_pd["mean"][4]

            max_cv = 0.25
            RNA_comp_pd["min"] = RNA_comp_pd["mean"][:4] - (max_cv * RNA_comp_pd["mean"][:4])
            RNA_comp_pd["max"] = RNA_comp_pd["mean"][:4] + (max_cv * RNA_comp_pd["mean"][:4])


            #ref
            RNA_comp_pd["mean"].at[4] = RNA_comp_pd["mean"][:4].sum()
            recal_RNA_ref = RNA_comp_pd.filter(like="mean", axis=1)

            for i in range(len(recal_RNA_ref.columns)):
                for coeff_i in range(len(RNAin)):
                    gpergP = recal_RNA_ref.iloc[coeff_i, i]

                    MW = RNAin.iloc[coeff_i, 1]
                    RNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                    coeff = RNAin.iloc[coeff_i, 4]
                    species = RNAin.iloc[coeff_i, 2]
                    compart = RNAin.iloc[coeff_i, 3]
                    if coeff_i == 0:
                        RNAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                    else:
                        RNAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                DNAout = RNApd.iloc[0, 5:8]
                coeff = DNAout.iloc[2]
                species = DNAout.iloc[0]
                compart = DNAout.iloc[1]
                RNApro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                RNAsyn_ref = RNAsub + " -> " + RNApro
                for ref_i in range(40):
                    RNAorDNA_Dict[ref_i] = RNAsyn_ref


            for RNA_i in range(4):
                RNA_num_range = range(4)
                RNA_num_list = [x for i, x in enumerate(RNA_num_range) if i != RNA_i]
                TOTAL_1 = RNA_comp_pd["mean"][RNA_num_list].sum()

                for min_or_max in ["min", "max"]:
                    RNA_norm = (RNA_comp_pd["mean"]) / TOTAL_1 * (1 - RNA_comp_pd[min_or_max][RNA_i])
                    RNA_norm[RNA_i] = RNA_comp_pd[min_or_max][RNA_i]
                    minmax_col = "RNAcoeff_" + str(RNA_comp_pd["g/g DNA"][RNA_i]) + "_" + str(min_or_max)
                    RNA_comp_pd[minmax_col] = RNA_norm
                    RNA_comp_pd.at[4, minmax_col] = RNA_comp_pd[minmax_col][:4].sum()

            recal_RNA = RNA_comp_pd.filter(like="RNAcoeff", axis=1)

            for i in range(len(recal_RNA.columns)):
                for coeff_i in range(len(RNAin)):
                    gpergP = recal_RNA.iloc[coeff_i, i]
                    MW = RNAin.iloc[coeff_i, 1]
                    RNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                    coeff = RNAin.iloc[coeff_i, 4]
                    species = RNAin.iloc[coeff_i, 2]
                    compart = RNAin.iloc[coeff_i, 3]
                    if coeff_i == 0:
                        RNAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                    else:
                        RNAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                RNAout = RNApd.iloc[0,5:8]
                coeff = RNAout.iloc[2]
                species = RNAout.iloc[0]
                compart = RNAout.iloc[1]
                RNApro = str(round(coeff,6)) + " " + species + "[" + compart + "]"
                RNAsyn = RNAsub + " -> " + RNApro
                RNAorDNA_Dict[recal_RNA.columns[i]] = RNAsyn
            for ref_i in range(50,86):
                RNAorDNA_Dict[ref_i] = RNAsyn_ref

            return RNAorDNA_Dict, recal_RNA


        elif DNA_or_RNA == "DNA":
            DNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            DNAin = DNApd.iloc[:4, 0:5]

            DNA_comp_pd = pd.read_excel(self.filename, sheet_name='DNAsyn', usecols="A:D", skiprows=10 ,convert_float=True)
            DNA_comp_pd = DNA_comp_pd.iloc[:5]
            print (DNA_comp_pd)
            DNA_comp_pd["mean"] = DNA_comp_pd["A1.2"] # g/g value
            max_cv = 0.25 #RNA_yesast_cv
            DNA_comp_pd["min"] = DNA_comp_pd["mean"][:4] - (max_cv * DNA_comp_pd["mean"][:4])
            DNA_comp_pd["max"] = DNA_comp_pd["mean"][:4] + (max_cv * DNA_comp_pd["mean"][:4])

            DNA_norm= (DNA_comp_pd["mean"])/DNA_comp_pd["mean"][0:4].sum()
            mean_col="DNA_coeff_" + "mean"
            DNA_comp_pd[mean_col]=DNA_norm
            DNA_comp_pd.at[4,mean_col] = DNA_comp_pd[mean_col][:4].sum()
            #ref
            DNA_comp_pd["mean"].at[4] = DNA_comp_pd["mean"][:4].sum()
            recal_DNA_ref = DNA_comp_pd.filter(like="mean", axis=1)

            for i in range(len(recal_DNA_ref.columns)):
                for coeff_i in range(len(DNAin)):
                    gpergP = recal_DNA_ref.iloc[coeff_i, i]

                    MW = DNAin.iloc[coeff_i, 1]
                    DNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                    coeff = DNAin.iloc[coeff_i, 4]
                    species = DNAin.iloc[coeff_i, 2]
                    compart = DNAin.iloc[coeff_i, 3]
                    if coeff_i == 0:
                        DNAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                    else:
                        DNAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                DNAout = DNApd.iloc[0, 5:8]
                coeff = DNAout.iloc[2]
                species = DNAout.iloc[0]
                compart = DNAout.iloc[1]
                DNApro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                DNAsyn_ref = DNAsub + " -> " + DNApro
                for ref_i in range(48):
                    RNAorDNA_Dict[ref_i] = DNAsyn_ref
                # return DNAsyn_ref

            for DNA_i in range(4):
                DNA_num_range = range(4)
                DNA_num_list = [x for i, x in enumerate(DNA_num_range) if i != DNA_i]
                TOTAL_1 = DNA_comp_pd["mean"][DNA_num_list].sum()

                for min_or_max in ["min", "max"]:
                    DNA_norm = (DNA_comp_pd["mean"]) / TOTAL_1 * (1 - DNA_comp_pd[min_or_max][DNA_i])
                    DNA_norm[DNA_i] = DNA_comp_pd[min_or_max][DNA_i]
                    minmax_col = "DNAcoeff_" + str(DNAin["Reactant"][DNA_i]) + "_" + str(min_or_max)
                    DNA_comp_pd[minmax_col] = DNA_norm
                    DNA_comp_pd.at[4, minmax_col] = DNA_comp_pd[minmax_col][:4].sum()

            recal_DNA = DNA_comp_pd.filter(like="DNAcoeff", axis=1)
            recal_DNA.at[4] = recal_DNA[:4].sum()

            for i in range(len(recal_DNA.columns)):
                for coeff_i in range(len(DNAin)):
                    gpergP = recal_DNA.iloc[coeff_i, i]
                    MW = DNAin.iloc[coeff_i, 1]
                    DNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                    coeff = DNAin.iloc[coeff_i, 4]
                    species = DNAin.iloc[coeff_i, 2]
                    compart = DNAin.iloc[coeff_i, 3]
                    if coeff_i == 0:
                        DNAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                    else:
                        DNAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                DNAout = DNApd.iloc[0,5:8]
                coeff = DNAout.iloc[2]
                species = DNAout.iloc[0]
                compart = DNAout.iloc[1]
                DNApro = str(round(coeff,6)) + " " + species + "[" + compart + "]"
                DNAsyn = DNAsub + " -> " + DNApro
                RNAorDNA_Dict[recal_DNA.columns[i]] = DNAsyn



            for ref_i in range(58,86):
                RNAorDNA_Dict[ref_i] = DNAsyn_ref

            return RNAorDNA_Dict, recal_DNA

    def CARB(self,Carb_per_DCW):
        CARBpd = pd.read_excel(self.filename, sheet_name='CARBsyn', usecols="A:H", convert_float=True)
        CARBin = CARBpd.iloc[:5, 0:5]


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

        # PROTEIN BIOMASS
        CARBsyn = CARBsub + " -> 1.0 carbohydrate[c]"
        return CARBsyn

    # 0.835736
    def LIPID(self, LIPID_wt_per_DCW):
        LIPIDDict = {}
        FA_dict = {}
        LIPIDpd = pd.read_excel(self.filename, sheet_name='LIPIDsyn', usecols="A:H", convert_float=True)
        LIPIDin = LIPIDpd.iloc[:16, 0:5]


        LIPIDin = LIPIDin[:-1]

        for coeff_i in range(len(LIPIDin)):
            molpermol = LIPIDin.iloc[coeff_i, 0]
            # MW = LIPIDin.iloc[coeff_i, 1]
            LIPIDin.at[coeff_i, "mmol/gDCW"] = LIPID_wt_per_DCW * molpermol / 3100 * 1000
            coeff = LIPIDin.iloc[coeff_i, 4]
            species = LIPIDin.iloc[coeff_i, 2]
            compart = LIPIDin.iloc[coeff_i, 3]
            if coeff_i == 0:
                LIPIDsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            else:
                LIPIDsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

        LIPIDsyn_ref = LIPIDsub + " -> 1.0 lipid[c]"

        for ref_i in range(56):  #
            LIPIDDict["ref" + str(ref_i)] = LIPIDsyn_ref
        ####################################
        FApd = pd.read_excel(self.filename, sheet_name='FattyAcidSyn', usecols="A:H", convert_float=True)
        FAin = FApd[:4]
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
        max_cv = 0.25
        FA_minmax_pd["min"] = FA_minmax_pd["norm_mean"] * (1 - max_cv)  ##추가
        FA_minmax_pd["max"] = FA_minmax_pd["norm_mean"] * (1 + max_cv)
        FA_list.remove("Sum")
        FA_component_list = FA_list
        print(FA_component_list)
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

        return LIPIDDict, FA_dict, FA_minmax_pd



    def BIOMASS(self,norm_fac):
        BIOMASSpd = pd.read_excel(self.filename, sheet_name='Biomass', usecols="A:H", convert_float=True)
        BIOMASSin = BIOMASSpd.iloc[:10, 0:5]

        for i in range(len(BIOMASSin)):
            if i > 6:
                coeff = BIOMASSin.iloc[i, 4] * norm_fac
            else:
                coeff = BIOMASSin.iloc[i, 4]
            species = BIOMASSin.iloc[i, 2]
            compart = BIOMASSin.iloc[i, 3]
            if i == 0:
                BIOMASSsub = str(round(coeff,6)) + " " + species + "[" + compart + "]"
            else:
                BIOMASSsub += " + " + str(round(coeff,6)) + " " + species + "[" + compart + "]"

        # out
        BIOMASSout = BIOMASSpd.iloc[:4, 5:8]
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
    # filename='&&& for python.xlsx'
    # file2save='name.xlsx'
    a = makeBiomass(filename=filenmae)
    b = a.RandomBiomass()

