#-*-coding:utf-8-*-
import pandas as pd
from pandas import ExcelWriter
from datetime import date, datetime
import random
import numpy as np
import os.path as path



class makeBiomass():

    def __init__(self, filename):
        date = datetime.now().strftime("%b%d %H;%M")
        self.filename = filename
        self.file2save = 'Ecoli_mono_+-25% biomass_{}.xlsx'.format(date)
        if path.isfile(file2save) == True :
            self.file2save=file2save + "(1)"
        # self.filetoSave=filetoSave

    def processBiomassFrame(self,df,in_or_out):
        self.df=df
        self.in_or_out=in_or_out

        if in_or_out == 0 :
            if len(self.df["Reactant"]) >=  len(self.df["Product"]) :
                df_processed=self.df[:len(df["Reactant"])]
            else :
                df_processed= self.df[:len(df["Product"])]

        elif in_or_out == -1 :
            df_processed = self.df.iloc[:,:list(self.df.columns).index("Product")]
            for i in range(len(self.df)) :
                if pd.isna(self.df["Reactant"][i]) == True :
                    df_processed=df_processed[:i]
                    break

        elif in_or_out == 1:
            df_processed = self.df.iloc[:, list(self.df.columns).index("Product"):]
            for i in range(len(self.df)):
                if pd.isna(self.df["Product"][i]) == True:
                    df_processed = df_processed[:i]
                    break
        return df_processed

    def RandomBiomass(self):
        pd.set_option("display.max_columns",999)
        BIOMASSlist=list()

        #Overall composition
        all = pd.read_excel(self.filename, sheet_name='Overall',usecols="A:W", convert_float=True,index_col=0)


        #USE average and stdev
        all["mean"]=all.mean(axis=1,skipna=True)
        all["stdev"]=all.iloc[:,:-1].std(axis=1,skipna=True)


        randomcoeffpd=pd.DataFrame({"Component":["Protein","DNA/Protein(col:(Norm_)Data Average,Stdev) or DNA", "RNA/Protein(col:Data Average,Stdev) or RNA",
                                                 "Lipid","LPS", "Murein", "Inorganic", "Soluble pool","Sum"],
                                    "Data Average":list(all["mean"]),"Data Stdev":list(all["stdev"]),
                                    "Norm_Data Average" : [np.nan]*9,
                                    "Random Average": [np.nan]*9,"Random Stdev" : [np.nan]*9 })
        #DNA or RNA * protein

        randomcoeffpd["Data Average"][-1] = randomcoeffpd["Data Average"][0]*(1+randomcoeffpd["Data Average"][1]+randomcoeffpd["Data Average"][2]) \
                                                  + randomcoeffpd["Data Average"][3:].sum(axis=0)

        randomcoeffpd["Data Average"][1]= randomcoeffpd["Data Average"][1] * randomcoeffpd["Data Average"][0]
        randomcoeffpd["Data Average"][2] = randomcoeffpd["Data Average"][2] * randomcoeffpd["Data Average"][0]

        varing_coeff=randomcoeffpd["Data Average"][:-3].sum(axis=0)
        fixed_coeff=randomcoeffpd[["Data Average"]][6:-1].sum(axis=0)

        #Normalization of average value
        for i in range(len(randomcoeffpd)):
            if i < 6 :
                randomcoeffpd["Norm_Data Average"][i]= randomcoeffpd["Data Average"][i] / varing_coeff * (1-fixed_coeff)
            else :
                randomcoeffpd["Norm_Data Average"][i]= randomcoeffpd["Data Average"][i]


        randomcoeffpd["Norm_Data Average"][-1]= randomcoeffpd["Norm_Data Average"][:-1].sum(axis=0)


        PROT_sensitivity = self.PROTEIN(randomcoeffpd["Norm_Data Average"][0])
        LIPID_sensitivity = self.LIPID(randomcoeffpd["Norm_Data Average"][3])
        RNA_sensitivity = self.DNA_or_RNA(randomcoeffpd["Norm_Data Average"][2],DNA_or_RNA="RNA")
        DNA_sensitivity = self.DNA_or_RNA(randomcoeffpd["Norm_Data Average"][1],DNA_or_RNA="DNA")

        LPSpDCW= randomcoeffpd.iloc[4][3]
        MUREINpDCW= randomcoeffpd.iloc[5][3]
        BIOMASSsyn= self.BIOMASS(LPSpDCW,MUREINpDCW)

        BIOMASSlist.append(BIOMASSsyn)



        with ExcelWriter(self.file2save) as writer:

            randomcoeffpd.to_excel(writer,sheet_name="Random coefficient")

            pd.DataFrame(PROT_sensitivity[0].items()).to_excel(writer, sheet_name='PROTsyn')
            pd.DataFrame(LIPID_sensitivity[0].items()).to_excel(writer, sheet_name='LIPIDsyn')
            LIPID_sensitivity[1].to_excel(writer, sheet_name='FATTYACID_info')
            pd.DataFrame(RNA_sensitivity[0].items()).to_excel(writer, sheet_name='RNAsyn')
            pd.DataFrame(DNA_sensitivity[0].items()).to_excel(writer, sheet_name='DNAsyn')

            pd.DataFrame(BIOMASSlist).to_excel(writer, sheet_name='biomass_ecoli')
            writer.save()

    def PROTEIN(self,wt_per_DCW,):

        PROTDict={}
        PROTpd= pd.read_excel(self.filename, sheet_name='PROTsyn', usecols="A:H" ,convert_float=True)
        PROTin=PROTpd.iloc[:20,0:5]

        AApd = pd.read_excel(self.filename, sheet_name="PROTsyn", usecols="J:W",convert_float=True)
        AApd["mean"] = AApd.mean(axis=1)
        AApd.at[20,"mean"]=AApd["mean"][:20].sum(axis=0)
        #ecoli_aa cv max
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
                minmax_col="aacoeff_" + str(AApd['g/g'][aa_i]) + "_" + str(min_or_max)
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
            PROTout = PROTpd.iloc[:2, 5:8]

            PROTout.at[0, "mmol/gDCW.1"] = PROTin.iloc[:20, 4].sum(axis=0)

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


        ##mean_ref식 만들기


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
        PROTout = PROTpd.iloc[:2, 5:8]

        PROTout.at[0, "mmol/gDCW.1"] = PROTin.iloc[:20, 4].sum(axis=0)

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
        for ref_i in range(34):
            PROTDict["ref" + str(ref_i)] = PROTsyn_ref

        return PROTDict, recal_AA

    def DNA_or_RNA(self,wt_per_DCW, DNA_or_RNA):
        RNAorDNA_Dict = {}

        # Nucleotide in
        #DNAcase
        if DNA_or_RNA == "RNA":
            RNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            RNAin = RNApd.iloc[:4, 0:5]

            RNA_comp_pd = RNApd.iloc[:5, 0:2]
            RNA_comp_pd["mean"] = RNA_comp_pd["g/g RNA"]
            RNA_comp_pd.at[4, "mean"] = RNA_comp_pd["mean"][:4].sum(axis=0)
            RNA_comp_pd["mean"] = RNA_comp_pd["mean"] / RNA_comp_pd["mean"][4]
            # ecoli_max_rna_cv
            max_cv = 0.25
            RNA_comp_pd["min"] = RNA_comp_pd["mean"][:4] - (max_cv * RNA_comp_pd["mean"][:4])
            RNA_comp_pd["max"] = RNA_comp_pd["mean"][:4] + (max_cv * RNA_comp_pd["mean"][:4])

            RNA_norm_ref = (RNA_comp_pd["mean"])/RNA_comp_pd["mean"][0:4].sum()
            mean_col="RNA_coeff_" + "mean"
            RNA_comp_pd[mean_col]=RNA_norm_ref
            RNA_comp_pd.at[4,mean_col] = RNA_comp_pd[mean_col][:4].sum()
            recal_RNA_ref=RNA_comp_pd.filter(like="RNA_coeff_mean",axis=1)

            for coeff_i in range(len(RNAin)):
                gpergP = recal_RNA_ref.iloc[coeff_i, 0]

                MW = RNAin.iloc[coeff_i, 1]
                RNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                coeff = RNAin.iloc[coeff_i, 4]
                species = RNAin.iloc[coeff_i, 2]
                compart = RNAin.iloc[coeff_i, 3]
                if coeff_i == 0:
                    RNAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                else:
                    RNAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

            # RNA out
            RNAout = RNApd.iloc[:2, 5:8]

            # ppi
            RNAout.at[0, "mmol/gDCW.1"] = RNAin.iloc[:, 4].sum(axis=0)

            for j in range(len(RNAout)):
                coeff = RNAout.iloc[j, 2]
                species = RNAout.iloc[j, 0]
                compart = RNAout.iloc[j, 1]
                if j == 0:
                    RNApro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                else:
                    RNApro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

            # PROTEIN BIOMASS
            RNAsyn_ref = RNAsub + " -> " + RNApro

            for ref_i in range(40):
                RNAorDNA_Dict["ref" + str(ref_i)] = RNAsyn_ref


            for RNA_i in range(len(RNA_comp_pd) - 1):
                lp_num_range = range(len(RNA_comp_pd) - 1)
                lp_num_list = [x for i, x in enumerate(lp_num_range) if i != RNA_i]
                TOTAL_1 = RNA_comp_pd["mean"][lp_num_list].sum()

                for min_or_max in ["min", "max"]:
                    RNA_norm = (RNA_comp_pd["mean"]) / TOTAL_1 * (1 - RNA_comp_pd[min_or_max][RNA_i])
                    RNA_norm[RNA_i] = RNA_comp_pd[min_or_max][RNA_i]
                    minmax_col = "RNAcoeff_" + str(RNApd["Reactant"][RNA_i]) + "_" + str(min_or_max)
                    RNA_comp_pd[minmax_col] = RNA_norm
                    RNA_comp_pd.at[4, minmax_col] = RNA_comp_pd[minmax_col][:4].sum()
            recal_RNA = RNA_comp_pd.filter(like="RNAcoeff_", axis=1)

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

                # RNA out
                RNAout = RNApd.iloc[:2, 5:8]

                # ppi
                RNAout.at[0, "mmol/gDCW.1"] = RNAin.iloc[:, 4].sum(axis=0)

                for j in range(len(RNAout)):
                    coeff = RNAout.iloc[j, 2]
                    species = RNAout.iloc[j, 0]
                    compart = RNAout.iloc[j, 1]
                    if j == 0:
                        RNApro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                    else:
                        RNApro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                # PROTEIN BIOMASS
                RNAsyn = RNAsub + " -> " + RNApro
                RNAorDNA_Dict[recal_RNA.columns[i]] = RNAsyn

            for ref_i in range(48,74):
                RNAorDNA_Dict["ref" + str(ref_i)] = RNAsyn_ref
            return RNAorDNA_Dict, recal_RNA


        elif DNA_or_RNA == "DNA":
            DNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            DNAin = DNApd.iloc[:4, 0:5]

            DNA_comp_pd = DNApd.iloc[:5, 0:2]
            DNA_comp_pd["mean"] = DNA_comp_pd["g/g DNA"]
            DNA_comp_pd.at[4, "mean"] = DNA_comp_pd["mean"][:4].sum(axis=0)
            DNA_comp_pd["mean"] = DNA_comp_pd["mean"] / DNA_comp_pd["mean"][4]
            # ecoli_max_RNA_cv
            max_cv = 0.25
            DNA_comp_pd["min"] = DNA_comp_pd["mean"][:4] - (max_cv * DNA_comp_pd["mean"][:4])
            DNA_comp_pd["max"] = DNA_comp_pd["mean"][:4] + (max_cv * DNA_comp_pd["mean"][:4])


            DNA_norm_ref = (DNA_comp_pd["mean"]) / DNA_comp_pd["mean"][0:4].sum()
            mean_col = "DNA_coeff_" + "mean"
            DNA_comp_pd[mean_col] = DNA_norm_ref
            DNA_comp_pd.at[4, mean_col] = DNA_comp_pd[mean_col][:4].sum()
            recal_DNA_ref = DNA_comp_pd.filter(like="DNA_coeff_mean", axis=1)

            for coeff_i in range(len(DNAin)):
                gpergP = recal_DNA_ref.iloc[coeff_i, 0]

                MW = DNAin.iloc[coeff_i, 1]
                DNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                coeff = DNAin.iloc[coeff_i, 4]
                species = DNAin.iloc[coeff_i, 2]
                compart = DNAin.iloc[coeff_i, 3]
                if coeff_i == 0:
                    DNAsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                else:
                    DNAsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                DNAout = DNApd.iloc[:2, 5:8]

                # ppi
                DNAout.at[0, "mmol/gDCW.1"] = DNAin.iloc[:, 4].sum(axis=0)

                for j in range(len(DNAout)):
                    coeff = DNAout.iloc[j, 2]
                    species = DNAout.iloc[j, 0]
                    compart = DNAout.iloc[j, 1]
                    if j == 0:
                        DNApro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                    else:
                        DNApro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                # PROTEIN BIOMASS
                DNAsyn_ref = DNAsub + " -> " + DNApro



            for ref_i in range(48):
                RNAorDNA_Dict[str(ref_i) + "mean"] = DNAsyn_ref

            for DNA_i in range(len(DNA_comp_pd) - 1):
                lp_num_range = range(len(DNA_comp_pd) - 1)
                lp_num_list = [x for i, x in enumerate(lp_num_range) if i != DNA_i]
                TOTAL_1 = DNA_comp_pd["mean"][lp_num_list].sum()

                for min_or_max in ["min", "max"]:
                    DNA_norm = (DNA_comp_pd["mean"]) / TOTAL_1 * (1 - DNA_comp_pd[min_or_max][DNA_i])
                    DNA_norm[DNA_i] = DNA_comp_pd[min_or_max][DNA_i]
                    minmax_col = "DNA_coeff_" + str(DNApd["Reactant"][DNA_i]) + "_" + str(min_or_max)
                    DNA_comp_pd[minmax_col] = DNA_norm
                    DNA_comp_pd.at[4, minmax_col] = DNA_comp_pd[minmax_col][:4].sum()
            recal_DNA = DNA_comp_pd.filter(like="DNA_coeff_0", axis=1)

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

                # DNA out
                DNAout = DNApd.iloc[:2, 5:8]

                # ppi
                DNAout.at[0, "mmol/gDCW.1"] = DNAin.iloc[:, 4].sum(axis=0)

                for j in range(len(DNAout)):
                    coeff = DNAout.iloc[j, 2]
                    species = DNAout.iloc[j, 0]
                    compart = DNAout.iloc[j, 1]
                    if j == 0:
                        DNApro = str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                    else:
                        DNApro += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"

                # PROTEIN BIOMASS
                DNAsyn = DNAsub + " -> " + DNApro
                RNAorDNA_Dict[recal_DNA.columns[i]] = DNAsyn



            for ref_i in range(56,74):
                RNAorDNA_Dict["mean" + str(ref_i)] = DNAsyn_ref

            return RNAorDNA_Dict, recal_DNA


    def LIPID(self,LIPID_wt_per_DCW):
        LIPIDDict = {}
        LIPIDpd= pd.read_excel(self.filename, sheet_name='LIPIDsyn', usecols="A:H" ,convert_float=True)

        #Lipid in

        LIPIDin = self.processBiomassFrame(df=LIPIDpd,in_or_out=-1)
        index_list=list(LIPIDin.index.values)+["Sum"]
        LIPIDminmax_df=pd.DataFrame({"Component":index_list})
        LIPIDminmax_df["norm_mean"] = LIPIDin["g/g Lipid"]/LIPIDin["g/g Lipid"].sum()
        LIPIDminmax_df.loc["Sum", "norm_mean"] = LIPIDminmax_df["norm_mean"][:-1].sum(axis=0)


        ##mean_ref

        for i in range(len(LIPIDin)):
            gpergLIPID=LIPIDin.iloc[i,0]
            MW=LIPIDin.iloc[i,1]
            LIPIDin.at[i,"mmol/gDCW"]= gpergLIPID/MW*1000*LIPID_wt_per_DCW

            coeff= LIPIDin.iloc[i,4]
            species= LIPIDin.iloc[i,2]
            compart= LIPIDin.iloc[i,3]

            if i==0 :
                LIPIDsub= str(round(coeff,6)) + " "+ species + "["+compart+"]"
            else :
                LIPIDsub += " + " + str(round(coeff,6)) + " " + species +"["+compart+"]"

        #Lipid out

        # LIPID BIOMASS
        LIPIDsyn_ref = LIPIDsub + " -> 1.0 lipid[c]"
        for ref_i in range(56):
            LIPIDDict[ref_i] = LIPIDsyn_ref


#FATTY ACID
        mono_norm_fac= LIPIDin["g/g Lipid"].sum()
        Lipid_mono_ratio_df=pd.DataFrame({"pe": np.nan*3,"pg": np.nan*3,"clpn": np.nan*3},index=["160","161","181"])

        for l_i,l in enumerate(LIPIDin["Reactant"]):
            if "pe" in l :
                col_loc= 0
            elif "pg" in l :
                col_loc =1
            elif "clpn" in l :
                col_loc= 2

            if "160" in l :
                row_loc=0
            elif "161" in l :
                row_loc=1
            elif "181" in l :
                row_loc=2
            Lipid_mono_ratio_df.iloc[row_loc,col_loc]= LIPIDin["g/g Lipid"][l_i] / mono_norm_fac

        print ("mono_nomr", mono_norm_fac)
        print (Lipid_mono_ratio_df)

        Lipid_mono_ratio_df["FA_Sum"]=Lipid_mono_ratio_df.sum(axis=1)
        Lipid_mono_ratio_df.loc["Sum"]=Lipid_mono_ratio_df.sum(axis=0)
        print(Lipid_mono_ratio_df)
        Fatty_acid_ratio=pd.DataFrame({"ref": list(Lipid_mono_ratio_df["FA_Sum"][:-1])},index=["160","161","181"])
        print (Fatty_acid_ratio)
        Lipid_kind_ratio= pd.DataFrame({"pe":[Lipid_mono_ratio_df.loc["Sum","pe"]],"pg":[Lipid_mono_ratio_df.loc["Sum","pg"]],
                                        "clpn":[Lipid_mono_ratio_df.loc["Sum","clpn"]]},index=["ref"])
        print (Lipid_kind_ratio)

        for fa in ["160","161","181"] :
            for min_or_max in ["min", "max"] :
                if min_or_max == "min" :
                    fa_min_max_val=Fatty_acid_ratio.loc[fa,"ref"] * 0.75
                elif min_or_max == "max":
                    fa_min_max_val = Fatty_acid_ratio.loc[fa,"ref"] * 1.25

                Fatty_acid_ratio[fa + "_" + min_or_max]=np.nan*len(Fatty_acid_ratio.index)
                for r in range(len(Fatty_acid_ratio.index)) :
                    Fatty_acid_ratio[fa + "_" + min_or_max][r]= Fatty_acid_ratio["ref"][r]/(1-Fatty_acid_ratio.loc[fa,"ref"])*(1-fa_min_max_val)
                Fatty_acid_ratio.loc[fa,fa+"_"+min_or_max] = fa_min_max_val

        print (Fatty_acid_ratio)

        #make biomass
        LIPIDDict= dict()
        for fa_mm_i in range(len(Fatty_acid_ratio.columns)):

            for l_i,l_kind in enumerate(["pe","pg","clpn"]):

                Lipid_mono_ratio_df.iloc[0,l_i] = Fatty_acid_ratio.iloc[0,fa_mm_i] * Lipid_kind_ratio[l_kind][0]
                Lipid_mono_ratio_df.iloc[1,l_i] =Fatty_acid_ratio.iloc[1,fa_mm_i] * Lipid_kind_ratio[l_kind][0]
                Lipid_mono_ratio_df.iloc[2,l_i] =Fatty_acid_ratio.iloc[2,fa_mm_i] * Lipid_kind_ratio[l_kind][0]

            lipid_coeff_dict=dict()
            for c in ["pe","pg","clpn"]:
                for f in ["160","161","181"]:
                    lipid_mono=c+f
                    monomer_g=Lipid_mono_ratio_df.loc[f,c]
                    lipid_coeff_dict[lipid_mono]=monomer_g


            for i in range(len(LIPIDin)):
                lipid_comp=LIPIDin["Reactant"][i]
                gpergP = lipid_coeff_dict[lipid_comp]
                MW = LIPIDin.iloc[i, 1]
                LIPIDin.at[i, "mmol/gDCW"] =  gpergP / MW * 1000 * LIPID_wt_per_DCW
                coeff = LIPIDin.iloc[i, 4]
                species = LIPIDin.iloc[i, 2]
                compart = LIPIDin.iloc[i, 3]
                if i == 0:
                    LIPIDsub = str(round(coeff, 6)) + " " + species + "[" + compart + "]"
                else:
                    LIPIDsub += " + " + str(round(coeff, 6)) + " " + species + "[" + compart + "]"
            # LIPID out
            LIPIDsyn = LIPIDsub + " -> 1.0 lipid[c]"

            LIPIDDict[Fatty_acid_ratio.columns[fa_mm_i]] = LIPIDsyn

        return LIPIDDict, Fatty_acid_ratio


    def BIOMASS(self,LPS_per_DCW,Murein_per_DCW):
        BIOMASSpd = pd.read_excel(self.filename, sheet_name='Biomass', usecols="A:H", convert_float=True)
        BIOMASSin = BIOMASSpd.iloc[:, 0:5]


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
        BIOMASSout = BIOMASSpd.iloc[:3, 5:8]

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

    #filename='&&& for python.xlsx'
    #file2save='name.xlsx'
    a=makeBiomass(filename=filename)
    b=a.RandomBiomass()



