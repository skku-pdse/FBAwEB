#-*-coding:utf-8-*-
import pandas as pd
from pandas import ExcelWriter
from datetime import datetime
import numpy as np
from pathlib import Path
import os

class makeBiomass():

    def __init__(self, filename):
        date = datetime.now().strftime("%b%d %H;%M")
        self.filename = filename
        saveDir="RandomBiomassExcel"
        Path(saveDir).mkdir(parents=True,exist_ok=True)
        savefile = 'Ecoli_mono_+-25% biomass_{}.xlsx'.format(date)
        self.file2save=os.path.join(saveDir,savefile)
        # self.filetoSave=filetoSave

    def _processBiomassFrame(self,df,in_or_out):
        if in_or_out == 0 :
            if len(df["Reactant"]) >=  len(df["Product"]) :
                df_processed=df[:len(df["Reactant"])]
            else :
                df_processed= df[:len(df["Product"])]

        elif in_or_out == -1 :
            df_processed = df.iloc[:,:list(df.columns).index("Product")]
            for i in range(len(df)) :
                if pd.isna(df["Reactant"][i]) == True :
                    df_processed=df_processed[:i]
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

    def Formulate_min_max(self):
        pd.set_option("display.max_columns",999)

        #Overall composition
        all = pd.read_excel(self.filename, sheet_name='Overall',usecols="A:W", convert_float=True,index_col=0)


        #USE average and stdev
        all["mean"]=all.mean(axis=1,skipna=True)
        all["stdev"]=all.iloc[:,:-1].std(axis=1,skipna=True)

        componentList =["Protein","DNA/Protein(col:(Norm_)Data Average,Stdev) or DNA", "RNA/Protein(col:Data Average,Stdev) or RNA",
                                                 "Carbohydrate", "Lipid","LPS", "Murein", "Inorganic", "Soluble pool","Sum"]
        randomcoeffpd=pd.DataFrame({"Component":componentList,
                                    "Data Average":list(all["mean"]),"Data Stdev":list(all["stdev"]),
                                    "Norm_Data Average" : [np.nan]*len(componentList),
                                    "Random Average": [np.nan]*len(componentList),"Random Stdev" : [np.nan]*len(componentList)})
        #DNA or RNA * protein

        randomcoeffpd["Data Average"][-1] = randomcoeffpd["Data Average"][0]*(1+randomcoeffpd["Data Average"][1]+randomcoeffpd["Data Average"][2]) \
                                                  + randomcoeffpd["Data Average"][3:].sum(axis=0)

        randomcoeffpd["Data Average"][1]= randomcoeffpd["Data Average"][1] * randomcoeffpd["Data Average"][0]
        randomcoeffpd["Data Average"][2] = randomcoeffpd["Data Average"][2] * randomcoeffpd["Data Average"][0]

        varing_coeff=randomcoeffpd["Data Average"][:-3].sum(axis=0)
        fixed_coeff=randomcoeffpd[["Data Average"]][7:-1].sum(axis=0)

        #Normalization of average value
        for i in range(len(randomcoeffpd)):
            if i < 7 :
                randomcoeffpd["Norm_Data Average"][i]= randomcoeffpd["Data Average"][i] / varing_coeff * (1-fixed_coeff)
            else :
                randomcoeffpd["Norm_Data Average"][i]= randomcoeffpd["Data Average"][i]


        randomcoeffpd["Norm_Data Average"][-1]= randomcoeffpd["Norm_Data Average"][:-1].sum(axis=0)


        PROT_sensitivity = self.PROTEIN(randomcoeffpd["Norm_Data Average"][0])
        LIPID_sensitivity = self.LIPID(randomcoeffpd["Norm_Data Average"][4])
        CARB_sensitivity = self.CARB(randomcoeffpd["Norm_Data Average"][3])
        RNA_sensitivity = self.DNA_or_RNA(randomcoeffpd["Norm_Data Average"][2],DNA_or_RNA="RNA")
        DNA_sensitivity = self.DNA_or_RNA(randomcoeffpd["Norm_Data Average"][1],DNA_or_RNA="DNA")

        LPSpDCW= randomcoeffpd.iloc[5][3]
        MUREINpDCW= randomcoeffpd.iloc[6][3]
        BIOMASSsyn= self.BIOMASS(LPSpDCW,MUREINpDCW)

        with ExcelWriter(self.file2save) as writer:

            randomcoeffpd.to_excel(writer,sheet_name="Random coefficient")
            pd.DataFrame(PROT_sensitivity[0].items()).to_excel(writer, sheet_name='PROTsyn')
            pd.DataFrame(LIPID_sensitivity[0].items()).to_excel(writer, sheet_name='LIPIDsyn')
            LIPID_sensitivity[1].to_excel(writer, sheet_name='FATTYACID_info')
            pd.DataFrame(CARB_sensitivity.items()).to_excel(writer, sheet_name='CARBsyn')
            pd.DataFrame(RNA_sensitivity[0].items()).to_excel(writer, sheet_name='RNAsyn')
            pd.DataFrame(DNA_sensitivity[0].items()).to_excel(writer, sheet_name='DNAsyn')

            pd.DataFrame(BIOMASSsyn.items()).to_excel(writer, sheet_name='biomass_ecoli')
            writer.save()

    def PROTEIN(self,wt_per_DCW):
        PROTDict={}
        PROTpd= pd.read_excel(self.filename, sheet_name='PROTsyn', usecols="A:H" ,convert_float=True)
        PROTin=self._processBiomassFrame(df=PROTpd,in_or_out=-1)

        AApd = pd.DataFrame(PROTin[['g/g Protein','Reactant']])
        AApd["mean"] = AApd.mean(axis=1)
        AApd.at[20,"mean"]=AApd["mean"][:20].sum(axis=0)

        for coeff_i in range(len(PROTin)):
            gpergP = PROTin.iloc[coeff_i, 0]
            MW = PROTin.iloc[coeff_i, 1]
            PROTin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000

        # PROT out
        PROTout=self._processBiomassFrame(df=PROTpd,in_or_out=1)
        for k_i, k in enumerate(PROTout["Product"]):
            if k.lower() == 'h2o':
                PROTout.at[k_i, "mmol/gDCW.1"] = PROTin.iloc[:20, 4].sum(axis=0)
            elif k.lower() in ("protein", "prot"):
                PROTout.at[k_i, "mmol/gDCW.1"] = 1
            else:  # amino acids
                for ami_i, ami in enumerate(PROTin["Reactant"]):
                    if k.lower() in ami.lower():
                        PROTout.at[k_i, "mmol/gDCW.1"] = PROTin.iloc[ami_i, 4]
                        break
        # REF_BIOMASS
        PROTsyn_ref = self._returnSynthesisEquation(df_in=PROTin, df_out=PROTout)
        PROTDict["ref"] = PROTsyn_ref

        ##
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
                minmax_col="aacoeff_" + str(AApd.iloc[aa_i,1]) + "_" + str(min_or_max)
                AApd[minmax_col]=AA_norm
                AApd.at[20,minmax_col] = AApd[minmax_col][:20].sum()
        recal_AA=AApd.filter(like="aacoeff",axis=1)


        for i in range(len(recal_AA.columns)): # aa min or max case
            for coeff_i in range(len(PROTin)):
                gpergP = recal_AA.iloc[coeff_i,i]
                MW = PROTin.iloc[coeff_i, 1]
                PROTin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000

        # PROT out
            # PROTout coeff cal
            for k_i, k in enumerate(PROTout["Product"]):
                if k.lower() == 'h2o':
                    PROTout.at[k_i, "mmol/gDCW.1"] = PROTin.iloc[:20, 4].sum(axis=0)
                elif k.lower() in ("protein", "prot"):
                    PROTout.at[k_i, "mmol/gDCW.1"] = 1
                else:  # amino acids
                    for ami_i, ami in enumerate(PROTin["Reactant"]):
                        if k.lower() in ami.lower():
                            PROTout.at[k_i, "mmol/gDCW.1"] = PROTin.iloc[ami_i, 4]
                            break

            PROTsyn = self._returnSynthesisEquation(df_in=PROTin,df_out=PROTout)
            PROTDict[recal_AA.columns[i]]=PROTsyn


        return PROTDict, recal_AA


    def DNA_or_RNA(self,wt_per_DCW, DNA_or_RNA):
        RNAorDNA_Dict = {}

        # Nucleotide in
        #DNAcase
        if DNA_or_RNA == "RNA":
            RNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            RNAin = self._processBiomassFrame(df=RNApd,in_or_out=-1)

            RNA_comp_pd = RNAin.iloc[:, :-1]
            RNA_comp_pd["mean"] = RNA_comp_pd['g/g RNA']

            for coeff_i in range(len(RNA_comp_pd)):
                gpergP = RNAin.iloc[coeff_i, 0]
                MW = RNAin.iloc[coeff_i, 1]
                RNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000

            RNAout = self._processBiomassFrame(df=RNApd, in_or_out=1)
            for k_i, k in enumerate(RNAout['Product']):
                if k.lower() in ['ppi', 'pi']:
                    RNAout.at[k_i, "mmol/gDCW.1"] = RNAin['mmol/gDCW'].sum()
            RNAsyn_ref = self._returnSynthesisEquation(df_in=RNAin, df_out=RNAout)
            RNAorDNA_Dict['ref'] = RNAsyn_ref

            max_cv = 0.25
            RNA_comp_pd["min"] = RNA_comp_pd["mean"][:4] - (max_cv * RNA_comp_pd["mean"][:4])
            RNA_comp_pd["max"] = RNA_comp_pd["mean"][:4] + (max_cv * RNA_comp_pd["mean"][:4])

            for RNA_i in range(4):
                RNA_num_range = range(4)
                RNA_num_list = [x for i, x in enumerate(RNA_num_range) if i != RNA_i]
                TOTAL_1 = RNA_comp_pd["mean"][RNA_num_list].sum()

                for min_or_max in ["min", "max"]:
                    RNA_norm = (RNA_comp_pd["mean"]) / TOTAL_1 * (1 - RNA_comp_pd[min_or_max][RNA_i])
                    RNA_norm[RNA_i] = RNA_comp_pd[min_or_max][RNA_i]
                    minmax_col = "RNAcoeff_" + str(RNA_comp_pd["Reactant"][RNA_i]) + "_" + str(min_or_max)
                    RNA_comp_pd[minmax_col] = RNA_norm
                    RNA_comp_pd.at[4, minmax_col] = RNA_comp_pd[minmax_col][:4].sum()
            recal_RNA = RNA_comp_pd.filter(like="RNAcoeff", axis=1)

            for i in range(len(recal_RNA.columns)):
                for coeff_i in range(len(RNAin)):
                    gpergP = recal_RNA.iloc[coeff_i, i]
                    MW = RNAin.iloc[coeff_i, 1]
                    RNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                for k_i, k in enumerate(RNAout['Product']):
                    if k.lower() in ['ppi', 'pi']:
                        RNAout.at[k_i, "mmol/gDCW.1"] = RNAin['mmol/gDCW'].sum()
                RNAsyn = self._returnSynthesisEquation(df_in=RNAin, df_out=RNAout)
                RNAorDNA_Dict[recal_RNA.columns[i]] = RNAsyn

            return RNAorDNA_Dict, recal_RNA



        elif DNA_or_RNA == "DNA":
            DNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            DNAin = self._processBiomassFrame(df=DNApd,in_or_out=-1)

            DNA_comp_pd = DNAin.iloc[:,:-1]
            DNA_comp_pd["mean"] = DNA_comp_pd["g/g DNA"]

            for coeff_i in range(len(DNA_comp_pd)):
                gpergP = DNAin.iloc[coeff_i, 0]
                MW = DNAin.iloc[coeff_i, 1]
                DNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000

            DNAout = self._processBiomassFrame(df=DNApd, in_or_out=1)
            for k_i, k in enumerate(DNAout['Product']):
                if k.lower() in ['ppi', 'pi']:
                    DNAout.at[k_i, "mmol/gDCW.1"] = DNAin['mmol/gDCW'].sum()
            DNAsyn_ref = self._returnSynthesisEquation(df_in=DNAin, df_out=DNAout)
            RNAorDNA_Dict['ref'] = DNAsyn_ref

            # ecoli_max_RNA_cv
            max_cv = 0.25
            DNA_comp_pd["min"] = DNA_comp_pd["mean"][:4] - (max_cv * DNA_comp_pd["mean"][:4])
            DNA_comp_pd["max"] = DNA_comp_pd["mean"][:4] + (max_cv * DNA_comp_pd["mean"][:4])

            for DNA_i in range(4):
                lp_num_range = range(4)
                lp_num_list = [x for i, x in enumerate(lp_num_range) if i != DNA_i]
                TOTAL_1 = DNA_comp_pd["mean"][lp_num_list].sum()

                for min_or_max in ["min", "max"]:
                    DNA_norm = (DNA_comp_pd["mean"]) / TOTAL_1 * (1 - DNA_comp_pd[min_or_max][DNA_i])
                    DNA_norm[DNA_i] = DNA_comp_pd[min_or_max][DNA_i]
                    minmax_col = "DNA_coeff_" + str(DNApd["Reactant"][DNA_i]) + "_" + str(min_or_max)
                    DNA_comp_pd[minmax_col] = DNA_norm
                    DNA_comp_pd.at[4, minmax_col] = DNA_comp_pd[minmax_col][:4].sum()
            recal_DNA = DNA_comp_pd.filter(like="DNA_coeff_", axis=1)

            for i in range(len(recal_DNA.columns)):
                for coeff_i in range(len(DNAin)):
                    gpergP = recal_DNA.iloc[coeff_i, i]
                    MW = DNAin.iloc[coeff_i, 1]
                    DNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000
                for k_i, k in enumerate(DNAout['Product']):
                    if k.lower() in ['ppi', 'pi']:
                        DNAout.at[k_i, "mmol/gDCW.1"] = DNAin['mmol/gDCW'].sum()
                DNAsyn = self._returnSynthesisEquation(df_in=DNAin, df_out=DNAout)
                RNAorDNA_Dict[recal_DNA.columns[i]] = DNAsyn

            return RNAorDNA_Dict, recal_DNA

    def CARB(self,Carb_per_DCW):
        CARB_Dict={}
        CARBpd = pd.read_excel(self.filename, sheet_name='CARBsyn', usecols="A:H", convert_float=True)
        CARBin = self._processBiomassFrame(df=CARBpd, in_or_out=-1)

        for i in range(len(CARBin)):
            gpergCARB = CARBin.iloc[i, 0]
            MW = CARBin.iloc[i, 1]
            CARBin.at[i, "mmol/gDCW"] = Carb_per_DCW * gpergCARB / MW * 1000

        CARBout = self._processBiomassFrame(df=CARBpd, in_or_out=1)
        for p_ind, p in enumerate(CARBout["Product"]):
            if any(c == p.lower() for c in ("carbohydrate", "carb", "carbo")):
                CARBout.at[p_ind, "mmol/gDCW.1"] = 1
        print(CARBout)

        # PROTEIN BIOMASS
        CARBsyn = self._returnSynthesisEquation(df_in=CARBin, df_out=CARBout)
        CARB_Dict['ref'] = CARBsyn
        return CARB_Dict


    def LIPID(self,LIPID_wt_per_DCW):
        LIPIDDict = {}
        LIPIDpd= pd.read_excel(self.filename, sheet_name='LIPIDsyn', usecols="A:H" ,convert_float=True)

        #Lipid in

        LIPIDin = self._processBiomassFrame(df=LIPIDpd,in_or_out=-1)
        index_list=list(LIPIDin.index.values)+["Sum"]
        LIPIDminmax_df=pd.DataFrame({"Component":index_list})
        LIPIDminmax_df["norm_mean"] = LIPIDin["g/g Lipid"]/LIPIDin["g/g Lipid"].sum()
        LIPIDminmax_df.loc["Sum", "norm_mean"] = LIPIDminmax_df["norm_mean"][:-1].sum(axis=0)


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

            LIPIDout=self._processBiomassFrame(df=LIPIDpd,in_or_out=1)
            LIPIDsyn=self._returnSynthesisEquation(df_in=LIPIDin,df_out=LIPIDout)

            LIPIDDict[Fatty_acid_ratio.columns[fa_mm_i]] = LIPIDsyn

        return LIPIDDict, Fatty_acid_ratio


    def BIOMASS(self,LPS_per_DCW,Murein_per_DCW):
        BIOMASSpd = pd.read_excel(self.filename, sheet_name='Biomass', usecols="A:H", convert_float=True)
        BIOMASSin = self._processBiomassFrame(df=BIOMASSpd,in_or_out=-1)

        #LPS
        LPSMW=BIOMASSin.iloc[0,1]
        BIOMASSin.at[0,"mmol/gDCW"]= LPS_per_DCW / LPSMW * 1000

        #Murein
        MureinMW=BIOMASSin.iloc[1,1]
        BIOMASSin.at[1,"mmol/gDCW"]= Murein_per_DCW / MureinMW *1000\

       # out
        BIOMASSout = self._processBiomassFrame(df=BIOMASSpd,in_or_out=1)

        # PROTEIN BIOMASS
        BIOMASSsyn = self._returnSynthesisEquation(df_in=BIOMASSin,df_out=BIOMASSout)
        BIOMASSdict={'ref':BIOMASSsyn}
        return BIOMASSdict


if __name__ == '__main__':
    # fileloc = ''
    # filename = fileloc + '\Ecoli test1_sensitivitiyOnly.xlsx'
    filename='Ecoli test1_sensitivityOnly.xlsx'
    a=makeBiomass(filename=filename)
    b=a.Formulate_min_max()



