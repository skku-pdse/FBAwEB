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
        savefile = 'Yeast_mono_+-25% biomass_{}.xlsx'.format(date)
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

    def Formulate_min_max(self):
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
                                    "Norm_Data Average": [np.nan] * len(componentList),
                                    "Random Average": [np.nan] * len(componentList), "Random Stdev": [np.nan] * len(componentList)})

        randomcoeffpd.at[len(componentList)-1, "Data Average"]= randomcoeffpd["Data Average"][:len(componentList)-1].sum(axis=0)

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

        with ExcelWriter(self.file2save) as writer:

            randomcoeffpd.to_excel(writer,sheet_name="Random coefficient")

            pd.DataFrame(PROT_sensitivity[0].items()).to_excel(writer, sheet_name='PROTsyn')
            pd.DataFrame(LIPID_sensitivity[0].items()).to_excel(writer, sheet_name='LIPIDsyn')
            pd.DataFrame(LIPID_sensitivity[1].items()).to_excel(writer, sheet_name='FATTYACIDsyn')
            pd.DataFrame(LIPID_sensitivity[2]).to_excel(writer, sheet_name="Fatty acid data")
            pd.DataFrame(DNA_sensitivity[0].items()).to_excel(writer, sheet_name='DNAsyn')
            pd.DataFrame(RNA_sensitivity[0].items()).to_excel(writer, sheet_name='RNAsyn')
            pd.DataFrame(CARBsyn.items()).to_excel(writer, sheet_name='CARBsyn')
            pd.DataFrame(BIOMASSsyn.items()).to_excel(writer, sheet_name='biomass_yeast')
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

        # PROT_REF out
        # PROTout coeff cal
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
        #RNA
        if DNA_or_RNA == "RNA":
            RNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            RNAin = self._processBiomassFrame(df=RNApd,in_or_out=-1)

            RNA_comp_pd = RNAin.iloc[:,:-1]
            RNA_comp_pd["mean"] = RNA_comp_pd['g/g RNA']


            for coeff_i in range(len(RNA_comp_pd)):
                gpergP = RNAin.iloc[coeff_i,0]
                MW = RNAin.iloc[coeff_i, 1]
                RNAin.at[coeff_i, "mmol/gDCW"] = wt_per_DCW * gpergP / MW * 1000

            RNAout = self._processBiomassFrame(df=RNApd,in_or_out=1)
            for k_i,k in enumerate(RNAout['Product']):
                if k.lower() in ['ppi', 'pi'] :
                    RNAout.at[k_i,"mmol/gDCW.1"]= RNAin['mmol/gDCW'].sum()
            RNAsyn_ref = self._returnSynthesisEquation(df_in=RNAin,df_out=RNAout)
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
                RNAsyn = self._returnSynthesisEquation(df_in=RNAin,df_out=RNAout)
                RNAorDNA_Dict[recal_RNA.columns[i]] = RNAsyn

            return RNAorDNA_Dict, recal_RNA


        elif DNA_or_RNA == "DNA":
            DNApd = pd.read_excel(self.filename, sheet_name=DNA_or_RNA + "syn", usecols="A:H", convert_float=True)
            DNAin = self._processBiomassFrame(df=DNApd,in_or_out=-1)

            DNA_comp_pd = DNAin.iloc[:,:-1]
            DNA_comp_pd["mean"] = DNA_comp_pd["g/g DNA"] # g/g value

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
                # return DNAsyn_ref

            max_cv = 0.25  # RNA_yesast_cv
            DNA_comp_pd["min"] = DNA_comp_pd["mean"][:4] - (max_cv * DNA_comp_pd["mean"][:4])
            DNA_comp_pd["max"] = DNA_comp_pd["mean"][:4] + (max_cv * DNA_comp_pd["mean"][:4])

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

                DNAsyn = self._returnSynthesisEquation(df_in=DNAin,df_out=DNAout)
                RNAorDNA_Dict[recal_DNA.columns[i]] = DNAsyn

            return RNAorDNA_Dict, recal_DNA

    def CARB(self,Carb_per_DCW):
        CARBpd = pd.read_excel(self.filename, sheet_name='CARBsyn', usecols="A:H", convert_float=True)
        CARBin = self._processBiomassFrame(df=CARBpd,in_or_out=-1)


        for i in range(len(CARBin)):
            gpergCARB = CARBin.iloc[i, 0]
            MW = CARBin.iloc[i, 1]
            CARBin.at[i, "mmol/gDCW"] = Carb_per_DCW * gpergCARB / MW * 1000

        CARBout =self._processBiomassFrame(df=CARBpd,in_or_out=1)
        # PROTEIN BIOMASS
        CARBsyn = self._returnSynthesisEquation(df_in=CARBin,df_out=CARBout)
        CARBdict={'ref': CARBsyn}
        return CARBdict

    # 0.835736
    def LIPID(self, LIPID_wt_per_DCW):
        LIPIDDict = {}
        FA_dict = {}
        LIPIDpd = pd.read_excel(self.filename, sheet_name='LIPIDsyn', usecols="A:H", convert_float=True)
        LIPIDin = self._processBiomassFrame(df=LIPIDpd,in_or_out=-1)

        for coeff_i in range(len(LIPIDin)):
            molpermol = LIPIDin.iloc[coeff_i, 0]
            # MW = LIPIDin.iloc[coeff_i, 1]
            LIPIDin.at[coeff_i, "mmol/gDCW"] = LIPID_wt_per_DCW * molpermol / 3100 * 1000

        LIPIDout=self._processBiomassFrame(df=LIPIDpd,in_or_out=1)
        LIPIDsyn_ref = self._returnSynthesisEquation(df_in=LIPIDin,df_out=LIPIDout)
        LIPIDDict["ref"] = LIPIDsyn_ref


        ####################################
        FApd = pd.read_excel(self.filename, sheet_name='FattyAcidSyn', usecols="A:H", convert_float=True)
        FAin = self._processBiomassFrame(df=FApd,in_or_out=-1)
        FA_list_default = FAin["Reactant"]
        FA_list = [x for x in FA_list_default if
                   x not in ("Sum", "sum", "Total", "total", "Summation", "summation")]
        FA_list.append("Sum")
        mean_list = FAin.iloc[:, 0].tolist()
        mean_list.append(0)
        FA_minmax_pd = pd.DataFrame({"Component": FA_list, "mean": mean_list})
        FA_minmax_pd = FA_minmax_pd.set_index(FA_minmax_pd.columns[0])

        FA_minmax_pd.loc["Sum", "mean"] = FA_minmax_pd["mean"][:-1].sum(axis=0)
        FA_minmax_pd["ref"] = FA_minmax_pd["mean"] / FA_minmax_pd.loc["Sum", "mean"]
        max_cv = 0.25
        FA_minmax_pd["min"] = FA_minmax_pd["ref"] * (1 - max_cv)  ##추가
        FA_minmax_pd["max"] = FA_minmax_pd["ref"] * (1 + max_cv)
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
                        FA_minmax_pd.loc[fa_name, new_col_name] = FA_minmax_pd.loc[fa_name, "ref"] / \
                                                                  (FA_minmax_pd.loc["Sum", "ref"] -
                                                                   FA_minmax_pd.loc[c, "ref"]) * \
                                                                  (1 - FA_minmax_pd.loc[c, min_or_max])
                FA_minmax_pd.loc["Sum", new_col_name] = FA_minmax_pd[new_col_name][:-1].sum(axis=0)

        min_max_col_name = [x for x in FA_minmax_pd.columns.values if any(m in x for m in ("_min", "_max"))]
        min_max_col_name = ["ref"] + min_max_col_name

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
        BIOMASSin = self._processBiomassFrame(df=BIOMASSpd,in_or_out=-1)
        for i in range(len(BIOMASSin)):
            if i > 6:
                BIOMASSin.iloc[i, 4] = BIOMASSin.iloc[i, 4] * norm_fac

        # out
        BIOMASSout = self._processBiomassFrame(df=BIOMASSpd,in_or_out=1)

        # PROTEIN BIOMASS
        BIOMASSsyn = self._returnSynthesisEquation(df_in=BIOMASSin,df_out=BIOMASSout)
        BIOMASSdict = {'ref':BIOMASSsyn}
        return BIOMASSdict


#
# if __name__ == '__main__':
#     fileloc = ''
#     filename=fileloc+'\Scerevisiae test1.xlsx'
#     a = makeBiomass(filename=filename)
#     b = a.Formulate_min_max()
#
