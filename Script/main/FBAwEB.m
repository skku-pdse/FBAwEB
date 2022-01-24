% % function [model_modi1] = FBAwEB(species,model_name,biomass_excel_name,number_of_biomass,starting_point,filename2save,FVAoption);
%%%%%%%%%%%%% species = applicable to  "ECOLI", "YEAST", "CHO" %%%%%%%%%%%%%
% The script should be modified when applied to other models.
% % INPUT
% %  species                 name of utilizing model's species, 
% %                          in this version CHO GEM(iCHO1766/iCHO2291), ECOLI(iML1515) and YEAST(yeastGEM_v8.0.0) are available 
% % 
% %  model_name              name of '.mat' format model          
% %                    
% %  biomass_excel_name      name of biomass excel file;
% %                          import biomass (ensemble biomass) excel file with '.xlsx'
% %                          (i.e., Ecoli_mono_+-25% biomass_Jan12% 02;04.xlsx)
% %  
% %  number_of_biomass      the number of sampled biomass equations that you want to implement FBAwEB 
% %	 
% %  starting_point         if implementation failed during MATLAB run, you can restart from the failure point
% %  
% %  filename2save          name of exporting excel file 
% %                         ("FBAwEB_datetime+filenae2save.xlsx" i.e., "FBAwEB_2020-01-01-00:00filename2save.xlsx")
% %           
% % OUTPUT
% %  excel file 
% %       Sheetname
% %             flux_table      "rxn_id", "rxn_formula", "subsystems", "resultant flux distributions of each biomass equation"
% %             summary_table   "rxn_id", "rxn_formula", "subsystems", "flux_median", "flux_standard deviation", "flux_max", "flux_min"
% %
% %
% %
%%
function [model_modi1, flux_table, summary_table] = FBAwEB(species,model_name,biomass_excel_name,number_of_biomass,starting_point,filename2save)
model_modi = model_name;

protsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','PROTsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
dnasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','DNAsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
rnasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','RNAsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
lipidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','LIPIDsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);


if species == "CHO"
        carbsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','CARBsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'biomass_cho','DNAsyn','RNAsyn','PROTsyn','biomass_cho_prod','DNAsyn_prod','RNAsyn_prod','PROTsyn_prod'});
        fattyacidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FAsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        fattyacudcoasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FACOAsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);

elseif species == "ECOLI"
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'Ec_biomass_iML1515_core_75p37M','Ec_biomass_iML1515_WT_75p37M'});     

elseif species == "YEAST"
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'r_2108','r_4041','r_4047','r_4048','r_4049','r_4050'});
        fattyacidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FAsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        carbsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','CARBsyn','Range','A:B','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
end

lowerBound_biomass = 0;
upperBound_biomass = 1000;

for p=starting_point+1:starting_point+number_of_biomass
    protsyn_table = protsyn_eqn_raw{p,1};
    lipidsyn_table = lipidsyn_eqn_raw{p,1};
    rnasyn_table = rnasyn_eqn_raw{p,1};
    dnasyn_table = dnasyn_eqn_raw{p,1};
    biomass_table = biomass_eqn_raw{p,1};
    
    model_modi1 = addReaction(model_modi1,'DNAsyn_modi','reactionFormula',dnasyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'RNAsyn_modi','reactionFormula',rnasyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'PROTsyn_modi','reactionFormula',protsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'LIPIDsyn_modi','reactionFormula',lipidsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    
  
          
    if species == "YEAST"
        carbsyn_table = carbsyn_eqn_raw{p,1};
        model_modi1 = addReaction(model_modi1,'CARBsyn_modi','reactionFormula',carbsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        fattyacidsyn_table = fattyacidsyn_eqn_raw{p,1};
        model_modi1 = addReaction(model_modi1,'FATTYACIDSyn','reactionFormula',fattyacidsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'biomass_yeast','reactionFormula',biomass_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = changeObjective(model_modi1, 'biomass_yeast');
        
    elseif species == "CHO"
        carbsyn_table = carbsyn_eqn_raw{p,1};
        model_modi1 = addReaction(model_modi1,'CARBsyn_modi','reactionFormula',carbsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        fattyacidsyn_table = fattyacidsyn_eqn_raw{p,1};
        fattyacudcoasyn_table = fattyacudcoasyn_eqn_raw{p,1};
        fattyacidsyn_m_table = strrep(fattyacidsyn_table,"[c]","[m]");
        fattyacidsyn_e_table = strrep(fattyacidsyn_table,"[c]","[e]");
        fattyacidsyn_x_table = strrep(fattyacidsyn_table,"[c]","[x]");
        fattyacidcoasyn_m_table = strrep(fattyacudcoasyn_table,"[c]","[m]");
        fattyacidcoasyn_x_table = strrep(fattyacudcoasyn_table,"[c]","[x]");
        model_modi1 = addReaction(model_modi1,'FATTYACIDSyn','reactionFormula',fattyacidsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'FATTYACUDCOASyn','reactionFormula',fattyacudcoasyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'FATTYACIDSyn_m','reactionFormula',fattyacidsyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'FATTYACIDSyn_x','reactionFormula',fattyacidsyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'FATTYACIDSyn_e','reactionFormula',fattyacidsyn_e_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'FATTYACUDCOASyn_m','reactionFormula',fattyacidcoasyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'FATTYACUDCOASyn_x','reactionFormula',fattyacidcoasyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = addReaction(model_modi1,'biomass_cho','reactionFormula',biomass_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = changeObjective(model_modi1, 'biomass_cho');
    elseif species == "ECOLI"
        model_modi1 = addReaction(model_modi1,'biomass_ecoli','reactionFormula',biomass_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
        model_modi1 = changeObjective(model_modi1, 'biomass_ecoli');
    end

    % check objective function
    checkObjective(model_modi1);
    if p ==1
        flux_list = zeros(size(model_modi1.rxns,1),1);
        eqn_list = cell(1,number_of_biomass);
    else
    end

    solutions_pfba = optimizeCbModel(model_modi1,'max','one');

    if isempty(solutions_pfba.x) == 0 
        FBA_MIN_overall_name_flux(:,2*p-1) = convertCharsToStrings(model_modi1.rxns);
        FBA_MIN_overall_name_flux(:,2*p) = num2cell(solutions_pfba.x);
        flux_list(:,p) = solutions_pfba.x;
        eqn_list {1,p} = biomass_table{1,1};
    else
        p
        FBA_MIN_overall_name_flux(:,2*p-1) = convertCharsToStrings(model_modi1.rxns);
        FBA_MIN_overall_name_flux(:,2*p) = zeros();
        eqn_list {1,p} = biomass_table{1,1};     
    end



    solutions_pfba.x = zeros();
end

flux_table = array2table(flux_list);
flux_table.Properties.VariableNames = biomass_eqn_raw.Properties.RowNames(starting_point+1:number_of_biomass,1);


flux_table.Properties.RowNames = model_modi1.rxns;
rxns_name_cell = model_modi1.rxns;
SUBSYS = model_modi1.subSystems;
ReactionFormula = printRxnFormula(model_modi1,'rxnAbbrList',rxns_name_cell,'printFlag',false);
FORFORMU = ReactionFormula;
flux_median = median(flux_table{:,:},2);
flux_std = std(flux_table{:,:},[],2);
flux_max = max(flux_table{:,:},[],2);
flux_min = min(flux_table{:,:},[],2);
summary_list = [flux_median flux_std flux_max flux_min];
summary_list = array2table(summary_list);
summary_list.Properties.RowNames = model_modi1.rxns;
summary_list.Properties.VariableNames = ["flux_median","flux_std","flux_max","flux_min"];
Flux_Distribution = [FORFORMU SUBSYS];
Flux_Distribution_table = cell2table(Flux_Distribution);
Flux_Distribution_table.Properties.VariableNames = ["RxnFormula","subSystems"];
flux_table = horzcat(Flux_Distribution_table,flux_table);
summary_table = horzcat(Flux_Distribution_table,summary_list);
t = string(datetime('now','TimeZone','Asia/Seoul','Format','y-MM-d_HH_mm_ss'));
formatSpec = "%s%s%s";

excel_name = sprintf(formatSpec,"FBAwEB_ ",t,filename2save);


writetable(flux_table,excel_name,'FileType','spreadsheet','Sheet','flux_table','WriteRowNames',true,'WriteVariableNames',true,'WriteMode','append','UseExcel',true);
writetable(summary_table,excel_name,'FileType','spreadsheet','Sheet','summary_table','WriteRowNames',true,'WriteVariableNames',true,'WriteMode','append','UseExcel',true);


end
