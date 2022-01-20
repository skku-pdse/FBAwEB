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
% %  filename2save          
% % 
% %           
% % OUTPUT
% %  excel file 
% %       Sheetname
% %             flux_table      "rxn_id", "rxn_formula", "subsystems", "resultant flux distributions of each biomass equation"
% %             summary_table   "rxn_id", "rxn_formula", "subsystems", "flux_median", "flux_standard deviation", "flux_max", "flux_min"
% %
% %
% %
function [model_modi1] = FBAwEB(species,model_name,biomass_excel_name,number_of_biomass,starting_point,filename2save);
model_modi = model_name;
% 
%
%
protsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','PROTsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
dnasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','DNAsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
rnasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','RNAsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
lipidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','LIPIDsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);


if species == "CHO"
        carbsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','CARBsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass_cho','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'biomass_cho','DNAsyn','RNAsyn','PROTsyn','biomass_cho_prod','DNAsyn_prod','RNAsyn_prod','PROTsyn_prod'});
        fattyacidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FATTYACIDsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        fattyacudcoasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FATTYACIDCOAsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);

elseif species == "ECOLI"
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass_ecoli','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'Ec_biomass_iML1515_core_75p37M','Ec_biomass_iML1515_WT_75p37M'});     

elseif species == "YEAST"
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass_yeast','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'r_2108','r_4041','r_4047','r_4048','r_4049','r_4050'});
        fattyacidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FATTYACIDsyn','Range',range,'ReadVariableNames',false);
        carbsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','CARBsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
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

    solutions_fba = optimizeCbModel(model_modi1,'max');
    solutions_pfba = optimizeCbModel(model_modi1,'max','one');


%     if isempty(solutions_fba.x) == 0 
%         FBA_overall_name_flux(:,2*p-1) = convertCharsToStrings(model_modi1.rxns);
%         FBA_overall_name_flux(:,2*p) = num2cell(solutions_fba.x);
%     else
%         p
%         FBA_overall_name_flux(:,2*p-1) = convertCharsToStrings(model_modi1.rxns);
%         FBA_overall_name_flux(:,2*p) = zeros();
%     end
%     
    if isempty(solutions_pfba.x) == 0 
        FBA_MIN_overall_name_flux(:,2*p-1) = convertCharsToStrings(model_modi1.rxns);
        FBA_MIN_overall_name_flux(:,2*p) = num2cell(solutions_pfba.x);
        flux_list(:,p) = solutions_pfba.x;
    else
        p
        FBA_MIN_overall_name_flux(:,2*p-1) = convertCharsToStrings(model_modi1.rxns);
        FBA_MIN_overall_name_flux(:,2*p) = zeros();
    end


%     solutions_fba.x = zeros();
    solutions_pfba.x = zeros();
end
save('modelIrrev_cho_etc.mat')
% mean_FBA = num2cell(FBA_overall_name_flux(1:end, 1:end));0H% mean_pFBA = num2cell(FBA_MIN_overall_name_flux(1:end, 1:end));
rxns_name = cellstr(convertCharsToStrings(model_modi1.rxns));
rxns_name_cell = model_modi1.rxns;
SUBSYS = model_modi1.subSystems;
ReactionFormula = printRxnFormula(model_modi1,rxns_name_cell);
FORFORMU = ReactionFormula;
flux_median = median(flux_list,2);
flux_max = max(flux_list,[],2);
flux_min = min(flux_list,[],2);
summary_list = [flux_median flux_max flux_min];

Flux_Distribution = [rxns_name FORFORMU SUBSYS];
t = string(datetime('now','TimeZone','Asia/Seoul','Format','y-MM-d_HH_mm_ss'));
formatSpec = '%s%s%s';
% formatSpec_fba = 'biomass_fba%s.xlsx';
excel_name = sprintf(formatSpec,"FBAwEB_",t,filename2save);
% excel_name_fba = sprintf(formatSpec_fba,t);
sheet_name_formula = 'Flux_Distribution';
% sheet_name1 = 'mean_pFBA';
% sheet_name2 = 'mean_FBA';
sheet_name4 = 'flux_list';
sheet_name_summary = 'SUMMARY';

xlswrite(excel_name,Flux_Distribution,sheet_name_formula);
% xlswrite(excel_name,mean_pFBA,sheet_name1);
% xlswrite(excel_name,mean_FBA,sheet_name2);
xlswrite(excel_name,flux_list,sheet_name4);
xlswrite(excel_name,summary_list,sheet_name_summary);



end
