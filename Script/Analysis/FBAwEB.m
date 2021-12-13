function [model_modi1] = FBAwEB(species,model_name,biomass_excel_name,number_of_biomass,starting_point,name,FVAoption);
%%%%%%%%%%%%% species = available to "CHO", "ECOLI", "YEAST" %%%%%%%%%%%%%

model_modi = model_name;
range = 'B2:B5001';
protsyn_eqn = readtable(biomass_excel_name,'Sheet','PROTsyn','Range',range,'ReadVariableNames',false);
dnasyn_eqn = readtable(biomass_excel_name,'Sheet','DNAsyn','Range',range,'ReadVariableNames',false);
rnasyn_eqn = readtable(biomass_excel_name,'Sheet','RNAsyn','Range',range,'ReadVariableNames',false);
lipidsyn_eqn = readtable(biomass_excel_name,'Sheet','LIPIDsyn','Range',range,'ReadVariableNames',false);
carbsyn_eqn = readtable(biomass_excel_name,'Sheet','CARBsyn','Range',range,'ReadVariableNames',false);



if species == "CHO"
        biomass_eqn = readtable(biomass_excel_name,'Sheet','biomass','Range',range,'ReadVariableNames',false);
        model_modi1 = removeRxns(model_modi, {'biomass_cho','DNAsyn','RNAsyn','PROTsyn','biomass_cho_prod','DNAsyn_prod','RNAsyn_prod','PROTsyn_prod'});
        fattyacidsyn_eqn = readtable(biomass_excel_name,'Sheet','FAsyn','Range',range,'ReadVariableNames',false);
        fattyacudcoasyn_eqn = readtable(biomass_excel_name,'Sheet','FACOAsyn','Range',range,'ReadVariableNames',false);
elseif species == "ECOLI"
        biomass_eqn = readtable(biomass_excel_name,'Sheet','biomass','Range',range,'ReadVariableNames',false);
        model_modi1 = removeRxns(model_modi, {'Ec_biomass_iML1515_core_75p37M','Ec_biomass_iML1515_WT_75p37M'});
        
elseif species == "YEAST"
        biomass_eqn = readtable(biomass_excel_name,'Sheet','biomass','Range',range,'ReadVariableNames',false);
        model_modi1 = removeRxns(model_modi, {'r_2108','r_4041','r_4047','r_4048','r_4049','r_4050'});
        fattyacidsyn_eqn = readtable(biomass_excel_name,'Sheet','FAsyn','Range',range,'ReadVariableNames',false);

end

lowerBound_biomass = 0;
upperBound_biomass = 1000;

for p=starting_point:number_of_biomass
    protsyn_table = protsyn_eqn{p,1};
    lipidsyn_table = lipidsyn_eqn{p,1};
    rnasyn_table = rnasyn_eqn{p,1};
    dnasyn_table = dnasyn_eqn{p,1};
    biomass_table = biomass_eqn{p,1};
    carbsyn_table = carbsyn_eqn{p,1};

   
    model_modi1 = addReaction(model_modi1,'DNAsyn_modi','reactionFormula',dnasyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'RNAsyn_modi','reactionFormula',rnasyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'PROTsyn_modi','reactionFormula',protsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'LIPIDsyn_modi','reactionFormula',lipidsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'biomass_modi','reactionFormula',biomass_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
    model_modi1 = addReaction(model_modi1,'CARBsyn_modi','reactionFormula',carbsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});

    if species == "YEAST"
            fattyacidsyn_table = fattyacidsyn_eqn{p,1};
            model_modi1 = addReaction(model_modi1,'FATTYACIDSyn','reactionFormula',fattyacidsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});

    elseif species == "CHO"
            fattyacidsyn_table = fattyacidsyn_eqn{p,1};
            fattyacudcoasyn_table = fattyacudcoasyn_eqn{p,1};
            fattyacidsyn_m_table = strrep(fattyacidsyn_table,"[c]","[m]");
            fattyacidsyn_e_table = strrep(fattyacidsyn_table,"[c]","[e]");
            fattyacidsyn_x_table = strrep(fattyacidsyn_table,"[c]","[x]");
            fattyacudcoasyn_m_table = strrep(fattyacudcoasyn_table,"[c]","[m]");
            fattyacudcoasyn_x_table = strrep(fattyacudcoasyn_table,"[c]","[x]");
            model_modi1 = addReaction(model_modi1,'FATTYACIDSyn','reactionFormula',fattyacidsyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
            model_modi1 = addReaction(model_modi1,'FATTYACUDCOASyn','reactionFormula',fattyacudcoasyn_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDSyn_m','reactionFormula',fattyacidsyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDSyn_x','reactionFormula',fattyacidsyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDSyn_e','reactionFormula',fattyacidsyn_e_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
            model_modi1 = addReaction(model_modi1,'FATTYACUDCOASyn_m','reactionFormula',fattyacudcoasyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});
            model_modi1 = addReaction(model_modi1,'FATTYACUDCOASyn_x','reactionFormula',fattyacudcoasyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'Biomass synthesis'});


    end

    
    model_modi1 = changeObjective(model_modi1, 'biomass_modi');
    checkObjective(model_modi1);

    solutions_fba = optimizeCbModel(model_modi1,'max');
    solutions_pfba = optimizeCbModel(model_modi1,'max','one');
    if FVAoption == true
        [minFlux, maxFlux] = fluxVariability(model_modi1, 'method', '1-norm');
        min_max = minus(maxFlux,minFlux);
        FVA_max_min(:,2*p-1) = convertCharsToStrings(model_modi1.rxns);
 
        FVA_max_min(:,2*p) = num2cell(min_max);
        max_FVA(:,p) = maxFlux;
        min_FVA(:,p) = minFlux;
        FOR_PCA(:,p) = min_max;
    end

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
formatSpec = string(name) + 'cho%s.xlsx';
% formatSpec_fba = 'biomass_fba%s.xlsx';
excel_name = sprintf(formatSpec,t);
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



    if FVAoption == true
        FVA_max_min = num2cell(FVA_max_min(1:end, 1:end));
        sheet_name3 = 'FVA_max--min';
        sheet_name5 = 'maxflux';
        sheet_name6 = 'minflux';
        sheet_name7 = 'for_PCA';
        xlswrite(excel_name,FVA_max_min,sheet_name3);
        xlswrite(excel_name,max_FVA,sheet_name5);
        xlswrite(excel_name,min_FVA,sheet_name6);
        xlswrite(excel_name,FOR_PCA,sheet_name7);
    end
end
