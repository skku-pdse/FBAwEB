% % function [model_modi1, flux_table, summary_table] = monoSensitivity(species,model_name,biomass_excel_name,filename2save,FVAoption)


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
function [model_modi1, flux_table, summary_table] = monoSensitivity(species,model_name,biomass_excel_name,filename2save,FVAoption)
model_modi = model_name;
% checkup list 
% ecoli AA mean 
%
%create emtpy table of equations
protsyn_eqn_raw = array2table(zeros(0,0));
dnasyn_eqn_raw = array2table(zeros(0,0));
rnasyn_eqn_raw = array2table(zeros(0,0));
lipidsyn_eqn_raw = array2table(zeros(0,0));
carbsyn_eqn_raw = array2table(zeros(0,0));
fattyacidsyn_eqn_raw = array2table(zeros(0,0));
fattyacidcoasyn_eqn_raw = array2table(zeros(0,0));

protsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','PROTsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
dnasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','DNAsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
rnasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','RNAsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
lipidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','LIPIDsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);


if species == "CHO"
        carbsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','CARBsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass_cho','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'biomass_cho','DNAsyn','RNAsyn','PROTsyn','biomass_cho_prod','DNAsyn_prod','RNAsyn_prod','PROTsyn_prod'});
        fattyacidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FATTYACIDsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        fattyacidcoasyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FATTYACIDCOAsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);

elseif species == "ECOLI"
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass_ecoli','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'Ec_biomass_iML1515_core_75p37M','Ec_biomass_iML1515_WT_75p37M'});     

elseif species == "YEAST"
        biomass_eqn_raw = readtable(biomass_excel_name,'Sheet','biomass_yeast','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        model_modi1 = removeRxns(model_modi, {'r_2108','r_4041','r_4047','r_4048','r_4049','r_4050'});
        fattyacidsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','FATTYACIDsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
        carbsyn_eqn_raw = readtable(biomass_excel_name,'Sheet','CARBsyn','Range','B:C','HeaderLines',1,'ReadVariableNames',false,'ReadRowNames',1);
end
%get size of equations 
[row_PROT, ] = size(protsyn_eqn_raw);
[row_DNA, ] = size(dnasyn_eqn_raw);
[row_RNA, ] = size(rnasyn_eqn_raw);
[row_LIPID, ] = size(lipidsyn_eqn_raw);
size_overall = vertcat(row_PROT,row_DNA,row_RNA,row_LIPID);

%allocate reactions
overall_cell_wo_ref = {protsyn_eqn_raw(2:end,1);dnasyn_eqn_raw(2:end,1);rnasyn_eqn_raw(2:end,1);lipidsyn_eqn_raw(2:end,1);[];[];[]}; 
overall_cell_ref = {protsyn_eqn_raw(1,1);dnasyn_eqn_raw(1,1);rnasyn_eqn_raw(1,1);lipidsyn_eqn_raw(1,1);[];[];[]};
% 1: prot;2: dna;3: rna;4:lipid;5: carb;6:FA;7: FACOA
if isempty(carbsyn_eqn_raw)==0
    [row_carb, ] = size(carbsyn_eqn_raw);
    size_overall = vertcat(size_overall,row_carb);
    overall_cell_wo_ref{5,1} = carbsyn_eqn_raw(2:end,1);
    overall_cell_ref{5,1} = carbsyn_eqn_raw(1,1);
end
if isempty(fattyacidsyn_eqn_raw)==0
    [row_FA, ] = size(fattyacidsyn_eqn_raw);
    size_overall = vertcat(size_overall,row_FA);
    overall_cell_wo_ref{6,1} = fattyacidsyn_eqn_raw(2:end,1);
    overall_cell_ref{6,1} = fattyacidsyn_eqn_raw(1,1);
end
if isempty(fattyacidcoasyn_eqn_raw)==0
    [row_FACOA, ] = size(fattyacidcoasyn_eqn_raw);
    size_overall = vertcat(size_overall,row_FACOA);
    overall_cell_wo_ref{7,1} = fattyacidcoasyn_eqn_raw(2:end,1);
    overall_cell_ref{7,1} = fattyacidcoasyn_eqn_raw(1,1);
end
biomass_table = biomass_eqn_raw{1,1}; 
%default upper/lower bound
lowerBound_biomass = 0;
upperBound_biomass = 1000;
% 1: prot;  2: dna;     3: rna; 4:lipid;    5: carb;    6:FA;   7: FACOA
%create empty flux list

[model_modi1] = rxn_ref_add(biomass_table,species,model_modi1,overall_cell_ref,lowerBound_biomass,upperBound_biomass);

rxn_name_ref = {'PROTsyn_modi','DNAsyn_modi','RNAsyn_modi','LIPIDsyn_modi','CARBsyn_modi','FATTYACIDsyn_modi','FATTYACIDCOAsyn_modi'};
flux_list = zeros(size(model_modi1.rxns,1),1);
for p=1:size(overall_cell_wo_ref,1) 
    if p ==7
        
    elseif isempty(overall_cell_wo_ref{p,1}) == 0
        for n = 1:size(overall_cell_wo_ref{p,1},1)

            if species  =="CHO" && (p == 6)
                fattyacidsyn_m_table = strrep(overall_cell_wo_ref{p,1}{n,1}{1,1},"[c]","[m]");
                fattyacidsyn_e_table = strrep(overall_cell_wo_ref{p,1}{n,1}{1,1},"[c]","[e]");
                fattyacidsyn_x_table = strrep(overall_cell_wo_ref{p,1}{n,1}{1,1},"[c]","[x]");
                model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_modi','reactionFormula',overall_cell_wo_ref{p,1}{n,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
                model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_m_modi','reactionFormula',fattyacidsyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
                model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_x_modi','reactionFormula',fattyacidsyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
                model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_e_modi','reactionFormula',fattyacidsyn_e_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});

                fattyacidcoasyn_m_table = strrep(overall_cell_wo_ref{p,1}{n,1}{1,1},"[c]","[m]");
                fattyacidcoasyn_x_table = strrep(overall_cell_wo_ref{p,1}{n,1}{1,1},"[c]","[x]");
                model_modi1 = addReaction(model_modi1,'FATTYACIDCOAsyn_modi','reactionFormula',overall_cell_wo_ref{p,1}{n,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
                model_modi1 = addReaction(model_modi1,'FATTYACIDCOAsyn_m_modi','reactionFormula',fattyacidcoasyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
                model_modi1 = addReaction(model_modi1,'FATTYACIDCOAsyn_x_modi','reactionFormula',fattyacidcoasyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});               
            else
                model_modi1 = addReaction(model_modi1,rxn_name_ref{1,p},'reactionFormula',overall_cell_wo_ref{p,1}{n,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
    %         checkObjective(model_modi1);
            end    
            [model_modi1 ,flux_list] = doFBA(flux_list,n,model_modi1,FVAoption);
            
        end
        
        if p==1
            
            var_name_ref = overall_cell_wo_ref{p,1}.Properties.RowNames;
            flux_table = array2table(flux_list);
            flux_table.Properties.VariableNames = var_name_ref;
        else
            var_name_ref = overall_cell_wo_ref{p,1}.Properties.RowNames;
            tmp_table = array2table(flux_list);
            tmp_table.Properties.VariableNames = var_name_ref;
            flux_table = horzcat(flux_table,tmp_table);

        end
    else
        
    end
    flux_list = zeros(size(model_modi1.rxns,1),1);
end
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

excel_name = sprintf(formatSpec,"MONO ",t,filename2save);


writetable(flux_table,excel_name,'FileType','spreadsheet','Sheet','flux_table','WriteRowNames',true,'WriteVariableNames',true,'WriteMode','append','UseExcel',true);
writetable(summary_table,excel_name,'FileType','spreadsheet','Sheet','summary_table','WriteRowNames',true,'WriteVariableNames',true,'WriteMode','append','UseExcel',true);


end

function [model_modi1] = rxn_ref_add(biomass_table,species,model_modi1,overall_cell_ref,lowerBound_biomass,upperBound_biomass)
 % add reference reaction
model_modi1 = addReaction(model_modi1,'PROTsyn_modi','reactionFormula',overall_cell_ref{1,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
model_modi1 = addReaction(model_modi1,'DNAsyn_modi','reactionFormula',overall_cell_ref{2,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
model_modi1 = addReaction(model_modi1,'RNAsyn_modi','reactionFormula',overall_cell_ref{3,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
model_modi1 = addReaction(model_modi1,'LIPIDsyn_modi','reactionFormula',overall_cell_ref{4,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
if isempty(overall_cell_ref{5,1})==0
    model_modi1 = addReaction(model_modi1,'CARBsyn_modi','reactionFormula',overall_cell_ref{5,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
end
if isempty(overall_cell_ref{6,1})==0
        if species == "YEAST"
            model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_modi','reactionFormula',overall_cell_ref{6,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
        elseif species == "CHO"
            fattyacidsyn_m_table = strrep(overall_cell_ref{6,1}{1,1}{1,1},"[c]","[m]");
            fattyacidsyn_e_table = strrep(overall_cell_ref{6,1}{1,1}{1,1},"[c]","[e]");
            fattyacidsyn_x_table = strrep(overall_cell_ref{6,1}{1,1}{1,1},"[c]","[x]");
            model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_modi','reactionFormula',overall_cell_ref{6,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_m_modi','reactionFormula',fattyacidsyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_x_modi','reactionFormula',fattyacidsyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDsyn_e_modi','reactionFormula',fattyacidsyn_e_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
        end
end
if isempty(overall_cell_ref{7,1})==0
        if species == "CHO"
            fattyacidcoasyn_m_table = strrep(overall_cell_ref{7,1}{1,1}{1,1},"[c]","[m]");
            fattyacidcoasyn_x_table = strrep(overall_cell_ref{7,1}{1,1}{1,1},"[c]","[x]");
            model_modi1 = addReaction(model_modi1,'FATTYACIDCOAsyn_modi','reactionFormula',overall_cell_ref{7,1}{1,1}{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDCOAsyn_m_modi','reactionFormula',fattyacidcoasyn_m_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
            model_modi1 = addReaction(model_modi1,'FATTYACIDCOAsyn_x_modi','reactionFormula',fattyacidcoasyn_x_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
        end
end
if species == "YEAST"
    model_modi1 = addReaction(model_modi1,'biomass_yeast','reactionFormula',biomass_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
    model_modi1 = changeObjective(model_modi1, 'biomass_yeast');

elseif species == "CHO"
    model_modi1 = addReaction(model_modi1,'biomass_cho','reactionFormula',biomass_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
    model_modi1 = changeObjective(model_modi1, 'biomass_cho');
elseif species == "ECOLI"
    model_modi1 = addReaction(model_modi1,'biomass_ecoli','reactionFormula',biomass_table{1,1},'reversible',false,'lowerBound',lowerBound_biomass,'upperBound',upperBound_biomass,'subSystem',{'FBAwEB'});
    model_modi1 = changeObjective(model_modi1, 'biomass_ecoli');
end
end

function [model_modi1 ,flux_list] = doFBA(flux_list,n,model_modi1,FVAoption)


solutions_pfba = optimizeCbModel(model_modi1,'max','one');
if FVAoption == true
    [minFlux, maxFlux] = fluxVariability(model_modi1, 'method', '1-norm');
    min_max = minus(maxFlux,minFlux);
    FVA_max_min(:,2*n-1) = convertCharsToStrings(model_modi1.rxns);

    FVA_max_min(:,2*n) = num2cell(min_max);
    max_FVA(:,n) = maxFlux;
    min_FVA(:,n) = minFlux;
    FOR_PCA(:,n) = min_max;
end


if isempty(solutions_pfba.x) == 0
    FBA_MIN_overall_name_flux(:,2*n-1) = convertCharsToStrings(model_modi1.rxns);
    FBA_MIN_overall_name_flux(:,2*n) = num2cell(solutions_pfba.x);
    flux_list(:,n) = solutions_pfba.x;
else
    flux_list(1,n) = "infeasible";

end

solutions_pfba.x = zeros();
end
