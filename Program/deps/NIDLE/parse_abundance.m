function [abundance, g_vect, report_tab]=parse_abundance(model_irrev,conditions,abundance_dat) 
%Function to select and map homomeric enzyme reactions to the protein
%abundance measurement
%INPUT: 
%   struc model_irrev:  irreversible iCre1355 GEM with GPR rules bearing the
%                       Proteins names found in abundance_dat
%   struc conditions:   cond.cond has be a cell array of character vectors
%                       that give the different conditions same as the
%                       columns in abundance_dat
%   tab abundance_dat:  table in which rows correspond to protein names
%                       present in the GPR rules of model_irrev and columns
%                       correspond to conditions given in conditions.cond
%OUTPUT:
%   stuc abundance
%       -abun:  array of cellular protein concentrations for normalized
%               with the cellular dry weight to mmol/gDW (size: n x m)
%       -genes: row names of abun given the proteins from protein_abun that
%               are present in the GPR rules of irrev_model (size: n x 1)
%       -cond:  column names of abun, giving the experimental condition in
%               which protein abundance was measured (size: m x 1 )
%       -GPR:   Homomeric and isoenzyme GPR rules in which all proteins
%               have abundance data in at least 1 condition (size: k x 1)
%       -reac:  reaction names linked to each entry in GPR size(size: k x
%               1)
%       -reacind:   index of reac in the vector of all reaction in
%                   model_irrev (size k x 1)
%   double vect:    A matrix (:D) of n_reaction x n_conditions vector containing 1
%                   for each reaction, for which all nescessary enzymes are found in
%                   abundance data in the respective condition
%   tab report_tab: A table containing the reaction name and subsystem
%                   together with the vect matrix
%Calculate cellular dw assuming constant cell density
meta_2022=readtable('Data/QconCAT_David20220124/QconCat_UPS.xlsx');
meta_2022=meta_2022(:, {'Characteristics_environmentalCondition_','Characteristics_cellVolume_'});

cond_map=containers.Map({'Control Type', 'dark', 'high cell density','no shaking','sodium chloride exposure', 'warm/hot temperature exposure'} ,...
    {'control', 'dark', 'highcell', 'noshaking', 'highsalt', 'hightemp'});
%aggregate cell volumen to mean value in place
for i=1:(size(meta_2022,1)/3)
    meta_2022.Characteristics_environmentalCondition_{i}=cond_map(meta_2022.Characteristics_environmentalCondition_{i});
    meta_2022.Characteristics_cellVolume_(i)=mean(meta_2022.Characteristics_cellVolume_(i:(i+2)));
    meta_2022((i+1):(i+2),:)=[];
end
%Using mass of stationary mixotrophic cell calculated from yang et al. ,
%Biotechnology for biofuels (2018)
dens=110000/meta_2022.Characteristics_cellVolume_(ismember(meta_2022.Characteristics_environmentalCondition_, 'control'));
meta_2022.cellDW=meta_2022.Characteristics_cellVolume_*dens;
%for mixotrophic growth data from 2021 assumen 4800 gDW
meta_2022=[meta_2022;
    {'UVM4', NaN, 110000;
    'Stop1', NaN, 110000;
    'Stop2', NaN, 110000}];
%sort 
%generate condition wise abundance 
    abundance.genes=abundance_dat.ProteinId;
    abundance.abun=zeros(size(abundance_dat, 1), length(conditions.cond));
    abundance.cond=cell(length(conditions.cond),1);
    for i=1:length(conditions.cond)
        abundance_con=table2array(abundance_dat(:,~cellfun(@isempty, strfind(abundance_dat.Properties.VariableNames, conditions.cond(i)))));
        abundance_con=median(abundance_con, 2);
        abundance_con=abundance_con/meta_2022.cellDW(ismember(meta_2022.Characteristics_environmentalCondition_, conditions.cond(i))); %devide by cell weight in fgDW/cell to get mmol/gDW
        abundance.abun(:,i)=abundance_con;
        abundance.cond(i)=conditions.cond(i);
    end

    [rules_type, g_vect, report_tab]=parse_rules(model_irrev, abundance);
    

%     %homomeric enzymes in model
%     homomeric_ind_model=find(rules_type==2);
%     homomeric_reac_model=model_irrev.rxns(homomeric_ind_model);
%     homomeric_en_model=cellfun(@(x) x{1}, findGenesFromRxns(model_irrev, homomeric_reac_model), 'UniformOutput', false); %this preserves the order of genes accoding to homomeric_reac_model
% 
%     %isomeric enzymes
%     homomeric_en_model=cellfun(@(x) x{1}, findGenesFromRxns(model_irrev, homomeric_reac_model), 'UniformOutput', false); %this preserves the order of genes accoding to homomeric_reac_model

    ind_model=find(ismember(rules_type, [2,3]));
    ind_reac_model=model_irrev.rxns(ind_model);
    homomeric_en_model=findGenesFromRxns(model_irrev, ind_reac_model); %this preserves the order of genes accoding to homomeric_reac_model
    
    %check if homomeric reactions indeed only have 1 enzyme assigned
    if any(cellfun(@length, homomeric_en_model(rules_type(ind_model)==2))>1)
        error('Parsing of GPR rules leads to homomeric catalyzed reactions with more than one gene assigned. Fix gene rule parsing')
    end
    
    keep_abun=zeros(length(abundance.genes),1);
    keep_reac=zeros(length(ind_model), 1);
    %find these enzymes in the abundance file
    for i=1:length(homomeric_en_model)
        if any(ismember(homomeric_en_model{i}, abundance.genes))
            %mark genes to retain them in abundance data
            keep_abun(ismember(abundance.genes, homomeric_en_model{i}))=1;
            %mark reaction to retained in the reaction list
            keep_reac(i)=1;
            %remove genes that are not in the abundance data
            homomeric_en_model{i}=homomeric_en_model{i}(ismember(homomeric_en_model{i}, abundance.genes));
        end
    end
    abundance.abun=abundance.abun(logical(keep_abun), :);
    abundance.genes=abundance.genes(logical(keep_abun));
    
    abundance.reac=ind_reac_model(logical(keep_reac));
    abundance.reacind=ind_model(logical(keep_reac));
    abundance.GPR=homomeric_en_model(logical(keep_reac));
    %check if GPR rules indeed only contain expressed enzymes
    if any(cellfun(@(x) ~all(ismember(x, abundance.genes)), abundance.GPR))
        error('Filtering of GPR rules for expressed proteins is broken...')
    end
end