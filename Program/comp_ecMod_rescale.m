%Code to test compare simulations between iCre1355 and sMOMENT smCre1355 and
%GECKO geCre1355
%% Inti Cobra toolbox
changeCobraSolver('gurobi', 'all')
changeCobraSolverParams('LP','feasTol', 1e-9)
if ~isfolder('Results/eccomp/abunvse')
    mkdir('Results/eccomp/abunvse')
end
if ~isfolder('Results/NIDLEecmods/')
    mkdir('Results/NIDLEecmods/')
end

%import protein biomass fraction from Boyle and Morgan 
Ptot={'auto', 0.261
    'mixo', 0.303
    'hetero', 0.222};
Ptot=cell2table(Ptot, 'VariableNames', {'name', 'Ptot'});
cstat_dat=readtable('Data/Cre1355/chemostat_dat.txt', 'TreatAsEmpty', {'NA'});
%add culture conditions id
cstat_dat.cond_id={'auto', 'auto', 'auto', 'auto', 'mixo', 'mixo', 'mixo', 'mixo', 'hetero', 'hetero'}';
%adapt reaction names as column names
cstat_dat.Properties.VariableNames(4:7)={'EX_co2_e', 'EX_nh4_e', 'EX_pi_e', 'EX_ac_e'};



%load models into cell structure, first two columns condition and model
%name 3rd column original model, 4th column unadapted geckomodel 5th column
%adapted gecko model, 6th column unadapted GECKO model parameterized with
%NIDLE results 7th column adapted GECKO model parameterized with NIDLE
gems=cell(5,7);
gems(:,1:2)=table2cell(cstat_dat(1:2:10, {'cond_id', 'Var1'}));
for i=1:size(gems,1)
    gems{i,3}=scaleBiomassCre(readCbModel(append('Data/Cre1355/iCre1355_', gems{i,1}, '_up.xml')), append('Biomass_Chlamy_', gems{i,1}), Ptot.Ptot(ismember(Ptot.name, gems{i,1})));
    %fix reactions in fba model the reactions in GECKO models have been
    %prefixed
    fix_r=cstat_dat.Properties.VariableNames(4:7);

    fix_r=fix_r(~isnan(table2array(cstat_dat(i*2-1, 4:7))));
    for r=fix_r
        if cstat_dat{i, r}<0 
            error('negative uptake rate detected, code is not adapted for secretion. Change setting of Exchange reaction bounds')
        end
        gems{i,3}=changeRxnBounds(gems{i,3}, r, -1*cstat_dat{i*2-1, r}, 'l');
        disp(append('Fixing reaction for', r, ' to ', num2str(-1*cstat_dat{i*2-1, r})))
    end
    gems{i,4}=readGKOmodel(append('Program/deps/GECKOcre/models/geCre1355', gems{i,2}, '/rawgeCre1355',gems{i,2}, '_batch.mat'));
    gems{i,5}=readGKOmodel(append('Program/deps/GECKOcre/models/geCre1355', gems{i,2}, '/adprawgeCre1355',gems{i,2}, '_batch.mat'));
    gems{i,6}=create_NIDLEmod(gems{i,4}, 'prot_', gems{i,1}, {'control', 'highcell', 'highsalt', 'hightemp', 'noshaking', 'UVM4', 'Stop1', 'Stop2', 'dark'}, {}, false);
    gems{i,7}=create_NIDLEmod(gems{i,5}, 'prot_', gems{i,1}, {'control', 'highcell', 'highsalt', 'hightemp', 'noshaking', 'UVM4', 'Stop1', 'Stop2', 'dark'}, {}, false);
    %Export the adapted GECKO model with NIDLE kcats for easier access
    writeCbModel(gems{i,7}, ['Results/NIDLEecmods/NDLadpraw_', gems{i,2} '.xml'])
end

%predict growth
res=nan(size(gems,1), size(gems, 2)-1);
res(:,1)=cstat_dat.mu_h_1_(1:2:10);
for i=1:size(res,1)
    for j=2:size(res,2)
        FBAsol=optimizeCbModel(gems{i,j+1});
        res(i,j)=FBAsol.f;
    end
end
res=[cstat_dat.Var1(1:2:10), array2table(res, 'VariableNames', {'exp_mu', 'FBA', 'GKOraw', 'GKOadp', 'NDLraw', 'NDLadp'})];
writetable(res, 'Results/eccomp/chemostat_comp.txt')
RMSE_FBA=sqrt(mean((res.exp_mu-res.FBA).^2))
RMSE_NDLadp=sqrt(mean((res.exp_mu-res.NDLadp).^2))
disp(['RMSE FBA is ', num2str(RMSE_FBA), '. RMSE NIDLE adapted is ' num2str(RMSE_NDLadp)])

%compare predicted enzyme usage
%import enzyme abundance 
abun_dat=readtable('Data/QconCAT_David20220124/abs_abundance/med_abun_all.tsv', 'FileType', 'text');
%remove non control conditions
abun_dat(:, {'highcell', 'highsalt', 'hightemp', 'noshaking'})=[];
%adapt Uniprot names
abun_dat.UPId=JGItoUP(abun_dat.ProteinId);
abun_dat(cellfun(@isempty, abun_dat.UPId), :)=[];
abun_dat=abun_dat(:, [1, size(abun_dat,2), 2:(size(abun_dat,2)-1)]);
%sum duplicated entries up and remove duplicates
rem_idx=[];
dups=duplicates(abun_dat.UPId);
for i=find(sum(dups,2))'
    js=find(dups(i,:));
    tmp=sum(table2array(abun_dat([i js], 3:end)), 1, 'omitnan');
    tmp(tmp==0)=NaN;
    abun_dat(i, 3:end)=array2table(tmp);
    rem_idx=[rem_idx, js];
end
abun_dat(rem_idx,:)=[];
%import flux info
fl_inf=readtable('Data/QconCAT_David20220124/fitted_acup.tsv', 'FileType', 'text');
[match, idx]=ismember(abun_dat.Properties.VariableNames(3:end), fl_inf.Condition);
%add additional info per condition %dry weight density for dark condition
%is calculated from cell volume in QCONCAT 2022 data as in
%parse_abundance.m function of nidle

mod_inf=table(abun_dat.Properties.VariableNames(3:end)', {'mixo', 'hetero', 'mixo', 'mixo', 'mixo'}', ...
    {48000, 28265, 48000, 48000, 48000}', fl_inf.mu(idx(match)),'VariableNames', {'cond', 'cond_id', 'DWdens', 'mu'});
%convert amol/cell to mmol/gDW using the specific dry weight of the
%condition
abun_dat(:, 3:end)=array2table(table2array(abun_dat(:, 3:end))./cell2mat(mod_inf.DWdens)');
%initialize result data objects
%fluxes and corresponding abundancies
res=cell(size(mod_inf, 1), 4);
%spearman correlation
spear=nan(size(mod_inf, 1), 4);
%RMSE on log scale
log_RMSE=nan(size(mod_inf, 1), 4);

%column names 
col_names={'rawGKO', 'adpGKO', 'rawGKONDL', 'adpGKONDL'};
for i=1:size(mod_inf, 1)
    
    mod_idx=find(ismember(gems(:,1), mod_inf.cond_id(i)));
    %get GECKO ecmodels of respective condition 
    %1st raw batch 2nd adapted batch 3rd NIDLE augmented raw batch 4th
    %NIDLE augmented adapted batch
    ecmodel=[gems(mod_idx(1),4:5), {create_NIDLEmod(gems{mod_idx(1),4}, 'prot_', mod_inf.cond_id{i}, ...
        {'control', 'highcell', 'highsalt', 'hightemp', 'noshaking', 'UVM4', 'Stop1', 'Stop2', 'dark'}, mod_inf.cond{i}, false), ...
        create_NIDLEmod(gems{mod_idx(1),5}, 'prot_', mod_inf.cond_id{i}, ...
        {'control', 'highcell', 'highsalt', 'hightemp', 'noshaking', 'UVM4', 'Stop1', 'Stop2', 'dark'}, mod_inf.cond{i}, false)}];
    if ~isequal(ecmodel{1}.rxns, ecmodel{2}.rxns)
        error('batch model and adapted model of same condition have different reactions. This is an internal error... :(')
    end
    %Since models have different reactions compile an extra table for each
    %condition
    dat=abun_dat(~isnan(table2array(abun_dat(:,(i+2)))), [2 (i+2)]);
    [match, idx]=ismember(dat.UPId, erase(ecmodel{1}.rxns, 'draw_prot_'));
    for j=1:length(ecmodel)
        ecmodel{j}=changeRxnBounds(ecmodel{j}, {['Biomass_Chlamy_' mod_inf.cond_id{i}]}, mod_inf.mu(i), 'b');
        ecmodel{j}=changeRxnBounds(ecmodel{j}, {'prot_pool_exchange'}, 1000, 'u');
        ecmodel{j}=changeObjective(ecmodel{j}, {'prot_pool_exchange'}, -1);
        pred=optimizeCbModel(ecmodel{j});
        tmpres=[dat(ismember(dat.UPId, erase(ecmodel{j}.rxns, 'draw_prot_')),:), table(pred.v(idx(match)))];
        res{i,j}=tmpres;
        plotdat=table2array(tmpres(:, 2:end));
        fidx=plotdat(:,2)>0;
        n=sum(fidx);
        spear(i,j)=corr(log10(plotdat(fidx,1)), log10(plotdat(fidx,2)), 'type', 'Spearman');
        log_RMSE(i,j)=sqrt(mean((log10(plotdat(fidx,1))-log10(plotdat(fidx,2))).^2));
        figure
        scatter(plotdat(fidx, 1), plotdat(fidx,2))
        title(['n = ' num2str(n)])
        set(gca, 'xscale', 'log')
        set(gca, 'yscale', 'log')
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0 0 8 8]);
        saveas(gcf, fullfile('Results', 'eccomp', 'abunvse', [mod_inf.cond{i} '_' col_names{j} '.svg']))
        close
        writetable(array2table(plotdat, "VariableNames", {'abs_prot', 'enz_usg'}), fullfile('Results', 'eccomp', 'abunvse', [mod_inf.cond{i} '_' col_names{j} '.tsv']), 'FileType', 'text', 'Delimiter','\t');
    end
end

spear=array2table(spear, 'VariableNames', col_names, 'RowNames', mod_inf.cond);
log_RMSE=array2table(log_RMSE, 'VariableNames', col_names, 'RowNames', mod_inf.cond);

writetable(spear, fullfile('Results', 'eccomp', 'abunvse', 'spearman.tsv'), 'FileType', 'text', 'Delimiter', 'tab', 'WriteRowNames',true);
writetable(log_RMSE, fullfile('Results', 'eccomp', 'abunvse', 'RMSE.tsv'), 'FileType', 'text', 'Delimiter', 'tab', 'WriteRowNames',true);
%match reactions to protein data
    

