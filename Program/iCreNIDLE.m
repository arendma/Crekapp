%% Adapted File from NIDLE manuscript
if exist('~/Software/cobratoolbox', 'dir')
    addpath '~/Software/cobratoolbox'
elseif ~exist('C:/cobratoolbox/', 'dir')
    error('No cobratoolbox path found')
end
addpath '../kcats/'
addpath '../Utilities/'
initCobraToolbox(false)
%Solver setup
changeCobraSolver('gurobi', 'all');
%Sicne gurobi is directly called in NIDLE() set feasibility tolerance
%as variable passed to the function 
feasTol=1e-6;
% Get a EC_number annotated irreversible model that splits all reactions
% reversible in mixo and heterotrohpic conditions, obtain conditions
% specific flux bounds for the proteomics conditions (acetate and photon
% flux are limited): numberic arrays lb ub (s. parse_bounds)
model_file.file={'../Data/Cre1355/iCre1355_mixo_upd.xml', '../Data/Cre1355/iCre1355_hetero_upd.xml'};
model_file.name={'mixo','hetero'};
cond_file='../Data/QconCAT_David20220124/fitted_acup.tsv';
[model_irrev, lb, ub, conditions]=parse_bounds(model_file, cond_file);

%load proteomics data from david in Amol per cell 
abundance_dat=readtable('../Data/QconCAT_David20220124/abs_abundance/med_abun_all.tsv', 'FileType', 'text');

%'Fit kcat information
[viri_Kcats, viri_Kcat_matches]=matchKcats(model_irrev, '../kcats/kcats-full-lineages.tsv', 'Chlamydomonas', 'Viridiplantae', false);
[euk_Kcats, euk_Kcat_matches]=matchKcats(model_irrev, '../kcats/kcats-full-lineages.tsv', 'Viridiplantae', 'Eukaryota', false);
[proc_Kcats, proc_Kcat_matches]=matchKcats(model_irrev, '../kcats/kcats-full-lineages.tsv', 'Cyanobacteria', 'Bacteria', false);
lit_kcats=table(viri_Kcats, viri_Kcat_matches, euk_Kcats, euk_Kcat_matches, proc_Kcats, proc_Kcat_matches);
%Export lit kcats indipendent of found Kmax values
lit_kcats_out=[model_irrev.rxns, model_irrev.subSystems, lit_kcats];
lit_kcats_out.Properties.VariableNames(1:2)= {'Reaction', 'SubSystem'};
writetable(lit_kcats_out, '../Data/Cre1355/Cre1355_lit_kcats.txt', 'Delimiter', '\t');

%abundance structure contains normalized abundance for all genes present in
%GPR rules of homomeric reaction w and w/o isozymes and a mapping between
%reactions, and GPR rules
%g_vect is a boolean matrix having a 1 if the corresponding reaction in the
%corresponding condition has enough associated proteins expressed to deem
%it active (dependent on wether it is a complex, isoenzyme etc)
[abundance, g_vect, report_tab] = parse_abundance(model_irrev, conditions,abundance_dat);


%report table on which genes have all enzymes expressed
writetable(report_tab, '../Results/Hammel_NIDLE/2022/found_enzymes.tab', 'fileType', 'Text', 'Delimiter', '\t')
%%
%Running pFBA and NIDLE-flux
scale=1e6; %the resulting flux will be 1e-6 times the original scale
conv=1e-6;
%threshold for pFBA to consider non-zero: 1e-10

% res_31_p=struct();
res_hamm_p.flux=pFBA(model_irrev,lb,ub,scale); 
res_hamm_n=NIDLE(model_irrev,g_vect,lb,ub,1e-4,scale, feasTol, 1);

%Computing Kapp

[Kapp_hamm_p,V_hamm_p,count_hamm_p]=getkapp(abundance,g_vect, res_hamm_p.flux,model_irrev, lit_kcats, conv,1e-4,(3600*1e-10), '../Results/Hammel_pFBA/Hammel_kcat_p');
[Kapp_hamm_n,V_hamm_n,count_hamm_n]=getkapp(abundance,g_vect,res_hamm_n.flux,model_irrev, lit_kcats,conv,1e-4,(3600*1e-10), '../Results/Hammel_NIDLE/2022/Hammel_kcat_n');

%print summary of kapp statistics and fluxes
printSummary(model_irrev, ub,count_hamm_p, V_hamm_p, conditions.cond, '../Results/Hammel_pFBA/pFBA_flux')
printSummary(model_irrev, ub,count_hamm_n, V_hamm_n, conditions.cond, '../Results/Hammel_NIDLE/2022/NIDLE_flux')

