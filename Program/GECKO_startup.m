%setup environment for gecko remove comments from code lines for de-novo
%run
function GECKO_startup(idx, assignkcat, logname)
%Function to create GECKO models form IMAM et al data and boyle & Morgen
%Ptot measures
%INPUT:
% - num start_idx:  either an integer giving the index from which to start
%                   builiding models, or an integer vector giving the
%                   model indices to loop over, default: 1
% - char assignkcat:    a character vector giving the path to the kcat.mat
%                       file in GECKOdir/models/ from which matched kcats 
%                       should be read, if empty kcats are matched to 
%                       reaction denvoe (runtime ~4h)
% - char logname:	a character vector giving the prefix of the log file
%			to which log file will be saved
if nargin<1
    idx=1:5;
elseif length(idx)==1
    idx=idx:5;
end

if nargin<2
    assignkcat=[];
end

if nargin<3
	logname='';
end

home=pwd();
%Adapt relative_proteomics
copyfile("Data/QconCAT_David20220124/abs_abundance/muns_relative_proteomics.txt", "Program/deps/GECKOcre/databases/relative_proteomics.txt")
%import protein biomass fraction from Boyle and Morgan switch between ptot
%by uncommenting/commenting
GKOinf={'auto', 0.261, {'EX_co2_e', 'EX_nh4_e', 'EX_pi_e'}, 
    'mixo', 0.303, {'EX_co2_e', 'EX_ac_e','EX_nh4_e', 'EX_pi_e' }
    'hetero', 0.222, {'EX_ac_e', 'EX_nh4_e', 'EX_pi_e'}};
GKOinf=cell2table(GKOinf, 'VariableNames', {'name', 'Ptot', 'C_upt'});
cstat_dat=readtable('Data/Cre1355/chemostat_dat.txt', 'TreatAsEmpty', {'NA'});
%add culture conditions id
cstat_dat.cond_id={'auto', 'auto', 'auto', 'auto', 'mixo', 'mixo', 'mixo', 'mixo', 'hetero', 'hetero'}';
%adapt reaction names as column names
cstat_dat.Properties.VariableNames(4:7)={'EX_co2_e', 'EX_nh4_e', 'EX_pi_e', 'EX_ac_e'};
%only keep every second row - data looks like every two replicates were
%taken fro m the sam starter culture
cstat_dat=cstat_dat(1:2:10,:);


%read in autotrophic model with JGI gene names
Cr55_a=readCbModel('Data/Cre1355/iCre1355_auto_up.xml');
Cr55_m=readCbModel('Data/Cre1355/iCre1355_mixo_up.xml');
Cr55_h=readCbModel('Data/Cre1355/iCre1355_hetero_up.xml');
models.mods={Cr55_a, Cr55_m, Cr55_h};
models.name={'auto', 'mixo', 'hetero'};
diary([logname date '.log'])
for i=idx
    base_mod=models.mods{ismember(models.name, cstat_dat.cond_id(i))}; 
    %adapt model for batch growth (set respective carbon sources to 1000)
    r=GKOinf.C_upt{ismember(GKOinf.name, cstat_dat.cond_id(i))};  
    disp('Old Flux bounds')
    printFluxBounds(base_mod, r')
    base_mod=changeRxnBounds(base_mod, r, -1000, 'l');
    disp('New Flux bounds')
    printFluxBounds(base_mod, r')
    cd Program/deps/GECKOcre/geckomat
    sav_parameters=getModelParameters();
    parameters=sav_parameters;
    parameters.Ptot=GKOinf.Ptot(ismember(GKOinf.name, cstat_dat.cond_id(i)));
    parameters.gR_exp=cstat_dat.mu_h_1_(i);
    save('../databases/parameters.mat', 'parameters')
    enhanceGEM(base_mod, 'COBRA', ['geCre1355' cstat_dat.Var1{i}], 'Cre1355', assignkcat)
    parameters=sav_parameters;
    save('../databases/parameters.mat', 'parameters')
end
cd(pwd)
    
