%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleBioMass(model,Ptot,GAM,scale_comp)
% 
% Benjamin Sanchez. Last update: 2018-10-23
% Ivan Domenzain.   Last update: 2019-09-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_model,GAM] = scaleBioMass(model,Ptot,GAM,scale_comp)
%Function to rescale Biomass according to proteomics measurements. It is
%supposed that the remaining cell composition remains unchanged
%Input: 
% - struct model: A iCre1355 model structured used in COBRA
% - double Ptot: measured relative cellular protein content in g/gDW 
%Output:
% - struct out_model:   iCre1355 model with rescaled coefficients in Biomass
%                                   pseudoreaction
% - GAM_coeff: GAM in mmol/gdW*h

if nargin>2
    warning('GAM refitting not implemented for scaleBioMass in C.reinhardtii')
    if nargin>3
        warning('scale_comp argument is not implemente for scaleBioMass in C. reinhardtii')
    end
end
    

if ~isfield(model, 'c') || sum(model.c)~=1
    error('vector c does not enable detection of biomass objective')
end
biomass_id=find(model.c);

% %correct for models with double compartment annotation - not present in
% RAVEN models --> deactivated
% model.mets=regexprep(model.mets, '\[\w\](\[\w\])', '$1');
% %for models with _c compartment annotation convert to novel compartment
% %annotation
% model.mets=regexprep(model.mets, '_\w([\w\])$', '$1');
% %import massweight info on metabolitesn from Ines
 mass_inf=readtable('../../databases/biomass_mw_Chlamy.csv');
% %adapt compartment abbrevitations for iCre1355
mass_inf.met_ID=regexprep(mass_inf.met_ID, '\[(\w)]$', '_$1');
bm_edukt_id= model.S(:, biomass_id)<0;
bm_prod_id=model.S(:, biomass_id)>0;

%Compile massweights for all biomass components
comp_ed=table(model.mets(bm_edukt_id), nan(sum(bm_edukt_id),1), repmat({'O'}, ...
    sum(bm_edukt_id),1), 'VariableNames', {'met_ID', 'MW', 'Type'});
%for all metabolites present in INes data take the MW from there
%for the sace of cell composition rescaling take the aa weights from ines
%as weight for the respective trnas in the pseudoreaction#
[match, ind] = ismember(regexprep(comp_ed.met_ID, 'trna_c', '__L_c'), mass_inf.met_ID);
comp_ed.MW(match)=mass_inf.MW(ind(match));
%set type info to protei for trna compounds
comp_ed.Type(contains(comp_ed.met_ID, 'trna_c'))={'P'};
%manually add glycine since it doesn't follow the __L[c] pattern
comp_ed.MW(contains(comp_ed.met_ID, 'glytrna_c'))=mass_inf.MW(contains(mass_inf.met_ID, 'gly_c'));
%Calculate the molecular weight for remaining unknowns from chemical
%formular (this is only problematic for trna ATP and h20 which are not
%drawn from the system and turning up in edukts
comp_ed.MW(isnan(comp_ed.MW))=computeMW(model, comp_ed.met_ID(isnan(comp_ed.MW)), true, true);

%convert to g/mmol
comp_ed.MW=comp_ed.MW/1000;
%disp('Rescale biomass reaction using the following massweights:')
%disp(comp_ed)

%calculate effective coeffecients for biomass (remove growht associate
%maintanance ATP demand (= H20 demand in pseudoreaction)
eff_coeff=model.S(bm_edukt_id, biomass_id);
GAM_coeff=eff_coeff(ismember(comp_ed.met_ID, {'atp_c'}));
GAM=abs(GAM_coeff);
%%ATTENTION this leads to positive effective coefficients since some water
%%is produced during biomass creation
eff_coeff(ismember(comp_ed.met_ID, {'atp_c'}))=eff_coeff(ismember(comp_ed.met_ID, {'atp_c'}))-GAM_coeff;
eff_coeff(ismember(comp_ed.met_ID, {'h2o_c'}))=eff_coeff(ismember(comp_ed.met_ID, {'h2o_c'}))-GAM_coeff;
%total mass of a cell (may be uequal to 1 because of differences in
%massweight of AA and protonation state
Tot=eff_coeff'*comp_ed.MW;

o_Ptot=eff_coeff(ismember(comp_ed.Type, {'P'}))'*comp_ed.MW(ismember(comp_ed.Type, {'P'}));

%Rescale the biomass
facP=-Ptot/o_Ptot;
facO=(Tot-(-Ptot))/(Tot-o_Ptot);

%new coefficients
n_coeff=nan(length(eff_coeff),1);
n_coeff(ismember(comp_ed.Type, {'P'}))=eff_coeff(ismember(comp_ed.Type, {'P'}))*facP;
n_coeff(ismember(comp_ed.Type, {'O'}))=eff_coeff(ismember(comp_ed.Type, {'O'}))*facO;
%check that standardization is kept
if n_coeff'*comp_ed.MW-Tot>0.0000001
    error('error when rescaling biomass. Total mass after rescaling is kept constant')
end

%Update biomass function 
model.S(bm_edukt_id, biomass_id)=n_coeff;

%add GAM to new coefficients
n_coeff(ismember(comp_ed.met_ID, {'atp_c', 'h2o_c'}))=n_coeff(ismember(comp_ed.met_ID, {'atp_c', 'h2o_c'}))+GAM_coeff;

%Update biomass function 
out_model=model;
out_model.S(bm_edukt_id, biomass_id)=n_coeff;

bm_prod_id=find(bm_prod_id);
[match, ind]=ismember(regexprep(model.mets(bm_prod_id), 'trna', ''), regexprep(comp_ed.met_ID, 'trna', ''));
out_model.S(bm_prod_id(match), biomass_id)=-n_coeff(ind(match));

%diplay new protein part
bm_edukt_id=find(bm_edukt_id);
P_n=out_model.S(bm_edukt_id(ismember(comp_ed.Type, {'P'})), biomass_id)'*comp_ed.MW(ismember(comp_ed.Type, {'P'}));
disp(append('New relative biomass content=', num2str(P_n)))


end