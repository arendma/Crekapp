function out_model =scaleBiomassCre(cre_model, biomass_id, Ptot)
%Function to rescale Biomass according to proteomics measurements. It is
%supposed that the remaining cell composition remains unchanged
%Input: 
% - struct cre_model: A iCre1355 model structured used in COBRA
% - int/char biomass_id:    Either character vector containing the rxn id of 
%                                       the biomass reaction or integer giving the Index 
%                                       of the biomass reaction in the model
% - double Ptot: measured relative cellular protein content in g/gDW 
%Output:
% - struct out_model:   iCre1355 model with rescaled coefficients in Biomass
%                                   pseudoreaction

% if biomass_id is given ass character vector convert ot index
if ischar(biomass_id)|| isa(biomass_id,'string')
    biomass_id=findRxnIDs(cre_model, biomass_id)
end

%correct for models with double compartment annotation
cre_model.mets=regexprep(cre_model.mets, '\[\w\](\[\w\])', '$1');
%for models with _c compartment annotation convert to novel compartment
%annotation
cre_model.mets=regexprep(cre_model.mets, '_\w([\w\])$', '$1');
%import massweight info on metabolitesn from Ines
mass_inf=readtable('Data/biomass_mw_Chlamy.csv');
%adapt compartment abbrevitations for iCre1355
mass_inf.met_ID=regexprep(mass_inf.met_ID, '_(\w)$', '[$1]');
bm_edukt_id= cre_model.S(:, biomass_id)<0;
bm_prod_id=cre_model.S(:, biomass_id)>0;

%Compile massweights for all biomass components
comp_ed=table(cre_model.mets(bm_edukt_id), nan(sum(bm_edukt_id),1), repmat({'O'}, ...
    sum(bm_edukt_id),1), 'VariableNames', {'met_ID', 'MW', 'Type'});
%for all metabolites present in INes data take the MW from there
%for the sace of cell composition rescaling take the aa weights from ines
%as weight for the respective trnas in the pseudoreaction#
[match, ind] = ismember(regexprep(comp_ed.met_ID, 'trna\[c\]', '__L[c]'), mass_inf.met_ID);
comp_ed.MW(match)=mass_inf.MW(ind(match));
%set type info to protei for trna compounds
comp_ed.Type(contains(comp_ed.met_ID, 'trna[c]'))={'P'};
%manually add glycine since it doesn't follow the __L[c] pattern
comp_ed.MW(contains(comp_ed.met_ID, 'glytrna[c]'))=mass_inf.MW(contains(mass_inf.met_ID, 'gly[c]'));
%Calculate the molecular weight for remaining unknowns from chemical
%formular (this is only problematic for trna ATP and h20 which are not
%drawn from the system and turning up in edukts
comp_ed.MW(isnan(comp_ed.MW))=computeMW(cre_model, comp_ed.met_ID(isnan(comp_ed.MW)), true, true);

%convert to g/mmol
comp_ed.MW=comp_ed.MW/1000;
disp('Rescale biomass reaction using the following massweights:')
disp(comp_ed)

%calculate effective coeffecients for biomass (remove growht associate
%maintanance ATP demand (= H20 demand in pseudoreaction)
eff_coeff=cre_model.S(bm_edukt_id, biomass_id);
GAM_coeff=eff_coeff(ismember(comp_ed.met_ID, {'atp[c]'}));
eff_coeff(ismember(comp_ed.met_ID, {'atp[c]'}))=eff_coeff(ismember(comp_ed.met_ID, {'atp[c]'}))-GAM_coeff;
eff_coeff(ismember(comp_ed.met_ID, {'h2o[c]'}))=eff_coeff(ismember(comp_ed.met_ID, {'h2o[c]'}))-GAM_coeff;
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
cre_model.S(bm_edukt_id, biomass_id)=n_coeff;

%add GAM to new coefficients
n_coeff(ismember(comp_ed.met_ID, {'atp[c]', 'h2o[c]'}))=n_coeff(ismember(comp_ed.met_ID, {'atp[c]', 'h2o[c]'}))+GAM_coeff;

%Update biomass function 
out_model=cre_model;
out_model.S(bm_edukt_id, biomass_id)=n_coeff;

bm_prod_id=find(bm_prod_id);
[match, ind]=ismember(regexprep(cre_model.mets(bm_prod_id), 'trna', ''), regexprep(comp_ed.met_ID, 'trna', ''));
out_model.S(bm_prod_id(match), biomass_id)=-n_coeff(ind(match));

%diplay new protein part
bm_edukt_id=find(bm_edukt_id);
P_n=out_model.S(bm_edukt_id(ismember(comp_ed.Type, {'P'})), biomass_id)'*comp_ed.MW(ismember(comp_ed.Type, {'P'}));
disp(append('New relative biomass content=', num2str(P_n)))

end
