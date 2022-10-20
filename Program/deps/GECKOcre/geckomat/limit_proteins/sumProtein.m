%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = sumProtein(model)
% Calculates protein content from the model.
%
% Benjamin Sanchez. Last update: 2018-10-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = sumProtein(model)
%Function to get the relative Protein Content according to biomass reaction
%Input: 
% - struct model: A iCre1355 model structured used in COBRA
%Output:
% - num P: relative protein content in g/gDW
if nargin>2
    warning('GAM refitting not implemented for scaleBioMass in C.reinhardtii')
    if nargin>3
        warning('scale_comp argument is not implemente for scaleBioMass in C. reinhardtii')
    end
end
    

if ~isfield(model, 'c') || sum(model.c)~=1
    error('vector c does not enable detection of biomass objective')
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

%calculate effective coeffecients for biomass (remove growht associate
%maintanance ATP demand (= H20 demand in pseudoreaction)
eff_coeff=model.S(bm_edukt_id, biomass_id);
P=abs(eff_coeff(ismember(comp_ed.Type, {'P'}))'*comp_ed.MW(ismember(comp_ed.Type, {'P'})));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
