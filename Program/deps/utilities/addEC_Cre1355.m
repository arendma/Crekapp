function [Cre1355model] = addEC_Cre1355(Cre1355model, supplfp)
%Function to map the EC reactionnumbers of Cre1355 from the supplementary
%xlsx to the COBRA modle
%INPUT:
% struct Cre1355model: A COBRA model structure of Cre1355 
% char supplfp: Path to the xlsx file from the Cre1355 publication
Cre1355_meta= readtable(supplfp);
replacePattern='[\,\&\;\-()\+''" _]';
%manually replace ECnumber annotations with invalid syntax
Cre1355_meta.ProteinClassification(ismember(Cre1355_meta.ReactionID, {'STARCH300S'}))={'2.4.1.21;2.4.1.18'};
Cre1355_meta.ProteinClassification(contains(Cre1355_meta.ProteinClassification, '+')) = ...
    regexprep(Cre1355_meta.ProteinClassification(contains(Cre1355_meta.ProteinClassification, '+')), ' \+ ', ';');
rxnsum=[setdiff(Cre1355model.rxns, Cre1355_meta.ReactionID), setdiff(Cre1355_meta.ReactionID, Cre1355model.rxns)];
rxnsum2=rxnsum;
for i=1:size(rxnsum2, 2)
    rxnsum2(:,i)=cellfun(@(x)regexprep(lower(x),replacePattern,''),rxnsum2(:,i),'un',0);
end
%get ID of corresponding entry in table
idx=double(size(rxnsum2,1));
for i=1:size(rxnsum2,1)
    tmp_ids=find(ismember(rxnsum2(:,2), rxnsum2(i,1)));
    if isempty(tmp_ids)
        idx(i)=0;
    elseif length(tmp_ids)==1
        idx(i)=tmp_ids;
    else
        error(append('Multiple matching reaction IDS in table for', rxnsum2{i,1}))
    end
end
%add 1 missing index by exclusion reasoning
idx(idx==0)=setdiff(1:length(idx), idx);
%adapt ordering
rxnsum(:,2)=rxnsum(idx,2);
%substitute xlsx reaction IDs with matched SMBL model IDs
[~, convert]=ismember(rxnsum(:,2), Cre1355_meta.ReactionID);
Cre1355_meta.ReactionID(convert)=rxnsum(:,1);
%add a EC Number information field to model
[~,convert]=ismember(Cre1355_meta.ReactionID, Cre1355model.rxns);
if any(convert==0)
    error('Not all reaction IDs match between model and xlsx')
end

tmpEC=Cre1355_meta.ProteinClassification;
% regex pattern to match E.C. numbers
ecPattern = '\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)';

% determine the total number of E.C. numbers in the model
tmpEC(cellfun(@isempty, regexp(tmpEC, ecPattern, 'match')))={''};

Cre1355model.rxnECNumbers=tmpEC(convert);
end