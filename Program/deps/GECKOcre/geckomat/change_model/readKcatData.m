%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = readKcatData(model_data,kcats, name)
% Reads the output of the getBRENDAdata module, saves a tab seperated
% overview of kcats and EC numbers matched to the model and creates the
% modified enzyme model, with additional metabolites (the enzymes) and
% reactions (for all isoenzymes and also the enzyme exchange reactions).
%
% INPUT:
% model_data        model and EC numbers and substrates/products from each
%                   reaction (output from "getECnumbers.m")
% kcats             kcats for each reaction/enzyme (output from
%                   "matchKcats.m")
%name      model name to give path to folder
%
% OUTPUT:
% eModel            modified model accounting for enzymes
% out_tab           A table where each row corresponds to a reaction in the
%                   irreversible model. column 1=rxns name; column 2= kcat
%                   entries for all isozymes; column 3: EC numbers for all
%                   isozymes: column 4: GECKO matching score for all
%                   isozyme kcats
% 
% Marius Arend             2021-06-05
% Cheng Zhang               2015-12-03
% Benjamin J. Sanchez       2018-08-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eModel, out_tab] = readKcatData(model_data,kcats)

n_r=size(kcats.forw.org_s, 1);
orig_tab=cell(n_r+sum(model_data.model.rev),1);
for i=1:size(orig_tab, 1)
    entry=cell(1,size(kcats.forw.org_s,2));
    for j=1:length(entry)

     %index of reverse reactions 
     rev_idx=find(model_data.model.rev);
     %for forward reactions
     if i<=n_r
         orig=find([kcats.forw.org_s(i,j)  kcats.forw.rest_s(i,j)  kcats.forw.org_ns(i,j) ...
             kcats.forw.org_sa(i,j) kcats.forw.rest_ns(i,j) kcats.forw.rest_sa(i,j)], 1, 'first');
         if ~isempty(orig)
            entry{j}=num2str(orig);
         end
     else %for reverse reactions
         orig=find([kcats.back.org_s(rev_idx(i-n_r),j)  kcats.back.rest_s(rev_idx(i-n_r),j) ...
             kcats.back.org_ns(rev_idx(i-n_r),j) kcats.back.org_sa(rev_idx(i-n_r),j) kcats.back.rest_ns(rev_idx(i-n_r),j) ...
             kcats.back.rest_sa(rev_idx(i-n_r),j)], 1, 'first');
         if ~isempty(orig)
             entry{j}=num2str(orig);
         end
     end
    end
     orig_tab{i}=strjoin(entry(~cellfun(@isempty, entry)), ';');
    
end
%Get kcat value for both directions:
Fkcat = kcats.forw.kcats;
Bkcat = kcats.back.kcats;
rev   = logical(model_data.model.rev);
kcats = [Fkcat;Bkcat(rev,:)];

%Update uniprots with both directions:
uniprots = [model_data.uniprots; model_data.uniprots(rev,:)];

%Update matched genes with both directions:
matchedGenes = [model_data.matchedGenes; model_data.matchedGenes(rev,:)];

%Update EC numbers
EC =[model_data.EC_numbers; model_data.EC_numbers(rev,:)];

%Convert to irreversible model with RAVEN function (will split in 2 any reversible rxn):
model = convertToIrrev(model_data.model);

%Save a table with an overview of kcats and EC numbers
out_tab=cell(size(orig_tab,1), 3);
for i=1:size(out_tab,1)
    kcat_entry={''};
    ec_entry={''};
    %check if EC Numbers and kcats have the same entries
    if ~isequal(find(kcats(i,:)>0), find(~cellfun(@isempty, EC(i,:))))
        error('Different number of entries and kcats detected. This is an internal error :(')
    end
    for k=find(kcats(i,:)>0)
        if k==1
            %convert h-1 to s-1 since BRENDA values are converted upon
            %import
            kcat_entry=num2str(kcats(i,k)/3600);
            ec_entry=EC(i,k);
        else
        kcat_entry=append(kcat_entry,';', num2str(kcats(i,k)/3600)); %convert h-1 to s-1
        ec_entry=append(ec_entry, ';', EC{i,k});
        end 
    end
    out_tab(i,:)=[model.rxns(i), kcat_entry, ec_entry ];
end
out_tab=[cell2table(out_tab, 'VariableNames', {'rxns', 'kcats', 'EC_Numbers'}), ...
    cell2table(orig_tab, 'VariableNames', {'MScores'})];
out_tab=[out_tab, array2table(orig_tab)];


%Convert original model to enzyme model according to uniprots and kcats:
eModel = convertToEnzymeModel(model,matchedGenes,uniprots,kcats);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
