function [o_CbModel] = fixGenes(i_CbModel, verbose)
% Function to remove duplicated entries in i_CmModel.genes and fix
% problematic gr Rules
% Input:
%     struc i_CbModel: A Cobratoolbox model >=v.2021
% Output:
%     struc o_CbModel: A Corbratoolbox model >=v.2021

%Create grRules fiels
if ~isfield(i_CbModel, 'grRules')
    i_CbModel=creategrRulesField(i_CbModel);
end

%find duplicate entries in genes
dup=find(sum(duplicates(i_CbModel.genes)));
o_CbModel=i_CbModel;
for dp=dup
    %For each duplicate find grRules with multiple occurences
    dp_gn=i_CbModel.genes{dp};
    matches=regexp(o_CbModel.grRules, dp_gn);
    if verbose
        disp(append('GPR rules including duplicates gene ',  dp_gn, ' before reduction:'))
        disp(o_CbModel.grRules(~cellfun(@isempty, matches)))
    end
    %replacement regular expression - don't declare in each loop again
    rep_regx=append('(or|and) ', dp_gn, ' ');
    
    %count number of occurences
    redund=cellfun(@length, matches);
    
    for j=find(redund>1)'
        %if there is more than 1 occurence split by brackets to keep the
        %very rare case of "(gup_gn and g1) or (dup_gn and g2)" - this
        %should be written as (dup_gn and (g1 or g2)"
        split_rule=strsplit(o_CbModel.grRules{j}, ')', 'CollapseDelimiters', false);
        rule_matches=regexp(split_rule, dp_gn);
        %only substitute if the gene turns up twice in the same pracket
        %         if any(cellfun(@length, rule_matches)>2)
        %             warning(append('Detected duplicated gene more than twice in the same bracket GPR', ...
        %                 'statement, This happens when when more then two transcripts are combined under', ...
        %                 'the same ID or there are errors in the GPR rules. The Code does not remove more', ...
        %                 'than 1 duplicated and will produce wrong rules. Aborting...'))
        %             pause
        %         end
        rule_redund=cellfun(@length, rule_matches);
        %from each bracket expression remove all duplicates
        for k=find(rule_redund>1)'
            n=rule_redund(k);
            while n>1
                split_rule(k)=regexprep(split_rule(k), rep_regx, '','once');
                n=n-1;
            end
        end
        o_CbModel.grRules{j}=strjoin(split_rule, ')');
    end
    if verbose
        disp(append('GPR rules including duplicates gene ',  dp_gn, 'after reduction:'))
        disp(o_CbModel.grRules(~cellfun(@isempty, matches)))
    end
end
%remove old gene related fields
ck_fields={'genes', 'geneNames', 'rules' , 'RxnGeneMat', 'proteins'};
for i=1:length(ck_fields)
    if isfield(o_CbModel, ck_fields(i))
        o_CbModel=rmfield(o_CbModel, ck_fields(i));
    end
end
o_CbModel=resolve_problemgrR(o_CbModel);
%create rules and genes field from updated gprRules
o_CbModel=generateRules(o_CbModel, 0);
%check fi all duplicates where removed
if any(duplicates(o_CbModel.genes), 'all')
    warning('Failed to remove all duplicates from genes')
end
%Check if gene names contain characters indicating corruption
if any(contains(o_CbModel.genes, {' ', ';'}))
    warning('Likely orrupted gene names detected')
    disp(o_CbModel.genes(contains(o_CbModel.genes, {' ', ';'})))
end
%
end