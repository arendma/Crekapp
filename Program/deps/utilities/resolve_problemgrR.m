function [o_model] =resolve_problemgrR(i_model)
%Function to resolve problematic occurences of ') and (', ') and' & 'and('
%by explicitly including all isoenzyme combinations of the complex in the
%grRules field
%INPUT:
% struc i_model: A GEM model compatible with COBRA
%OUTPUT:
%struc o_model: The GEM with resolved .grRules field
if ~isfield(i_model, 'grRules')
   i_model=creategrRulesField(i_model);
end 

problm=find(contains(i_model.grRules, {') and (', ') and', 'and ('}));
o_model=i_model;
for i=problm'
    %divide string according to brackets
    splitpoints=[];
    grRule=i_model.grRules{i};
    
    splitgrR=strsplit(grRule, ' and ');
    bracket_c=cellfun(@(x) count(x, '('), splitgrR);
    % if there is more than 1 opening bracket in the second to las element
    % this could be an indication of problematic strings
    if any(bracket_c(2:end)>1)
        warning("Possibly problematic bracketing detecting in grRule:")
        disp(grRule)
    end
    %remove brackets and split by or to get isozyme groups of complexes
    for j=1:length(splitgrR)
        splitgrR{j}= regexprep(splitgrR{j}, '\( ?| ?\)', '');
        splitgrR{j}=strsplit(splitgrR{j}, ' or ');
    end
    gRmat=splitgrR{1}';
    for j=2:length(splitgrR)
        new_col=repelem(splitgrR{j},size(gRmat,1))';
        gRmat=[repmat(gRmat, [length(splitgrR{j}),1]), new_col];
    end
    %Create all possible complex rules
    new_complex=cell(size(gRmat, 1), 1);
    for j=1:size(gRmat,1)
        new_complex{j}=append('( ', strjoin(gRmat(j, :), ' and '), ' )');
    end
    if any(duplicates(new_complex),'all')
        error('Duplicated new complex rules found. Aborting...')
    end
    
    %concatenate rules
    o_model.grRules{i}=append('(', strjoin(new_complex', ' or '), ')');
end
end