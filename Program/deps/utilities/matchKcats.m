function [kcatAll,match] = matchKcats(model, kcatFileName, taxTerm, taxLimit, pruneEC)
%% kcatAll = findKcatFromEC(model, substrates, kcatFileName, taxTerm, taxLimit(opt), pruneEC(opt))
% This function scans a file containing turnover numbers for matching
% E.C. numbers (assigned to reactions) in a metabolic model.
% The results from matching an E.C. number are, if possible, further
% filtered for matching substrates or given taxonomic classification (any level).
% If kcat(s) have been found, the overall maximum is returned chosen per
% reaction (per E.C. number and over all associated E.C. numbers).
% Compounds, which are often non-specific are excluded from the matching
% procedure.
% Input:
%   struct model:           COBRA metabolic model structure
%   char kcatFileName:      path to file containing turnover numbers
%                           columns: E.C.-Number;substrate names;taxonomic
%                           lineage;turnover value;temperature
%                           the file has no header line!
%   char taxTerm:           taxonomic classification, can be organism name
%                           but also higher-level terms
%   char taxLimit:          (optional) limit for taxonomic range, default:''
%   logical pruneEC:        (optional) whether or not pruning should be
%                           used to match E.C. numbers, default: false
% Output:
%   double kcatAll:         array containing one kcat per matched reaction;
%                           if unassigned: kcatAll(i)==0
%   double match:           quality of the match 1: taxon + substrate; 2:
%                           substrate, 3:taxon, 4:EC Number; 0:NO match

% initialize outpur variable and counters
kcatAll = nan(size(model.rxns));
match= zeros(size(model.rxns));
ecNotFoundCount = 0;
prunedECMatchCount = 0;
substrateMatchCount = 0;
fullTaxMatchCount = 0;
completeMatchCount = 0;

% get substrates and products for every reaction
substrates = cell(numel(model.rxns),1);
products = cell(numel(model.rxns),1);
for i = 1:numel(substrates)
    substrates{i} = model.metNames(model.S(:,i)<0);
    products{i} = model.metNames(model.S(:,i)>0);
end

% regex pattern to match E.C. numbers
ecPattern = '\d+\.\d+\.\d+\.\d+';

% determine the total number of E.C. numbers in the model
tmpEC = regexp(model.rxnECNumbers, ecPattern, 'match');
nEC  = sum(~cellfun(@isempty,[tmpEC{:}])); clear tmpEC

% loop over all reactions
n = numel(model.rxns);
fprintf('\n')
disp(['Obtaining maximum kcats for ', num2str(n), ' reactions.'])
fprintf('\n')
for i = 1:n
    if mod(i,100)==0
        disp(['Processed ', num2str(i), ' reactions ...'])
    end
    
    % get all E.C. numbers assigned to the current reaction
    ec = regexp(model.rxnECNumbers{i}, ecPattern, 'match');
    
    % only proceed if the reaction has at least one E.C. number assigned
    if ~isempty(ec)
        % initialize kcat array for the current reaction
        tmpResTable = cell2table(cell(0,2),'VariableNames',{'Var4','Var5'});
        for j = 1:numel(ec)
            % initial E.C. number level and default status
            level = 4;
            status = 1;
            % define the initial query
            query = ['^', ec{j}, '\>'];
            
            while status ~= 0 && level ~= 1
                
                % call unix grep to scan turnover value file
                if ispc
                    [status, res] = system(['findstr "', query, '" ', regexprep(kcatFileName, '\/', '\')]);    
                else 
                    [status, res] = system(['grep "', query, '" ', kcatFileName]);
                end
                % check whether we have a zero exit code
                if status == 0
                    
                    % process the result to obtain a table
                    res = strtrim(res);
                    res = strsplit(res, {'\t', '\n'});
                    
                    % a check if every row was five column entries
                    try
                        res = cell2table(reshape(res, 4, numel(res)/4)');
                    catch
                        disp(res)
                    end
                    
                    % if a taxonomic limit is definedm all considered
                    % entries must contain this term
                    if ~isempty(taxLimit)
                        res=res(contains(res.(3),taxLimit),:);
                    end

                    % find matching products
                    productMatchIdx = findSubstrateMatches(res.(2),products{i});
                    % exclude those entries for which the substrate matches
                    % the product of the current reaction
                    res = res(~productMatchIdx,:);
                    
                    if isempty(res)
                        status = 1;
                    end
                    
                    % find matching substrates
                    substrateMatchIdx = findSubstrateMatches(res.(2),substrates{i});
                    
                    % find entries for the given taxonomic term
                    taxMatchIdx = ~cellfun(@isempty, regexpi(res.(3), taxTerm));
                    
                    if any(substrateMatchIdx&taxMatchIdx)
                        substrateMatchCount = substrateMatchCount + 1;
                        fullTaxMatchCount = fullTaxMatchCount + 1;
                        completeMatchCount = completeMatchCount + 1;
                        Var5=repmat(1, sum(substrateMatchIdx&taxMatchIdx),1);
                        % substrate and fungal match
                        tmpResTable = [tmpResTable;res(substrateMatchIdx&taxMatchIdx,4), table(Var5)];
                    elseif any(substrateMatchIdx)
                        substrateMatchCount = substrateMatchCount + 1;
                        Var5=repmat(2,sum(substrateMatchIdx),1);
                        % only substrate match
                        tmpResTable = [tmpResTable;res(substrateMatchIdx,4),table(Var5)];
                    elseif any(taxMatchIdx)
                        fullTaxMatchCount = fullTaxMatchCount + 1;
                        Var5=repmat(3, sum(taxMatchIdx), 1);
                        % only taxonomic match
                        tmpResTable = [tmpResTable;res(taxMatchIdx,4),table(Var5)];
                    else
                        Var5=repmat(4,size(res,1),1);
                        % if no filter possible, take all results
                        tmpResTable = [tmpResTable;res(:,4),table(Var5)];
                    end
                    
                end
                
                if pruneEC && (status ~=0 || isempty(res))
                    % remove the lower EC level
                    query = regexprep(query, '\.[\d\-]+(\\>)*$', '');
                    level = level - 1;
                elseif level < 4
                    prunedECMatchCount = prunedECMatchCount + 1;
                else
                    level = 1;
                end
                
            end
            
            if status ~= 0 && ~isempty(regexp(ec(j),ecPattern))
                % only count as 'not found' if the reaction has an E.C.
                % number assigned
                ecNotFoundCount = ecNotFoundCount + 1;
            end
        end
        
        if status == 0 && ~isempty(tmpResTable)
            tmpResTable.(1) = str2double(tmpResTable.(1));
            idxMaxKcat = find(tmpResTable.(1)==max(tmpResTable.(1)));
            kcatAll(i) = tmpResTable.(1)(idxMaxKcat(1));
            match(i)= tmpResTable.(2)(idxMaxKcat(1));
            clear tmpResTable
        end
    end
    
end

fprintf('\nfinished!\n')
disp(['total number of E.C. numbers: ', num2str(nEC)])
disp(['complete matches: ', num2str(completeMatchCount)])
disp(['substrate matches: ', num2str(substrateMatchCount)])
disp(['taxonomic matches: ', num2str(fullTaxMatchCount)])
disp(['E.C. class matches: ', num2str(prunedECMatchCount)])
disp([num2str(ecNotFoundCount), ' E.C. numbers could not be matched'])
disp('')

end
