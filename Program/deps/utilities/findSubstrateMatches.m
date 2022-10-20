function [matchingIdx,matchingIDs] = findSubstrateMatches(substrates1,substrates2,exclude)
%% matchingIdx = findSubstrateMatches(substrates1,substrates2,exclude(opt))
% attempt to match substrate names after conversion to lower case and
% removing problematic symbols
% Input:
%     cellstr/cell substrates1:       array of substrate names, the
%                                     function will look for matches
%                                     of the first array in the second array
%     cellstr substrates2
%     cellstr exclude(optional):      array containing metabolite names,
%                                     which should be ignored when matching
%                                     the substrates; if contained, will
%                                     return no match for the given entry
% Output:
%     logical matchingIdx:            contains matching substrates, has
%                                     dimension of the first substrate array
%     cellstr matchingIDs:            contains the matched substrates of
%                                     the first array (in adapted form)
if nargin < 2
    error('USAGE: findSubstrateMatches(substrates1,substrates2,exclude(opt))')
elseif ~iscell(substrates1)
    error('first array is not of type cellstr')
elseif ~iscellstr(substrates2)
    error('second array is not of type cellstr')
elseif nargin < 3
        exclude = {'nadh', 'nadp', 'nadph', 'nad', 'atp', 'adp', 'h2o', 'coa', 'h',...
            'fadh2', 'fad'};
end
% define problematic characters which should be removed
replacePattern='[\,\&\;\-()\+''" _]';

% prepare arrays
substrates1=cellfun(@(x)regexprep(lower(x),replacePattern,''),substrates1,'un',0);

substrates2=regexprep(lower(substrates2),replacePattern,'');
for j=1:numel(exclude)
    substrates2=cellfun(@(x)regexprep(x,['^',exclude{j},'$'],''),substrates2,'un',0);
end

matchingIdx=cellfun(@(x)any(ismember(x,substrates2)),substrates1);
matchingIDs=substrates1(matchingIdx);
    
end
