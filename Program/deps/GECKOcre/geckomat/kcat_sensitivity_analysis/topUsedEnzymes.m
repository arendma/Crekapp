function T = topUsedEnzymes(fluxes,model,conditions,name,writeFile,n)
% topUsedEnzymes
%
%   Function that gets an ecModel_batch and a matrix of FBA solution vectors
%   and calculate the top ten enzyme usages for every solution (conditions)
%   in a mass-wise way. The results are written in the topUsedEnzymes.txt
%   file and stored in the container folder.
%
%   fluxes      matrix with FBA solution vectors 
%               [# Rxns in model, # simulated conditions]
%   model       ecGEM structure used to generate FBA solution vectors
%   conditions  cell array of strings, with names identifying the conditions
%               for each solution vector
%               solution vector
%   name        string with model abbreviation (opt, default 'ecModel')
%   writeFile   logical, whether a file should be written at the location
%               gecko/models/'name' (opt, default true)
%   n           number of top used enzymes to be listed. (opt, default 10)
%
%   T           Table including the top "n" used proteins (names and usages)
%               for each specified condition
%
%   Usage: T = topUsedEnzymes(fluxes,model,conditions,name,writeFile,n)
% 
%   Ivan Domenzain    Last edited 2018-12-05
%   Eduard Kerkhoven  Last edited 2018-12-05
%

if nargin < 6
    n = 10;
end
if nargin < 5
    writeFile = true;
end
%Find the enzyme usage reactions
usages = find(~cellfun(@isempty,strfind(model.rxns,'prot_')));
%Exclude protein pool exchange reaction
usages = usages(1:end-1);
%Avoid exceeding vector dimension
if n>length(usages)
    n=length(usages);
end
usages = fluxes(usages,:);
[~, nConds] = size(usages);
outputFile  = cell(n,2*nConds);
colNames    = cell(1,2*nConds);

for i=1:length(usages(1,:))
    %Units conversion [mmol/gDwh] -> [g prot/gDwh]
    enzUsage = usages(:,i).*model.MWs;
    %Get the mass fraction of the proteome for every protein in the mod
    enzUsage = enzUsage/sum(enzUsage);
    %Sort and keep the top ten used enzymes
    [enzUsage,Indexes] = sort(enzUsage,'descend');
    enzNames = model.enzymes(Indexes);
    enzUsage = enzUsage(1:n);
    enzNames = enzNames(1:n);
    outputFile(:,(2*i)-1) = enzNames;
    outputFile(:,(2*i))   = num2cell(enzUsage);
    %Create column names for output table/file
    colNames(1,(2*i)-1)   = strcat('prots_',conditions(i));
    colNames(1,(2*i))     = strcat('Usages_',conditions(i));
end
%Write the top-ten used enzyme names and their percentage usages for
%every condition on the output file
i=1:length(usages(1,:));
i=(2*i);
outputFile = truncateValues(outputFile,i);
T = cell2table(outputFile,'VariableNames',colNames);
if writeFile
    writetable(T,['../../models/' name '/' name '_topUsedEnzymes.txt'])
end
end