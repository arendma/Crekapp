function[report] = mapp_cre_NIDLE(ecModel, savname)
%Function to genearte mapping between reactionin raw GECKO model and 
%Reaction in NIDLE model of Chlamydomonas reinhardtii
%INPUT:
% - struct ecModel: A GECKO format enzyme constrained Cre model
% - savname: prefix for saving the mapping in
%            /Results/Hammel_NIDLE/2022/<savname>_mapping.tsv

NIDLE_kapp=readtable('Results/NIDLE/kcat_n.tsv', 'FileType','text', 'Delimiter', '\t');

%match reactions using perfect and fuzzy matching of substringst
match=zeros(length(ecModel.rxns), size(NIDLE_kapp, 1));
fzmatch=zeros(length(ecModel.rxns), size(NIDLE_kapp, 1));
for i=1:size(NIDLE_kapp,1)
    query=regexprep(NIDLE_kapp.Rxns{i}, {'_f$', '_b$'}, {'', '_REV'});
    match(:,i)=contains(ecModel.rxns, query, 'IgnoreCase', true);
    for j=1:length(ecModel.rxns)
        min_dist=fzsearch(ecModel.rxns{j},query, 0, 1);
        min_dist=min_dist{1};
        fzmatch(j,i)=min_dist(1);
    end 
end

%create report cell array first column davidi reaction name, second colunn
%perfect matches of substrings, third column model reaction abbreviations
%with shortes levensthein distance
report=cell(size(NIDLE_kapp,1),4);
report(:,1)=NIDLE_kapp.Rxns;
for i=1:size(NIDLE_kapp,1)
    perfmatch=ecModel.rxns(logical(match(:,i)));
   
    %if there is only one perfect match already assign final match
    if length(perfmatch)==1
        report(i,2)=perfmatch;
        report(i,4)=perfmatch;
    %if there are two perfect matches and they only differ in '_REV' take
    %the forward reaction as final match
    elseif length(perfmatch)==2 & strcmp(regexprep(perfmatch(1), '_REV', ''), ...
            regexprep(perfmatch(2), '_REV', ''))
         report{i,2}=strjoin(ecModel.rxns(logical(match(:,i))), '; ');
         report(i,4)=perfmatch(~contains(perfmatch, '_REV'));
    %if there are exact matchs if No<digit> and _f is removed or No<digit> is removed
    %and _b is substituted by _REV assign them as (concatenated) final
    %match
    elseif sum(ismember(regexprep(perfmatch, 'No\d{1,2}$', ''), regexprep(report{i,1}, {'_f$', '_b$'}, {'', '_REV'} )))>0
        report{i,2}=strjoin(ecModel.rxns(logical(match(:,i))), '; ');
        report{i,4}=strjoin(perfmatch(ismember(regexprep(perfmatch, 'No\d{1,2}$', ''), regexprep(report{i,1}, {'_f$', '_b$'}, {'', '_REV'} ))), '; ');
    else
        report{i,2}=strjoin(ecModel.rxns(logical(match(:,i))), '; ');
    end
    report{i,3}=strjoin(ecModel.rxns(fzmatch(:,i)==min(fzmatch(:,i))), '; ');
end
report=cell2table(report, 'VariableNames', {'reactionName', 'perfect_match', 'min_levensthein_match', 'final_match'});
writetable(report, fullfile('Results', 'NIDLE', [savname, '_mapping.tsv']), 'FileType','text', 'Delimiter','\t')
if any(cellfun(@isempty, report.final_match))
    error(['Exported ' savname '_mapping.tsv has missing final match entries, manual curation nescessary'])
end
end

