function [] = printSummary(model, ub, count, fluxes,con_names, outf)
%Function that takes bounds, and model used for flux estimation
% as well as resulting flux and kapp statistics
% and prints a human readable output of the mean result
% statistics over all conditions
% If a filename is given a tsv of upper bounds and predicted fluxes
% Is also saved
%INPUT:
%   struc   model:      The COBRA model used for flux prediction
%   double  ub:         The upper bounds of fluxes used in prediction for
%                       in each condition
%   char    con_names:  The names of the different conditions 
%   struc   count:      A count structure returned by getkapp()
%   double  fluxes:     Predicted fluxes for the Kapp estimation in each
%                       condition

n_fluxes=size(fluxes,1);
idx=[find(findExcRxns(model));findRxnIDs(model, 'Biomass_Chlamy_mixo')];
for i=1:size(ub,2)
    disp(append('bounds and fluxes in condition', string(i)))
    temp_tab=table(model.rxns(idx), ub(idx,i), fluxes(idx, i));
    temp_tab.Properties.VariableNames(2:3)=join([repmat(con_names(i),2,1),{'ub', 'flux'}']);
    disp(temp_tab)
    if i==1
        res_tab=temp_tab;
    else
        res_tab=[res_tab, temp_tab(:,2:3)];
    end
end

disp(append('ON AVERAGE ' , string(mean(count.nonzero)), ' of ', string(n_fluxes), ...
    ' reaction in the model carried flux larger then the threshhold 1e-10'))
disp(append(string(mean(count.nonzero2)), ' of ', string(mean(count.withabun)), ...
    ' reactions catalyzed by homomeric (iso) enzymes with known abundance carried flux'))
disp(append('Allowing for the calculation of the k_app values of ', string(mean(count.kapp)+count.lsqrgood), ...
    ' reactions corresponding to ' ,string(mean(count.homoenzyme2)), ...
    ' enzymes'))
disp(append('The relativ amount of reactions that carried flux for all reactions for which enzyme abundance was measured was', ...
    string(mean(count.ratio))))
disp(append('Additionally, for ', string(count.lsqrgood) , ' reactions with all isoenzymes expressed in standard condition a kapp could be determined'))
%Create a barplot stacked barplot [kapps(mutiple expressed isozyme
%reactions); kapps(single expressed isozyme); kapps(homomer); no_kapps(but abundande and flux))
%
plotmat=[repelem(count.lsqrgood, length(count.nonzero))', count.singisozymes, count.kapp-count.singisozymes, count.multisoenzyme-count.lsqrgood];
p=bar(plotmat, 'stacked');
set(p, {'DisplayName'}, {'Isozyme','Single Isozyme','Single Homomer', 'Not determined'}')
legend()
xticklabels(con_names)
org_ylim=ylim();
ylim([0 org_ylim(2)*1.2])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 4]);
saveas(gcf, fullfile(fileparts(outf), 'kappperrxn.svg'))
close
if exist('outf', 'var')
    if ~exist(fileparts(outf), 'dir')
        mkdir(fileparts(outf))
    end
    writetable(res_tab, append(outf, '.tsv'), 'FileType', 'text', 'Delimiter', 'tab')
end