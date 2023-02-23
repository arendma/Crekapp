function[NIDLE_mod]=create_NIDLEmod(raw_batchmod, enzMetPfx, modcond, cnames, leaveout,rxnbased )
%function to create a batchmodell from NIDLE data of C. reihardtii
%INPUT:
% - struc raw_batchmod: A unmodified ecModel obtained from GECKO
% - char enzMetPfx: Character vector marking enzyme pseudometabolites
% - char modcond: Character indicating the specific condition model
%                 ('auto', 'mixo', 'hetero')
% - cell cnames:    column names of kapp columns in NIDLE output (condition
%                   names)
% - leavout:    A character vector or cell array giving names of conditions to
%               be excluded when calculating the max kapp
% - logic rxnbased: Logical, if true uses reaction centric kapp output from
%                   NIDLE ('kcat_n.tsv'), else uses enzyme centric output
%                   from NIDLE ('kcat_n_genewise.tsv')
if nargin <6 
    rxnbased=true;
    if nargin<5
    leaveout={};
    end
end
enzMetIdx = find(contains(raw_batchmod.mets, enzMetPfx) &~ismember(raw_batchmod.mets, {'prot_pool'}));
RxnIdx=find(~contains(raw_batchmod.rxns, 'draw_prot_'));
if rxnbased
    kapp=readtable('Results/NIDLE/kcat_n.tsv', 'FileType','text');
else
    kapp=readtable('Results/NIDLE/kcat_n_genewise.tsv', 'FileType','text');
    kapp.UP_ID=JGItoUP(kapp.ID5_5);
    disp(append(string(size(kapp, 1)-length(unique(kapp.UP_ID))), ' duplicated UP_ID entries were detected and all entries with these gene ID removed.'))
    kapp(find(sum(duplicates(kapp.UP_ID),2)),:)=[];
end
%calculate maximum kapp
kapp.kappmax=max(table2array(kapp(:,setdiff(cnames, leaveout))), [],2);
kapp(isnan(kapp.kappmax),:)=[];
kapp(kapp.kappmax==0,:)=[];

enz_missing={};
NIDLE_mod=raw_batchmod;
if rxnbased
    %get reaction mapping information
    if ~isfile(['Results/NIDLE/' modcond '_mapping.tsv'])
        map=mapp_cre_NIDLE(raw_batchmod, modcond);
    else
        map=readtable(['Results/NIDLE/' modcond '_mapping.tsv'], 'FileType','text');
    end
    %filter unmatched
    map=map(~cellfun(@isempty, map.final_match),[1,4]);
    %duplicate isozyme entries.
    i=1;
    while i<=size(map,1)
        if contains(map.final_match(i), ';')
            isozms=strsplit(map.final_match{i}, '; ')';
            tmp_map=map(repelem(i, length(isozms)), :);
            tmp_map.final_match=isozms;
            map=[map(1:(i-1), :); tmp_map; map((i+1):end,:)];
            i=i+length(isozms);
        else
            i=i+1;
        end
    end
  

    for i=1:size(map,1)
        if ~any(ismember(kapp.Rxns, map.reactionName(i)))
            %due to leavout samples not all Rxns in map have to be in kapp
            %table
            continue
        end
%         %retrieve davidi kapp%                 disp('oldkcat:')
         new_kapp=-1/(kapp.kappmax(ismember(kapp.Rxns, map.reactionName(i)))*3600);
        %retrieve max presto kcat 
        enz_idx=find(raw_batchmod.S(enzMetIdx,ismember(raw_batchmod.rxns, map.final_match(i))));
        if length(enz_idx)~=1
            if length(unique(raw_batchmod.S(enzMetIdx(enz_idx),ismember(raw_batchmod.rxns, map.final_match(i)))))==1
                %complexes should not be found but accept them if kcat is
                %identical
                disp('Complex')
                NIDLE_mod.S(enzMetIdx(enz_idx), ismember(raw_batchmod.rxns, map.final_match(i)))=new_kapp;
            else
                if length(enz_idx)==0
                    enz_missing=[enz_missing, raw_batchmod.rxns{ismember(raw_batchmod.rxns, map.final_match(i))}];

                end
                warning(['reaction with missing enzyme or enzyme complex detected comparison impossible. skipping reaction ', ...
                    raw_batchmod.rxns{ismember(raw_batchmod.rxns, map.final_match(i))}])
            end
        else

            NIDLE_mod.S(enzMetIdx(enz_idx), ismember(raw_batchmod.rxns, map.final_match(i)))=new_kapp;
        end
    end
else
    %if enyzme centric values are used set the kcat for all reaction this
    %enzyme is involved in to the NIDLE value 
    warning("For enzymes with nidle kappmax the kcat for all reactions they participate in is changed. This leads to complexes with different stochiometries")
    for i=1:size(kapp,1)
        enz_idx=ismember(NIDLE_mod.mets, ['prot_', kapp.UP_ID{i}]);
        if sum(enz_idx)==0
           enz_missing=[enz_missing, kapp.ID5_5(i)];
           warning(['reaction with missing enzyme detected comparison impossible. skipping enzyme ', ...
                    kapp.ID5_5{i}])
        elseif sum(enz_idx)==1
           enzrxn_idx=intersect(find(NIDLE_mod.S(enz_idx,:)), RxnIdx);
           %Update kcats
%            disp('-1/Old kcats:')
%            disp(NIDLE_mod.S(enz_idx,enzrxn_idx))
           NIDLE_mod.S(enz_idx,enzrxn_idx)=-1/(kapp.kappmax(i)*3600);
%            disp('-1/New kapp:')
%            disp(NIDLE_mod.S(enz_idx,enzrxn_idx))
        elseif sum(enz_idx)>1
            error(append("duplicated protein metabolite name detected: prot_",kapp.UP_ID{i}))
        end
    end
end
end