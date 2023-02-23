function [rules_type,vect, report_tab]=parse_rules(model,abundance_u)
%Function to find if all nescessary  proteins to catalyze a reaction are
%expressed
%INPUT: 
% - struct model: An irreversible model used by NIDLE
% - struct abundance_u:
%Output: 
% - double rules_type: An integer giving the type of GPR rule for each
% reaciton
% - double vect: A matrix (:D) of n_reaction x n_conditions vector containing 1
% for each reaction, for which all nescessary enzymes are found in
% abundance data in the respective condition
% - tab report_tab: A table containing the reaction name and subsystem
% together with the vect matrix
n_r=size(model.S,2);
if ~isfield(model, 'grRules')
    model=creategrRulesField(model);
end
%1) Find the type of GPR rules_type
%first filter for multimers
ind1=find(contains(model.grRules,' or ('));
ind2=find(contains(model.grRules,') or '));
complex_iso=union(ind1,ind2); %if there is an or and a bracket present there is a complex with an isoenzyme

ind1=find(contains(model.grRules,' and ('));
ind2=find(contains(model.grRules,') and '));
complex_and=union(ind1,ind2); %if there is an and and a bracket present its a heteromer

only_and=setdiff(complex_and,complex_iso); %these are heteromers without isoenzymes involved

rules_type=zeros(n_r,1);
%rules_type==0: not defined yet
%rules_type==1: no GPR rules_type
%rules_type==2: single gene
%rules_type==3: isoenzyme
%rules_type==4: complex
%rules_type==5: complex isoenzyme
%rules_type==6: complex 'complex'
rules_type(complex_iso)=5;
rules_type(only_and)=6;
l=1:n_r;
rest=setdiff(l,union(complex_iso,only_and));
for i=1:n_r
    if find(rest==i)
        if isempty(model.grRules{i})
            rules_type(i)=1;
        elseif contains(model.grRules{i}, 'or')
            rules_type(i)=3;
        elseif contains(model.grRules{i}, 'and')
            rules_type(i)=4;
        else
            rules_type(i)=2;
        end
    end
end



%Part2: assign vector g according to the conditions

n_cond=size(abundance_u.abun,2);

g_vect=NaN(n_r,n_cond);

single=find(rules_type==2);
iso=find(rules_type==3);
compl=find(rules_type==4);
compl2=find(rules_type==6);
compl_iso=find(rules_type==5);

for cond=1:n_cond
    %1. single genes
    
    for i=1:length(single)
        reac_ind=single(i);
        [~,gene_ind]=ismember(model.grRules(reac_ind),abundance_u.genes);
        if gene_ind~=0     
           g_vect(reac_ind,cond)=abundance_u.abun(gene_ind,cond);
        end
    end 
    
    %2. isoenzymes
    for i=1:length(iso)
        reac_ind=iso(i);
        genes=split(model.grRules{reac_ind},' or ');
        genes=erase(genes,{'(', ' '});
        genes=erase(genes,')');
        [~,gene_ind]=ismember(genes,abundance_u.genes);
        gene_ind2=gene_ind(gene_ind~=0);
        if ~isempty(gene_ind2)
            g_vect(reac_ind,cond)=max(abundance_u.abun(gene_ind2,cond),[], "omitnan");
            %we only care if there is at least one enzyme with measured
            %abundance nanmax(NaN,NaN)=NaN
        end
    end

    %2. complex
    for i=1:length(compl)
        reac_ind=compl(i);
        genes=split(model.grRules{reac_ind},' and ');
        genes=erase(genes,{'(', ' '});
        genes=erase(genes,')');
        [~,gene_ind]=ismember(genes,abundance_u.genes);
        gene_ind2=gene_ind(gene_ind~=0);
        if length(gene_ind2)==length(gene_ind)
            ind=find(isnan(abundance_u.abun(gene_ind2,cond)));
            if ~isempty(ind)
                g_vect(reac_ind,cond)=NaN;
            else
                g_vect(reac_ind,cond)=min(abundance_u.abun(gene_ind2,cond));
            end
        end
    end
    
   %3. complex complex
    for i=1:length(compl2)
        reac_ind=compl2(i);
        genes=strsplit(model.grRules{reac_ind},' and ');
        genes=erase(genes,{'(', ' '});
        genes=erase(genes,')');
        [~,gene_ind]=ismember(genes,abundance_u.genes);
        gene_ind2=gene_ind(gene_ind~=0);
        if length(gene_ind2)==length(gene_ind)
            ind=find(isnan(abundance_u.abun(gene_ind2,cond)));
            if ~isempty(ind)
                g_vect(reac_ind,cond)=NaN;
            else
                g_vect(reac_ind,cond)=min(abundance_u.abun(gene_ind2,cond));
            end
        end
    end 
    
    % 4. complex isoenzymes
    for i=1:length(compl_iso)
        reac_ind=compl_iso(i);
        complex=split(model.grRules{reac_ind},' or ');
        count2=NaN(length(complex),1);
        for m=1:length(complex)
            genes=split(complex(m),' and ');
            genes=erase(genes,{'(', ' '});
            genes=erase(genes,')');
            [~,gene_ind]=ismember(genes,abundance_u.genes);
            gene_ind2=gene_ind(gene_ind~=0);
            if length(gene_ind2)==length(gene_ind)
                ind=find(isnan(abundance_u.abun(gene_ind2,cond)));
                if ~isempty(ind)
                    count2(m)=NaN;
                else
                    count2(m)=min(abundance_u.abun(gene_ind2,cond));
                end
            end
        end
        g_vect(reac_ind,cond)=max(count2,[], "omitnan");
    end
end

vect=zeros(n_r,n_cond);
vect(~isnan(g_vect))=1;
report_tab=table(model.rxns, model.rxnNames, model.subSystems, model.rxnECNumbers, rules_type);
report_tab=[report_tab, num2cell(vect)];
report_tab.Properties.VariableNames= [{'rxns', 'rxnNames', 'SubSystems', 'ECNumber', 'GPR_ruletype'}, abundance_u.cond'];
end