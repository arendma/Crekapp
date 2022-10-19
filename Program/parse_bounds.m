function [irrev_mod, lb, ub, cond] = parse_bounds(model_file, cond_file, obj_rf)
%Function that differently boudn version of the same model in model_dir
%and generates a model having all possible reversible reactions, 
%and a set of bounds combining the original modle bounds and 
%infor on conditions in the table stored at cond_dir
% INPUT:
% - cell model_file: A struc with two fields:
%                   - file: cell array of character vectors giving the relative
%                   paths to differently bounded model versions
%                   - name: cell array of character vectors giving the
%                   names of the different models

% - char cond_file: one character vector giving the path to a table with
%                   conditions infos having columns: Condition(condition
%                   name), mu(growth rate), model( model versions to use in
%                   this condition) and a variable number of columns that
%                   have the name of reactions in the irreversible model as
%                   headers and contain the upper bound for this reaction
%                   in the given condition
% -num obj_rf:  a single number giving the relative relaxation of the
%               biomass objective (default 0.05)
%OUTPUT:
% - struct irrev_mod:   COBRA model structrue of irreversible model reversing
%                       any reaction that is reversible in the differently
%                       bounded models in mod_file.file
%- num lb, ub:  numeric arrays where each column corresponds to the
%               lower(lb) or upper (ub) bounds of irrev_mod parse from the
%               input models and the cond_file
%- struct cond: A structure having one field 
%               - cond: a cell array of character vectors giving the
%               different conditions names for lb and ub (weird format is
%               due to NIDLE compatibility)
if nargin<3
    obj_rf=0.05;
end
models=cell(length(model_file.file),1);
%Read in models with different bounds
for i=1:length(model_file.file)
    models{i}=readCbModel(model_file.file{i});
    if i==1
        models{i}=addEC_Cre1355(models{i}, '../Data/Cre1355/tpj13059-sup-0011-tables5.xlsx');
        %Fix problematic GPR rules
        models{i}=resolve_problemgrR(models{i});
    end
    if sum(models{i}.c)==0
        error('Loaded models do not have reaction marked for optimization, fitting of bounds for nidle impossible')
    end
end
n_r=length(models{1}.rxns);
%generate model with maximum positive and minimum negative bound
model_lb=nan(n_r,length(model_file.file));
model_ub=model_lb;
for i=1:length(models)
    model_lb(:,i)=models{i}.lb;
    model_ub(:,i)=models{i}.ub;
end
max_mod=models{1};
max_mod.lb=min(model_lb,[], 2);
max_mod.ub=max(model_ub,[], 2);
%get irreversible model 
irrev_mod=convertToIrreversible(max_mod);

%read in condition data and set upper and lower bounds
cond_info=readtable(cond_file, 'FileType', 'text');
cond.cond=cond_info.Condition;
lb=zeros(length(irrev_mod.lb), size(cond_info,1));
ub=lb;
for i=1:size(cond_info,1)
    %get model id
    idx=find(ismember(model_file.name, cond_info.model(i)));
    %regenerate default model bounds
    %for postivie non reversible reactions keep model bounds
    forw_r=find(models{idx}.lb>=0 & models{idx}.ub>=0);
    lb(forw_r,i)=models{idx}.lb(forw_r);
    ub(forw_r,i)=models{idx}.ub(forw_r);
    %for negative non reversible reactions invert bounds
    backw_r=find(models{idx}.lb<0 & models{idx}.ub<=0);
    lb(backw_r,i)=-1*models{idx}.ub(backw_r);
    ub(backw_r,i)=-1*models{idx}.lb(backw_r);
    %adjust upper bound of reversible reaction
    rev_r=find(models{idx}.lb<0 & models{idx}.ub>0);
    ub(rev_r,i)=models{idx}.ub(rev_r);
    ub(irrev_mod.match(rev_r),i)=-1*models{idx}.lb(rev_r);
    
    % fit condition specicif paramters
    %mu 
    lb(find(models{idx}.c), i)=(1-obj_rf)*cond_info.mu(i);
    ub(find(models{idx}.c), i)=(1+obj_rf)*cond_info.mu(i);
    %Uptake reaction
    for j=find(~(cellfun(@isempty, regexp(cond_info.Properties.VariableNames, '^EX'))))
    ub(irrev_mod.match(ismember(models{idx}.rxns, cond_info.Properties.VariableNames(j))),i)=cond_info{i, j};
    end
end
end

