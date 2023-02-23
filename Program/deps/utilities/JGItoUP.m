function up_ID = JGItoUP(JGI_ID)
%Function to convert a JGIv.5.6 transcript ID to the corresponding UNIPROT
%gene ID
%INPUT:
% - cell/char JGI_ID =either a cell array of  character vector of JGI.v5.5 transcript ID
%OUTPUT:
% -char up_ID =character vector of correspoinding UNIPROT ID
%
%Dependencies:
% - duplicates.m from COBRA Toolbox
% - conversion table in tsv format from biomart
% - the GECKO generated ProtDatabase.mat database
%NOTE: This function is based on the for loop in adapt_uniprot.m
charout=false;
if ischar(JGI_ID)
    charout=true;
    JGI_ID={JGI_ID};
end
%import biomart information obtained from JGI
%Query:
%https://phytozome-next.jgi.doe.gov/biomart/martview?VIRTUALSCHEMANAME=zome_mart&ATTRIBUTES=phytozome.default.features.gene_name1|phytozome.default.features.uniprot_id|phytozome.default.features.peptide_name|phytozome.default.features.enzyme_id&FILTERS=phytozome.default.filters.organism_id."281"&VISIBLEPANEL=resultspanel
mart_inf=readtable('Data/mart_Cre_Uniprot.txt');
%
%To use it as a conversion from transcript and to Uniprot remove rows which
%are duplicates except the EC value
%duplicates functio is from COBRA toolbox
mart_inf(logical(sum(duplicates(mart_inf(:,1:3)))),:) = [];
%remove rows with missing uniprot info
mart_inf(cellfun(@isempty, mart_inf.UniProtID),:)=[];
%create a dictionary
JGI2UP=containers.Map(mart_inf.PeptideName, mart_inf.UniProtID);

%import Prot information database constructed by GECKO 06-2022
load('Data/GECKOCre_ProtDatabase.mat')
up_ID=cell(size(JGI_ID));
for i=1:length(JGI_ID)
    new_gn='';
try
    new_gn=JGI2UP(JGI_ID{i});
catch
    %% 2. try to find gene in kegg
    in_kegg=contains(kegg(:,3), JGI_ID{i});
    if sum(in_kegg)==1 %if one is found set the gene name accordingly
        idx=find(in_kegg);
        %Kegg contains 1JGI-to-manyUNiprot mappings
        %Check for these cases and only take the first ID
        if contains(kegg{idx, 1}, ' ')
            ss=strsplit(kegg{idx,1});
            new_gn=ss{1};
        else
            new_gn=kegg{idx, 1};
        end
    elseif sum(in_kegg)>1
        warning(append(JGI_ID{i}, ' is linked to several ids in kegg. Skipping...'))
    else
        %% 3. try to find gene in swissprot
        in_swissprot=contains(swissprot(:,3), JGI_ID{i}, 'IgnoreCase', true);
        if sum(in_swissprot)==1 %if one is found set the gene name accordingly
            idx=find(in_swissprot);
            new_gn=swissprot{idx, 1};
        elseif sum(in_swissprot)>1
            disp(append(JGI_ID{i}, ' is linked to several ids in swissprot. Skipping...'))
        else
            %% 4. try to find transcript unspecific ENSEMBLE ID in
            %swissprot - depends on JGItoENS.m
            ENSid=JGItoENS(JGI_ID{i});
            in_swissprot=contains(swissprot(:,3), ENSid);
            if sum(in_swissprot)==1 %if one is found set the gene name accordingly
                disp(append(JGI_ID{i}, ' is only found via transcript unspecific ENSEMBLE ID'))
                idx=find(in_swissprot);
                new_gn=swissprot{idx, 1};
            elseif sum(in_swissprot)>1
                warning(append(JGI_ID{i}, ' is linked to several ids in swissprot. Skipping...'))
            else
                warning(append(JGI_ID{i}, ' can not be mapped to UNIPROT. Skipping...'))
            end
        end
    end
end
if ~isempty(new_gn)
    up_ID{i}=new_gn;
end
end
if charout
    up_ID=up_ID{1};
end
end

function ENSEMBLEID = JGItoENS(JGI_CreID)
% Function to convert a JGIv.5.5. id semantically to a ENSEMBLE ID 
% If JGI transcript ID is given the transcript identifier will be removed
% INPUT: 
%  -  char JGI_CreID: A character vector of the JGI ID 
% OUTPUT: 
%  - char ENSEMBLEID: A character vector of the corresponding ENSEMBLE ID 
    ENSEMBLEID= regexprep(JGI_CreID, '\.t.*$', '');
    ENSEMBLEID= regexprep(ENSEMBLEID, '^Cre', 'CHLRE_');
    ENSEMBLEID= regexprep(ENSEMBLEID, '\.', '');
    ENSEMBLEID= append(ENSEMBLEID, 'v5');
end

