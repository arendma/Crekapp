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