
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MW = calculateMW(sequence)
% Calculates de molecular weight of a protein.
% 
% INPUT: Sequence (can include extra things such as spaces or numbers).
% OUTPUT: Molecular weight number
% 
% Benjam�n S�nchez. Last edited: 2015-04-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MW = calculateMW(sequence)

% A	Alanine
% B	Aspartic acid or Asparagine
% C	Cysteine
% D	Aspartic acid
% E	Glutamic acid
% F	Phenylalanine
% G	Glycine
% H	Histidine
% I	Isoleucine
% J	Leucine or Isoleucine
% K	Lysine
% L	Leucine
% M	Methionine
% N	Asparagine
% O	Pyrrolysine
% P	Proline
% Q	Glutamine
% R	Arginine
% S	Serine
% T	Threonine
% U	Selenocysteine
% V	Valine
% W	Tryptophan
% X	any
% Y	Tyrosine
% Z	Glutamic acid or Glutamine

aa_codes = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N', ...
            'O','P','Q','R','S','T','U','V','W','X','Y','Z'};
aa_MWs   = [71.08 114.60 103.14 115.09 129.11 147.17 57.05 137.14 ...
            113.16 113.16 128.17 113.16 131.20 114.10 255.31 97.12 ...
            128.13 156.19 87.08 101.10 150.04 99.13 186.21 126.50 ...
            163.17 128.62];

MW = 18;
for i = 1:length(aa_codes)
    count = length(strfind(sequence,aa_codes{i}));
    MW = MW + count*aa_MWs(i);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%