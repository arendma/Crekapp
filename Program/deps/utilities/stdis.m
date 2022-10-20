function dl = stdis(r,p,n)

% dl = stdis(R,P,1) builds matrix for calculation of Levenshtein distances
% between strings R and P;
% dl = stdis(R,P,0) builds matrix for calculation of Levenshtein distances
% for matching P to substrings of R. 
% Wagner-Fischer algorithm is used.

li = numel(r) + 1;
lu = numel(p) + 1;       
dl(lu,li) = 0; % memory allocation
dl(:,1) = 0:lu - 1;
if n
  dl(1,:) = 0:li - 1;
end
%Distance
for i = 2:lu
   bi = p(i-1);
   for j = 2:li
      k = 1;
      if strcmp(r(j-1),bi)
         k = 0;
      end
   dl(i,j) = min([dl(i-1,j-1) + k,dl(i-1,j) + 1,dl(i,j-1)+1]);
   end
end
