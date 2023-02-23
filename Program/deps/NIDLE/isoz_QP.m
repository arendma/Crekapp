function res = isoz_QP(v,E, Eps)
%Function that finds the minimum L2 estimates of kapps
%INPUT:
% - num v: column vector of fluxes in different conditions
% - mat E:  numeric matrix that gives the enzyme concentrations, rows give
%           the condition while columns give the different isoenzymes
%- num Eps: A real number giving the lower bound of kcat values (default:
%           1e-10 (lowest kcat in BRENDA+SABIORK 0722=5.8e-10
if nargin<3
    Eps=3600*1e-10;
end
%take the top n(=number of isozymes) conditions ranked by decreasing flux
% [~, top]=sort(v, 'descend');
n=size(E,1);
m=size(E,2);
% v_top=v(:,top(1:n));
% E_top=E(:, top(1:n));
%k1E11+k1E21 +d1=v1
Aflux=[E, eye(n,n)];
bflux=v;

%sum(d)-t<=0
%-sum(d)-t<=0
% Asum=[zeros(2*n,m), [eye(n,n);-eye(n,n)], [-eye(n,n);-eye(n,n)]];
% bsum=zeros(2*n,1);
%min(d'*Q*d)
model.Q=sparse([zeros(m+n,m),[zeros(m,n);1*eye(n,n)]]);
model.lb=[repmat(Eps, m,1);-inf(n,1)];
model.A=sparse(Aflux); 
model.sense=repelem('=', n);
model.rhs=bflux;
sol=gurobi(model);

if strcmp(sol.status,'OPTIMAL')
    res.kapp=sol.x(1:m);
    res.RMSD=sqrt(sum(sol.x((m+1):(m+n)).^2)/n);
else
    warning(['GUROBI has non optimal exit status: ' sol.status])
end
end
