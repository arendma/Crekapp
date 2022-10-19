function[V_solution_p]=pFBA(model_irrev,lb,ub,scale)
n_cond=size(lb,2);
[n_m,n_r]=size(model_irrev.S);

f=ones(n_r,1);
V_solution_p=zeros(n_r,n_cond);

biomass_idx=find(contains(model_irrev.rxns, 'Biomass'));


for cond=1:n_cond
    biomass=biomass_idx(ub(biomass_idx,cond)~=0);
    m=struct();
    m.obj = f;
    Aeq=model_irrev.S;
    beq=zeros(n_m,1);
    Aeq(:,biomass)=Aeq(:,biomass)*scale;
    
    m.A = sparse(Aeq);
    n = size(m.A, 2);
    m.vtype = repmat('C', n, 1);
    m.sense =  repmat('=',size(Aeq,1),1);
    m.rhs = full( beq(:));
    
    m.lb=lb(:,cond)*scale; 
    m.ub=ub(:,cond)*scale;
    m.lb(biomass)=lb(biomass,cond);
    m.ub(biomass)=ub(biomass,cond);
    
    params = struct();
    params.FeasibilityTol=1e-6;
    x = gurobi(m,params);
    
    if ~strcmp(x.status,'INFEASIBLE')
        sol_p=x.x;
        V=sol_p(1:n_r);
        V_solution_p(:,cond)=V;
    end
end
