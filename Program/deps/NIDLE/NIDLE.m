function [new]=NIDLE(model,g_vect,lb,ub,ep,scale, feasTol, Threads,  OutputFlag)
%FUNCTION to construct the NIDLE MILP
%INPUT:
% struct model: The irreversible model in COBRA format
% double g_vect:    A matrix (:D) of n_reaction x n_conditions vector containing 1
%                   for each reaction, for which all nescessary enzymes are found in
%                   abundance data in the respective condition
% double lb:    A matrix n_reactions x n_conditions bearing the lower bounds
%               of reactions
% double ub:    Same as lb but containing upper bounds of reactions
% double ep:    Number giving the lower limit of reaction flux to be 
%               considered an active reaction
% double scale: A factor to scale all reaction fluxes except the biomass
%               flux by (here stochiometric coefficients are scaled)
% double feasTol:  Gurobi solver feasibility tolerance
% double threads:   Param passed to gurobi solver (default 0; 
%                   takes the available number of virtual cores or 32 if 
%                   more available)
% double OutputFlag: Param passed to gurobi solver (default 0; no output)

if nargin<9
    OutputFlag=0;
    if nargin<8
        Threads=0;
    end
end
biomass_idx=find(contains(model.rxns, 'Biomass'));


[n_m,n_r]=size(model.S);
n_cond=size(lb,2);
new.y=zeros(n_r,n_cond);
new.flux=zeros(n_r,n_cond);


for cond=1:size(lb,2)
    biomass=biomass_idx(ub(biomass_idx,cond)~=0);
    if length(biomass)>1
        error('Detected more than 1 biomass reaction with non zero constraints')
    end
    ind_active=find(g_vect(:,cond)~=0);
    %Assign y variables to reactions with not nan abundance
    %Since the reversible modle contained reaction ending with _f and _b we
    %need match field as additional selector
    ind_b=find(~(cellfun(@isempty, regexp(model.rxns(ind_active),'_b$'))) & model.match(ind_active)~=0);
    ind_f=find(~(cellfun(@isempty, regexp(model.rxns(ind_active),'_f$'))) & model.match(ind_active)~=0);
   

    
    n_y=length(ind_active); 
    y0=n_r+1; %index where the y variable starts


    %%Maximizing sum of y

    %a) Inequality matrix
    %Only Y variable constraints only for reactions associated with genes of known abundance
    %Two inequalities of length(abun_genes) as rows, and n_r+n_y as columns
    %One inequality of length of reversible abundance reactions as rows
    I=eye(n_r);
    I_abun=I(ind_active,:);
    I_y=eye(n_y);
    

    Vmin=lb(:,cond)*scale;
    Vmax=ub(:,cond)*scale;
    Vmin(biomass)=lb(biomass,cond);
    Vmax(biomass)=ub(biomass,cond);

    A1=[ -I_abun, -(Vmin(ind_active)-ep).*I_y;
        I_abun, -(Vmax(ind_active)-ep).*I_y];

    A2=zeros(length(ind_f),n_y);
    for i=1:length(ind_f)
        A2(i,ind_f(i))=1;
        A2(i,ind_b(i))=1;
    end

    A=[A1;[zeros(length(ind_f),n_r),A2]];

    b=[-Vmin(ind_active);ep*ones(n_y,1);ones(length(ind_f),1)];

    Aeq=[model.S,zeros(n_m,n_y)];
    Aeq(:,biomass)=Aeq(:,biomass)*scale;
    beq=zeros(n_m,1);
   
    f=[zeros(1,n_r),-ones(1,n_y)]; %default gurobi is minimization thus the objective coefficient is -1 here instead of maximizin 1

    
    int=n_r+1:n_r+n_y;

    lower=[Vmin;zeros(n_y,1)];
    upper=[Vmax;ones(n_y,1)];

    
    m=struct();
    m.obj = f; % coefficienct for the objective function (c)
    m.A = [sparse(A); sparse(Aeq)]; % A must be sparse
    n = size(m.A, 2);
    m.vtype = repmat('C', n, 1);
    m.vtype(int) = 'I';
    m.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
    m.rhs = full([b(:); beq(:)]); % rhs must be dense
    m.lb=lower;
    m.ub=upper;
    
    params = struct();
    params.FeasibilityTol=feasTol;
    params.IntFeasTol=1e-3;
    params.OutputFlag=OutputFlag;
    params.Threads=Threads
    x = gurobi(m,params);
    
    %in case of infeasibility try supplying partial solution with higher
    %feasibility tolerance
    count=0;
    while count<50 && strcmp(x.status, 'INFEASIBLE')
        if count==0
            disp('Solution at feasTol was not found, Trying again supplying partial solution')
        end
        %set maximum feasibility tolerance to epsilon-0.1*epsilon to ensure
        %0 fluxes will have y=0
        partparams=params;
        partparams.FeasibilityTol=ep-0.1*ep;
        %solve
        x = gurobi(m, partparams);
        sol_y=find(x.x(y0:end)==1);
        %fix 10% of y==1 
        fix_idx=randsample(sol_y, floor(length(sol_y)/10), false);
        partsol=m;
        partsol.lb(y0-1+fix_idx)=1;
        partsol.ub(y0-1+fix_idx)=1;
        %run again with initial feasibility tolerance
        partparams.FeasibilityTol=feasTol;
        x=gurobi(partsol, partparams);
        count=count+1;
    end
%         
    
%     %formulate a glpk model - does not converge in condition 5 and 7 
%     if false
%         ctype= [repmat('U',size(A,1),1); repmat('S',size(Aeq,1),1)];
%         glpkparams.tolint=1e-3;
%         [xopt, fmin, status, extra]=glpk(m.obj, m.A, m.rhs, m.lb, m.ub, ctype, m.vtype, 1, glpkparams);
%     end

    if ~strcmp(x.status,'INFEASIBLE')
        sol=x.x;
        y=sol(y0:end);

        %Minimizing sum of fluxes
        Aeq2=[Aeq;zeros(1,n_r), ones(1,n_y)];
        beq2=[beq;sum(abs(y-1)<1e-3)];


        f2=zeros(1,n_r+n_y);
        f2(1:n_r)=1; %minimize the total sum of fluxes
        %%sol_new=intlinprog_adap(f2,int,A,b,Aeq2,beq2,lower,upper);
        m.obj = f2;
        m.A = [sparse(A); sparse(Aeq2)]; % A must be sparse

        m.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq2,1),1)];
        m.rhs = full([b(:); beq2(:)]); % rhs must be dense

        x2 = gurobi(m, params);
        count=0;
        while count<50 && strcmp(x2.status, 'INFEASIBLE')
            if count==0
                disp('Solution at feasTol was not found, Trying again supplying partial solution')
            end
            sol_y=find(x.x(y0:end)==1);
            %fix 10% of y==1
            fix_idx=randsample(sol_y, floor(length(sol_y)/10), false);
            partsol=m;
            partsol.lb(y0-1+fix_idx)=1;
            partsol.ub(y0-1+fix_idx)=1;
            
            x2=gurobi(partsol, params);
            count=count+1;
        end
        if ~strcmp(x2.status,'INFEASIBLE')
            sol_new=x2.x;
            y_new=sol_new(n_r+1:end);

            new.flux(:,cond)=sol_new(1:n_r);
            new.y(ind_active,cond)=y_new;
        end
    end
end
