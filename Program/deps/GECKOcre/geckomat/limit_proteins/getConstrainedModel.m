function [ecModel_batch,OptSigma, changes] = getConstrainedModel(ecModel,modifications,name, raw, parameters)
% getConstrainedModel
%
% Function that gets a GEM with kinetic data and returns an enzyme 
% constrained model, either with individual enzyme levels or with the total
% measured protein content.
% 
% The model parameters are automatically curated (querying the BRENDA
% files) for reaching the experimental maximal growth rate on glucose
% minimal media (if it was overconstrained). Then, the average saturation 
% factor for enzymes is fitted for the same growth conditions. Finally the 
% fitted model is simulated on the same conditions and the top ten used 
% enzymes (mass-wise) are saved in a file as an output, this can be used for 
% suggesting further parameters curation (enzyme usages > 10% of the total 
% proteome).
%
% INPUT:
%   ecModel         An enzyme constrained model.
%   modifications	List of manually changed kcats (will be skipped in the
%                   kcat flexibilization).
%   name            Name for the ecModel, should be the same as the
%                   directory in which the model output files should be stored.
%   raw    (optional)if true a raw batch model without automatic kcat correction
%             and adapted uptake rate is returend; alsosigma is not refitted 
%             (default = false)
%   parameters  (optional) A structure from which parameters can be supplied instead
%                        of getModelParameters 
%
% OUTPUT:
%   ecModel_batch	The enzyme constrained model under batch conditions.
%   OptSigma        Optimized saturation factor.
%   changes T   able  containing information on the changed kcats
%
% Usage: [ecModel_batch,OptSigma] = getConstrainedModel(ecModel,parameters,modifications,name)
% 
% Marius Arend             2021-10-13
% Benjamin J. Sanchez   2018-12-11
% Ivan Domenzain        2019-07-13
%
%

% default is to apply kcat correction
if nargin <4
    raw = false;
end
%Get model parameters values if no parameters argument is passed 
if nargin< 5
    cd ..
    parameters = getModelParameters;
end
    c_source   = parameters.c_source;
    sigma      = parameters.sigma;
    Ptot       = parameters.Ptot;
    gRate      = parameters.gR_exp;
    cd limit_proteins

%Get f (estimated mass fraction of enzymes in model)
[f,~] = measureAbundance(ecModel.enzymes);

disp('Fitted f:')
disp(f)
warning('Kcats of model are optimized using the input bounds')
% if ~raw
%     %Change media to batch conditions
%     cd ../kcat_sensitivity_analysis
%     ecModel = changeMedia_batch(ecModel,c_source);
%     cd ../limit_proteins
% end

%Get a preliminary enzyme constrained model for performing the Kcats
%sensitivity analysis
%here protein content of biomass is rescaled and GAM is fitted
[ecModel_batch,~,~] = constrainEnzymes(ecModel,f, raw);

if raw 
    OptSigma=sigma;
    %stop function here since model is not adapted for optimization
    changes=[];
    return
end
solution            = solveLP(ecModel_batch,1);
if ~isempty(solution.f)
    %Set the media according to the experimental conditions
    cd ../kcat_sensitivity_analysis
    ObjIndex = find(ecModel_batch.c);
    % If the model is overconstrained
    if (gRate-solution.x(ObjIndex))>0
        disp('The ECmodel is overconstrained!')
        %Perform a sensitivity analysis on the objective function with
        %respect to the individual Kcat coefficients, the algorithm will
        %iterate replacing the top limiting value according to the maximum
        %value available in BRENDA for the same EC number until the objective
        %is no longer underpredicted
        [ecModel_batch, changes] = modifyKcats(ecModel_batch,gRate,modifications,name);
    else
        disp('The ECmodel is not overconstrained.')
        changes=[];
    end
    %The sigma factor is reffited for the specified conditions (constraints in the model)
    disp('***************************************************************')
    OptSigma          = sigmaFitter(ecModel_batch,Ptot,gRate,f);
    enzymePos         = strcmp(ecModel_batch.rxns,'prot_pool_exchange');
    currentEnzymeUB   = ecModel_batch.ub(enzymePos);
    newEnzymeUB       = currentEnzymeUB*OptSigma/sigma;
    ecModel_batch     = setParam(ecModel_batch,'ub','prot_pool_exchange',newEnzymeUB);
    %Simulate growth on minimal media and export to the output folder:
    % 1) the exchange fluxes to the file "exchangeFluxes.txt"
    % 2) the top ten used enzymes to the file "topUsedEnzymes.txt"
    solution = solveLP(ecModel_batch,1);
    if ~isempty(solution.x)
        disp('Saving simulation results files...')
        fluxFileName = ['../../models/' name '/' name '_exchangeFluxes.txt'];
        printFluxes(ecModel_batch,solution.x,true,10^-6,fluxFileName);
        topUsedEnzymes(solution.x,ecModel_batch,{'Min_glucose'},name);
    end 
    cd ../limit_proteins   
else
    disp('ecModel with enzymes pool constraint is not feasible')
end
end
