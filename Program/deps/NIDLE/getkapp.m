function [Kapp_matrix,V_matrix,count, Kapp_genes]=getkapp(abundance,g_vect,sol_flux, model_irrev,lit_kcats, conv,thre, eps, outf)
%Function to estimape Kapp values for non-complex catalysed reactions
%from condition wise fluxes and abundancies, if only one catalysing isomer
% is quantified the kapp is calculated by v(i)/E(i) 
%if all isozymes are expressed a L1 estimate is found for the linear system 
%of equations including all measurements of mixotrophic standard conditions in which 
% all isozymes are quantified condition (n_condition>= n_isozymes), while enforcing
%kapps to be >=0: 
%A(i,j)=Concentration of isomere j, in condition i
%B(i) = flux of reaction in condtition i
%argmin||A*x-B||
% Kapp_genes is only updated with the obtained x (=kapp) if there is no
% condition specific estimate from practicall homomeric reactions
% for Kapp_matrix kapp average weighted by the isoenzyme concentrations 
% is saved
%INPUT:
%  struc abundance: An abundance structure giving the abundancies and
%     reactions linked to homomeric enzymes
%  logical g_vect: a logical matrix giving the reactions linked to a enzyme 
%     dectected abundance in this condition
%  double sol_flux: matrix giving the fluxes returned from NIDLE for the 
%     given conditions
%  struct model_irrev: The COBRA modle on which simulations are based
%  double lit_kcats: A double vector giving literature kcats from Brenda and 
%     Sabio RK for each reaciton in the irreversible modle (if nothing
%     is found the element is NaN
%  double conv: The conversion vector which fluxes have to be scaled with to
%     get mmol/gDWh
%  double thre: Threshold below which converted fluxes are considered
%     to be 0
%  double eps: Minimum bound for kapps in L1 estimate of isozyme kapp
%  char outf: (optional) if given The filepath to which the results should be
%     saved as tsv
%OUTPUT:
% num Kapp_matrix:  A matrix with rows equal to number of reactions and
%                   columns equal to number of conditions bearing the
%                   determined kapp values in s^-1 (saved under <outf>.tsv)
% num V_matrix: A matrix of fluxes converted back to mmol/gDWh
% struct count: Summary statistics for each condition(see inline comments
%               for details on the fields
% num Kapp_genes:   A matrix where each row represents a protein and columns
%                   correspond to conditions, where entries give the maximum 
%                   kapp obtained for all catalysed reactions in the
%                   respective condition. This is saved to
%                   <outf>_genewise.tsv

%test the isozLP function
testE=[1,2,1;2,1,3;2,3,2];
testv=[18; 37; 33];
if  sum(abs(isoz_LP(testv, testE).kapp-[2;3;10]))>1e-8 || sum(abs(isoz_QP(testv, testE).kapp-[2;3;10]))>1e-8
    error('isoz_LP function does not return correct solution for soluble system')
end
%Generate output dir if not existing
if ~isfolder(fileparts(outf))
    mkdir(fileparts(outf))
end
%cleanup sol_flux 
%Result for fluxes reconverted to mmol/gDW*h
V_matrix=sol_flux;
V_matrix(V_matrix<=thre)=0;
V_matrix(~(contains(model_irrev.rxns, 'Biomass')),:)=V_matrix(~(contains(model_irrev.rxns, 'Biomass')),:)*conv;


%reactions with only 1 expressed catalysing enzymes
hom_idx=cellfun(@length, abundance.GPR)==1;
%isomeric reactions with several quantified enzyme
isz_idx=cellfun(@length, abundance.GPR)>1;

R_index=abundance.reacind(hom_idx);
[~,gidx]=ismember(vertcat(abundance.GPR{hom_idx}), abundance.genes);
R_ab_matrix=abundance.abun(gidx(gidx~=0),:);
R_names=abundance.genes(gidx(gidx~=0));
n=size(R_ab_matrix,2);
%all reactions catalyzed by a given gen
cat_rec=cell(length(abundance.genes), 1);
for i=1:length(abundance.genes)
    cat_rec{i}=abundance.reacind(cellfun(@(x) any(ismember(x, abundance.genes(i))), abundance.GPR));
end

%estimated Kapp values NaN for values without protein abundance 0 for
%values without flux
Kapp_matrix=nan(size(sol_flux,1),n);

%maximum kapp for each enzyme of all reactions it can catalyze in the given
%condition
Kapp_genes=nan(length(abundance.genes),n);

%number of reactions with flux larger than the threshold
count.nonzero=zeros(n,1); 

%number of reactions with at least one isoenzyme with non-nan abundance 
count.withabun=zeros(n,1);
%number of reactions with flux larger than the threshold with non-nan abundance 
%enzymes
count.nonzero2=zeros(n,1); 


%number of reactions for which only one enzyme had non-zero abundance -->
%determined kapp
count.kapp=zeros(n,1);

%numer of ENZMYES  with determined kapp
count.homoenzyme2=zeros(n,1);

%number of genes with isozyme GPR rules
count.isozymes=sum(isz_idx);
%number of isoenzyme reactions with only one expressed isoenzyme for which
%kapp can be determined
count.singisozymes=zeros(n,1);

%number of isoenzymes reactions with multiple expressed isoenzymes so not
%kapp can be easily obtained
count.multisoenzyme=zeros(n,1);

%number of isoenzyme reactions without flux
count.nofisoenzyme=zeros(n,1);
%number of active reactions with active gprRules versus total reactions
%with active rule
count.ratio=zeros(n,1);


for cond=1:n
    V=V_matrix(:,cond);
    const=find(g_vect(:,cond)==1);
    %calculate kapp for homomeric enzymes
    Kapp_matrix(R_index,cond)=V_matrix(R_index,cond)./(R_ab_matrix(:,cond)*3600);
    %get enzyme centered maximum kcat
    Kapp_genes(:,cond)=cellfun(@(x) max(Kapp_matrix(x,cond),[], 'omitnan'), cat_rec);
    count1=0;
    count2=0;
    count3=0;
    for i=find(isz_idx)'
        %if only one isoenzyme  is expressed calculate apparent catalytic
        %in the same way as for the homomeric enzymes
        [~, g_GPR_idx]=ismember(abundance.GPR{i}, abundance.genes);
        temp_abun=abundance.abun(g_GPR_idx,cond);
        if sum(~isnan(temp_abun))==1
            tmp_kapp=V_matrix(abundance.reacind(i), cond)/(temp_abun(~isnan(temp_abun))*3600);
            Kapp_matrix(abundance.reacind(i),cond)= tmp_kapp;
            Kapp_genes(ismember(abundance.genes, abundance.GPR{i}(~isnan(temp_abun))), cond)= ...
                max([Kapp_genes(ismember(abundance.genes, abundance.GPR{i}(~isnan(temp_abun))), cond), tmp_kapp], [], 'omitnan'); 
            count1=count1+1;
        elseif sum(~isnan(temp_abun))>1 && V(abundance.reacind(i))>0
            count2=count2+1;
%         elseif V(abundance.reacind(i))>0 %This heuristic is very
%         questionable since the total abundance of protein is assumed to
%         partake in both reactions... it also only leads to 1 determined
%         kapp + in control conditions so it is skipped
%             %if a kapp is already nown for the other isozemes in the
%             %reaction excep substract the vmax from vi and calculate the
%             %result form the difference
%             tmp_isz=abundance.GPR{i}(~isnan(temp_abun));
%             [~, g_tisz_idx]=ismember(tmp_isz, abundance.genes);
%             known_kapp=Kapp_genes(g_tisz_idx,cond);
%             if sum(isnan(known_kapp)|known_kapp==0)==1
%                 %calculate maximum capacity for other isozymes
%                 isz_cap=sum(known_kapp(~(isnan(known_kapp)|known_kapp==0))*3600.*...
%                     abundance.abun(g_tisz_idx(~(isnan(known_kapp)|known_kapp==0)), cond));
%                 %divide difference to total flux to obtain kapp 
%                 tmp_kapp=(V_matrix(abundance.reacind(i),cond)-isz_cap)/(abundance.abun(g_tisz_idx(~(isnan(known_kapp)|known_kapp==0)), cond)*3600);
%                 if tmp_kapp>0
%                     count2=count2+1;
%                     Kapp_genes(g_tisz_idx(~(isnan(known_kapp)|known_kapp==0)), cond)=tmp_kapp;
%                 end
%            end
        else
            count3=count3+1;
        end
    end
    count.nonzero(cond)=sum(V>0);
    count.nonzero2(cond)=sum(V(unique(vertcat(cat_rec{~isnan(abundance.abun(:,cond))})))>0);
    count.kapp(cond)=length(intersect(find(Kapp_matrix(:,cond)~=0),find(~isnan(Kapp_matrix(:,cond)))));
    count.homoenzyme2(cond)=sum(~isnan(Kapp_genes(:,cond)) & Kapp_genes(:,cond)~=0);
    count.withabun(cond)=length(unique(vertcat(cat_rec{~isnan(abundance.abun(:,cond))})));
    count.ratio(cond)=sum(V(const)>0)/sum(g_vect(:,cond)~=0);
    count.singisozymes(cond)=count1;
    count.multisoenzyme(cond)=count2;
    count.nofisoenzyme(cond)=count3;
end
%check if number of obtained kcats + multi isozyme reactions sum up to
%number of reacitons for which flux and abundance is available
if any((count.kapp+count.multisoenzyme-count.nonzero2)~=0, 'all')
    error('number of determined kapps is not equal to the number of reactions for which abundance and non-zero flux is available minus multi-isoenzyme reactions. This is an internal error')
end

%to Obtain kapp estimates for reaction with multiple expressed isozymes we
%assume the kapp between the 4 mixotrophic control conditions is identical
%since this is an approximation we only update kapps for isozymes
%that were not assigned any kapp in any condition yet
upd_idx=all(isnan(Kapp_genes)|Kapp_genes==0, 2);
stand_mix=ismember(abundance.cond, {'control', 'UVM4', 'Stop2'});
stand_mix_idx=find(stand_mix);
R2_dist=[];
RMSD_dist=[];
fit_isz_idx=[];
QP_kapp=[];
LP_kapp=[];
count4=0;
for i=find(isz_idx)'
    if all(V_matrix(abundance.reacind(i), stand_mix)>0, 'all') ... in case reaction carries flux in all condtions
            &&sum(all(~isnan(abundance.abun(ismember(abundance.genes, abundance.GPR{i}), stand_mix)),1))>=length(abundance.GPR{i}) % and all isoenzymes are quantified in as many conditions as there are isoenzymes
        count4=count4+1;
        complete_cond=all(~isnan(abundance.abun(ismember(abundance.genes, abundance.GPR{i}), stand_mix)),1);
        %calculate the LSQR estimate of kapps from the
        %resulting linear system of equations
        tmp_g_idx=ismember(abundance.genes, abundance.GPR{i});
        Etemp=abundance.abun(tmp_g_idx, stand_mix_idx(complete_cond))';
        vtemp=V_matrix(abundance.reacind(i), stand_mix_idx(complete_cond))';
        sol=isoz_QP(vtemp, Etemp, eps);
        LP_sol=isoz_LP(vtemp, Etemp, eps);
        QP_kapp=[QP_kapp;sol.kapp];
        LP_kapp=[LP_kapp;LP_sol.kapp];
        %kapp=lsqr(A, aprod(1), B)
        kapp=sol.kapp;

        %And the quality of fit
        R2=1-sum((Etemp*kapp-vtemp).^2)/sum((vtemp-mean(vtemp)).^2);
        R2_dist=[R2_dist, R2];
        %and root mean squared deviation
        RMSD=sqrt(sum((Etemp*kapp-vtemp).^2)/length(vtemp));
        RMSD_dist=[RMSD_dist RMSD];
        %If kapps are larger than eps
        if all(kapp>eps)
            fit_isz_idx=[fit_isz_idx, i];
            %set the gene specific kapp for updatable isozymes to the maximum
            %of the obtained and stored value
            tmp_g_idx=find(tmp_g_idx);
            for j=find(upd_idx(tmp_g_idx))'
                Kapp_genes(tmp_g_idx(j),stand_mix_idx(complete_cond))=repmat(max([kapp(j)/3600,  Kapp_genes(tmp_g_idx(j),find(stand_mix_idx,1))], [], 'omitnan'), 1, sum(complete_cond));
                %check again if all values are identical 
                if ~(all(isnan(Kapp_genes(tmp_g_idx(j),stand_mix_idx(complete_cond))), 'all') || length(unique(Kapp_genes(tmp_g_idx(j),stand_mix_idx(complete_cond))))==1)

                    error('Gene specific kapp storage array contains different values for isozymes that are determined over multiple conditions. This is an internal error')
                end
            end
            complete_cond_idx=find(complete_cond);
            for j=1:sum(complete_cond)
                %set the kapp for each reaction to the average of the obtained kapp
                %weighted by the protein abundances in the given condition
                Kapp_matrix(abundance.reacind(i),stand_mix_idx(complete_cond_idx(j)))=(Etemp(j,:)*kapp)/sum(Etemp(j,:))/3600;
            end
        else
            disp('Excluded solution')
            disp(kapp)
        end
    end
end
        %Number of isomeric reactions for which data in all mixotrophic
        %standard conditions was available so that a model could be fit
        count.lsqrfit=count4;
        %Number of isomeric reactions for which a valuable solution was
        %obtained
        count.lsqrgood=length(fit_isz_idx);
        %Export a scatterplot of L1 vs L2 values
        scatter(log10(LP_kapp), log10(QP_kapp), 'fill')
        ylabel('log_{10} min L2 kapp')
        xlabel("log_{10} min L1 kapp")
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0 0 4 4]);
        saveas(gcf, fullfile(fileparts(outf), 'LPvsQP.svg'))
        close


%         %save a RMSD scatter plot marking the solutions that have
%         negative - REMOVED WITH Epsilon constrain no R2 filter
%         %R2
%         %plot a poxplot of log10 values 
%         boxplot(log10(RMSD_dist))
%         hold on 
%         %plot the removed RMSD of solutions with negative R2
%         scatter(ones(sum(R2_dist<=0),1).*(1+(rand(sum(R2_dist<=0),1)-0.5)/10), log10(RMSD_dist(R2_dist<=0)), 'r', 'filled')
%         title('Removed isoenzyme kapp estimates')
%         ylabel('RMSD kapp*E-v')
%         set(gcf, 'PaperUnits', 'inches');
%         set(gcf, 'PaperPosition', [0 0 4 4]);
%         saveas(gcf, fullfile(fileparts(outf), 'excluded_isozymekapp.svg'))
%         close
        %plot density plots of kapp values for homomeric/1-isozyme enzmes
        %vs isozyme enzymes
        
        ksdensity(log10(Kapp_matrix(setdiff(find(~isnan(Kapp_matrix(:,1))& Kapp_matrix(:,1)~=0), abundance.reacind(fit_isz_idx)),1)))
        hold on
        ksdensity(log10(Kapp_matrix(abundance.reacind(fit_isz_idx),1)))
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0 0 4 4]);
        saveas(gcf, fullfile(fileparts(outf), 'cond1_homvsiszkappdist.svg'))
        close
 %Export the estimated Kapps
 if exist('outf', 'var')
     %export nidle flux solution
     %array indicating if reaction kapp was calculated using the isozyme QP
     isz_calc=zeros(size(V_matrix,1),1);
     isz_calc(abundance.reacind(fit_isz_idx))=1;
     export_V=[table(model_irrev.rxns, isz_calc),array2table(V_matrix, 'VariableNames', abundance.cond)];
     writetable(export_V, append(outf, '_flux.tsv'), 'FileType', 'text', 'Delimiter', 'tab')
     if ~exist(fileparts(outf), 'dir')
         mkdir(fileparts(outf))
     end
     Res_idx=any(~isnan(Kapp_matrix) & ~(Kapp_matrix==0),2);
     res_rxns= model_irrev.rxns(Res_idx);
     model_irrev=creategrRulesField(model_irrev);
     res_JGIg= model_irrev.grRules(Res_idx);
     up_mod=model_irrev;
     up_mod=rmfield(up_mod, 'grRules');
     up_mod=convertToUniprot(up_mod);
     res_UPg= up_mod.grRules(Res_idx);
     res_EC= model_irrev.rxnECNumbers(Res_idx);
     res_lkcat= lit_kcats(Res_idx,:);
     res_Kcats= Kapp_matrix(Res_idx,:);
     res_iszcalc=isz_calc(Res_idx);
     %export a overview table
     res_tab=[table(res_rxns, res_iszcalc, res_JGIg, res_UPg, res_EC),res_lkcat, array2table(res_Kcats)];
     res_tab.Properties.VariableNames=[{'Rxns', 'IsozymeQP', 'JGIgeneID', 'UniprotgeneID', 'ECNumber'}, lit_kcats.Properties.VariableNames, abundance.cond'];
     writetable(res_tab, append(outf, '.tsv'), 'FileType', 'text', 'Delimiter', 'tab')
     %export gene overview table
     gene_tab=[table(abundance.genes(~all(isnan(Kapp_genes) | Kapp_genes==0, 2)), 'VariableNames', {'ID5_5'}), array2table(Kapp_genes(~all(isnan(Kapp_genes) | Kapp_genes==0, 2),:), 'VariableNames', abundance.cond)]; 
     writetable(gene_tab, append(outf, '_genewise.tsv'), 'FileType', 'text', 'Delimiter', 'tab')
     %Old Code for deprecated output
     %      %export a table for integration as manual kcat modification into GECKO
%      gko_tab=table(res_rxns, res_UPg, res_EC, res_UPg, max(res_Kcats, [], 2));
%      writetable(gko_tab, fullfile(fileparts(outf), 'manual_data.txt'), 'Delimiter', 'tab', 'WriteVariableNames', false)
%      %take he maximum kcat for proteins
%     for prot=unique(gko_tab.res_UPg)'
%         gko_tab.Var5(ismember(gko_tab.res_UPg, prot))=max(gko_tab.Var5(ismember(gko_tab.res_UPg, prot)));
%     end
%     writetable(gko_tab, fullfile(fileparts(outf), 'manual_data.txt'), 'Delimiter', 'tab', 'WriteVariableNames', false)

 end
 disp('Returned Kapp values are in the unit [s⁻1] to be used in the context of models they have to be reconverted to [h^-1]')
 
end