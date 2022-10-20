function out_model = convertToUniprot(model)
    out_model=model;
    out_model.genes=JGItoUP(model.genes);
    if isequal(model.genes, model.geneNames)
        out_model.geneNames=out_model.genes
    end
     %check for whitespace or empty gene names
    ws=find(contains(out_model.genes, {' ', ';'}));
    emt=find(cellfun(@isempty, out_model.genes));
    if  ~(isempty(ws) && isempty(emt))
        disp('Corrupted gene names detected:')
        disp(out_model.genes([ws, emt]))
        error('Aborting...')
    end
    if any(duplicates(out_model.genes), 'all')
        dup=[find(sum(duplicates(out_model.genes))), find(sum(duplicates(out_model.genes), 2))'];
        disp('Duplicated gene names introduced')
        dup_gn=table(out_model.genes(dup), model.genes(dup), 'VariableNames', {'New_ID', 'Old_ID'});
        [~,sidx]=sort(dup_gn.New_ID);
        dup_gn=dup_gn(sidx,:);
        disp(dup_gn)
        %remove duplicates 
        out_model=fixGenes(out_model, true);
    end
    %copy the content of genes field to geneUniprotID field
    out_model.geneUniprotID=out_model.genes;
end
    