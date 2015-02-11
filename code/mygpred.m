%function [to_write, ind_order] = mygpred(altdat, rr, sel, add_con, ml, rem_blank, fname)
function [to_write, ind_order] = mygpred(altdat, rr, sel, add_con, pred_struct, rem_blank, fname)

netfind = {'NC'};%,'SP'};%,'SP_SS'}; %,'dgidb_NC','dgidb_SP'};
%coexpr = ;
mysets = [{'gene_intersection','gene_enriched','pathway_genes'}, netfind, 'C_coexpr'];
gene_tab = {'MD','C'};
isel = find(sel);
alts = fieldnames(altdat);
alts = pred_struct.alts;
pred_mat = pred_struct.pred_mat;
gres = [altdat.(alts{1}).MD(isel), altdat.(alts{1}).C(isel), cell(length(isel),length(mysets)*(length(alts)))];
ind = 3;
ind_order = [];
csel = zeros(1,size(gres,2));

for alt = 1:length(alts)
  my_pred = pred_mat(:,alt) == 1;
  %my_pred = altdat.(alts{alt}).pred==1;
  %if strcmp(ml,'LR')  my_pred = altdat.(alts{alt}).logistic_regression.prediction; end
  altsel = isel(my_pred(sel) & rr) ;% & altdat.(alts{alt}).relative_risk(isel) > 1);
  put_at = find(my_pred(sel) & rr);
  ind_order = [ind_order; put_at(~ismember(put_at, ind_order))];
  fprintf('%s %d\n',alts{alt},length(altsel));
  

  %%%%% now add in gene info
  glist = altdat.(alts{alt}).gene_intersection(altsel);
  rsel = ~ismember(glist,{'','NA'});
  gres(put_at(rsel),ind) = glist(rsel);

  glist = altdat.(alts{alt}).gene_intersection(altsel);
  rsel = ~ismember(glist,{'','NA'}) & altdat.(alts{alt}).gene_enrichment(altsel) < .05;
  gres(put_at(rsel),ind + 1) = glist(rsel);
  %if ismember('pathway_overlap', altdat.(alts{alt}).logistic_regression.predictors)
  %if ismember('pathway_correlation', altdat.(alts{alt}).logistic_regression.predictors)
    glist = altdat.(alts{alt}).pathway_genes(altsel);
    gcon = altdat.(alts{alt}).pathway_sel(altsel);
    gcon2 = altdat.(alts{alt}).pathway_cancer_gene(altsel);
    %rsel = ~ismember(glist,{'','NA'}) & altdat.(alts{alt}).pathway_correlation(altsel) < .05;
    rsel = altdat.(alts{alt}).pathway_correlation(altsel) < .05;
    gadd = glist(rsel);
    if add_con 
      idx = find(rsel); 
      con_txt = cell(length(idx), 1);

      gadd2 = '';
      for pidx = 1:length(idx)
        pathways = strsplit(gcon{idx(pidx)},';');
        pathway_cancer_gn = strsplit(gcon2{idx(pidx)},';');
        pdiddy = strcat(pathways , '(', pathway_cancer_gn, ')');
        con_txt{pidx} = join_str(pdiddy,';');
        
        % yes?
        gadd_split = strsplit(gadd{pidx},';');
        gadd{pidx} = join_str(strcat(gadd_split ,'->', pathways, '(', pathway_cancer_gn, ')'),';');
        if strcmp(gadd{pidx},'NA->NA(NA)') gadd{pidx} = 'correlation only'; end
      end
      %gadd = strcat(gadd, '//', con_txt);
      
    end
    gres(put_at(rsel),ind+2) = gadd;
    %else csel(ind+1) = 1; end

  for nf = 1:length(netfind)
    %if ~ismember([netfind{nf} '_set'], altdat.(alts{alt}).logistic_regression.predictors) 
      fprintf('NO %s\n',[netfind{nf} '_set']);
      %csel(ind+1) =1 ; continue; 
    %end
    glist = altdat.(alts{alt}).([netfind{nf} '_genes'])(altsel);
    gcon = altdat.(alts{alt}).([netfind{nf} '_ginfo'])(altsel);
    rsel = ~ismember(glist,{'','NA'}) & ...
           (altdat.(alts{alt}).([netfind{nf} '_set'])(altsel) < .05 | ...
            repmat(strcmp(netfind(nf), 'dgidb_NC'), length(glist),1));
    gadd = glist(rsel);
    if add_con 
      idx = find(rsel);
      for pidx = 1:length(idx)
        gadd_split = strsplit(gadd{pidx},',');
        gcon_split = strsplit(gcon{idx(pidx)},';');
        gadd{pidx} = join_str(strcat(gadd_split, '->',gcon_split),';');
      end
      %gadd = strcat(gadd, '//', gcon(rsel)); 
    end
    gres(put_at(rsel),ind +nf+2) = gadd;
  end
  
  column = ind + length(netfind) + 3;
  %{
  noncomexpr = altdat.(alts{alt}).max_coexpression(altdat.(alts{alt}).relative_risk > -1 & ...
                                                   altdat.(alts{alt}).relative_risk < 1);
  glist = altdat.(alts{alt}).coexpr_max_M(altsel);
  gcon = altdat.(alts{alt}).coexpr_max_C(altsel);
  %lr = altdat.(alts{alt}).logistic_regression;
  rsel = ~ismember(glist,{''}) & altdat.(alts{alt}).max_coexpression(altsel) > (mean(noncomexpr) + std(noncomexpr)); %lr.covar_mat_binary(altsel,ismember(lr.covar_mat_names, 'max_coexpression'));
  gadd = glist(rsel);
  if add_con gadd = strcat(gadd, '//', gcon(rsel)); end
  gres(put_at(rsel), column) = gadd;
  %}
  glist = altdat.(alts{alt}).coexpr_cg_M(altsel);
  gcon = altdat.(alts{alt}).coexpr_cg_C(altsel);
  rsel = ~ismember(glist,{''}) & altdat.(alts{alt}).coex_CG(altsel) < .05 & altdat.(alts{alt}).coex_CG(altsel)  > -1; 
  gadd = glist(rsel);
  %if add_con gadd = strcat(gadd, '//', gcon(rsel)); end
  if add_con 
      idx = find(rsel);
      for pidx = 1:length(idx)
        gadd_split = strsplit(gadd{pidx},';');
        gcon_split = strsplit(gcon{idx(pidx)},';');
        gadd{pidx} = join_str(strcat(gadd_split, '->',gcon_split),';');
      end
      %gadd = strcat(gadd, '//', gcon(rsel)); 
  end

  gres(put_at(rsel), column) = gadd;
  %{
  glist = altdat.(alts{alt}).coexpr_mg_M(altsel);
  gcon = altdat.(alts{alt}).coexpr_mg_C(altsel);
  rsel = ~ismember(glist,{''}) & altdat.(alts{alt}).coex_MG(altsel) < .05; %lr.covar_mat_binary(altsel,ismember(lr.covar_mat_names, 'mean_max_coexpression'));
  gadd = glist(rsel);
  if add_con gadd = strcat(gadd, '//', gcon(rsel)); end
  gres(put_at(rsel), column + 2) = gadd;
  %}
  %%% finish iteration
  gene_tab = [gene_tab, strcat(alts{alt}, '.', mysets)];
  ind = ind + length(mysets);
end

to_write = [gene_tab(csel==0);gres(ind_order,csel==0)];
ind_order_noblank = ind_order;
if rem_blank
  keep = sum(cellfun('size',to_write(:,~strfind_bool(to_write(1,:),'dgidb')),2) > 0,2) > 2;
  keep_col = sum(cellfun('size',to_write,2) > 0,1) > 1;
  to_write= to_write(keep,keep_col);
  ind_order_noblank = ind_order_noblank(keep(2:end));
end
writeCellTxt(to_write, ['figwork/mygpred.' fname '.' num2str(add_con) '.' num2str(rem_blank) '..xls'])
