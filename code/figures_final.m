addpath('code')


%% 2a
[altdat, sel, rr] = load_altdat('comorbidity_analysis_mendelian');
fig_2a(altdat, sel)

%% 2b
x = dlmread('figwork/skcm_rt.xls','\t',1,1);
md_c = struct('text',{x});
md = 'Rubinstein-Taybi'; cancer = 'SKCM'
keeps = md_c.text(:,4) > .05; % removing "pancancer"
sk_q = log_cap(md_c.text(keeps,3));
[x,ix] = sort(sk_q,'descend');
plot(1:length(ix), sk_q(ix),'b','linewidth',2); hold on
yl = ylim;
rt_b = md_c.text(keeps,1);
rt_b = rt_b(ix);
for i = 1:length(ix)
  if rt_b(i)==1
    plot([i,i], yl, 'r')
  end
end
plot(1:length(ix), sk_q(ix),'b','linewidth',2); hold on
legend({cancer, md},'Location','Northeast','FontSize',14)
xlabel(['Pathway index (ordered by ' cancer ' enrichment)'],'FontSize',14)
ylabel({'Cancer enrichment (-log10)'},'FontSize',14)
set(gca,'FontSize',14)

xlim([0,length(ix)])
ylim(yl)
h = gcf
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(h,'figwork/fig2b','-dpdf','-r0')


%% 2c
gene = 'PTK6';
nm = 'Ectodermal';
allm = {'EDARADD','LAMA3','ITGB4','KRT14','WNT10A','ITGA6','GJB6','EDA','PLEC','KRT5','LAMC2','LAMB3','COL7A1','EDAR','COL17A1'};
        
gene_coex = parseText(['figwork/' gene '.txt']);
gene_coex.coex = cell2num(gene_coex.text(:,2));

x = ismember(gene_coex.text(:,1),allm);
boxplot(gene_coex.coex, x, 'labels',{'All genes',nm},'symbol','')
%ylabel({'Pearson coexpression'; ['with ' gene]},'FontSize',14)        
ylabel([gene ' Pearson coexpression'],'FontSize',14)        
set(findall(gcf,'type','text'),'fontSize',14)

%% 3c
%%%%%%%%%%%%
rpl_mdm = parseText('figwork/rpl_mdm.txt','rowname',true);
rpl_mdm.status = cell2num(rpl_mdm.text);
rpls = rpl_mdm.status(1:3,:);
rpls(rpls > 0) = 0;
mdm = rpl_mdm.status(4,:);
mdm(mdm < 0) = 0;
rpl_stat = sum(rpls < 0,1);
both= sum(mdm > 0 & rpl_stat > 0);
mdm_sum = sum(mdm > 0);
any_rpl = sum(rpl_stat > 0);
for i = 1:size(rpls,1)
  fprintf('%s-MDM2 = %1.3f\n',rpl_mdm.rowname{i}, hygecdf(sum(rpls(i,:) < 0 & mdm > 0),length(mdm),mdm_sum,sum(rpls(i,:) < 0)));
end
fprintf('Any RPL = %1.3f\n',hygecdf(both,length(mdm),mdm_sum,any_rpl))

mdm(mdm > 1) = 1;
rpl_mdm_clean = [mdm; rpls([1,3,2],:)];
srt = sortrows(rpl_mdm_clean',1:4)';
HeatMap(srt, 'rowlabels',{'MDM2'; 'RPL5';'RPL11';'RPS7'},'standardize',3,'colormap',genColorMap('rwn',64))
pdf('figwork/3c')

%% figure 4
%4a

[altdat, sel, rr] = load_altdat('comorbidity_analysis_mendelian');

comorb = parseText('data_source/comorbidities','colname',true)
code_comorb = struct('rowname',{unique(altdat.peak_mut.MD(sel))},'colname',{unique(comorb.text(~ismember(comorb.text(:,1),{'ALCL','ALL','DLBCL','PTCL'}),1))});
code_comorb.text = zeros(length(code_comorb.rowname),length(code_comorb.colname));
for r = 1:length(code_comorb.rowname)
  for c = 1:length(code_comorb.colname)
    rsel = ismember(comorb.text(:,2),code_comorb.rowname(r)) & ismember(comorb.text(:,1), code_comorb.colname(c));
    if sum(rsel) == 1
      code_comorb.text(r,c) = cell2num(comorb.text(rsel,3));
    end
  end
end
code_comorb.clinical_score = sum(code_comorb.text > 1, 2);
bval = 0:max(code_comorb.clinical_score);
[ct, xv] = hist(code_comorb.clinical_score,bval);
pbval = binopdf(bval, length(code_comorb.colname), sum(sum(code_comorb.text > 1))/length(code_comorb.text(:)));

harray = bar(xv, [pbval', ct'/sum(ct)])
set(harray(1),'facecolor','blue')
set(harray(2),'facecolor','red')
legend({'expected at random','observed'},'Location','Northeast','FontSize',14)
xlabel('Number of cancer comorbidities per Mendelian Disease','FontSize',14)
ylabel('Fraction of Mendelian diseases','FontSize',14)
ylim([0,.3])
xlim([-.5,13.5])
set(gca,'FontSize',14,'TickLength',[ 0 .2 ],'XTick',[0:2:max(xv)])
pdf('figwork/fig_4a_distribution_comorb')

%4b
fig_4b('comorbidity_analysis_mendelian_germline');

%% S1
% s1a
ng_m = parseText('figwork/numg_per_icd.txt');
ng_m.ct = cell2num(ng_m.text(:,2));
hist(ng_m.ct)
xlabel('Num genes per Mendelian disease group','FontSize',14)
ylabel('Num Mendelian disease groups','FontSize',14)
set(gca,'FontSize',14)
pdf('figwork/s1a_numg_per_MD')

% s1b
ng_c = parseText('figwork/numg_per_c_peakmut.txt');
ng_c.ct = cell2num(ng_c.text(:,2));
csel = ~ismember(ng_c.text(:,1), {'CENSUS','MY_PAN_CAN'});
hist(ng_c.ct(csel))
xlabel('Num genes per cancer','FontSize',14)
ylabel('Num cancers','FontSize',14)
set(gca,'FontSize',14)
pdf('figwork/s1b_numg_per_canc_focal')

%% S2
[altdat, sel, rr] = load_altdat('comorbidity_analysis_mendelian');

%%% 
subplot(2,2,1)
randct = dlmread('comorb_aggregate/rand_ge.100000.txt');
rand_vs_true(randct, 41)
title('Genes shared')

subplot(2,2,2)
randct = dlmread('comorb_aggregate/rand_path_b.100000.txt');
rand_vs_true(randct, 136)
title('Pathways shared')

subplot(2,2,3)
randct = dlmread('data_processed/biogrid/randcts');
rand_vs_true(randct, 797)
title('BioGRID edges')

subplot(2,2,4)
randct = dlmread('data_processed/humannet.9.unwt/randcts');
rand_vs_true(randct, 296)
title('HumanNet edges')

h = gcf
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(h,'rand_vs_true_comorbid','-dpdf','-r0')

h = gcf
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(h,'biogrid_distr','-dpdf','-r0')



%% s3
[altdat, sel, rr] = load_altdat('comorbidity_analysis_mendelian');
manual_predictors = [altdat.peak_mut.gene_enrichment < .05, altdat.peak_mut.pathway_correlation < .05, altdat.peak_mut.NC_set < .05,altdat.peak_mut.coex_CG < .05 & altdat.peak_mut.coex_CG > -1];
pred_sel = manual_predictors(sel,:);

pnames = {'Gene enrichment','Pathway correlation','Network neighbors','Coexpression'}
fracs = sum(pred_sel(rr,:))/sum(rr);
fracs_nocomorb = sum(pred_sel(~rr,:))/sum(~rr);
colors = ['r','k','m','b'];

% s3a
for i = 1:length(pnames)
  subplot(2,2,i)
  %harray = bar([[fracs(i),1-fracs(i)]; [fracs_nocomorb(i),1-fracs_nocomorb(i)]],.2,'stacked')
  %set(harray(1),'facecolor','y')
  %set(harray(2),'facecolor','k')
  harray = bar([fracs(i); fracs_nocomorb(i)],.2,colors(i))
  set(gca,'XTickLabel',{'Comorbid','Not comorbid'},'FontSize',10)
  p = 1 - hygecdf(sum(pred_sel(rr,i))-1, size(pred_sel,1), sum(rr), sum(pred_sel(:,i)));
  xlim([.5,2.5])
  title({pnames{i}; [' p = ' num2str(p,'%1.4f')]})
end
pdf('figwork/s3a'); close
% s3b
harray = bar([fracs; fracs_nocomorb]);
p = zeros(1,4)
for i = 1:length(p)
  set(harray(i),'facecolor',colors(i))
  p(i) = 1 - hygecdf(sum(pred_sel(rr,i))-1, size(pred_sel,1), sum(rr), sum(pred_sel(:,i)));
end
legend(strcat(pnames', ' p = ',num2cellstr(p')),'location','eastoutside','FontSize',14)
set(gca,'XTickLabel',{'Comorbid','Not comorbid'},'FontSize',14)
ylabel('Fraction of disease pairs','FontSize',14)
pdf('figwork/s3b'); close