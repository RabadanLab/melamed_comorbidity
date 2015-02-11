function fig_4b(path)
%
%'generemove_binarypath/report_pairs.xls'
[altdat, sel, rr] = load_altdat(path);
unk_canc = {'THCA','LAML','OV','HNSC'};
pairs = parseText([path '/report_pairs.xls'],'colname',true);
pairs.gene_enrichment = cell2num(pairs.text(:,ismember(pairs.colname, 'gene_enrichment')));
pairs.pathway_correlation = cell2num(pairs.text(:,ismember(pairs.colname, 'pathway_correlation')));
pairs.coex_CG = cell2num(pairs.text(:,ismember(pairs.colname, 'coex_CG')));
pairs.biogrid = cell2num(pairs.text(:,ismember(pairs.colname, 'biogrid_set')));
pairs.humannet_set = cell2num(pairs.text(:,ismember(pairs.colname, 'humannet_set')));
pairs.MD = pairs.text(:,ismember(pairs.colname, 'MD'));
pairs.C = pairs.text(:,ismember(pairs.colname, 'C'));
[u,c] = uniq_count_cell(pairs.MD);

md3 = u(c > 3);
[altdat, sel, rr] = load_altdat(path);
cancers = [unique(pairs.C); unk_canc'];
mat_rr = zeros(length(md3), length(cancers));
for m = 1:length(md3)
  for ci = 1:length(cancers)
    mat_rr(m,ci) = altdat.peak_mut.relative_risk(ismember(altdat.peak_mut.MD,md3(m)) & ismember(altdat.peak_mut.C,cancers(ci)));
    if mat_rr(m,ci) > 0
      mat_rr(m,ci) = +(mat_rr(m,ci) > 1);
    end
  end
end

cg = clustergram(mat_rr,'rowlabels',md3,'columnlabels',cancers,'symmetric',false,'standardize',3,'rowpdist','euclidean','linkage','ward','columnpdist','euclidean')

[x, cluster_c] = ismember(cg.ColumnLabels, cancers);
[x, cluster_r] = ismember(cg.RowLabels, md3);
close all hidden;

sel = sel | (ismember(altdat.peak_mut.MD, md3) & ismember(altdat.peak_mut.C, unk_canc));
pairs = [altdat.peak_mut.MD, altdat.peak_mut.C];
metric_values = [altdat.peak_mut.gene_enrichment, altdat.peak_mut.pathway_correlation, altdat.peak_mut.coex_CG, altdat.peak_mut.biogrid_set, altdat.peak_mut.humannet_set];
metric_values = metric_values(sel,:);
pairs = pairs(sel,:);
metric_values(metric_values == -1) = 1;
for i = 1:size(metric_values,2)
  metric_values(:,i) = mafdr(metric_values(:,i), 'BHFDR',true); %manual_predictors(rr==1,i),'BHFDR',true);
end

figure
set(gcf,'Position',[0,0,40,40],'units','centimeters') % 30 60
set(gcf,'Position',[0,0,40,40],'units','centimeters') % 30 60
set(gca,'units','centimeters','Position',[15,2,20,35])

axis equal
hold on
md3_draw = md3(cluster_r);
cancers_draw = cancers(cluster_c);
xlim([0,length(cg.ColumnLabels)])
ylim([0,length(cg.RowLabels)])
xnames = cancers_draw;

gray_purple = [0.8373    0.9000    0.9745]; %0.8039    0.8784    0.9686]; %0.8706    0.9216    0.9804];  %[ 0.7098    0.6392    0.7608] ; % [0.5804    0.5804    0.7333];
light_gray = [0.8000    0.8000    0.8];%0.8784    0.8706    0.8706]; %0.8314    0.8157    0.7843]; %

for i = 1:size(mat_rr,1)
  for j = 1:size(mat_rr,2)
    psel = ismember(altdat.peak_mut.MD,md3_draw(i)) & ismember(altdat.peak_mut.C,cancers_draw(j));
    rr = altdat.peak_mut.relative_risk(psel);
    color = [1,1,1];
    if rr > 1
      color = gray_purple;
      rectangle('Position',[j-1,i-1,1,1],'FaceColor',gray_purple); 
    end
    if rr < 0
      rectangle('Position',[j-1,i-1,1,1],'FaceColor',light_gray); 
    end
    if rr == 0
      rectangle('Position',[j-1,i-1,1,1],'FaceColor',[1 1 1]); 
    end
    ind = find(ismember(pairs(:,1),md3_draw(i)) & ismember(pairs(:,2),cancers_draw(j)));
    if metric_values(ind,3) < .1 % coex
      plot(j-.5,i-.5,'ob','MarkerSize',5)
    end
    if metric_values(ind,2) < .1 %pathway
      plot(j-.5,i-.5,'om','MarkerSize',10);%,'MarkerEdgeColor',[0.8157    0.5843    0.0471])
    end
    if metric_values(ind,4) < .1 % biogrid
      plot(j-.5,i-.5,'oc','MarkerSize',18,'LineWidth',1);%,'MarkerEdgeColor',[0.2000    0.8000         0]) %[1.0000    0.8000    0.2000] = orange
    end
    if metric_values(ind,5) < .1 % biogrid
      plot(j-.5,i-.5,'ok','MarkerSize',20,'LineWidth',1);%,'MarkerEdgeColor',[0.2000    0.8000         0]) %[1.0000    0.8000    0.2000] = orange
    end
    if metric_values(ind, 1) < .1 %GE
      plot(j-.5,i-.5,'or','MarkerSize',25,'LineWidth',1)
    end
  end
end
xticklabel_rotate(.5:1:length(cancers_draw),90,cancers_draw,'FontSize',13)
set(gca,'Ytick',.5:1:(length(md3_draw) - .5),'Yticklabel',md3_draw,'FontSize',13,...
        'TickLength',[ 0 0 ]);
box on            
h = gcf;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(h,'figwork/fig_4b','-dpdf','-r0')
close

%% legend
figure
axis equal
xlim([0,15]); ylim([0,15])
hold on
rectangle('Position',[1,5,1,1],'FaceColor',gray_purple); 
text(2.1,5.5,'comorbid','HorizontalAlignment','left','FontSize',14)
rectangle('Position',[1,3.5,1,1],'FaceColor',[1 1 1]); 
text(2.1,4,'no comorbidity','FontSize',14)
rectangle('Position',[1,2,1,1],'FaceColor',light_gray); 
text(2.1,2.5,'unmeasured','FontSize',14)
%
plot(7.5,6,'ob','MarkerSize',5)
text(8.5,6,'Coexpression metric','FontSize',14)
plot(7.5,4.8,'om','MarkerSize',10)
text(8.5,4.8,'Pathway metric','FontSize',14)
plot(7.5,3.6,'oc','MarkerSize',18,'LineWidth',1);
text(8.5,3.6,'BioGRID metric','FontSize',14)
plot(7.5,2.4,'ok','MarkerSize',20,'LineWidth',1);
text(8.5,2.4,'HumanNet metric','FontSize',14)
plot(7.5,1.2,'or','MarkerSize',25,'LineWidth',1);
text(8.5,1.2,'Gene enrichment','FontSize',14)
pdf('figwork/fig_4b_legend')
