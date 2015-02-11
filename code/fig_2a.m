function fig_2a(altdat, sel,comorbid_convo_out)

%%% HeatMap fig for comorb/shared genetic

md = unique(altdat.peak_mut.MD(sel));
c = unique(altdat.peak_mut.C(sel));
mat_rr = zeros(length(md), length(c));
mat_ge = zeros(length(md), length(c));
mat2 = mat_rr;
for m = 1:length(md)
  for ci = 1:length(c)
    rr = altdat.peak_mut.relative_risk(ismember(altdat.peak_mut.MD,md(m)) & ismember(altdat.peak_mut.C,c(ci)));
    ge = altdat.peak_mut.gene_enrichment(ismember(altdat.peak_mut.MD,md(m)) & ismember(altdat.peak_mut.C,c(ci)));
    %ge = altdat.peak_mut.NC_best(ismember(altdat.peak_mut.MD,md(m)) & ismember(altdat.peak_mut.C,c(ci)));
    mat_rr(m,ci) = rr;
    mat_ge(m,ci) = +(ge < .05);
    mat2(m,ci) = +(rr > 1);
    if rr > 1 & ge < .05 mat2(m,ci) = 2; 
    elseif rr <= 1 & ge < .05 mat2(m,ci) = 3; end
  end
end
to_rm = sum(mat_rr,2)==0;
mat_rr = mat_rr(~to_rm,:);
md = md(~to_rm);
mat_ge = mat_ge(~to_rm,:);
%cg = clustergram(mat,'rowlabels', md,'columnlabels',c,'symmetric',false,'standardize',3,'Colormap',genColorMap('rw',64),'DisplayRange',8, 'rowpdist','correlation');
cg = clustergram(mat_rr,'rowlabels',md,'columnlabels',c,'symmetric',false,'standardize',3,'Colormap',genColorMap('rw',64),'DisplayRange',8,'rowpdist','euclidean','linkage','ward','columnpdist','euclidean')

[x, cluster_c] = ismember(cg.ColumnLabels, c);
[x, cluster_md] = ismember(cg.RowLabels, md);
close all hidden

%% x & box version
figure
set(gcf,'Position',[0,0,30,45],'units','centimeters')
set(gcf,'Position',[0,0,30,45],'units','centimeters')
set(gca,'units','centimeters','Position',[10,5,25,30])
axis equal
hold on
%ismember(md_names(i),'Glucose-6-Phosphate Dehydrogenase Deficiency') & ismember(cancer_names(j),'LUAD')
[x, cluster_md] = sort(sum(mat_rr > 1,2));
rr_to_draw = mat_rr(cluster_md, cluster_c);
ge_to_draw = mat_ge(cluster_md, cluster_c);
xlim([0,length(cg.ColumnLabels)])
ylim([0,length(cg.RowLabels)])
cancer_names = cg.ColumnLabels;
md_names = md(cluster_md); %cg.RowLabels;
%comorbid_convo_out = 'jan28out100000.xls';
comorbid_info = parseText(comorbid_convo_out);
comorbid_info = comorbid_info.text;
for i = 1:size(rr_to_draw,1)
  for j = 1:size(rr_to_draw,2)
    if rr_to_draw(i,j) > 1
      %rectangle('Position',[j-1,i-1,1,1],'FaceColor',[ 0.0392    0.1412    0.4157]); %,'EdgeColor',[ 0.0392    0.1412    0.4157])
      col = [0.6549    0.5765    0.7294];
      col = [0.5608    0.5608    0.8157];
      col = [0.5804    0.5804    0.7333];
      col = [0.8373    0.9000    0.9745]; % same color as figure 4B
      rectangle('Position',[j-1,i-1,1,1],'FaceColor',col) %,'EdgeColor',col)
    else
      rectangle('Position',[j-1,i-1,1,1],'FaceColor',[1 1 1]) %,'EdgeColor',col)
    end
    comorbid_ind = ismember(comorbid_info(:,1), md_names{i}) & ismember(comorbid_info(:,2), cancer_names{j});
    if sum(comorbid_ind) > 0 %ge_to_draw(i,j) 
      geneo = comorbid_info{comorbid_ind,3};
      color = [1.0000    0.8000         0];
      edgecolor = 'k';
      if length(strsplit(geneo,',')) == 2 color = 'r'; edgecolor = 'r'; end
      if length(geneo) > 0
        plot(j-.5,i-.5,'o','MarkerSize',6,'MarkerFaceColor',color,'MarkerEdgeColor',edgecolor) %,'LineWidth',1.5)
      end
      %plot([j-.85,j-.15],[i-.85,i-.15],'Color',[0 1 0],'LineWidth',2)
      %plot([j-.85,j-.15],[i-.15,i-.85],'Color',[0 1 0],'LineWidth',2)
      %%plot([j-1,j],[i-1,i],'Color',[0 1 0],'LineWidth',2)
      %%plot([j-1,j],[i,i-1],'Color',[0 1 0],'LineWidth',2)
    end
  end
end

xticklabel_rotate(.5:1:length(cancer_names),90,cancer_names,'FontSize',12)
set(gca,'Ytick',.5:1:(length(md_names) - .5),'Yticklabel',md_names,'FontSize',12,...
        'TickLength',[ 0 0 ]);
box on            
h = gcf;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(h,'figwork/2A_blue_red','-dpdf','-r0')
print(h,'2A_blue_red','-dpdf','-r0')


figure
axis equal
xlim([0,15]); ylim([0,15])
hold on
rectangle('Position',[1,5,1,1],'FaceColor',[0.8373    0.9000    0.9745]); % same color as figure 4Bgray_purple); 
text(2.1,5.5,'comorbid','HorizontalAlignment','left','FontSize',16)
rectangle('Position',[1,3.5,1,1],'FaceColor',[1 1 1]); 
text(2.1,4,{'comorbidity'; 'undetected'},'FontSize',16)
%
color = [1.0000    0.8000         0];
plot(6.4,5.5,'o','MarkerSize',18,'MarkerFaceColor',color,'MarkerEdgeColor','k')
text(6.9,5.5,'One gene shared','FontSize',16)
color = 'r'
plot(6.4,4,'o','MarkerSize',18,'MarkerFaceColor',color,'MarkerEdgeColor',color);
text(6.9,4,'Two genes shared','FontSize',16)
h = gcf;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
%print(h,'figwork/2A_blue_red','-dpdf','-r0')
print(h,'2A_legend','-dpdf','-r0')

