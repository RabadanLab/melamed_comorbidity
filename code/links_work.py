import mendelian_code
import comorbidity_genetics
import pickle
import gene_info
import mendelian_mutation
ncbi_genes= "data_source/Homo_sapiens.gene_info.protein_coding"
symbol_info, alias_to_symbol = gene_info.process_ncbi(ncbi_genes)

dat_dir = 'md_data_jul15'
#save_file = dat_dir + '_germline' + '/dat.pkl', 
icd_gene_clinical, cancer_info, mendelian_genes = mendelian_code.load_associations_annotations(conservative=True)
alterations = mendelian_mutation.load_alterations('data_source/cancer_alterations', alias_to_symbol)
f = open(dat_dir + '_germline' + '/dat.pkl','w')
pickle.dump((icd_gene_clinical, cancer_info, mendelian_genes, alterations, set([])), f)
f.close()
f = open('data_source/census/germline_to_remove')
cens_germline = f.read().split("\n")

icd_gene_clinical, cancer_info, mendelian_genes = mendelian_code.open_associations_nogermline(dat_dir + '/dat_w_germline.pkl', cens_germline)

pkl = open(dat_dir + '/dat.pkl','w')
pickle.dump((icd_gene_clinical, cancer_info, mendelian_genes, alterations, set(cens_germline)),pkl)
pkl.close()


icd_gene_clinical, cancer_info, mendelian_genes, alterations, remove_known =  pickle.load(open(dat_dir + '/dat.pkl'))

pway = mendelian_mutation.load_pathways('data_source/CPDB_pathways_genessymbol.tab', symbol_info, alias_to_symbol);
alt_en = mendelian_mutation.alteration_enrichments(alterations, 'june30_reduc', remove_genes = set([]))
md_en = mendelian_mutation.md_enrichments(icd_gene_clinical, 'jul15MD', remove_known, dat_dir)

mendelian_mutation.write_pathway_enrichments(md_en, pway, 'p', 'md_pway_jul15')

qsub -N network -V -cwd -l mem=3G,time=6:: -b y "python /ifs/home/c2b2/rr_lab/rdm2114/mendelian/code/run_net.py NEIGHBOR_COUNT 1000 hprd dat.pkl"

g = 'FAM58A'
coex_mat.loc[:,g].to_csv(g + 'ecto.txt',sep='\t')
coexpression_c.loc[notm,g].to_csv(g + 'notm.txt',sep='\t')

alt_scores = comorbidity_genetics.comorbidity_scores(dat_dir,'comorb_jul15')

counts = zip(*tuple([[path[0] < .05, path[1] < .05, path[0] < .05 and path[1] < .05, path[0] < .05 or path[1] < .05] for path in zip(p_value, alt_en['peak_mut']['SKCM']['p'])]))

alt_scores = comorbidity_genetics.comorbidity_scores('md_data_june30_conservative_coding','comorb_june30')

icd_gene_clinical, cancer_info, mendelian_genes = mendelian_code.open_associations_annotations()
icd_gene_clinical, cancer_info, mendelian_genes = mendelian_code.open_associations_nogermline()

alt_enrich = mendelian_mutation.alteration_enrichments(alterations, 'focal_only')
pkl = open('mendelian_disease_alterations_nogermline.pkl','w')
pickle.dump((icd_gene_clinical, cancer_info, mendelian_genes, alterations),pkl)
pkl.close()

icd_gene_clinical, cancer_info, mendelian_genes, alterations = pickle.load(open('mendelian_disease_alterations_nogermline'))
#qsub -N nc_hprd -V -cwd -l mem=3G,time=2:: -b y "python /ifs/home/c2b2/rr_lab/rdm2114/mendelian/code/run_net.py NEIGHBOR_COUNT 1000 hprd mendelian_disease_alterations_nogermline.pkl"    

f = open('peak_mut_p.txt','wb')
writer = csv.writer(f, delimiter = '\t', quoting=csv.QUOTE_NONE, quotechar="|")
#writer.writerow(['canc'] + pathways.keys()) 
#writer.writerow([tumor_dat] + enrichments[alt][tumor_dat].tolist())
writer.writerow(['p'] + enrichments[alt].keys())
x = zip(*tuple([enrichments[alt][t]['p'] for t in enrichments[alt]]))
pk = pathways.keys()
writer.writerows([[pk[i]] + list(x[i]) for i in range(len(x))])
f.close()

import scipy
import networkx
from itertools import chain
import scipy.stats as stats
import pandas as pd
import numpy as np
import mendelian_mutation


alt_scores = comorbidity_genetics.comorbidity_scores('mendelian_disease_alterations_nogermline.pkl','comorb_focal', set(cens_germline))

cg_score = (1 - stats.binom.cdf((coex_mat > rho_cut).sum(axis = 0) - 1, len(msel), (coexpression_c.loc[notm, csel] > rho_cut).sum(axis = 0)/float(len(notm))))
                    
sc = pd.DataFrame(np.array(cg_score),index = coex_mat.columns,columns=['p'])
v = sc.sort(['p'])

v = pd.DataFrame(scipy.vstack(((gene_scores < .05*float(nrand)).sum(axis=0),gscore_vs_background,gene_scores.min(axis=0), background_probability)),columns = rand_score[alt]['cancers'])
                                                       1 - stats.binom.cdf((gene_scores < .05*float(nrand)).sum(axis=0),gene_scores.shape[0],
                                                   background_probability),

gs = pd.DataFrame(gene_scores, columns = rand_score[alt]['cancers'], index = mend_gn)
                                background_probability)),
gscore_vs_background = [0]*len(rand_score[alt]['cancers'])
            for (c_i, c) in enumerate(rand_score[alt]['cancers']): 
                background_genes = background_set - alterations[alt][c][1]
                p = (background[alt].loc[background_genes, c] < .05).sum()/float(len(background_genes))
                background_probability[c_i] =p

net_scores = mendelian_mutation.open_network_scores()
alt_enrichments = mendelian_mutation.open_alteration_enrichments()

xx = zip(*tuple([stats.ranksums(coexpression_c.loc[notm, csel_g], coex_mat.loc[:,csel_g]) for csel_g in csel]))
bh_cg_score = my_bh_fdr(np.array(xx[1]))
bh_cg_score[np.array(xx[0]) > 0] = 1
ix = pd.DataFrame(np.zeros([len(bh_cg_score),1]), index=csel, columns=['p'])
ix.iloc[:,0] = bh_cg_score
srt = ix.sort(['p'])
srt[:5]




for alt in alt_scores:
    expr_guys = []
    for expr_guy in expr_stats:
        expr_uncom = [alt_scores[alt][expr_guy][i] for i in inds 
                      if alt_scores[alt]['relative_risk'][i] <= 1]
        expr_com = [alt_scores[alt][expr_guy][i] for i in inds 
                    if alt_scores[alt]['relative_risk'][i] > 1]
        tstat, ttest_Score = stats.ttest_ind(expr_uncom, expr_com, equal_var = False)
        tstat, rs_score = stats.ranksums(expr_uncom, expr_com)
        expr_guys.append(ttest_Score)
        print alt + '\t' + expr_guy + '\t' + '%1.3f' % scipy.array(expr_com).mean() + '(' + '%1.3f' % scipy.array(expr_com).std() + ')'+ '\t' +   '%1.3f' % scipy.array(expr_uncom).mean() + '(' + '%1.3f' % scipy.array(expr_uncom).std() + ')'+ '\t' + str(ttest_Score) + "\t" + str(rs_score)



network_findings = mendelian_mutation.open_network_scores()

get_mendelian_diseases(icd_gene_clinical)    

x =  scipy.nonzero((score['amp_peak']['genewise_disease_score'] < .05).sum( axis = 1)  > 1)[0]
 mdx = [score['amp_peak']['MD'][i] for i in x]

mut_dir = 'data_source/cancer_alterations';
altfile = os.listdir(mut_dir)
cn_file = [ alt_file for alt_file in altfile if 'all_thresholded.by_genes.txt' in alt_file ]
for cn in cn_file:
    cancer = cn.split(".")[0]
    if len([ci for ci in cancer_info if cancer in cancer_info[ci]['TCGA']]) == 0:
        continue
    cn_thres = open(mut_dir + '/' + cn)
    genes = []
    amped = []
    deled = []
    for gene_cn in cn_thres:
        if gene_cn.startswith('Gene'):
            continue
        gene_line = gene_cn.strip().split()
        gene = gene_line[0]
        if gene in alias_to_symbol and \
          not gene in alias_to_symbol[gene] and len(alias_to_symbol[gene]) == 1:
            gene = alias_to_symbol[gene][0]
        if gene in symbol_info and not gene in genes:
            genes.append(gene)
            amped.append(str(len([x for x in gene_line[3:] if int(x) > 0])))
            deled.append(str(len([x for x in gene_line[3:] if int(x) < 0])))
    muted = ['0']*len(genes)
    muts = open(mut_dir + '/' + cancer + '.mut')
    for gene_mut in muts:
        gene_line = gene_mut.strip().split('\t')
        gene = gene_line[1]
        if not gene in genes:
            if gene in alias_to_symbol and len(alias_to_symbol[gene]) == 1:
                gene = alias_to_symbol[gene]
        if gene in genes:
            muted[genes.index(gene)] = gene_line[5]
    numalt = open(cancer + '.gene_alt.txt','w')
    numalt.write('genes\tamped\tdeled\tmuted\n')
    for gene in zip(*tuple([genes, amped, deled, muted])):
        numalt.write('\t'.join(gene) + '\n')
    numalt.close()

mut_dir = 'data_source/cancer_alterations'
cn = 'GBM.all_thresholded.by_genes.txt'
glist = ['RPL5','RPS7','RPL11','MDM2','TP53']
cn_thres = open(mut_dir + '/' + cn)
header = cn_thres.readline()

arr = pd.DataFrame(np.zeros([len(glist),len(header.strip().split("\t")) - 3]), index=glist)

for gene_cn in cn_thres:
    gene_line = gene_cn.strip().split()
    gene = gene_line[0]
    if  gene in glist:
        print gene
        arr.loc[gene,:] = [int(i) for i in gene_line[3:]]
arr.to_csv('rpl_mdm.txt',sep='\t',header=False)
    
ct = open('lesion_ct.xls','w')
for a in alteration:
    for c in alteration[a]:
        ct.write( c + '\t' + a + '\t' + str(len(alteration[a][c][1])) + '\n')
ct.close()



    scipy.savetxt('randres.txt', rand_res, header='\t'.join([code_genes[code_gene_ids.index(score_g_id[g])] for g in range(len(score_g_id))]),delimiter='\t')
    scipy.savetxt('randres.txt', scipy.vstack((code_scores.transpose(), rand_res)), header='\t'.join([code_genes[code_gene_ids.index(score_g_id[g])] for g in range(len(score_g_id))]),delimiter='\t')

    hist = scipy.histogram(code_scores, bins=[0,5,10,15,20,25,30])
    h, c = scipy.histogram(code_scores)
    g2 = humannet.copy()

    
import scipy.cluster.hierarchy as sch
#copy_number, mutation, cn_mut, all_cancer_types = mendelian_mutation.load_mutations(cancer_info)
## now have:
# icd_gene_clinical:  gene_omim, omim_clinical, cancer_assoc
# gene_cancer_assoc:  orpha, omim 2 structs of cancer assoc with scores for each 
# cancer_info: which also has its associated mendelian diseases...
disease_assoc = mendelian_mutation.gene_assoc_scores(gene_cancer_assoc,  cancer_info, mendelian_genes, )
    import os
mendelian_mutation.gene_assoc_scores(gene_cancer_assoc, cancer_info, mendelian_genes, 'per_set')


######### plots
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

Y = sch.linkage(disease_assoc, method='centroid')
Z = sch.dendrogram(Y, orientation='right')
index = Z['leaves']    

Y = sch.linkage(disease_assoc.T, method='centroid')
Z = sch.dendrogram(Y, orientation='right')
index2 = Z['leaves']

D1 = squareform(pdist(disease_assoc, metric='euclidean'))
D2 = squareform(pdist(disease_assoc.T, metric='euclidean'))

    Y = sch.linkage(disease_assoc, method='centroid')
    Z = sch.dendrogram(Y, orientation='right')
    index = Z['leaves']    

f = figure(figsize=(8, 8))



