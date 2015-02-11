
import mendelian_code
import comorbidity_genetics
import pickle
import gene_info
import mendelian_mutation
import comorbid_convo
import os
from itertools import chain
'''
import sys
sys.path.append('code')
'''
ncbi_genes= "data_source/Homo_sapiens.gene_info.protein_coding"


def load_cancer_alteration_info():
    symbol_info, alias_to_symbol = gene_info.process_ncbi(ncbi_genes)
    alterations = mendelian_mutation.load_alterations('data_source/cancer_alterations', alias_to_symbol)
    pway = mendelian_mutation.load_pathways('data_source/CPDB_pathways_genessymbol.tab', symbol_info, alias_to_symbol);
    alt_en = mendelian_mutation.alteration_enrichments(alterations, 'data_processed/cancer_pathways', remove_genes = set([]))
    return alterations

def comorbidity_pairs(dat_dir):
    tab = comorbidity_genetics.comorbidity_scores(dat_dir,'comorbidity_analysis_' + dat_dir)

def comorbidity_aggregation(dat_dir, outname, rnd, nets):
    return comorbid_convo.comorbidity_convo(dat_dir, outname, rnd, nets)
    
def load_mendelian_disease_info(dat_dir, alterations):
    # set up the data, link in the comorbidities
    # save data with germline and alterations
    icd_gene_clinical, cancer_info, mendelian_genes = \
      mendelian_code.load_associations_annotations(conservative=True)
    germline_dir = dat_dir + '_germline'
    if not os.path.exists(germline_dir): os.mkdir(germline_dir)
    f = open(germline_dir + '/dat.pkl','w')
    pickle.dump((icd_gene_clinical, cancer_info, mendelian_genes, alterations, set([])), f)
    f.close()
    md_en = mendelian_mutation.md_enrichments(icd_gene_clinical, germline_dir + '/pathway', set([]), germline_dir)
    
    # set up the germline census genes to remove from the analysis
    f = open('data_source/germline_to_remove')
    cens_germline = f.read().split("\n")

    # then remove germline and save the no-germline version to another directory
    icd_gene_clinical, cancer_info, mendelian_genes, removed_genes = \
      mendelian_code.open_associations_nogermline(germline_dir + '/dat.pkl', cens_germline)
    if not os.path.exists(dat_dir): os.mkdir(dat_dir)
    pkl = open(dat_dir + '/dat.pkl','w')
    pickle.dump((icd_gene_clinical, cancer_info, mendelian_genes, alterations, removed_genes),pkl)
    pkl.close()
    md_en = mendelian_mutation.md_enrichments(icd_gene_clinical, dat_dir + '/pathway', removed_genes, dat_dir) 

def setup_networks(biogrid_source_file, humannet_source_file, mendelian_dir, mendelian_germ_dir, nrand=1000):
    ncbi_genes= "data_source/Homo_sapiens.gene_info.protein_coding"
    # hprd_source_file = 'data_source/HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'
    # hprd_proc_file = data_processed/hprd/"
    import subprocess
    import os
    os.mkdir(hprd_dir)
    import humannet
    humannet_dir = 'data_processed/humannet.9.unwt'
    humannet.humannet(humannet_source_file, ncbi_genes, humannet_dir)
    
    import biogrid
    biogrid_dir = 'data_processed/biogrid'
    biogrid.biogrid(biogrid_source_file, ncbi_genes, biogrid_dir)
    prepare_randomization(biogrid_dir, mendelian_dir, nrand=nrand)
    prepare_randomization(humannet_dir, mendelian_dir, nrand=nrand)
    prepare_randomization(biogrid_dir, mendelian_germ_dir, nrand=nrand)
    prepare_randomization(humannet_dir, mendelian_germ_dir, nrand=nrand)


def fantom():
    import fantom
    symbol_info, alias_to_symbol = gene_info.process_ncbi(ncbi_genes)
    alterations = mendelian_mutation.open_alterations()
    alt = 'peak_mut'
    cancer_gn = set(list(chain.from_iterable([alterations[alt][x][1] for x in alterations[alt]])))

    to_cor = fantom.fantom_tpm_to_expr('data_source/hg19.cage_peak_tpm_ann.osc.txt', symbol_info, alias_to_symbol, cancer_gn)
    fantom.do_cor(to_cor)

def figure_files(dat_dir):
    icd_gene_clinical, cancer_info, mendelian_genes, alterations, remove_known = \
      pickle.load(open(dat_dir + '/dat.pkl'))

    #
    if not os.path.exists('figwork'): os.mkdir('figwork')
    import pandas as pd
    import numpy as np

    # Fig 2b
    alt_enrichments = mendelian_mutation.open_alteration_enrichments()
    md_enrichments = mendelian_mutation.open_md_enrichments(dat_dir)
    pd.DataFrame({'RTb':md_enrichments['Combined Heart and Skeletal Defects']['b'],
                  'SKp':alt_enrichments['peak_mut']['SKCM']['p'],
                  'SKq':alt_enrichments['peak_mut']['SKCM']['q'],
                  'pancanq':alt_enrichments['peak_mut']['MY_PAN_CAN']['q']}).to_csv('figwork/skcm_rt.xls',sep='\t')

    # Fig 2d
    coexpression_c, NUM_SAMPLES, rho_cut = mendelian_mutation.open_coexpression()
    coexpression_c.loc[:,'PTK6'].to_csv('figwork/PTK6.txt',sep='\t')

    # Fig 3c
    mut_dir = 'data_source/cancer_alterations'
    cn = 'GBM.all_thresholded.by_genes.txt'
    glist = ['RPL5','RPS7','RPL11','MDM2']
    cn_thres = open(mut_dir + '/' + cn)
    header = cn_thres.readline()

    arr = pd.DataFrame(np.zeros([len(glist),len(header.strip().split("\t")) - 3]), index=glist)

    for gene_cn in cn_thres:
        gene_line = gene_cn.strip().split()
        gene = gene_line[0]
        if  gene in glist:
            arr.loc[gene,:] = [int(i) for i in gene_line[3:]]
    arr.to_csv('figwork/rpl_mdm.txt',sep='\t',header=False)

    ## S1a
    ct = open('figwork/numg_per_icd.txt','w')
    for icd in icd_gene_clinical:
        ct.write( icd + '\t' + str(len(icd_gene_clinical[icd]['gene_omim'].keys())) + '\n')
    ct.close()
    
    ## S1b
    ct = open('figwork/numg_per_c_peakmut.txt','w')
    alt = 'peak_mut'
    for c in alterations[alt]:
        ct.write( c + '\t'  + str(len(alterations[alt][c][1])) + '\n')
    ct.close()
    
