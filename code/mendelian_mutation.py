import scipy
import os
import pdb
import scipy.stats as stats
import pandas as pd
import numpy as np
import math
import csv
from itertools import chain
import pickle
import time

        
def enrich_comorb(alt_scores, covariate, inds, pcut=.05):
        covar = zip(alt_scores['relative_risk'], alt_scores[covariate])
        #inds2 = [i for i in inds if alt_scores[covariate] > -1]
        tab = [sum(x) for x in 
               zip(*tuple([[covar[i][0] > 1, covar[i][1] < pcut, covar[i][0] > 1 and covar[i][1] < pcut] 
                           for i in inds if covar[i][1] > -1]))]
        score = 1
        if len(tab) > 0: 
                score = 1 - stats.hypergeom.cdf(tab[2] - 1,
                                                (pd.DataFrame(alt_scores[covariate]).iloc[inds] > -1).sum()[0],
                                                tab[0], tab[1])
        #print covariate + '\t' + str(score) #+ '\t' + '\t'.join([str(x) for x in tab])
        #pdb.set_trace()
        return score

        
def pathway_enrich(gset, pathways, remove_genes = set([])):
    all_pathway_genes = set(chain.from_iterable([list(xx) for xx in pathways.values()])) - remove_genes
    gset_ixn =( gset & all_pathway_genes ) - remove_genes
    pathway_scores = [0]*len(pathways)
    for (p_i, p) in enumerate(pathways):
        #pathway_scores[p_i] = max(1 - stats.hypergeom.cdf(len(gset_ixn & pathways[p])-1, len(all_pathway_genes), len(gset_ixn), len(pathways[p] - remove_genes)),0)
        pathway_scores[p_i] = 1 - stats.binom.cdf(len(gset_ixn & pathways[p])-1, len(pathways[p] - remove_genes),
                                                  len(gset_ixn)/float(len(all_pathway_genes)))
    return scipy.array(pathway_scores)

def my_bh_fdr(p_val_vec):
    index = scipy.argsort(p_val_vec)
    exp_err = scipy.vstack((float(len(p_val_vec))/scipy.arange(1,len(p_val_vec) + 1)*p_val_vec[index],
                                      scipy.tile(1, [1, len(p_val_vec)]))).min(axis = 0)
    exp_err = scipy.vstack((exp_err,exp_err[scipy.r_[0,scipy.arange(len(exp_err)-1)]])).max(axis=0)
    #scipy.r_[index[0], index[range(len(index)-1)]
    resort_index = scipy.argsort(index)                 
    return exp_err[resort_index]



def write_pathway_enrichments(enrichments, pathways, p_or_q, out_name):
    f = open(out_name + '.' + p_or_q + '.txt','wb')
    writer = csv.writer(f, delimiter = '\t')
    writer.writerow([''] + enrichments.keys())
    x = zip(*tuple([enrichments[t][p_or_q] for t in enrichments]))
    pk = pathways.keys()
    writer.writerows([[pk[i]] + list(x[i]) for i in range(len(x))])
    f.close()
                                   
def alteration_enrichments(alterations, out_name, remove_genes = set([])):
    pathways = open_pathways()              
    enrichments = dict()
    for alt in alterations:
        enrichments[alt] = dict()
        for tumor_dat in alterations[alt]:
            p_value = pathway_enrich(alterations[alt][tumor_dat][1], pathways, set(remove_genes))
            q_value = my_bh_fdr(p_value)
            enrichments[alt][tumor_dat] = {'p':p_value,'q':q_value}
            write_pathway_enrichments(enrichments[alt], pathways, 'q', out_name + '.'  + alt)
        # also put to txt file
    pkl = open('data_processed/alteration_enrichments.pkl','w')
    pickle.dump(enrichments,pkl)
    pkl.close()
    return enrichments

def md_enrichments(icd_gene_clinical, out_name, remove_genes = set([]), save_dir='./'):
    pathways = open_pathways()              
    enrich = dict()
    for (m_i, icd) in enumerate(icd_gene_clinical):
        gene_list = set(icd_gene_clinical[icd]['gene_omim'].keys())
        p_value = pathway_enrich(gene_list, pathways, set(remove_genes))
        q_value = my_bh_fdr(p_value)
        b_value = np.array([1 if len(gene_list & pathways[p])>0 else 0 for p in pathways])
        enrich[icd] = {'p':p_value, 'q':q_value, 'b':b_value}
        #print icd + ' num q < .1 = ' , len([q for q in q_value if q < .1])
    write_pathway_enrichments(enrich, pathways, 'q', out_name + '_MD')
    write_pathway_enrichments(enrich, pathways, 'p', out_name + '_MD')
    pkl = open(save_dir + '/md_enrichments.pkl','w')
    pickle.dump(enrich,pkl)
    pkl.close()
    return enrich



def get_tumor_list(icd_gene_clinical, alterations, cancer_info):
    tumor_data_list = alterations['peak_mut'].keys()
    list2 = set()
    for icd in icd_gene_clinical:
        for cancer in tumor_data_list:
            rr = [icd_gene_clinical[icd]['cancer_assoc'][cancer_icd] 
                  for cancer_icd in cancer_info 
                  if cancer in cancer_info[cancer_icd]['TCGA']]
            if len(rr)>0 and rr[0] > 1:
                list2 |= set([cancer])
    tumor_data_list = list2
    return tumor_data_list
    
def open_network_scores(dat_dir):
    import os
    to_ret = dict()
    net = dat_dir + '/NEIGHBOR_COUNT.biogrid.2000.pkl'
    if os.path.exists(net):
        nc, neighbor_count = pickle.load(open(net))
        to_ret['biogrid'] = neighbor_count
    net = dat_dir + '/NEIGHBOR_COUNT.humannet.9.unwt.1000.pkl'
    if os.path.exists(net):
        nch, neighbor_counth_unwt = pickle.load(open(net))
        to_ret['humannet'] = neighbor_counth_unwt
    # 'hnetUnwt':neighbor_counth_unwt,
    #to_ret = {'biogrid':neighbor_count,  'humannet':neighbor_counth_unwt} 

    return to_ret
    
def open_alteration_enrichments():
    pkl = open('data_processed/alteration_enrichments.pkl')  
    enrichments = pickle.load(pkl);
    pkl.close()
    return enrichments

def open_md_enrichments(dat_dir):
    enrichments = pickle.load(open(dat_dir + '/md_enrichments.pkl'))
    return enrichments

def open_coexpression():
    NUM_SAMPLES = 889;
    cg_corr_file = 'data_processed/fantom/fantom_corrcoef.' + str(NUM_SAMPLES) + '.pkl'
    #mg_corr_file = 'fantom/fantom_corrcoef_MG.' + str(NUM_SAMPLES) + '.pkl'
    version = '.'.join(pd.__version__.split('.')[1:])
    coex_cg = 1
    if float(version) < 11.1:
        #print 'version old'
        coex_cg = pd.DataFrame.load(cg_corr_file)

    else:
        #print 'version new'
        coex_cg = pd.read_pickle(cg_corr_file)        
    #coex_mg = pd.DataFrame.load(mg_corr_file)
    t_inv = stats.t.isf(.05/(coex_cg.shape[0]*coex_cg.shape[1]), NUM_SAMPLES)
    rho_cut = t_inv/((NUM_SAMPLES - 2 + t_inv**2)**.5)

    return coex_cg, NUM_SAMPLES, rho_cut #, coex_mg
    
def load_pathways(pathway_tab, symbol_info, alias_to_symbol):
    pathways = dict()
    ptab = open(pathway_tab)
    ptab.readline()
    for pline in ptab:
        pathway_info = pline.strip().split('\t')
        if not pathway_info[1] in ['PID','PharmGKB']:
            continue
        pgenes = pathway_info[2].split(',')
        for (g_i, g) in enumerate(pgenes):
            if not g in symbol_info and g in alias_to_symbol and len(alias_to_symbol[g]) == 1:
                pgenes[g_i] = alias_to_symbol[g][0]
        pathways[pathway_info[0]] = set(pgenes)
    import pickle
    pkl = open('data_processed/pathways.pkl','w')
    pickle.dump(pathways,pkl)
    pkl.close()
    return pathways

def open_pathways():
    import pickle
    pkl = open('data_processed/pathways.pkl')
    pathways = pickle.load(pkl);
    pkl.close()
    return pathways

############################################
def open_alterations():
    pkl = open('data_processed/alterations.pkl')  
    import pickle
    alterations = pickle.load(pkl);
    pkl.close()
    return alterations

def load_alterations(mut_dir, alias_to_symbol):
    alteration = {'amp_peak':dict(), 'del_peak':dict()} #, 'amp_region':dict(), , 'del_region':dict()}

    altfile = os.listdir(mut_dir)
    cn_file = [ alt_file for alt_file in altfile if 'all_data_by_genes.txt' in alt_file ]
    for cn in cn_file:
        cancer = cn.split(".")[0]
        cn_reads = open(mut_dir + '/' + cn)
        cn_genes = set()
        for gene_cn in cn_reads:
            if not gene_cn.startswith('Gene'):
                cn_genes.add(gene_cn.split('\t')[0])
        for cn_type in ['amp','del']:
            table = open(mut_dir + '/' + cancer + '.table_' + cn_type + '.conf_99.txt')
            #region = set()
            peak = set()
            for sig_amp in table:
                processed = sig_amp.strip().split("\t")
                #region |= set(processed[9].split(","))
                cn_gene_list = set(processed[11].split(","))
                if len(cn_gene_list) < 50:
                    peak |= cn_gene_list
            #alteration[cn_type + '_region' ][cancer] = (cn_genes, region)
            alteration[cn_type + '_peak' ][cancer] = (cn_genes, peak)
    mut_file = [ altf for altf in altfile if altf.endswith('.mut') ]                        
    mutation = dict()
    for mut in mut_file: 
        cancer = mut.split(".")[0]
        mutsig = open(mut_dir + '/' + cancer + '.mut')
        mut_tested = set()
        mut_significant = set()
        for gene_line in mutsig:
            if gene_line.startswith('rank'): continue
            gene_info = gene_line.strip().split("\t")
            gene = gene_info[1]
            mut_tested.add(gene)
            if float(gene_info[19]) < .25:
                mut_significant.add(gene)
        mutation[cancer] = (mut_tested, mut_significant)

    cn_mut = {'peak_mut':dict()} #'amp_mut':dict(),'del_mut':dict(),
    for cancer_name in mutation:
        if cancer_name in alteration['amp_peak']:
            #cn_mut['amp_mut'][cancer_name] = (alteration['amp_peak'][cancer_name][0] | mutation[cancer_name][0],
            #                                  alteration['amp_peak'][cancer_name][1] | mutation[cancer_name][1])
            #cn_mut['del_mut'][cancer_name] = (alteration['del_peak'][cancer_name][0] | mutation[cancer_name][0],
            #                                  alteration['del_peak'][cancer_name][1] | mutation[cancer_name][1])
            cn_mut['peak_mut'][cancer_name] = (alteration['del_peak'][cancer_name][0] | mutation[cancer_name][0],
                                              alteration['del_peak'][cancer_name][1] | alteration['amp_peak'][cancer_name][1] | mutation[cancer_name][1])

    alteration['mutation'] = mutation
    for cn_mut_type in cn_mut:
        alteration[cn_mut_type] = cn_mut[cn_mut_type]

    from itertools import chain
    for alt in alteration:
        tested = set(list(chain.from_iterable([alteration[alt][cx][0] for cx in alteration[alt]])))
        altered = set(list(chain.from_iterable([alteration[alt][cx][1] for cx in alteration[alt]])))
        if len(altered) > 0:
            alteration[alt]['MY_PAN_CAN'] = (tested, altered)

    ## load and add cancer census
    """
    with open(mut_dir + '/census_deletion') as f:
        alteration['del_peak']['CENSUS'] = (alteration['del_peak']['MY_PAN_CAN'][0], set(f.read().splitlines()))
    with open(mut_dir + '/census_amplification') as f:
        alteration['amp_peak']['CENSUS'] = (alteration['amp_peak']['MY_PAN_CAN'][0], set(f.read().splitlines()))
    with open(mut_dir + '/census_mutation') as f:
        alteration['mutation']['CENSUS'] = (alteration['mutation']['MY_PAN_CAN'][0], set(f.read().splitlines()) )
    alteration['peak_mut']['CENSUS'] = (alteration['mutation']['MY_PAN_CAN'][0], 
                                        set(list(chain.from_iterable([alteration[alt]['CENSUS'][1] for alt in ['mutation','amp_peak','del_peak']]))))
    """                                        

    for alt_type in alteration:
        for alt_dat in alteration[alt_type]:
            tested = set()
            for gene in alteration[alt_type][alt_dat][0]:
                if gene in alias_to_symbol:
                    if not gene in alias_to_symbol[gene]:
                        tested |= set(alias_to_symbol[gene])
                    else:
                        tested |= set([gene])
            sig = set()
            for gene in alteration[alt_type][alt_dat][1]:
                if gene in alias_to_symbol:
                    if not gene in alias_to_symbol[gene]:
                        sig |= set(alias_to_symbol[gene])
                    else:
                        sig |= set([gene])
            alteration[alt_type][alt_dat] = (tested, sig)

    import pickle
    pkl = open('data_processed/alterations.pkl','w')
    pickle.dump(alteration,pkl)
    pkl.close()
    return alteration


def ranksummat(mat, msel, notm):
    return stats.ranksums(mat[notm], mat[msel])

def expression_test(coexpression_c, cancer_alt, gene_list, rho_cut=0):
    t = [i for (i,g) in enumerate(coexpression_c.index) if g in gene_list];
    f = [i for (i,g) in enumerate(coexpression_c.index) if g not in (set(gene_list) | cancer_alt)];
    csel = list(set(coexpression_c.columns) & cancer_alt)
    score = -1
    coex_cg_C = ''
    coex_cg_M = ''
    if len(t) > 0 and len(csel) > 0:

        tmp = coexpression_c.loc[:, csel]
        xx = np.apply_along_axis(ranksummat, 0, tmp, t, f)
        bh_cg_score = my_bh_fdr(np.array(xx[1]))
        bh_cg_score[np.array(xx[0]) > 0] = 1
        score = bh_cg_score.min()
        if bh_cg_score.min() < .05 and rho_cut > 0:
            coex_cg_C = \
              ';'.join(list(tmp.columns[np.nonzero(bh_cg_score < .05)[0]]))
            coex_cg_M = \
              ';'.join([','.join(tmp.index[tmp.iloc[t,i] > rho_cut])
                        for i in np.nonzero(bh_cg_score < .05)[0]])
    return score, coex_cg_M, coex_cg_C

