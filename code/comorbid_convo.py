import mendelian_mutation
import mendelian_code
import pickle
import network_ops
import scipy.stats as stats
import numpy as np
import pandas as pd
import pdb

alt = 'peak_mut'
def comorbidity_convo(data_in, out_name, rnd, nets):  #, remove_known = set([])
    #if not os.path.exists(out_name): os.mkdir(out_name)
    #out_name = out_name + "/"

    pathways = mendelian_mutation.open_pathways()              
    pathways_keys = pathways.keys()
    alt_enrichments = mendelian_mutation.open_alteration_enrichments()
    md_pathways = mendelian_mutation.open_md_enrichments(data_in)

    icd_gene_clinical, cancer_info, mendelian_genes, alterations, remove_known = \
      pickle.load(open(data_in + '/dat.pkl'))

    tumor_data_list = mendelian_mutation.get_tumor_list(icd_gene_clinical, alterations, cancer_info)
    print out_name
    out_name = out_name + '/'
    import os
    if not os.path.exists(out_name):
        print 'making' + out_name 
        os.mkdir(out_name)
    #coexpression_c, NUM_SAMPLES, rho_cut = mendelian_mutation.open_coexpression()
    out = open(out_name + 'pairs' + str(rnd) + '.xls','w')    
    sum_overlaps = 0
    rands = np.zeros([0,rnd])

    sum_path_bin = 0
    path_rands_bin = np.zeros([0,rnd])

    sum_coexpressions = 0
    rand_coexpressions = np.zeros([0,rnd])

    mypancan_not_enrich = alt_enrichments[alt]['MY_PAN_CAN']['q'] > .05
    i = 0
    for (m_i, icd) in enumerate(icd_gene_clinical):

        gene_list = set(icd_gene_clinical[icd]['gene_omim'].keys())
        for cancer in tumor_data_list:
            if len(gene_list) == 0:
                continue

            rr = [icd_gene_clinical[icd]['cancer_assoc'][cancer_icd] 
                  for cancer_icd in cancer_info 
                  if cancer in cancer_info[cancer_icd]['TCGA']]
            if len(rr) == 0: continue
            if rr[0] <= 1: continue
            #print icd
            i += 1
            #if cancer in alterations[alt]:
            cancer_alt = alterations[alt][cancer][1]
            cancer_alt_test_overlap = alterations[alt][cancer][1] - remove_known
            overlap = gene_list & cancer_alt_test_overlap
            sum_overlaps += len(overlap)
            
            #alt_scores[alt]['mendelian_gene_overlap'][ind] = len(overlap)/float(len(gene_list))
            genome_tested = alterations[alt][cancer][0]
            rgen = stats.binom.rvs(len(cancer_alt_test_overlap),
                                len(genome_tested & set(gene_list) - remove_known)/float(len(genome_tested)),
                                size = rnd)
            rands = np.vstack((rands, rgen)).sum(axis=0)

            ### binary
            cancpath = mypancan_not_enrich & (alt_enrichments[alt][cancer]['q'] < .1) 
            bin_md_path = [mypancan_not_enrich[i] & (len(gene_list & pathways[p]) > 0)
                        for (i,p) in enumerate(pathways)]
            if not icd in ['poop']: #'"Pervasive, Specified Congenital Anomalies"']: #, 'Combined Heart and Skeletal Defects']:
                #pdb.set_trace()
                sum_path_bin += sum(bin_md_path & cancpath)
            if sum(bin_md_path) > 0 and sum(cancpath) > 0 and not icd in ['poop']: #'"Pervasive, Specified Congenital Anomalies"']: #, 'Combined Heart and Skeletal Defects']:
                rgen = stats.hypergeom.rvs(sum(mypancan_not_enrich), sum(cancpath), sum(bin_md_path), size = rnd)
                path_rands_bin = np.vstack((path_rands_bin, rgen)).sum(axis = 0)
            p_pathway_string = ''
            psel = list(np.nonzero(cancpath & bin_md_path)[0])
            if len(psel) > 0:
                for p_both in psel:
                    p_pathway_string += ';' + ','.join(pathways[pathways_keys[p_both]] & gene_list) + ' -> ' + \
                      pathways_keys[p_both].replace(';',' ') + '(' + ','.join(pathways[pathways_keys[p_both]] & cancer_alt) + ')'

            '''
            csel = set(coexpression_c.columns) & cancer_alt
            msel = (set(coexpression_c.index) & set(gene_list)) - set(csel)
            if len(csel) > 0 and len(msel) > 0:
                coex_mat = coexpression_c.loc[msel, csel]
                notm = set(coexpression_c.index) - msel - csel
                background_prob = (coexpression_c.loc[notm, csel] > rho_cut).sum(axis = 0)/float(len(notm))
                successes = (coex_mat > rho_cut).sum(axis = 0)
                cg_score = 1 - stats.binom.cdf(successes - 1, len(msel), background_prob)
                #if sum(successes) > 0: 
                p_use = background_prob[cg_score == min(cg_score)]
                rgen = np.zeros((len(p_use),rnd))
                for p in range(len(p_use)):
                    rgen[p,:] = stats.binom.rvs(len(msel), p_use[p], size = rnd)
                rand_coexpressions = np.vstack((rand_coexpressions, rgen.mean(axis=0))).sum(axis=0)
                sum_coexpressions += (successes[cg_score == cg_score.min()]).mean()
            '''
            towrite = [icd, cancer, ','.join(overlap), p_pathway_string] #q_pathway_string, 
            out.write('\t'.join([str(x) for x in towrite]) + '\n')

    out.close()
    scores = pd.DataFrame(np.zeros([4,2]), index =['Gene_Overlap','Shared_pathways'] +nets, columns=['count','p'])
    #print('GE: tot=' + str(sum_overlaps) + '\tp=' +str((rands >= sum_overlaps).sum()/float(rnd)))
    #print('Path_b tot=' + str(sum_path_bin) + '\tp=' + str((path_rands_bin >= sum_path_bin).sum()/float(rnd) ))
    scores.loc['Gene_Overlap',:] = [sum_overlaps, (rands >= sum_overlaps).sum()/float(rnd)]
    scores.loc['Shared_pathways',:] = [sum_path_bin, (path_rands_bin >= sum_path_bin).sum()/float(rnd)]
    #print('Coex tot=' + str(sum_coexpressions) + '\tp=' + str((rand_coexpressions >= sum_coexpressions).sum()/float(rnd) ))
    #+ '\nCoex_p=' + str((rand_coexpressions >= sum_coexpressions).sum()/float(rnd))+ '\n')
    #pkl = open(out_name + 'rand' + str(rnd) + '.pkl','w')
    #pickle.dump((sum_overlaps, rands, sum_path_p, path_rands_p, sum_coexpressions, rand_coexpressions, sum_path_bin, path_rands_bin), pkl)
    #pkl.close()
    for net in nets:
        netdir = 'data_processed/' + net
        [ct, randcts] = network_ops.comorbid_count_compare(netdir, icd_gene_clinical, cancer_info, alterations)
        x = open(netdir +'/randcts','w')
        x.write('\n'.join([str(i) for i in randcts]) + '\n')
        x.close()
        scores.loc[net,:] = [ct, sum(randcts >= ct)/float(len(randcts))]

    np.savetxt(out_name + 'rand_ge' + str(rnd)  + '.txt', rands, delimiter ='\t')
    np.savetxt(out_name +  'rand_path_b.' + str(rnd) +'.txt', path_rands_bin, delimiter ='\t')
    return scores
    #np.savetxt(out_name + '.' + str(rnd) + '.rand_coex.txt', rand_coexpressions, delimiter ='\t')

