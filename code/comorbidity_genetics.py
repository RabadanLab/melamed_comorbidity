from mendelian_mutation import *
alterations_to_analyze = ['peak_mut']

def comorbidity_scores(data_in, out_name, report = True):  #, remove_known = set([])
    pathways = open_pathways()              
    pathways_keys = pathways.keys()
    if not os.path.exists(out_name): os.mkdir(out_name)
    out_name = out_name + "/"
    
    icd_gene_clinical, cancer_info, mendelian_genes, alterations, remove_known = \
      pickle.load(open(data_in + '/dat.pkl'))

    tumor_data_list = alterations[alterations_to_analyze[0]].keys()
    #for alt_type in alterations:
    #    tumor_data_list = tumor_data_list.union(set(alterations[alt_type].keys()))
    #tumor_data_list = list(tumor_data_list)
    alt_enrichments = open_alteration_enrichments()
    md_enrichments = open_md_enrichments(data_in)

    coexpression_c, NUM_SAMPLES, rho_cut = open_coexpression()
    # pearson is distr by t-distr w/ num-samples df --> choose p cutoff and then get 
    #  t-inverse -> transform to get rho cutoff
    t_inv = stats.t.isf(.05/(coexpression_c.shape[0]*coexpression_c.shape[1]), NUM_SAMPLES)
    rho_cut = t_inv/((NUM_SAMPLES - 2 + t_inv**2)**.5)
    #pcut = .05/ (coexpression['p'].shape[0] * coexpression['p'].shape[1])
    #coexpression_p = coexpression_p < 
    #pdb.set_trace()
    to_ret = 0
    #(neighbor_count, shortest_path, sp_seed) = open_network_scores()
    import pdb
    #pdb.set_trace()
    network_findings = open_network_scores(data_in)
    alt_scores = dict()
    num_comorbidities = len(icd_gene_clinical)*len(tumor_data_list)
    #expr_stats = ['max_coexpression','coex_CG', 'coex_fisher'] #,'coex_MG_CG'] #, 'mean_max_coexpression','coex_MG']
    expr_stats = ['coex_CG'] #,'coex_binom'] #,'coex_MG_CG'] #, 'mean_max_coexpression','coex_MG']
    alt = 'peak_mut'

    alt_scores = {'MD':list(chain.from_iterable([[icd]*len(tumor_data_list) 
                                                      for icd in icd_gene_clinical])),
                       'C':tumor_data_list*len(icd_gene_clinical),
                       'relative_risk':[-1]*num_comorbidities,
                       'gene_enrichment':[-1]*num_comorbidities,
                       'gene_intersection':['NA']*num_comorbidities,
                       'pathway_overlap':[-1]*num_comorbidities,
                       'pathway_correlation':[-1]*num_comorbidities,
                       'pathway_correlation_rho':[-1]*num_comorbidities,
                       'pathway_genes':['NA']*num_comorbidities,
                       'pathway_sel':['NA']*num_comorbidities,
                       'pathway_cancer_gene':['NA']*num_comorbidities,
                       'coexpr_cg_M':['NA']*num_comorbidities,
                       'coexpr_cg_C':['NA']*num_comorbidities,
                       'coexpr_p':['NA']*num_comorbidities,
                       'coex_z':[-1]*num_comorbidities}
                       #'coexpr_mg_M':['NA']*num_comorbidities,
                       #'coexpr_mg_C':['NA']*num_comorbidities,
                       #'coexpr_max_M':['NA']*num_comorbidities,
                       #'coexpr_max_C':['NA']*num_comorbidities,

    #'mendelian_gene_overlap':[-1]*num_comorbidities,
    for score in network_findings:
        print score + '_set'
        alt_scores[score + '_set'] = [-1]*num_comorbidities
        alt_scores[score + '_genes'] = ['NA']*num_comorbidities
        alt_scores[score + '_ginfo'] = ['NA']*num_comorbidities
    for score in expr_stats:
        alt_scores[score] = [-1]*num_comorbidities

    header = ['MD','C',
              'gene_intersection','pathway_genes','pathway_sel','pathway_cancer_gene',
              'coexpr_cg_M','coexpr_cg_C'] + \
              list(chain.from_iterable([[score + '_genes', score  + '_ginfo'] for score in network_findings])) #'coexpr_max_M','coexpr_max_C','coexpr_mg_M','coexpr_mg_C']
    to_crosstab = ['gene_enrichment', 'pathway_correlation','pathway_overlap'] + expr_stats + [score + '_set' for score in network_findings] #,'coex_CG'] #,'coex_fisher'] #,'coex_MG_CG'] #,'coex_MG','mean_max_coexpression']        

    header = header + ['relative_risk', 'pathway_correlation_rho'] + to_crosstab
                                                    
    stat_tab_h = open(out_name + 'stat_tab.xls','w')
    stat_tab_h.write('\t'.join(['alt','pathway_correlation'] + to_crosstab + expr_stats) + '\n')
    

    ind = -1
    mypancan_not_enrich = alt_enrichments[alt]['MY_PAN_CAN']['q'] > .05
    pathways_keys_nopancan = [pathways_keys[i] for (i,el) in enumerate(mypancan_not_enrich) if el==True]
    for (m_i, icd) in enumerate(icd_gene_clinical):

        gene_list = set(icd_gene_clinical[icd]['gene_omim'].keys())
        for cancer in tumor_data_list:
            ind += 1
            if len(gene_list) == 0:
                continue

            rr = [icd_gene_clinical[icd]['cancer_assoc'][cancer_icd] 
                  for cancer_icd in cancer_info 
                  if cancer in cancer_info[cancer_icd]['TCGA']]
            if len(rr) == 1:
                alt_scores['relative_risk'][ind] = rr[0]
            #if cancer in alterations[alt]:
            cancer_alt = alterations[alt][cancer][1]
            cancer_alt_test_overlap = alterations[alt][cancer][1] - remove_known
            overlap = gene_list & cancer_alt_test_overlap
            #alt_scores[alt]['mendelian_gene_overlap'][ind] = len(overlap)/float(len(gene_list))
            genome_tested = alterations[alt][cancer][0]
            binom = 1 - stats.binom.cdf(len(overlap) - 1, len(cancer_alt_test_overlap), 
                                        len(genome_tested & set(gene_list) - remove_known)/float(len(genome_tested)))
            alt_scores['gene_enrichment'][ind] = binom
            alt_scores['gene_intersection'][ind] = ','.join(overlap)

            #### pathways
            pathway_cor, corp = stats.spearmanr(md_enrichments[icd]['b'][mypancan_not_enrich], 
                                                alt_enrichments[alt][cancer]['q'][mypancan_not_enrich]) #alt_enrichments[alt][cancer]['q'][mypancan_not_enrich])
            alt_scores['pathway_correlation'][ind] = corp if not math.isnan(pathway_cor) else 1
            alt_scores['pathway_correlation_rho'][ind] = pathway_cor
            canc_enrich = alt_enrichments[alt][cancer]['p'][mypancan_not_enrich] < .05
            md_enrich = md_enrichments[icd]['b'][mypancan_not_enrich] #md_enrichments[icd]['p'][mypancan_not_enrich] < .05
            res = np.histogram2d(canc_enrich, md_enrich,2)[0]
            #counts = zip(*tuple([[path[0] < .05, path[1] < .05, path[0] < .05 and path[1] < .05, path[0] < .05 or path[1] < .05] for path in zip(md_enrichments[icd]['p'], alt_enrichments[alt][cancer]['p'])]))

            alt_scores['pathway_overlap'][ind] = 1
            if res[1,1] > 0:
                #pdb.set_trace()
                alt_scores['pathway_overlap'][ind] = max(0, 1 - stats.hypergeom.cdf(res[1,1] - 1, 
                                                                                    sum(mypancan_not_enrich), 
                                                                                    sum(res[1,:]), sum(res[:,1])))
            psel = list(scipy.nonzero(np.array(md_enrich) & np.array(canc_enrich))[0])

            if len(psel) > 0:
                path_genes = []
                path = []
                path_cancer_genes = []
                for p_both in psel:
                    path_genes.append(','.join(pathways[pathways_keys_nopancan[p_both]] & gene_list))
                    path.append(pathways_keys_nopancan[p_both].replace(';',' '))
                    path_cancer_genes.append(','.join(pathways[pathways_keys_nopancan[p_both]] &
                                                      cancer_alt))
                alt_scores['pathway_genes'][ind] = ';'.join(path_genes)
                alt_scores['pathway_sel'][ind] = ';'.join(path)
                alt_scores['pathway_cancer_gene'][ind] = ';'.join(path_cancer_genes)
            for nscore in network_findings:
                alt_scores[nscore + '_set'][ind] = \
                  network_findings[nscore][alt]['disease_score'].loc[icd, cancer]
                n_m_i = list(network_findings[nscore][alt]['disease_score'].index).index(icd)
                c_i = list(network_findings[nscore][alt]['disease_score'].columns).index(cancer)
                alt_scores[nscore + '_genes'][ind] = \
                  network_findings[nscore][alt]['set_connected'][n_m_i][c_i]
                alt_scores[nscore + '_ginfo'][ind] = \
                  network_findings[nscore][alt]['set_connection'][n_m_i][c_i]
            csel = set(coexpression_c.columns) & cancer_alt
            msel = (set(coexpression_c.index) & set(gene_list)) - set(csel)
            if len(csel) > 0 and len(msel) > 0:
                '''
                notm = set(coexpression_c.index) - msel - csel
                cg_score = (1 - stats.binom.cdf((coex_mat > rho_cut).sum(axis = 0) - 1, 
                                                len(msel), 
                                                (coexpression_c.loc[notm, csel] > rho_cut).sum(axis = 0)/float(len(notm))))
                alt_scores['coex_binom'][ind] = my_bh_fdr(np.array(cg_score)).min()
                '''
                #cg_score = np.array([stats.ranksums(coexpression_c.loc[notm, csel_g], coex_mat.loc[:,csel_g])[1]
                #                    for csel_g in csel])
                #notm = set(coexpression_c.index) - msel - csel
                t = [i for (i,g) in enumerate(coexpression_c.index) if g in msel];
                f = [i for (i,g) in enumerate(coexpression_c.index) if g not in (msel | csel)];
                tmp = coexpression_c.loc[:, csel]
                xx = np.apply_along_axis(ranksummat, 0, tmp, t, f)
                bh_cg_score = my_bh_fdr(np.array(xx[1]))
                bh_cg_score[np.array(xx[0]) > 0] = 1 ## only look for significant positive 
                alt_scores['coex_CG'][ind] = bh_cg_score.min()
                alt_scores['coex_z'][ind] = np.array(xx[0])[bh_cg_score.argmin()]
                if bh_cg_score.min() < .05:
                    #pdb.set_trace()
                    alt_scores['coexpr_cg_C'][ind] = \
                      ';'.join(list(tmp.columns[np.nonzero(bh_cg_score < .05)[0]]))
                    #alt_scores['coexpr_cg_M'][ind] = \
                    #  ';'.join([','.join(tmp.index[tmp.iloc[t,i] > rho_cut])
                    #            for i in np.nonzero(bh_cg_score < .05)[0]])
                    mdlab = tmp.index[t]
                    alt_scores['coexpr_cg_M'][ind] = \
                      ';'.join([','.join(mdlab[tmp.iloc[t,i] > rho_cut].values)
                                for i in np.nonzero(bh_cg_score < .05)[0]])
                    alt_scores['coexpr_p'][ind] = ';'.join([str('%1.2e' % x)
                                                           for x in bh_cg_score[np.nonzero(bh_cg_score < .05)]])
    h = open(out_name + 'comorb_scores.pkl','w')
    pickle.dump(alt_scores, h)
    h.close()
    if not report: return
    ## get scores and also write file
    ## make contingency table (custom for each score type.  etc.
    inds = [i for i in range(len(alt_scores['relative_risk'])) if alt_scores['relative_risk'][i] > -1]
    stat_tab = [enrich_comorb(alt_scores, covar, inds) for covar in to_crosstab] #+ expr_guys
    stat_tab_h.write(alt + '\t' + '\t'.join([str(x) for x in stat_tab]) + '\n')
    tab = np.zeros((2,2))
    for ind in inds:
        similar = int(min([abs(alt_scores[covar][ind]) for covar in to_crosstab]) < .05)
        comorb = int(alt_scores['relative_risk'][ind] > 1)
        tab[similar, comorb] += 1
    score = 1 - stats.hypergeom.cdf(tab[1,1] - 1, len(inds), sum(tab[1,:]), sum(tab[:,1]))
    to_ret = pd.DataFrame(np.array(stat_tab + [score]), index = to_crosstab + ['Any Metric'],
                          columns = ['comorb. association'])
    alt_report = open(out_name + alt + '.xls','w')
    alt_report.write('\t'.join(header) + '\n')
    for comorb in zip(*tuple([[str(x) for x in alt_scores[var]] for var in header])):
        alt_report.write('\t'.join(comorb) + '\n')
    alt_report.close()
    #pdb.set_trace()
    columns = ['MD','C','gene_enrichment','gene_intersection',
               'pathway_correlation','pathway_genes','pathway_sel','pathway_cancer_gene',
                'coex_CG','coexpr_cg_M','coexpr_cg_C','coexpr_p'] + \
        list(chain.from_iterable([[score + s for s in ['_set','_genes','_ginfo']] for score in network_findings]))
    dfpairs = pd.DataFrame(alt_scores)
    dfpairs = dfpairs.loc[np.array(alt_scores['relative_risk']) > 1, columns]
    to_correct = ['gene_enrichment','pathway_correlation','coex_CG'] + \
      [score + '_set' for score in network_findings]
    sig = np.zeros([dfpairs.shape[0],1])
    towrite = dfpairs.loc[:,['MD','C']]
    for tc in to_correct:
        tocorrect = np.array(dfpairs.loc[:,tc])
        tocorrect[tocorrect > -1] = my_bh_fdr(tocorrect[tocorrect > -1])
        sig[:,0] = (sig[:,0]==1) | ((tocorrect < .1) & (tocorrect > -1))
        towrite.loc[:,tc] = tocorrect
    scores = network_findings.keys() + ['pathway','coexpression','geneIntersection']
    for score in scores:
        towrite.loc[:,score] = ['']*towrite.shape[0]
    coldict = dict()
    for score in scores:
        coldict[score] = towrite.columns.get_loc(score)

    for i in range(towrite.shape[0]):
        if (towrite.iloc[i,:].loc['pathway_correlation'] < .1) & (towrite.iloc[i,:].loc['pathway_correlation'] > -1):
            psets = zip(*tuple((dfpairs.iloc[i,:].loc['pathway_genes'].split(';'),
                                dfpairs.iloc[i,:].loc['pathway_sel'].split(';'),
                                dfpairs.iloc[i,:].loc['pathway_cancer_gene'].split(';'))))
            towrite.iloc[i,coldict['pathway']] = ';'.join([res[0] + ' -> ' + res[1] + '(' + res[2] + ')' for res in psets])
        if (towrite.iloc[i,:].loc['coex_CG'] < .1) & (towrite.iloc[i,:].loc['coex_CG'] > -1):
            coex = zip(*tuple((dfpairs.iloc[i,:].loc['coexpr_cg_M'].split(';'),
                               dfpairs.iloc[i,:].loc['coexpr_cg_C'].split(';'),
                               dfpairs.iloc[i,:].loc['coexpr_p'].split(';'))))
            towrite.iloc[i,coldict['coexpression']] = ';'.join([res[0] + ' -> ' + res[1] + '(' + res[2] +')' for res in coex])
        #if (towrite.iloc[i,:].loc['gene_enrichment'] < .1) & (towrite.iloc[i,:].loc['gene_enrichment'] > -1):
        towrite.iloc[i,coldict['geneIntersection']] = dfpairs.iloc[i,:,].loc['gene_intersection']
        for score in network_findings:
            if (towrite.iloc[i,:].loc[score + '_set'] < .25) & (towrite.iloc[i,:].loc[score + '_set'] > -1):
                conn = zip(*tuple((dfpairs.iloc[i,:].loc[score + '_genes'].split(','),
                                   dfpairs.iloc[i,:].loc[score + '_ginfo'].split(';'))))
                #conn = zip(*tuple((dfpairs.iloc[i,:].loc[score + '_genes'].split(','),
                #                   [x.replace(';',',') for x in dfpairs.iloc[i,:].loc[score + '_ginfo'].split(',')])))
                towrite.iloc[i,coldict[score]] = ';'.join([res[0] + ' -> ' + res[1] for res in conn])
    import csv
    reord = ['MD','C','gene_enrichment','geneIntersection','pathway_correlation','pathway','coex_CG','coexpression'] + \
      list(chain.from_iterable([[network + '_set', network] for network in network_findings]))
    towrite.iloc[sig[:,0]==True,:].loc[:,reord].to_csv(out_name + 'report_pairs.xls',sep='\t',quoting=csv.QUOTE_NONE, index=False)

    stat_tab_h.close()
    return to_ret #alt_scores
