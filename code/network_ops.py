import mendelian_code
import mendelian_mutation
import networkx
import scipy
import scipy.stats as stats
from itertools import chain
import pdb
import random
import pandas as pd
import numpy as np
import pickle

alterations_to_analyze = ['peak_mut'] # 'mutation','amp_peak','del_peak', 'peak_mut']
def comp(measure, rand, score): 
    if measure=='SHORTEST_PATH':
        return rand <= score
    elif measure=='NEIGHBOR_COUNT':
        return rand >= score

def randomized_score(measure, nrand, netdir, data_in):
    from itertools import chain
    
    icd_gene_clinical, cancer_info, mendelian_genes, alterations, germline_list = pickle.load(open(data_in  + '/dat.pkl'))
    
    if measure == 'SHORTEST_PATH_SEED_SHUFFLE':
        return shortest_path_seed_shuffle(hprd, alterations, icd_gene_clinical, nrand)

    scores = network_measure(measure,   netdir + '/network', alterations, icd_gene_clinical, cancer_info)
    rand_score = dict()
    for alt in scores:
        rand_score[alt] = {'disease_score':scores[alt]['disease_score'].copy(deep=True),
                           'gene_score':scores[alt]['gene_score'].copy(deep=True)}
        #'cancers':scores[alt]['cancers'],'MD':scores[alt]['MD']}
        rand_score[alt]['disease_score'].iloc[:,:] = 0
        rand_score[alt]['gene_score'].iloc[:,:] = 0

    # do random iterations!
    for i in range(1,nrand + 1):
        print str(i),
        #net = networkx.read_edgelist(netdir + '/rand/' + str(i) + '.net', nodetype = str)
        rand_res = network_measure(measure, netdir + '/rand/' + str(i) + '.net', alterations, icd_gene_clinical,cancer_info, scores)
        for alt in scores:
            rand_score[alt]['disease_score'] += comp(measure, rand_res[alt]['disease_score'],
                                                     scores[alt]['disease_score'])
            rand_score[alt]['gene_score'] += comp(measure, rand_res[alt]['gene_score'],
                                                     scores[alt]['gene_score'])
    for alt in rand_score:
        rand_score[alt]['disease_score'] *= 1/float(nrand)

    #bg_nc, bg_sp = pickle.load(open(netdir + '/background.pkl'))
    #background = bg_nc
    #if measure=='SHORTEST_PATH': background = bg_sp
    add_genewise_score(icd_gene_clinical, rand_score, scores, nrand, alterations)

    pkl = open(data_in + '/' + measure + '.' +netdir + '.'+ str(nrand) + '.pkl','w')
    pickle.dump((scores, rand_score),pkl)
    pkl.close()

    return scores, rand_score


def add_genewise_score(icd_gene_clinical, rand_score, measure, nrand, alterations):
    all_mend_gn = mendelian_code.get_mendelian_genes(icd_gene_clinical)
    for alt in rand_score:
        
        ### add genewise scores
        genewise_score = rand_score[alt]['disease_score'].copy(deep=True)
        genewise_score.iloc[:,:] = 1
        #disease_score = rand_score[alt]['disease_score'].copy(deep=True)
        #disease_score.iloc[:,:] = 1
        gene_winner = [['']*rand_score[alt]['disease_score'].shape[1] for i in range(len(icd_gene_clinical))]
        gene_connect = [['']*rand_score[alt]['disease_score'].shape[1] for i in range(len(icd_gene_clinical))]        
        set_connected = [['']*rand_score[alt]['disease_score'].shape[1] for i in range(len(icd_gene_clinical))]
        set_connection = [['']*rand_score[alt]['disease_score'].shape[1] for i in range(len(icd_gene_clinical))]        

        for row, md in enumerate(rand_score[alt]['disease_score'].index):
            #mend_gn = icd_gene_clinical[md]['gene_omim'].keys()
            #if alt == 'del_peak' and md == 'Specified Anomalies of the Musculoskeletal System':
            #    pdb.set_trace()
            mend_gn = [g for g in all_mend_gn if g in icd_gene_clinical[md]['gene_omim']]

            if len(mend_gn) == 0: continue
            gene_scores = rand_score[alt]['gene_score'].loc[mend_gn,:]
            

            #B background_set =  set(background[alt].index) - set(mend_gn)
            #B background_probability = [0]*gene_scores.shape[1]
            #B gscore_vs_background = [0]*gene_scores.shape[1]
            gene_sp = ['']*gene_scores.shape[1]
            connect = ['']*gene_scores.shape[1]
            set_gene_con = ['']*gene_scores.shape[1]
            set_connect = ['']*gene_scores.shape[1]
            genewise_score.iloc[row,:] = gene_scores.min(axis=0)    
            gsa = np.array(gene_scores)
            gsel = np.nonzero(scipy.logical_and(scipy.logical_or(gsa == scipy.tile(gsa.min(axis=0), (gsa.shape[0],1)),
                                                                gsa <= .05*float(nrand)),
                                                gsa < nrand))

            for (c_i, c) in enumerate(gene_scores.columns):
                #B background_genes = background_set - alterations[alt][c][1]
                #B p = (background[alt].loc[background_genes, c] < .05).sum()/float(len(background_genes))
                #B background_probability[c_i] =p
                #B n_siggene = (gene_scores.loc[:,c] < .05*float(nrand)).sum(axis=0)
                #B gscore_vs_background[c_i] = 1 - stats.binom.cdf(n_siggene - 1, gene_scores.shape[0], p)
                #B if n_siggene == 0:
                #B    gscore_vs_background[c_i] = 1.0
                genes =[]
                conns = []
                for i in gsel[0][gsel[1]==c_i]:
                    g = gene_scores.index[i]
                    genes += [g]
                    conns += [measure[alt]['connection'][all_mend_gn.index(g)][c_i].strip(',')]
                gene_sp[c_i] = ','.join(genes)
                connect[c_i] = ';'.join(conns)
                genes = []
                conns = []
                #if (md =='"Pervasive, Specified Congenital Anomalies"' ) and (c == 'GBM'): pdb.set_trace()
                for g in gene_scores.index:
                    c_conn = measure[alt]['connection'][all_mend_gn.index(g)][c_i].strip(",")
                    if not c_conn == '':
                        genes += [g]
                        conns += [c_conn]
                set_gene_con[c_i] = ','.join(genes)
                set_connect[c_i] = ';'.join(conns)

            #B genewise_score.iloc[row,:] = gscore_vs_background
            gene_winner[row] = gene_sp
            gene_connect[row] = connect
            set_connected[row] = set_gene_con
            set_connection[row] = set_connect

        rand_score[alt]['genewise_disease_score'] = genewise_score*1/float(nrand)
        rand_score[alt]['genewise_disease_best'] = gene_winner
        rand_score[alt]['genewise_disease_connect'] = gene_connect
        rand_score[alt]['set_connected'] = set_connected
        rand_score[alt]['set_connection'] = set_connection
        


def network_measure(measure, network, alterations, icd_gene_clinical, cancer_info, scores=[]):
    if measure=='SHORTEST_PATH':
        # if compare to REAL scores, then this is random run, dont get connector
        get_connector = True
        if len(scores) > 0: get_connector = False
        return shortest_path(network, alterations, icd_gene_clinical, get_connector, scores)
    if measure=='NEIGHBOR_COUNT':
        alt_scores = dict()
        net = networkx.read_edgelist(network, nodetype = str)
        for alt in alterations_to_analyze:
            alt_scores[alt] = neighbor_count_md_c(net, alterations[alt], icd_gene_clinical, cancer_info)
        return alt_scores



def neighbor_count_md_c(network, alteration_type, icd_gene_clinical, cancer_info, comorbid_only = False, comorb_perm = False):
    an_edge = network.edges()[0]
    weighted = len(network[an_edge[0]][an_edge[1]].keys()) > 0
    mend_disease = mendelian_code.get_mendelian_diseases(icd_gene_clinical)
    cancers = alteration_type.keys()
    #neighbor_count_mat = scipy.zeros([len(mend_disease), len(cancers)])
    all_mend_gn = mendelian_code.get_mendelian_genes(icd_gene_clinical)
    #gn_count_mat = scipy.zeros([len(all_mend_gn), len(cancers)])

    disease_score = pd.DataFrame(scipy.zeros([len(mend_disease), len(cancers)]), index = mend_disease, columns = cancers)
    gene_score = pd.DataFrame(scipy.zeros([len(all_mend_gn), len(cancers)]), index = all_mend_gn, columns = cancers)

    gene_connection = [['']*len(cancers) for i in range(len(all_mend_gn))]
    v = zip(*tuple([(c, canc) for canc in cancers for c in cancer_info if canc in cancer_info[c]['TCGA'] ]))
    canc_tab = dict(zip(v[1],v[0]))
    if comorb_perm:
        for md in mend_disease:
            canc_tab_in =[k for k in icd_gene_clinical[md]['cancer_assoc'].keys() if k in canc_tab.values()]
            canc_tab_val = scipy.random.permutation([ icd_gene_clinical[md]['cancer_assoc'][c] for c in canc_tab_in])
            icd_gene_clinical[md]['cancer_assoc'] = dict(zip(canc_tab_in, canc_tab_val))
        
    for (m, md) in enumerate(mend_disease):
        mend_gn = icd_gene_clinical[md]['gene_omim'].keys()
        md_score = [0]*len(cancers)
        if len(mend_gn)==0:  
            continue
        for (c, canc) in enumerate(cancers):
            #if not canc in canc_tab: continue
            #if (md =='"Pervasive, Specified Congenital Anomalies"' ) and (canc == 'GBM'): pdb.set_trace()
            if comorbid_only and (not canc in canc_tab or icd_gene_clinical[md]['cancer_assoc'][canc_tab[canc]] < 1): continue
            seed_gene = list(alteration_type[canc][1])
            md_subgr = networkx.subgraph(network, set(mend_gn) | set(seed_gene))
            for gn in mend_gn:
                if gn in md_subgr and len(md_subgr[gn]) > 0:
                    res = zip(*tuple([[neighbor, md_subgr[gn][neighbor]]
                                                for neighbor in networkx.all_neighbors(md_subgr,gn) 
                                                if neighbor in seed_gene]))
                    if len(res) > 0:
                        gene_result, weights = res
                        g_i = all_mend_gn.index(gn)
                        if len(gene_result) > 0:
                            gene_connection[g_i][c] = ','.join(gene_result)
                        connection = len(gene_result)
                        if weighted:
                            connection = sum([w['weight'] for w in weights])
                        gene_score.loc[gn, canc] = connection
                        md_score[c] += connection
                        
        disease_score.loc[md,:] = md_score
    return {'disease_score':disease_score, 'gene_score':gene_score, 
            'connection':gene_connection} # 'cancers':cancers, 'MD':mend_disease}


def neighbor_count_comorbid(network, alteration_type, icd_gene_clinical, cancer_info,
                            comorbid_only = False, comorb_perm = False, weighted=False):
    tumor_data_list = alteration_type.keys()
    list2 = set()
    for icd in icd_gene_clinical:
        for cancer in tumor_data_list:
            rr = [icd_gene_clinical[icd]['cancer_assoc'][cancer_icd] 
                  for cancer_icd in cancer_info 
                  if cancer in cancer_info[cancer_icd]['TCGA']]
            if len(rr)>0 and rr[0] > 1:
                list2 |= set([cancer])
    cancers = list2

    mend_disease = mendelian_code.get_mendelian_diseases(icd_gene_clinical)
    #neighbor_count_mat = scipy.zeros([len(mend_disease), len(cancers)])
    all_mend_gn = mendelian_code.get_mendelian_genes(icd_gene_clinical)
    #gn_count_mat = scipy.zeros([len(all_mend_gn), len(cancers)])

    disease_score = pd.DataFrame(scipy.zeros([len(mend_disease), len(cancers)]), index = mend_disease, columns = cancers)
    gene_score = pd.DataFrame(scipy.zeros([len(all_mend_gn), len(cancers)]), index = all_mend_gn, columns = cancers)

    gene_connection = [['']*len(cancers) for i in range(len(all_mend_gn))]
    v = zip(*tuple([(c, canc) for canc in cancers for c in cancer_info if canc in cancer_info[c]['TCGA'] ]))
    canc_tab = dict(zip(v[1],v[0]))
    BIG_COUNT = 0
    i = 0
    for (m, md) in enumerate(mend_disease):
        mend_gn = icd_gene_clinical[md]['gene_omim'].keys()
        md_score = [0]*len(cancers)
        if len(mend_gn)==0:  
            continue
        canc_gn = set()
        for (c, canc) in enumerate(cancers):
            if not canc in canc_tab: continue
            if comorbid_only and (not canc in canc_tab or icd_gene_clinical[md]['cancer_assoc'][canc_tab[canc]] < 1): continue
            i += 1
            canc_gn |= alteration_type[canc][1]
        #pdb.set_trace()            
        md_subgr = networkx.subgraph(network, set(mend_gn) | set(canc_gn))
        for gn in mend_gn:
            if gn in md_subgr:
                if not weighted:
                    gene_result = [neighbor for neighbor in networkx.all_neighbors(md_subgr,gn) 
                                    if neighbor in canc_gn]
                    BIG_COUNT += len(gene_result)
                else:
                    BIG_COUNT += sum([md_subgr[gn][neighbor]['weight']
                                      for neighbor in networkx.all_neighbors(md_subgr,gn)
                                      if neighbor in canc_gn])
    #print 'tot pair: ' + str(i) + ' tot count: ' + str(BIG_COUNT)
    return BIG_COUNT

def random_network(graph, dirwrite, start=0, end=0, weighted=False):
    import os
    #dirn = 'rand'
    if not os.path.exists(dirwrite): os.mkdir(dirwrite)
    todo = range(1000)
    if start > 0:
        todo = range(start, end)
    for i in todo:
        towrite = dirwrite  + '/' + str(i) + '.net'
        if os.path.exists(towrite): continue
        rnd = graph
        if not weighted:
            rnd = networkx.double_edge_swap(graph, nswap=graph.number_of_edges()*4,max_tries = graph.number_of_edges()*8)
        else:
            rnd = double_edge_swap_weighted(graph, nswap=graph.number_of_edges()*4,max_tries = graph.number_of_edges()*8)
        networkx.write_edgelist(rnd, towrite, data=weighted,delimiter="\t")
        
def comorbid_count_compare(net_dir, icd_gene_clinical, cancer_info, alterations, weighted=False):
    # = 'humannet.9'
    graph = networkx.read_edgelist(net_dir + '/network',nodetype=str)
    ct = neighbor_count_comorbid(graph, alterations['peak_mut'], icd_gene_clinical, cancer_info, comorbid_only = True, weighted=weighted)
    import os
    randdir = net_dir + '/rand/'
    randnets = os.listdir(randdir)
    x = scipy.zeros([len(randnets)])
    for i,f in enumerate(randnets):
        net = networkx.read_edgelist(randdir + f, nodetype = str, data=weighted)
        x[i] = neighbor_count_comorbid(net, alterations['peak_mut'], icd_gene_clinical, cancer_info, comorbid_only = True, weighted = weighted)    
    print 'comorbid_edges= ' + str(ct) + "\tngreater=" +str(sum(x >= ct)) + '\tp=' + str(sum(x >= ct)/float(len(randnets)))
    return ct, x

def prepare_randomization(netdir, dat_dir, nrand=1000, weighted=False):

    print netdir + " " + dat_dir + " " + str(nrand) + " weighted=" + str(weighted)

    net = networkx.read_edgelist(netdir + '/network',nodetype=str)
    random_network(net, netdir + '/rand', start=1, end=nrand + 1, weighted=weighted)
    import os
        
    score, rand_score = randomized_score('NEIGHBOR_COUNT', nrand, netdir, dat_dir)
    #pkl = open(dat_dir + '/NEIGHBOR_COUNT.' +netdir + str(rand) + '.pkl','w')

    # qsub -N humannet.7.unwt -V -cwd -l mem=4G,time=8:: -b y 'python -c "import network_ops; network_ops.prepare_randomization(\"humannet.7.unwt\",\"mendelian\",1000)"'
    # qsub -N humannet.7.wt -V -cwd -l mem=4G,time=8:: -b y 'python -c "import network_ops; network_ops.prepare_randomization(\"humannet.7.wt\",\"mendelian\",1000,weighted=True)"'



