import numpy as np
import pandas as pd
import pickle
outdir = 'data_processed/fantom'

def fantom_tpm_to_expr(tpm_file, symbol_info, alias_to_symbol, cancer_genes):
    import os
    try:
        os.mkdir(outdir)
    except OSError:
        pass
    import subprocess
    get_gns = 'cut -f 2 ' + tpm_file + ' | grep -v chr\|ColumnVariables | tr "," "\n" | cut -d "@" -f 2 | sort -u > ' + outdir + '/gns'
    subprocess.call(get_gns,shell=True)
    # cut -f 2 hg19.cage_peak_tpm_ann.osc.txt | grep -v chr\|ColumnVariables | tr "," "\n" | cut -d "@" -f 2 | sort -u > gns
    gns  = open(outdir + '/gns')
    fantom_gns = gns.read().split("\n")

    import csv
    al_in_fant = (set(fantom_gns) & set(alias_to_symbol.keys())) - set(symbol_info.keys())
    x=  [al for al in al_in_fant if len(alias_to_symbol[al]) > 1]
    al_in_fant = al_in_fant - set(x)
    deal = [alias_to_symbol[al][0] for al in al_in_fant]
    aliased = (set(symbol_info.keys()) - set(fantom_gns)) & set(deal)
    fantom_gene_names = list((set(fantom_gns) & set(symbol_info.keys())) | aliased)

    tpm = csv.reader(open(tpm_file),delimiter='\t')
    hdr = tpm.next()
    while hdr[0].startswith("##"):
        hdr = tpm.next()
    logfile = open(outdir + '/log','w')
    expr_dat_all = pd.DataFrame(np.zeros([len(fantom_gene_names), len(hdr) - 7]), index=fantom_gene_names,
                                columns = hdr[7:])
    for peakline in tpm:
        for peak in peakline[1].split(","):
            peak_gene = peak.split("@")
            if len(peak_gene) > 1 and not peak_gene[1].startswith('chr'):
                exprs = np.array([float(xi) for xi in peakline[7:]])
                gene_name = peak_gene[1]
                try: 
                    expr_dat_all.loc[gene_name,:] += exprs
                except KeyError:
                    try:
                        expr_dat_all.loc[alias_to_symbol[gene_name],:] += exprs
                    except KeyError:
                        if gene_name in alias_to_symbol:
                            logfile.write('Aliased: fant name = ' + gene_name + ' al to ' + ','.join(alias_to_symbol[gene_name]) + '\n')
    logfile.close()
    #import pdb
    #pdb.set_trace()
    save_to_pickle(expr_dat_all, outdir + '/fantom_all_expr.pkl')
    to_cor_gns = outdir + '/to_cor.pkl'
    pkl = open(to_cor_gns,'w')
    pickle.dump((list(cancer_genes & set(expr_dat_all.index)), expr_dat_all.index), pkl)
    return to_cor_gns

def save_to_pickle(data, name):
    version = '.'.join(pd.__version__.split('.')[1:])
    if float(version) < 11.1:
        #print 'version old'
        data.save(name)
    else:
        #print 'version new'
        data.to_pickle(name)
            
def do_cor(cor_gene_file):
    expr_dat_all = pd.load(outdir + '/fantom_all_expr.pkl')
    (col_gn, row_gn) = pickle.load(open(cor_gene_file))
    col_gn = col_gn[:10]
    cor_mat = pd.DataFrame(np.zeros([len(row_gn), len(col_gn)]), index = row_gn, columns = col_gn)
    import time
    sta = time.clock()
    stds = expr_dat_all.T.std()
    for (c_i, c) in enumerate(col_gn):
        for (m_i, m) in enumerate(row_gn):
            cor_mat.iloc[m_i, c_i] = np.corrcoef(expr_dat_all.loc[m,:], expr_dat_all.loc[c,:])[0,1] 
        if c_i % 500 == 0:
            print 'at c = ' + str(c) + ' # ' + str(c_i)
    sto = time.clock()
    save_to_pickle(cor_mat, outdir + '/fantom_corrcoef.' + str(expr_dat_all.shape[1]) + '.pkl') # + cor_gene_file.split(".")[0] 
    print 'done ' + str(sto - sta)

if __name__ == '__main__':
    import sys
    do_cor(sys.argv[1])

