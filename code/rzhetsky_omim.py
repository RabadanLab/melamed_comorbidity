import annotation_work

import subprocess
def write_orphas(orpha_hits, fname):
    x = open(fname,'w')
    x.write("\n".join(orpha_hits.keys())+"\norpha\n")
    x.close()
    do='grep -wf ' + fname + ' orpha/orpha_omim.txt > '+fname+ '.xls'    
    subprocess.call(do, shell=True)



from sets import Set


import link_rzhetsky_orpha_omim
#reload ( link_rzhetsky_orpha_omim    )
morbidmap="OMIM/morbidmap.txt"
orpha_parsed='orph2'
ncbi_genes= "Homo_sapiens.gene_info.protein_coding"
mmc3= "mmc3_2"
link_rzhetsky_orpha_omim.load_n_link(morbidmap, orpha_parsed, ncbi_genes, mmc3)


#for o in omims: print o + list(omim_dict[o]['disorder'])[0][0]
            
