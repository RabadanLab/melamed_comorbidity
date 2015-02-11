import gene_info
import networkx
ncbi_genes = "data_source/Homo_sapiens.gene_info.protein_coding"
def humannet(humannet_file, ncbi_genes, netdir):
    ncbi = open(ncbi_genes)
    ncbi_id = dict()
    for gn in ncbi:
        gline = gn.strip().split("\t")
        ncbi_id[gline[0]] = gline[1]
    symbol_info, alias_to_symbol = gene_info.process_ncbi(ncbi_genes)

    hnet = open(humannet_file) #'/ifs/scratch/c2b2/rr_lab/rdm2114/mendelian_work/humannet/HumanNet.v1.join.txt') #)    
    hnet.readline() # header
    humannet = networkx.Graph()
    use_weight = True
    rep = open('humannet_rep','w')
    i = 0
    j = 0
    k = 0
    u = 0
    for net_line in hnet:
        edge = net_line.strip().split("\t")
        if not edge[0] in ncbi_id or not edge[1] in ncbi_id:
            u += 1
        if not edge[0] in ncbi_id:
            i += 1
            rep.write(edge[0] + '\n')
            continue
        if not edge[1] in ncbi_id:
            j += 1
            rep.write(edge[1] + '\n')
            continue
        weight = float(edge[len(edge)-1])
        if weight < 2.169:
            continue
        k += 1
        g1 = ncbi_id[edge[0]]
        g2 = ncbi_id[edge[1]]
        if use_weight:
            humannet.add_edge(g1, g2, weight=weight)
        else:
            humannet.add_edge(g1, g2)
    rep.close()            


    networkx.write_edgelist(humannet, netdir + '/network', data=use_weight,delimiter="\t")
            x = list([list(chain.from_iterable([graph[mdg][cg]['weight'] for cg in cancer_alt if cg in graph[mdg]])
                           for mdg in gene_list if mdg in graph]))

