import networkx
import os

def biogrid(mitab, ncbi_genes, netdir):
    #mitab = 'BIOGRID-ORGANISM-Homo_sapiens-3.2.119.mitab.txt'
    graph = parseMitab(mitab, ncbi_genes)
    graph.remove_node('UBC')
    for gene in graph:
        if gene in graph[gene]:
            graph.remove_edge(gene, gene)
    os.mkdir('biogrid')
    networkx.write_edgelist(graph, netdir + '/network', data=False,delimiter="\t")
    #adj = networkx.to_scipy_sparse_matrix(graph) #adj = networkx.to_numpy_matrix(graph)
    return graph

def parseMitab(mitab, ncbi_gene_annot):
    
    ncbi = open(ncbi_gene_annot)
    ncbi_id = dict()
    for gn in ncbi:
        gline = gn.strip().split("\t")
        ncbi_id[gline[0]] = gline[1]

    biogrid = open(mitab) #)    
    biogrid.readline() # header
    bg_graph = networkx.Graph()
    rep = open('rep','w')
    other_org = open('other_org','w')
    for ixn_line in biogrid:
        ixn = ixn_line.strip().split("\t")
        if check_tax(ixn[9], rep) or check_tax(ixn[10], rep):
            other_org.write(ixn_line)
            continue

        g1 = get_ncbi_id(ixn[0], ncbi_id, str(i), rep)
        g2 = get_ncbi_id(ixn[1], ncbi_id, str(i), rep)
        #print 'g1:' + g1 + '\tg2:' + g2
        if g1 in ncbi_id and g2 in ncbi_id:
            bg_graph.add_edge(ncbi_id[g1], ncbi_id[g2])
        #g1 = get_gene(ixn[2], ixn[4],symbol_info, alias_to_symbol, str(i))
        #g2 = get_gene(ixn[3], ixn[5], symbol_info, alias_to_symbol, str(i))
    rep.close()
    return bg_graph

def check_tax(tax_chunk, rep):
    tax = tax_chunk.split(":")
    try:
        if tax[0]=="taxid":
            if int(tax[1])==9606:
                return False
            else:
                return True
    except:
        rep.write('failed:' + tax_chunk + '\n')

def get_ncbi_id(gene_chunk, ncbi_id, lineno, rep):
    p1 = gene_chunk.split("|")
    found = '0'
    for possible in p1:
        annot = possible.split(":")
        if annot[0] == 'entrez gene/locuslink':
            if annot[1] in ncbi_id:
                if not found == '0':
                    rep.write(lineno + "--2 found\t" + str(found) + '\n')
                found = annot[1]
    if found == '0':
        rep.write( lineno + "--missed!\t--" + gene_chunk + '\n')
    return found

def get_gene(gene_chunk, alt_chunk, symbol_info, alias_to_symbol, lineno):
    p1 = gene_chunk.split("|")
    p1_aliased = ''
    p1_direct = ''
    for possible in p1:
        gn = possible.split(":")[1]
        if alt_chunk == '': # then this IS the alt-chunk, so split off the paren
            gn = gn.split("(")[0]
        if gn in alias_to_symbol:
            if gn in symbol_info:
                if not p1_direct == '':
                    print lineno + "--2 found\t" + p1_found + "\t" + gn
                p1_direct = gn
            else:
                p1_aliased = alias_to_symbol[gn][0]
                if len(alias_to_symbol[gn]) > 1:
                    print lineno + 'multi-alias:' + gn
    if p1_direct == '' and p1_aliased == '':
        next_try = get_gene(alt_chunk, '', symbol_info, alias_to_symbol, lineno)
        if next_try == '':
            print lineno + "--missed!\t--" + gene_chunk
        else:
            return next_try
    if not p1_direct == '':
        return p1_direct
    else:
        return p1_aliased

