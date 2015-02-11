chromlist = [str(x) for x in range(1,23) ]
chromlist.extend(['X','Y'])
arm_list = ['p','q']

def process_ncbi(ncbi_gene_info_path):
    alias_to_symbol = dict()
    gene_info = dict()
    genes = open(ncbi_gene_info_path) #../Homo_sapiens.gene_info.protein_coding")
    errfile = open('data_processed/log_process_ncbi.err','w')
    for line in genes:
        gene_line = line.strip().split('\t')
        symbol = gene_line[1]
        ncbi_id = gene_line[0]
        chrom = gene_line[5]
        cytoband = gene_line[6]
        chrom2, arm, band = maploc_to_band(cytoband)
        if not chrom2 in chromlist:
            errfile.write("MAPLOC_NOT_PARSED:" + "gene " + symbol + " " +  cytoband + "\n") 
        if not chrom2 == chrom:
            errfile.write("Chrom_misID: " + chrom + " CB" + cytoband + " -->" + chrom2 + "\n")
        gene_info[symbol] = {'ncbi_id':ncbi_id, 'chr':chrom, 'map_location':cytoband, 'arm':arm, 'band':band}
        if symbol in alias_to_symbol:
            errfile.write( "Symbol-Was-Alias\t" + alias +"\t" + ",".join(alias_to_symbol[alias]) + "\t" + ncbi_id + "\t OldChr: " + gene_info[alias_to_symbol[symbol][0]]['chr'] + "\t NewChr:" + chrom + "\n")
            alias_to_symbol[symbol].append(symbol)
        else:
            alias_to_symbol[symbol] = [symbol]
        if not gene_line[3] == "-":
            aliases = gene_line[3].split("|")
            for alias in aliases:
                if alias in alias_to_symbol:
                    # print "NCBI:" + ncbi_id
                    #print str(len(alias_to_symbol[alias])) + " GRR " +  alias_to_symbol[alias][0]
                    errfile.write( "DoubleAlias\t" + alias +"\t" + ",".join(alias_to_symbol[alias]) + "\t" + ncbi_id + "\n" )
                    alias_to_symbol[alias].append(symbol)
                else:
                    alias_to_symbol[alias] = [symbol]
    genes.close()
    errfile.close()
    return gene_info, alias_to_symbol 


def maploc_to_band(map_loc):
    map_loc_chrom = map_loc.replace("p","q").split("q")[0]
    if not map_loc_chrom in chromlist:
        try2 = map_loc.split(".")
        if len(try2) > 1 and try2[1] in chromlist:
            map_loc_chrom = try2[1]
    map_loc_arm = 'n'
    map_loc_band = '0'
    if len(map_loc) > len(map_loc_chrom) and map_loc[len(map_loc_chrom):(len(map_loc_chrom) +1)] in arm_list:
        map_loc_arm = map_loc[len(map_loc_chrom):(len(map_loc_chrom) +1)]
        if len(map_loc) > (len(map_loc_chrom) + len(map_loc_arm)):
            map_loc_band = map_loc[(len(map_loc_chrom) + len(map_loc_arm)):len(map_loc)]
    return map_loc_chrom, map_loc_arm, map_loc_band


from sets import Set
def unaliased_set(genes, query_chr, query_arm, alias_to_symbol, symbol_info, omim_err):
    # print 'Getting arm:' + query_arm + " chr:"+ query_chr + " genes = "+ ",".join(genes)
    genes_dealiased = Set()
    for g in genes:
        if not g in alias_to_symbol:
            omim_err.write( "NOT_ANNOT:" + g + "\n")
            continue
        dealiased_same_chrom = []
        for proper_symbol in alias_to_symbol[g]:
            if symbol_info[proper_symbol]['chr'] == query_chr or not query_chr in chromlist:
                omim_err.write( 'adding ' + proper_symbol + ' is in chr ' + symbol_info[proper_symbol]['chr'] + " as with " + query_chr +"\n")
                dealiased_same_chrom.append(proper_symbol)
        if len(dealiased_same_chrom) > 1 and query_arm in arm_list:
            dealiased_same_chrom_arm = []
            for proper_symbol in dealiased_same_chrom:
                if not symbol_info[proper_symbol]['arm'] in arm_list or symbol_info[proper_symbol]['arm'] == query_arm:
                    dealiased_same_chrom_arm.append(proper_symbol)
            if len(dealiased_same_chrom_arm) > 1:
                omim_err.write( "MULTI_ALIASED_CHR_ARM: gene\t"+ g + "\tmatches = "+ ",".join(dealiased_same_chrom_arm) + "\n")
            if len(dealiased_same_chrom) == 0:
                omim_err.write("CHR_ARM_NO_MATCH: gene = " + g + "\n")
            else:
                dealiased_same_chrom = dealiased_same_chrom_arm
        if len(dealiased_same_chrom) == 0:
            omim_err.write( "MISSING_ALIAS_SAMECHR: gene = " + g + "\n")
        else:
            genes_dealiased |= Set(dealiased_same_chrom)
    return genes_dealiased
