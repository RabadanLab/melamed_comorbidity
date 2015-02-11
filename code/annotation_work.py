from sets import Set
import gene_info
def omim_annotation(omim_file, gene_annot): # "morbidmap.txt",'../Homo_sapiens.gene_info.protein_coding')
    #print " annot is " + gene_annot
    symbol_info, alias_to_symbol = gene_info.process_ncbi(gene_annot)
    omimmap = open(omim_file)
    omim_dict = dict()
    gene_dict = dict()
    locus_dict = dict()
    ddesc_dict = dict()
    omim_err = open("data_processed/log_morbidmap.err",'w')
    for line in omimmap:
        linkage = line.strip().split("|")
        #print line.strip()
        genes = [item.strip() for item in linkage[1].split(",")]
        map_loc_chrom, map_loc_arm, map_loc_band = gene_info.maploc_to_band(linkage[3])
        genes_dealiased = gene_info.unaliased_set(genes, map_loc_chrom, map_loc_arm, alias_to_symbol, symbol_info, omim_err)
        omim_locus = linkage[2]
        disease_info_method = linkage[0].split("(") # first split off the "method" ie: (3)
        omim = "NA"
        link_map="0"
        disease = disease_info_method[0]
        if len(disease_info_method) > 1:
            link_map = disease_info_method[len(disease_info_method)-1].strip().strip(")")
        if len(disease_info_method) > 2:  # there's a paren in the descripton...
            disease_info_method = ['('.join(disease_info_method[0:(len(disease_info_method)-1)])]
        disease_omim = disease_info_method[0].strip().split(" ")
        if len(disease_omim[len(disease_omim) - 1])  < 6:
            disease_omim = [" ".join(disease_omim)]
        try:
            omim = str(int(disease_omim[len(disease_omim)-1]))
            disease = " ".join(disease_omim[0:(len(disease_omim)-1)])
        except ValueError:
            ### last try: see if can split by comma
            disease_omim = disease_info_method[0].split(",")
            if len(disease_omim[len(disease_omim) - 1]) < 6:
                disease_omim = [",".join(disease_omim)]
            try:
                omim = str(int(disease_omim[len(disease_omim)-1]))
                disease = ",".join(disease_omim[0:(len(disease_omim)-1)])
            except ValueError:
                omim = omim_locus
        if omim == "NA":
            omim_err.write("NO_OMIM_FOUND " + line.strip()+ "\n")
        if omim in omim_dict:
            omim_dict[omim]['genes'] |= genes_dealiased
            omim_dict[omim]['omim_locus'] |= Set([omim_locus])
        else:
            omim_dict[omim] = {'disorder':Set(),'genes':genes_dealiased,'omim_locus':Set([omim_locus])}
        omim_dict[omim]['disorder'].add((disease, link_map))
        for gene in genes_dealiased:
            if not 'omim' in symbol_info[gene]:
                symbol_info[gene]['omim'] = Set([omim])
                symbol_info[gene]['omim_locus'] = Set([omim_locus])
            else:
                symbol_info[gene]['omim'].add(omim)
                symbol_info[gene]['omim_locus'] |= Set([omim_locus])
        if omim_locus in locus_dict:
            locus_dict[omim_locus]['omim'].add(omim)
            if len(locus_dict[omim_locus]['genes'] ^ genes_dealiased) > 0:
                omim_err.write("New_locus_genes " + omim_locus + " for omim " + omim+ "\n")
                locus_dict[omim_locus]['genes'] |= genes_dealiased
        else:
            locus_dict[omim_locus]={'omim':Set([omim]), 'genes':genes_dealiased}
        if disease in ddesc_dict:
            omim_err.write("multiple_disease_entry for " + disease+ "\n")
            ddesc_dict[disease]['omim'].add(omim_locus)
            ddesc_dict[disease]['omim_locus'] |= Set([omim_locus])
        elif not disease == "":
            ddesc_dict[disease]={'omim':Set([omim]),'omim_locus':Set([omim_locus])}

    outomim = open("data_processed/log_morbidmap_by_omim.txt","w")
    for omim in omim_dict:
        omim_disorders = [dis[0] + "(" +dis[1] + ")" for dis in omim_dict[omim]['disorder']]
        omim_disorders_no_linkinfo = [dis[0] for dis in omim_dict[omim]['disorder']]
        if len(omim_dict[omim]['omim_locus']) > 1:
            omim_err.write("MULTI_LOCUS_PER_OMIM " + omim + " loci:" + ",".join(omim_dict[omim]['omim_locus'])+ "\n")
        if len(omim_dict[omim]['disorder']) > 1:
            omim_err.write(" MULTI_DISORDER_PER_OMIM " + omim + " disorder:" + ",".join(omim_disorders_no_linkinfo)+ "\n")
        outomim.write("\t".join([omim,";".join(omim_disorders_no_linkinfo)]) + "\t")
        outomim.write("\t".join([";".join(omim_disorders), ",".join(omim_dict[omim]['genes']),",".join(omim_dict[omim]['omim_locus'])]))
        outomim.write("\n")
    outomim.close()


    outgenes = open("data_processed/log_morbidmap_genes.txt","w")
    for gene in symbol_info:
        if not 'omim' in symbol_info[gene]:
            continue
        if len(symbol_info[gene]['omim_locus']) > 1:
            omim_err.write(gene + " MULTI-LOCUS " + ",".join(symbol_info[gene]['omim_locus'])+ "\n")
        if len(symbol_info[gene]['omim']) > 1:
            omim_err.write(gene + " MULTI-DISEASE " + ",".join(symbol_info[gene]['omim'])+ "\n")
        outgenes.write("\t".join([gene, ",".join(symbol_info[gene]['omim']), ",".join(symbol_info[gene]['omim_locus'])])+"\n")
    outgenes.close()    

    outdisorder = open("data_processed/log_morbidmap_disorders.txt","w")
    for dis in ddesc_dict:
        if len(ddesc_dict[dis]['omim_locus']) > 1:
            omim_err.write(dis + " DISEASE-MULTI-LOCUS " + ",".join(ddesc_dict[dis]['omim_locus']))
        if len(ddesc_dict[dis]['omim']) > 1:
            omim_err.write(dis + " DISEASE-MULTI-OMIM " + ",".join(ddesc_dict[dis]['omim']))
        outdisorder.write("\t".join([dis, ",".join(ddesc_dict[dis]['omim']), ",".join(ddesc_dict[dis]['omim_locus'])])+"\n")
    outdisorder.close()    


    omim_err.close()
    return omim_dict, symbol_info, alias_to_symbol, locus_dict
