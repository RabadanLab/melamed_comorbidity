import re
import annotation_work
import pdb
def print_orpha_omim_gene(alias, orpha_entry=[], omim="", omim_dict={}):
    if alias == "GET_HEADER":
        return "MMC3_Summary_Name\tMMC3_Assoc_Disorder\talias\tOMIM\tOMIM_disorder\torpha\torpha_disorder\torpha_genes\tdeath\tprevalence\torpha_synonyms"
    else:
        omim_disorders_no_linkinfo = [dis[0] for dis in omim_dict[omim]['disorder']]
        return "\t".join([alias, omim, ";".join(omim_disorders_no_linkinfo), orpha_entry['ICD10'], orpha_entry['orpha'], orpha_entry['name'],";".join(orpha_entry['genes']), orpha_entry['death'], orpha_entry['prevalence'], orpha_entry['synonyms']])

def get_icd10_orpha(icd10number): 
    ic = icd10number.split(".")
    icd10 = ic[0]
    if len(ic) > 1:
        icd10 = ic[0] + "." + ic[1][0]
    return icd10

def rem_RDM(omim_number):
    return omim_number[0:6]


def link_up(gene_2_omim, omim_2_clinical, omim, omim_entry, orpha_entry, alias):
    if not alias in gene_2_omim: gene_2_omim[alias] = {'omim':set()}
    gene_2_omim[alias]['omim'].add(omim)
    if not omim in omim_2_clinical: 
        omim_2_clinical[omim] = omim_entry
        omim_2_clinical[omim]['found_genes'] = set()
        omim_2_clinical[omim]['orpha_links'] = dict()
    omim_2_clinical[omim]['found_genes'].add(alias)
    if not orpha_entry == "NA":
        omim_2_clinical[omim]['orpha_links'][orpha_entry['orpha']] = orpha_entry
    return gene_2_omim, omim_2_clinical

def link_with_orpha(mmc3_line, orpha_by_icd10, orpha_by_omim, alias_to_symbol, symbol_info, omim_dict, locus_dict, outputs, remove_unconfirmed):
    icd9_line = mmc3_line.strip().split("\t")
    icd10_numbers = set(map(get_icd10_orpha, icd9_line[2].strip('"').split(",")))
    cCHRorfXX = re.compile("C([0-9]+|X|Y)ORF([0-9]+)\Z")  ## check for Cxorfx
    ## according ot orpha, get list of omims for this ICD9 and their genes
    ### look for orpha to support the gene --> omim mappings
    # First, find the genes in orpha *with this ICD10* and map them to their orpha info
    orpha_gene_info = dict()  ### look up info about orpha entries with this ICD10
    for icd10 in icd10_numbers:
        for icd10_use in [icd10, icd10.split(".")[0]]:
            #icd10_use = icd10
            #if not icd10 in orpha_by_icd10 and icd10.split(".")[0] in orpha_by_icd10:
            #icd10_use = icd10.split(".")[0]
            if not icd10_use in orpha_by_icd10:
                continue
            for orpha_line in orpha_by_icd10[icd10_use]:
                gene_list = orpha_line['gene']
                # corresponding_synonyms = orpha_line['synonyms'].split(";")
                for g in gene_list: #range(len(gene_list)):
                    matches_orf = cCHRorfXX.match(g)
                    if matches_orf:
                        #print "match to " + matches_orf.group(0) 
                        g = "C" + matches_orf.group(1) + "orf" + matches_orf.group(2)
                    if not g in alias_to_symbol:
                        #print "UNFOUND_ORPHA_GENE " + g + " in " + icd9_line[0]
                        continue
                    for g_mainname in alias_to_symbol[g]:
                        ## if it's not already linked to an OMIM then forget it...
                        if True: # !!!! taking this out! 'omim' in symbol_info[g_mainname]:
                            if not g_mainname in orpha_gene_info:
                                orpha_gene_info[g_mainname] = dict()
                            orpha_to_add = orpha_line
                            orpha_to_add['ICD10_match'] = icd10_use
                            orpha_gene_info[g_mainname][orpha_to_add['orpha']] = orpha_to_add

    #return these 2 dicts
    gene_2_omim = dict()
    omim_2_clinical = dict()
    genes = set([ g.strip() for g in icd9_line[4].strip('"').split(",") ])
    omims = set()
    omims_in_same_orpha_entries = set()
    orpha_hits = dict()
    gene_to_link_up = []
    ## Next, for every gene in MMC3, use MMC3 ICD10 info to link it to orpha by gene name & omim
    for g in genes:
        if not g in alias_to_symbol:
            continue
        got_gene_omim = False
        ## 1) go through aliases, see if  they are in Orpha lines with this ICD10 code.  If the code matches, the gene/alias matches AND the 
        ## omim of the gene matches the omim of the orpha line, then use all this info together
        for alias in alias_to_symbol[g]:
            if not 'omim' in symbol_info[alias]:
                #print "Alias has no disease: gene " + g +  " -->main-name:" + alias
                continue
            ## using a double confirmation:  if the orpha & rzhetsky ICD10 codes agree on the same genes 
            # *and* OMIM and orpha agree on the same categories, it's probably good
            for omim in symbol_info[alias]['omim']:
                if alias in orpha_gene_info:
                    for orpha_entry in orpha_gene_info[alias].values():
                        ## check if the gene's omim is on the Orpha omims
                        #print "match res should be " + str(omim in orpha_entry['omim'])
                        omim_disorders_no_linkinfo = [dis[0] for dis in omim_dict[omim]['disorder']]
                        out_line = [g, alias, omim, ";".join(omim_disorders_no_linkinfo), orpha_entry['name']]
                        if omim in orpha_entry['omim']:
                            #print " res TRUE (should = " + str(omim in orpha_entry['omim'])
                            #print "orpha match is " + orpha_entry['name']
                            omims.add(omim)
                            omims_in_same_orpha_entries |= set(orpha_entry['omim']) # (orpha_entry['omim'] - Set(omim))
                            orpha_hits[orpha_entry['orpha']] = orpha_entry
                            outputs[0].write("\t".join(out_line) + "\n")
                            got_gene_omim = True
                            gene_2_omim, omim_2_clinical = link_up(gene_2_omim, omim_2_clinical, omim, omim_dict[omim], orpha_entry, alias)
                            to_write = [icd9_line[0], icd9_line[2], icd9_line[3], alias, 
                                        omim, ";".join(omim_disorders_no_linkinfo), 
                                        orpha_entry['ICD10_match'], orpha_entry['orpha'], 
                                        orpha_entry['name'],";".join(orpha_entry['gene']), 
                                        orpha_entry['death'], orpha_entry['prevalence'], orpha_entry['synonyms']]
                            outputs[2].write("\t".join(to_write) + "\n") # + print_orpha_omim_gene(alias, orpha_entry, omim, omim_dict) + "\n")
                        else:
                            outputs[1].write("\t".join(out_line) + "\n")
        ## experiment: NO omim for some
        if not got_gene_omim:
            for alias in alias_to_symbol[g]:
                if alias in orpha_gene_info:
                    #pdb.set_trace()
                    #if len(orpha_gene_info[alias]) > 1:
                    #    pdb.set_trace()
                    for orpha_entry in orpha_gene_info[alias].values():
                        #print orpha_entry['name']
                        out_line = [g, alias, 'NO_OMIM', "NO_!_OMIM", orpha_entry['name']]
                        orpha_hits[orpha_entry['orpha']] = orpha_entry
                        outputs[0].write("\t".join(out_line) + "\n")
                        got_gene_omim = True
                        gene_2_omim, omim_2_clinical = link_up(gene_2_omim, omim_2_clinical, 'NA', {}, orpha_entry, alias)
                        to_write = [icd9_line[0], icd9_line[2], icd9_line[3], alias, 
                                    'NA', 'NA', 
                                    orpha_entry['ICD10_match'], orpha_entry['orpha'], 
                                    orpha_entry['name'],";".join(orpha_entry['gene']), 
                                    orpha_entry['death'], orpha_entry['prevalence'], orpha_entry['synonyms']]
                        outputs[2].write("\t".join(to_write) + "\n") # + print_orpha_omim_gene(alias, orpha_entry, omim, omim_dict) + "\n")

        ## 2) Try to get OMIM for unmatched by ICD10.  IF can get it then try to get Orpha.  Assume that if there's only 1 OMIM for that 
        ## gene, and the gene is not an alias/only an alias for one thing, then htat's the right OMIM
        if not got_gene_omim and ((g in alias_to_symbol[g] and 'omim' in symbol_info[g]) or 
                                  (len(alias_to_symbol[g])==1 and 'omim' in symbol_info[alias_to_symbol[g][0]])):
            symbol_use = alias_to_symbol[g][0]
            if len(alias_to_symbol[g]) > 1: symbol_use = g
            omim_set = list(symbol_info[symbol_use]['omim'])
            filter_omim_not_locus = [om for om in omim_set if not om in locus_dict ]
            found_omim = ''
            if len(omim_set) == 1 and not remove_unconfirmed:
                omims.add(omim_set[0])
                found_omim = omim_set[0] 
                got_gene_omim = True
                outputs[3].write("DEFAULT " + g  + " -->OMIM:" + list(symbol_info[symbol_use]['omim'])[0] + '\n')
            elif len(filter_omim_not_locus) == 1 and not remove_unconfirmed:
                omims.add(filter_omim_not_locus[0])
                found_omim = filter_omim_not_locus[0] 
                got_gene_omim = True
                outputs[3].write( "DEFAULT_NONLOCUS " + g  + " -->OMIM:" + list(symbol_info[symbol_use]['omim'])[0] + '\n')
            # By default, this is most likely the real OMIM linked to the gene and the same one the mmc3_line wants
            # From OMIM, try to see if there's a matching entry in orpha:
            gene_omim_match = False
            if found_omim in orpha_by_omim:
                ## got an OMIM: get the Orpha's that have this omim and see if any have my gene = real match
                for orpha_matching_omim in orpha_by_omim[found_omim]:
                    for orpha_match_omim_gene in set(orpha_matching_omim['gene']):
                        if orpha_match_omim_gene in alias_to_symbol and symbol_use in alias_to_symbol[orpha_match_omim_gene]:
                            gene_omim_match = True
                            omims_in_same_orpha_entries |= set(orpha_matching_omim['omim']) # (orpha_match_omim_gene['omim'] - Set(filter_omim_not_locus))
                            orpha_hits[orpha_matching_omim['orpha']] = orpha_matching_omim
                            omim_disorders_no_linkinfo = [dis[0] for dis in omim_dict[found_omim]['disorder']]
                            to_write = [icd9_line[0], icd9_line[2], icd9_line[3], symbol_use, 
                                        found_omim, ";".join(omim_disorders_no_linkinfo), 
                                        ",".join(orpha_matching_omim['icd10']), orpha_matching_omim['orpha'], 
                                        orpha_matching_omim['name'],";".join(orpha_matching_omim['gene']), 
                                        orpha_matching_omim['death'], orpha_matching_omim['prevalence'], orpha_matching_omim['synonyms']]
                            outputs[2].write("\t".join(to_write) + "\n") # + print_orpha_omim_gene(alias, orpha_entry, omim, omim_dict) + "\n")
                            gene_2_omim, omim_2_clinical = link_up(gene_2_omim, omim_2_clinical, found_omim, omim_dict[found_omim], orpha_matching_omim, symbol_use)
            if got_gene_omim and not gene_omim_match:  ## this si the best we're going to get for this gene, so write that
                omim_disorders_no_linkinfo = [dis[0] for dis in omim_dict[omim]['disorder']]
                to_write = [icd9_line[0], icd9_line[2], icd9_line[3], symbol_use, 
                            found_omim, ";".join(omim_disorders_no_linkinfo)]
                to_write.extend(['NA']*7)
                outputs[2].write("\t".join(to_write) + "\n") # + print_orpha_omim_gene(alias, orpha_entry, omim, omim_dict) + "\n")
                gene_2_omim, omim_2_clinical = link_up(gene_2_omim, omim_2_clinical, found_omim, omim_dict[found_omim], 'NA', symbol_use)
        if not got_gene_omim:
            gene_to_link_up.append(g)
    #[ om for om in omims if len(set(omim_dict[om]['genes']) -  set(gene_2_omim.keys()) ) > 0]
    # As with the "Default" ones Giving up on Orpha for these, but still want them in our list
    for missing_gene in gene_to_link_up:
        symbol_use = missing_gene
        if len(alias_to_symbol[missing_gene]) == 1:
            symbol_use = alias_to_symbol[missing_gene][0]
        if not missing_gene in alias_to_symbol[missing_gene] or not 'omim' in symbol_info[symbol_use]:
            outputs[3].write( "GENE_UNLINKED_NO_OMIM: gene" + missing_gene + " symbol_use=" + symbol_use  + " for " + ",".join(icd10_numbers) + '\n')
            continue
        ## intersect the SYMBOL's linked omims with our CANDIDATE omims
        
        intersection = set(symbol_info[symbol_use]['omim']) & (omims_in_same_orpha_entries | omims)
        if len(intersection) > 0:
            omims |= intersection
            for found_omim in intersection:
                for orpha_matching_omim in orpha_by_omim[found_omim]:
                    for orpha_match_omim_gene in orpha_matching_omim['gene']:
                        if orpha_match_omim_gene in alias_to_symbol and symbol_use in alias_to_symbol[orpha_match_omim_gene]:
                            gene_omim_match = True
                            #omims_in_same_orpha_entries |= Set(orpha_matching_omim['omim']) # (orpha_match_omim_gene['omim'] - Set(filter_omim_not_locus))
                            orpha_hits[orpha_matching_omim['orpha']] = orpha_matching_omim
                            omim_disorders_no_linkinfo = [dis[0] for dis in omim_dict[found_omim]['disorder']]
                            to_write = [icd9_line[0], icd9_line[2],icd9_line[3], symbol_use, 
                                        found_omim, ";".join(omim_disorders_no_linkinfo), 
                                        ",".join(orpha_matching_omim['icd10']), orpha_matching_omim['orpha'], 
                                        orpha_matching_omim['name'],";".join(orpha_matching_omim['gene']), 
                                        orpha_matching_omim['death'], orpha_matching_omim['prevalence'], orpha_matching_omim['synonyms']]
                            gene_2_omim, omim_2_clinical = link_up(gene_2_omim, omim_2_clinical, found_omim, omim_dict[found_omim], orpha_matching_omim, symbol_use)
                            outputs[2].write("\t".join(to_write) + "\n") # + print_orpha_omim_gene(alias, orpha_entry, omim, omim_dict) + "\n")
        else:
            outputs[3].write("GENE_NO_OMIM_MATCH: gene" + missing_gene + " symbol_use=" + symbol_use  + " for " + ",".join(icd10_numbers) + " omims = " + ",".join(symbol_info[symbol_use]['omim']) + '\n')
    return gene_2_omim, omim_2_clinical

def load_n_link(morbidmap, orpha_parsed, ncbi_genes, mmc3, conservative_setting=False):
    omim_dict, symbol_info, alias_to_symbol, locus_dict = annotation_work.omim_annotation(morbidmap,ncbi_genes)

    orpha_info = open(orpha_parsed)
    orpha_by_icd10 = dict()
    orpha_by_omim = dict()
    for orpha_disease in orpha_info:
        orpha = orpha_disease.strip().split("\t")
        if orpha[0] == "orpha":  # first line
            continue
        genes = orpha[8].split(";")
        icd10_list = map(get_icd10_orpha, orpha[2].split(","))
        omim_list = map(rem_RDM, orpha[1].split(","))
        line_info = {'orpha':orpha[0], 'omim':omim_list, 'name':orpha[3], 'name2':orpha[4], 'prevalence':orpha[5], 'death':orpha[7], 'onset':orpha[4], 'gene':genes, 'synonyms':orpha[9], 'icd10':icd10_list}
        for icd10_code in icd10_list:
            if not icd10_code in orpha_by_icd10:
                orpha_by_icd10[icd10_code] = []
            orpha_by_icd10[icd10_code].append(line_info)
        for omim in omim_list:
            if not omim in orpha_by_omim:
                orpha_by_omim[omim] = []
            orpha_by_omim[omim].append(line_info)

    accepted = open('data_processed/log_link_accepted.txt','w')  ## ICD10 & gene & omim match of Orpha, Omim, MMC3
    rejected = open('data_processed/log_link_rejected.txt','w')  ## cases where ICD10 matches, gene matches but omims do not match.  One line per Gene & OMIM
    loglog = open('data_processed/log_link','w')
    output = open('data_processed/supp_tab_1.txt','w')
    output.write("MMC3_Summary_Name\tMMC3_ICD10\tMMC3_Assoc_Disorder\talias\tOMIM\tOMIM_disorder\tICD10\t")
    output.write("orpha\torpha_disorder\torpha_genes\tdeath\tprevalence\torpha_synonyms\n")
    mend_icd9 = open(mmc3)
    output_handles = [accepted, rejected, output, loglog]
    # output_handles = [sys.stdout]*3
    icd9_genetic_clinical = dict()
    for mmc3_line in mend_icd9:
        if mmc3_line.startswith('Summary Name'):
            continue       
        icd9_gene_2_omim, icd9_omim_2_clinical = link_with_orpha(mmc3_line, orpha_by_icd10, orpha_by_omim, alias_to_symbol, symbol_info, omim_dict, locus_dict, output_handles, conservative_setting)
        icd9_genetic_clinical[mmc3_line.strip().split("\t")[0]] = {'gene_omim':icd9_gene_2_omim, 'omim_clinical':icd9_omim_2_clinical}
    for h in output_handles:
        h.close()
    # get all genes with a mendelian disease...
    mendelian_genes = set()
    for gene in symbol_info:
        if 'omim' in symbol_info[gene]:
            mendelian_genes.add(gene)
    return icd9_genetic_clinical, mendelian_genes





