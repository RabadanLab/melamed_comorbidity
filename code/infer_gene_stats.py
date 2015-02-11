from collections import Counter 
UNIFORM_PREV = .0001
def prevalence_number(prevalence_code):
    prev = 0
    if prevalence_code == ">1 / 1 000":
        prev = float(1)/1000
    elif prevalence_code == "1-5 / 10 000":
        prev = float(5)/10000
    elif prevalence_code == "1-9 / 100 000":
        prev = float(9)/100000
    elif prevalence_code == "1-9 / 1 000 000":
        prev = float(9)/1000000
    elif prevalence_code == "<1 / 1 000 000":
        prev = float(1)/1000000
    return prev


def average_survival(survival_code):
    survival = 100
    if survival_code == "Adult":
        survival = 45
    elif survival_code == "Child / adolescent":
        survival = 16
    elif survival_code == "Before age 5":
        survival = 4
    elif survival_code == "Young adult":
        survival = 30
    return survival

def infer_omim_clinical(code_info, considered_age=18, nm='NO_WRITE'):
    # for each omim entry, go throug the linked clinical entries and infer freq, death
    rep = open('data_processed/log_report.txt','a')

    ### omim is the basic unit of disease... GENE stats are derived from OMIM stats
    for omim in code_info['omim_clinical']:
        # average the prevalence of hte non-too-young sub-diagnoses
        prevalences = []
        num_links = 0
        est_age = 0
        young_prevalences = []
        young_num_links = 0
        young_est_age = 0
        for orpha_link in code_info['omim_clinical'][omim]['orpha_links'].values():
            orpha_age = average_survival(orpha_link['death'])
            if orpha_age > considered_age:
                est_age = est_age + orpha_age
                prevalences.append(prevalence_number(orpha_link['prevalence']))
                num_links += 1
            else:
                young_est_age = young_est_age + orpha_age
                young_prevalences.append(prevalence_number(orpha_link['prevalence']))
                young_num_links += 1
        est_prevalence = 0
        if num_links > 0:
            est_prevalence = max(prevalences)
            est_age = float(est_age)/num_links
        elif young_num_links > 0:  # if the only assoc are to young-death diseases...
            est_prevalence = max(young_prevalences)
            est_age = float(young_est_age)/young_num_links
        else:
            # else = there are NO links, of any age! --> assume this is NOT lethal
            est_age = 100
        code_info['omim_clinical'][omim]['est_age'] = est_age
        code_info['omim_clinical'][omim]['est_prevalence'] = est_prevalence
    # for each gene go through the linked OMIMs and infer from there the 
    #   - age, 
    #   - # genes co-linked here
    #   - frequency
    prevalence_list = []
    for gene in code_info['gene_omim']:
        # average the prevalence of hte non-too-young sub-diagnoses
        # again if the only diagnosis is a young-death one then use that
        prevalence = []
        num_links = 0
        est_age = 0
        num_linked_genes = 0
        young_prevalence = []
        young_num_links = 0
        young_est_age = 0
        young_num_linked_genes = 0
        for omim in code_info['gene_omim'][gene]['omim']:
            survival = code_info['omim_clinical'][omim]['est_age']
            #print "GENE: " + gene  + " omim:"+omim + " age:" +str( survival) + " prev:" + str(code_info['omim_clinical'][omim]['est_prevalence'] ) + " ng:" + str(len(code_info['omim_clinical'][omim]['found_genes']))
            if survival > considered_age:
                prevalence.append(code_info['omim_clinical'][omim]['est_prevalence'])
                est_age += survival
                num_linked_genes += len(code_info['omim_clinical'][omim]['found_genes'])
                num_links += 1
            else:
                young_prevalence.append(code_info['omim_clinical'][omim]['est_prevalence'])
                young_est_age += survival
                young_num_linked_genes += len(code_info['omim_clinical'][omim]['found_genes'])
                young_num_links += 1
        est_prevalence = 0
        if num_links > 0:
            est_prevalence = max(prevalence)
            est_age = float(est_age)/num_links
            num_linked_genes = float(num_linked_genes)/num_links
        elif young_num_links > 0:  # if the only assoc are to young-death diseases...
            est_prevalence = max(young_prevalence)
            est_age = float(young_est_age)/young_num_links
            num_linked_genes = float(young_num_linked_genes)/young_num_links
        code_info['gene_omim'][gene]['est_age'] = est_age
        code_info['gene_omim'][gene]['est_prevalence'] = est_prevalence
        code_info['gene_omim'][gene]['num_genes'] = num_linked_genes
        if est_prevalence > 0:
            prevalence_list.append(est_prevalence)
        if not nm == 'NO_WRITE':
            rep.write('\t'.join([nm, gene, str(est_prevalence), str(est_age), str(num_linked_genes), ','.join(code_info['gene_omim'][gene]['omim'])])+ "\n")
    rep.close()
    code_info['mode_mendelian_prevalence'] = 0
    code_info['max_mendelian_prevalence'] = 0
    if len(prevalence_list) > 0:
        counted_prev = Counter(prevalence_list)
        code_info['mode_mendelian_prevalence'] = counted_prev.most_common(1)[0][0]
        code_info['max_mendelian_prevalence'] = max(prevalence_list)
        import pdb
        #pdb.set_trace()
    return code_info


def infer_p_I(code_info, considered_age, default_max, prev_mode='orpha'):
    code_info['p_I'] = 0
    default_prev_use = code_info['mode_mendelian_prevalence']
    if default_prev_use == 0:
        default_prev_use = default_max
    for gene in code_info['gene_omim']:
        if code_info['gene_omim'][gene]['est_age'] > considered_age:
            prev_use = default_prev_use #code_info['mode_mendelian_prevalence']
            if code_info['gene_omim'][gene]['est_prevalence'] > 0:
                prev_use = code_info['gene_omim'][gene]['est_prevalence']
            if prev_mode == 'uniform':
                prev_use = UNIFORM_PREV
            gene_contribution = 1/code_info['gene_omim'][gene]['num_genes']*prev_use
            code_info['p_I'] += gene_contribution
    return code_info

def infer_gene_score(code_info, gene, relative_risk, low_age, default_max, prev_mode='orpha', clinical_code='NO_WRITE'):
    import pdb
    #pdb.set_trace()
    if code_info['gene_omim'][gene]['est_age'] <= low_age:
        return 0

    # p(M | I) = p(I | M)*p(M)/p(I)
    # p(I) = sum of prevalences of the mendelian diseases, where *Code-wide MODE* is substituted for missing vals
    # p(I | M) = 1 if age is > low_age (if not, return 0, above)
    prev_use = code_info['mode_mendelian_prevalence']
    if code_info['gene_omim'][gene]['est_prevalence'] > 0:
        prev_use = code_info['gene_omim'][gene]['est_prevalence']
    if prev_use == 0:
        prev_use = default_max
    if prev_mode == 'uniform':
        prev_use = UNIFORM_PREV
    p_g_given_m = 1/code_info['gene_omim'][gene]['num_genes']
    p_g_given_i = p_g_given_m*prev_use/code_info['p_I']

    ## then, p(G | C) = sum(relevant I_j) P(G | I_j)*P(I_j | C)
    #     p(I_j | C ) = P(C | I_j)/P(C) * P(I_j)
    #                 = (Relative risk ) * max-prevalence
    # *note* different definition of P(I_j) from above: compares P(I_j) to other P(I)
    max_prevalence = code_info['max_mendelian_prevalence']
    if max_prevalence == 0:
        max_prevalence = default_max
    if prev_mode == 'uniform':
        max_prevalence = UNIFORM_PREV
    p_i_given_c = relative_risk*max_prevalence
    if not clinical_code == 'NO_WRITE':
        score_out = open('data_processed/logs_' + prev_mode + '.gene_scores.txt','a')
        score_out.write('\t'.join([clinical_code, gene, str(p_g_given_i*p_i_given_c), str(relative_risk), str(code_info['p_I']), str(code_info['gene_omim'][gene]['num_genes']), str(p_g_given_i), str(p_i_given_c)]) + "\n")
        score_out.close()
    return p_g_given_i*p_i_given_c
