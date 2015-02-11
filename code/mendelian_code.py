import link_rzhetsky_orpha_omim
import pickle
from itertools import chain
import pdb

def load_associations_annotations(conservative=False): #save_file = 'mendelian_code_associations.pkl', 
    morbidmap="data_source/OMIM/morbidmap.txt"
    orpha_parsed='data_source/orph2'
    ncbi_genes= "data_source/Homo_sapiens.gene_info.protein_coding"
    mmc3= "data_source/mmc3_2"
    comorbidity_file = "data_source/comorbidities"
    cancer_tab = 'data_source/cancer_table.txt'

    icd_gene_clinical, mendelian_genes = link_rzhetsky_orpha_omim.load_n_link(morbidmap, orpha_parsed, ncbi_genes, mmc3, conservative)
    cancer_info= disease_info(cancer_tab)
    icd_gene_clinical, gene_cancer_assoc = clinical_link(icd_gene_clinical, cancer_info, comorbidity_file)
    '''
    pkl = open(save_file,'w')
    pickle.dump((icd_gene_clinical, cancer_info, mendelian_genes),pkl)
    pkl.close()
    '''
    return icd_gene_clinical, cancer_info, mendelian_genes

def open_associations_annotations():
    pkl = open('data_processed/mendelian_code_associations.pkl')
    icd_gene_clinical, cancer_info, mendelian_genes = pickle.load(pkl)
    return icd_gene_clinical, cancer_info, mendelian_genes

def open_associations_nogermline(pkl, cens_germline):
    icd_gene_clinical, cancer_info, mendelian_genes, alterations, gn = pickle.load(open(pkl))
    for icd in icd_gene_clinical:
        new_genes = dict()
        for g in icd_gene_clinical[icd]['gene_omim']:
            if not g in cens_germline:
                new_genes[g] = icd_gene_clinical[icd]['gene_omim'][g]
        icd_gene_clinical[icd]['gene_omim'] = new_genes
    mendelian_genes = list(set(mendelian_genes) - set(cens_germline))
    
    return icd_gene_clinical, cancer_info, mendelian_genes, cens_germline

def get_mendelian_genes(icd_gene_clinical):
    all_mend_gns = list(set(chain.from_iterable([icd_gene_clinical[icd]['gene_omim'].keys() 
                                                 for icd in icd_gene_clinical])))
    all_mend_gns.sort()
    return all_mend_gns

def get_mendelian_diseases(icd_gene_clinical):
    all_md = icd_gene_clinical.keys()
    all_md.sort()
    return all_md

def disease_info(info_file):
    infos = open(info_file)
    cancer_info = dict()
    header = infos.readline().strip().split("\t")
    for cancer in infos:
        cancer = cancer.strip().split("\t")
        cancer_info[cancer[0]] = {'control_incidence':float(cancer[1]), 'young':int(cancer[2]),
                                  'TCGA':[] if len(cancer) < 4 else cancer[3].split(',')}
    return cancer_info


def clinical_link(icd_gene_clinical, cancer_info, comorbidity_file):
    cancer_assoc=dict()
    for c in cancer_info:
        cancer_assoc[c] = 0

    gene_cancer_assoc = dict()
    for icd in icd_gene_clinical:
        icd_gene_clinical[icd]['cancer_assoc'] = cancer_assoc.copy()
        for gene in icd_gene_clinical[icd]['gene_omim']:
            if not gene in gene_cancer_assoc:
                gene_cancer_assoc[gene] = {'uniform':cancer_assoc.copy(), 'orpha':cancer_assoc.copy()}

    #icd_gene_clinical['Anophthalmos/Micropthalmos']['cancer_assoc']['Melanoma']
    import infer_gene_stats
    from collections import Counter         
    for cancer in cancer_info:
        #cancer = "Melanoma"
        min_age = cancer_info[cancer]['young']
        comorbidities = open(comorbidity_file)
        comorbid_findings = []
        for mendelian_line in comorbidities:
            mendelian = mendelian_line.strip().split("\t")
            if not mendelian[0] == cancer or not mendelian[1] in icd_gene_clinical or float(mendelian[2]) < 1:
                continue
            comorbid_findings.append([mendelian[1], float(mendelian[2])])
            icd_gene_clinical[mendelian[1]]['cancer_assoc'][cancer] = float(mendelian[2])
            #print 'assigning ' + mendelian[1] + ' assoc with ' + cancer + ' as '+ mendelian[2]
        #cancer_info[cancer]['mendelian_code'] = comorbid_findings
        code_genes = set()
        rep = open('data_processed/log_report.txt','w')
        rep.write('Diagnosis\tGene\tPrev\tSurvival\tNum Gene\tOMIMs\n')
        rep.close()
        max_prevalences = []
        for mend_disease in range(len(comorbid_findings)):
            mend_name = comorbid_findings[mend_disease][0]
            icd_gene_clinical[mend_name] = infer_gene_stats.infer_omim_clinical(icd_gene_clinical[mend_name], min_age, mend_name)
            code_genes |= set(icd_gene_clinical[mend_name]['gene_omim'].keys())
            if not icd_gene_clinical[mend_name]['max_mendelian_prevalence'] == 0:
                max_prevalences.append(icd_gene_clinical[mend_name]['max_mendelian_prevalence'])
        for prev_mode in ['orpha','uniform']:
            #prev_mode = 'orpha'# 'uniform' # 'orpha'
            counted_prev = Counter(max_prevalences)
            default_prevalence = counted_prev.most_common(1)[0][0]
            for mend_disease in range(len(comorbid_findings)):
                mend_name = comorbid_findings[mend_disease][0]
                icd_gene_clinical[mend_name] = infer_gene_stats.infer_p_I(icd_gene_clinical[mend_name], min_age, default_prevalence, prev_mode)
            #score_out = open(cancer.replace("/","-") + '.' + prev_mode + '.gene_scores.txt','w')
            #score_out.write('Clinical_code\tGene\tP\tRel_risk\tp(I)\tnum_gene\tp(G | I)\tP(I | C)\n')
            #score_out.close()
            score_results = dict()
            #score2 = open(prev_mode + '.final_score.txt','w')
            #score2.write('Gene\tScore\n')
            for mend_gene in code_genes:
                score = 0
                for c_ind in range(len(comorbid_findings)):
                    comorbidity = comorbid_findings[c_ind]
                    if mend_gene in icd_gene_clinical[comorbidity[0]]['gene_omim']:
                        score += infer_gene_stats.infer_gene_score(icd_gene_clinical[comorbidity[0]], mend_gene, comorbidity[1], min_age, default_prevalence, prev_mode, comorbidity[0])
                score_results = score
                #score2.write('\t'.join([mend_gene, str(score)])+'\n')
                gene_cancer_assoc[mend_gene][prev_mode][cancer] = score
            #score2.close()
    return icd_gene_clinical, gene_cancer_assoc
