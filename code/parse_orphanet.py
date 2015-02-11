import xml.etree.ElementTree as ET

tree = ET.parse('en_product1.xml')   ### 6794 diseases
root = tree.getroot()

## parse diseases and get omim links
ddict = dict()
for disease in root[0]:    
    name = disease.findall("./Name")[0].text
    orpha = disease.findall("./OrphaNumber")[0].text
    references = disease.findall("./ExternalReferenceList/ExternalReference")
    omim_id=[]
    icd10=[]
    for refi in range(len(references)):
        if references[refi].findall("./Source")[0].text == "OMIM":
            omim_id.append(references[refi].findall("./Reference")[0].text)
        if references[refi].findall("./Source")[0].text == "ICD10":
            icd10.append(references[refi].findall("./Reference")[0].text)
    synonym_list = disease.findall("./SynonymList/Synonym")
    synonyms = []
    for s in range(len(synonym_list)):
        synonyms.append(synonym_list[s].text)
    # print name + " " + orpha + " " + omim_id + " " + ",".join(synonyms)
    ddict[orpha] = {'name':name, 'omim':",".join(omim_id), 'ICD10':",".join(icd10), 'synonym':";".join(synonyms)}

## for each disease, get prevalences etc
tree = ET.parse('en_product2.xml')  ## again 6794
root = tree.getroot()
for disease in root[0]:    
    name = disease.findall("./Name")[0].text
    orpha = disease.findall("./OrphaNumber")[0].text
    print orpha
    x = disease.findall("./ClassOfPrevalence")[0]
    if len(x) > 0:
        ddict[orpha]['prevalence'] = x[0].text
    x = disease.findall("./AverageAgeOfOnset")[0]
    if len(x) > 0:
        ddict[orpha]['onset'] = x[0].text
    x = disease.findall("./AverageAgeOfDeath")[0]
    if len(x) > 0:
        ddict[orpha]['death'] = x[0].text

## get genes
tree = ET.parse('en_product6.xml')  ## again 6794
root = tree.getroot()
for disease in root[0]:    
    orpha = disease.findall("./OrphaNumber")[0].text
    genes = disease.findall("./DisorderGeneAssociationList/DisorderGeneAssociation/Gene")
    gres = []
    sres = []
    for g in range(len(genes)):
        gres.append(genes[g].findall("./Symbol")[0].text)
        synonyms = []
        syn = genes[g].findall("./SynonymList/Synonym")
        for s in range(len(syn)):
            synonyms.append(syn[s].text)
        sres.append(",".join(synonyms))
    ddict[orpha]['gene'] = ";".join(gres)
    ddict[orpha]['gene_synonym'] = ";".join(sres)

fields = ["orpha","omim","ICD10","name","synonym","prevalence","onset","death","gene","gene_synonym"]
out_tab = open('orpha_omim.txt','w')
out_tab.write("\t".join(fields) + "\n")
i = 1
import re
for k in ddict:
    out_tab.write(k + "\t") 
    #print "####### " +  k
    i = i + 1
    for f in range(1,len(fields)):
        field = "NA"
        if fields[f] in ddict[k].keys():
            if not ddict[k][fields[f]] == None and not ddict[k][fields[f]] == "":
                field=ddict[k][fields[f]]
                field = re.sub(u"[^\x20-\x7f]+",u"",field)
        #print fields[f] + " " + field
        out_tab.write(field + "\t")
    out_tab.write("\n")
out_tab.close()

'''
query = "79399"
myd = []
for disease in root[0]:    
    name = disease.findall("./Name")[0].text
    orpha = disease.findall("./OrphaNumber")[0].text
    if orpha == query:
        myd = disease
        break

    ########
i= 1
for disorder in root[0]:
    print "########"
    print disorder.tag, disorder.attrib
    for disorder_info in disorder:
        print disorder_info.tag; disorder_info.attrib
    i = i + 1
    print 
    if i > 5:
        print i
        break
    
for i in range(len(root[0][1])):
    print "i is " + str(i) + " level:" + root[0][1][i].tag
    for j in range(len(root[0][1][i])):
        print "j is " + str(j) 
        for k in root[0][1][i].attrib.keys():
            print k + " " + root[0][1][i].attrib[k]
    if root[0][1][i].tag == "ExternalReferenceList":
    for k in root[0][1][i].attrib.keys():
        print k + " " + root[0][1][i].attrib[k]            



#tree = ET.parse('en_product1.xml')
root = tree.getroot()

for d in range(len(root[0][1])):

    print "i is " + str(i) + " level:" + root[0][1][i].tag  + " text:" + root[0][1][i].text
    disease = root[0][1]    
'''
