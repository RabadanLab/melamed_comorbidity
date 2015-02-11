cut -f 1 mmc3_tableS3_icd9_omim_names.txt  | sort -u > z
cut -f 2 mmc4_tableS4_comorbidity.txt | sort -u > x

	Bartters Syndrome
Bartter Syndrome
Congenital Hirschsprung Disease
	Congenital Hirschsprungs Disease
	Downs Syndrome
Down Syndrome
	Edwards Syndrome
Edward Syndrome
	Fragile X Syndrome
Fragile-X Syndrome
Friedreich Ataxia
	Friedreichs Ataxia
	Klinefelters Syndrome
Klinefelter Syndrome
Li Fraumeni and Related Syndromes
	Li-Fraumeni and Related Syndromes
Mendelian Disease
	Pataus Syndrome
Patau Syndrome
	Spinocerebellar ataxia
Spinocerebellar Ataxia
	Summary Name
	Turners Syndrome
Turner_Syndrome


cut -f 1 mmc3_tableS3_icd9_omim_names.txt  | sort -k1,1 -t $'\t' > mmc31
cut -f 1 Mendelian_Genes_fixedRM.txt | sort -k1,1 -t $'\t'  > mgrm1 ## then remove "Anomalies of Musculoskelatal System"
paste mmc31 mgrm1 > tmpcomb
 awk -F"\t" '$1!=$2' tmpcomb # verify
sed '1d' mmc3_tableS3_icd9_omim_names.txt | sort -k1,1 -t $'\t' > x
join -t $'\t' tmpcomb x > y
head -n 1 mmc3_tableS3_icd9_omim_names.txt > res
cut -f 2- y >> res
mv res icd9_mendelian_name.txt

head -n 1 Homo_sapiens.gene_info | tr " " "\t" | cut -f 3-16 > Homo_sapiens.gene_info.protein_coding
awk -F"\t" '$3=="TERC" || $3=="ATXN8OS" || $3=="ATXN8OS" || $3=="RMRP" || $3=="RNU4ATAC" || $3=="H19"|| $10=="protein-coding" || $10=="unknown" || $10=="other"|| $10=="tRNA"' Homo_sapiens.gene_info | cut -f 2- | awk 'BEGIN{FS=OFS="\t"}{if($10!=$2){if($4=="-"){$4==$10}else{$4=$4"|"$10};}print $0}' >> Homo_sapiens.gene_info.protein_coding

alias python=/ifs/home/c2b2/rr_lab/siz2102/apps/anaconda/bin/python
cut -f 5 ../from_rzhetsky/icd9_mendelian_name.txt | sed 's/"//g' | tr "," "\n" | sort -u > rg

cut -f 4 ../from_rzhetsky/icd9_mendelian_name.txt | sed 's/"//g' | tr ";" "\n" | sort -u > rm
cut -f 1 morbidmap_disorders.txt | sort > mm

awk '{print "ACC\t"$0}' accepted_links.txt > acc_rej
awk '{print "REJ\t"$0}' rejected_links.txt >> acc_rej

cut -f 1-3 rejected_links.txt | sort -u > rejected_gene_omim
cut -f 1-3 accepted_links.txt | sort -u > accepted_gene_omim
comm -23 rejected_gene_omim accepted_gene_omim > rejonly


grep Q79.8 problems | cut -d" " -f 2 | sed 's/gene//g'


    0.0016