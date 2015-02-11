# in gistic2 dir
for cn in `ls -d *0.0`; do
    cname=`echo $cn | cut  -d "_" -f 2 | cut -d "-" -f 1`
    if [ ! -e ../analyses__2013_09_23/$cname/20130923 ]; then mkdir -p ../analyses__2013_09_23/$cname/20130923; fi
    
    mv $cn ../analyses__2013_09_23/$cname/20130923
done

awk -F"\t" '{print $3}' ~/mendelian/complex_to_TCGA.txt  | tr "," "\n" | sort -u | sed '/^$/d' | grep -v Mutation > mut_data
awk -F"\t" '{print $2}' ~/mendelian/complex_to_TCGA.txt  | tr "," "\n" | sort -u | sed '/^$/d' | grep -v Copy > cn_data
./firehose_get -tasks mutsig  analyses latest `cat mut_data | tr "\n" " "`

cn_suffix=CopyNumber_Gistic2.Level_4.2013092300.0.0
mut_suffix="-T*.MutSigNozzleReportMerged.Level_4.2013092300.0.0"
mut_prefix="gdac.broadinstitute.org_"
# etc
rootdir=`pwd`
while read mut; do
    echo $mut_list
    cd analyses__2013_09_23/$mut/20130923
    gunzip $mut_prefix$mut$mut_suffix.tar.gz
    tar -xvf $mut_prefix$mut$mut_suffix.tar
    cd $rootdir
done < mut_data
    
output_dir=~/mendelian/mut_info
while read mut; do
    dirwrite=`grep $mut ~/mendelian/cancer_info.txt | cut -f 1`
    mkdir -p $output_dir/$dirwrite
    #awk -F "\t" '{if($20 < .25){print $2}}' analyses__2013_09_23/$mut/20130923/$mut_prefix$mut$mut_suffix/$mut*.sig_genes.txt | sort > $output_dir/$dirwrite/$mut.mutations
    ln -s `pwd`/analyses__2013_09_23/$mut/20130923/$mut_prefix$mut$mut_suffix/$mut*.sig_genes.txt $output_dir/$dirwrite/$mut.mutsig
done < mut_data

cn_suffix="-T*.CopyNumber_Gistic2.Level_4.2013092300.0.0"
cn_prefix="gdac.broadinstitute.org_"
while read cn; do
    dirwrite=`grep $cn ~/mendelian/cancer_info.txt | cut -f 1`
    mkdir -p $output_dir/$dirwrite
    echo $cn to $dirwrite
    gistic_res=analyses__2013_09_23/$cn/20130923/$cn_prefix$cn$cn_suffix

    deletions=$gistic_res/table_del.conf_99.txt
    rm $output_dir/$dirwrite/$cn.deletions
    ln -s `pwd`/$deletions $output_dir/$dirwrite/$cn.deletions
    #cut -f 10 $deletions | sed '1d' | tr "," "\n" | sort -u | sed '/^$/d' >  $output_dir/$dirwrite/$cn.deleted_region
    #cut -f 12 $deletions | sed '1d' | tr "," "\n" | sort -u | sed '/^$/d' >  $output_dir/$dirwrite/$cn.deleted_peak
    #cut -f 15 $deletions | sed '1d' | tr "," "\n" | sort -u | sed '/^$/d' >  $output_dir/$dirwrite/$cn.deleted_top

    amplification=$gistic_res/table_amp.conf_99.txt
    rm $output_dir/$dirwrite/$cn.amplifications
    ln -s `pwd`/$amplification $output_dir/$dirwrite/$cn.amplifications
    #cut -f 10 $amplification | sed '1d' | tr "," "\n" | sort -u | sed '/^$/d' >  $output_dir/$dirwrite/$cn.amplified_region
    #cut -f 12 $amplification | sed '1d' | tr "," "\n" | sort -u | sed '/^$/d' >  $output_dir/$dirwrite/$cn.amplified_peak
    #cut -f 15 $amplification | sed '1d' | tr "," "\n" | sort -u | sed '/^$/d' >  $output_dir/$dirwrite/$cn.amplified_top
    rm $output_dir/$dirwrite/$cn.genes
    ln -s `pwd`/$gistic_res/all_data_by_genes.txt $output_dir/$dirwrite/$cn.genes
done < cn_data


comm -12 mend_genes cn_genes | sort > to_test
comm -12 to_test mut_info/Melanoma/SKCM.deleted_peak | sort > delmend
comm -12 to_test mut_info/Melanoma/SKCM.amplified_peak | sort > ampmend
##################
#### all cancers!
./firehose_get -tasks mutsig analyses 2013_09_23 `cat mut_data | awk '{print "~"$0}' | tr "\n" " "`
./firehose_get -tasks mutsig gistic analyses 2013_09_23 OV  PAAD  PCPG  READ  SARC  SKCM  THCA

cn_suffix=CopyNumber_Gistic2.Level_4.2013092300.0.0
mut_suffix="-T*.MutSigNozzleReportMerged.Level_4.2013092300.0.0"
mut_prefix="gdac.broadinstitute.org_"
# etc
rootdir=`pwd`
for ziptared in `ls -d analyses__2013_09_23/*/20130923/$mut_prefix*{$mut_suffix,$cn_suffix}.tar.gz`; do
    dirn=`dirname $ziptared`
    bn=`basename $ziptared .tar.gz`
    if [ ! -e $dirn/$bn ]; then
        echo Look for $dirn/$bn
        echo gunzip $ziptared
        echo tar -xvf $dirn/$bn.tar -C $dirn
    fi
done

del_table=table_del.conf_99.txt
amp_table=table_amp.conf_99.txt
cn_genes=all_data_by_genes.txt
for cntype in {$del_table,$amp_table}; do
    for cn_res in `ls analyses__2013_09_23/*/20130923/$mut_prefix*$cn_suffix/$cntype`; do
        #cut -f 12 $cn_res | sed '1d' | tr "," "\n" | sort -u | sed '/^$/d' >>  $cntype
        #cut -f 1 `dirname $cn_res`/$cn_genes | sed '1d' >>  $cn_genes
        cancer_name=`echo $cn_res | cut -d "/" -f 2`
        ln -s `pwd`/$cn_res ~/mendelian/cancer_alterations/$cancer_name.$cntype
        ln -s `pwd`/$cn_genes ~/mendelian/cancer_alterations/$cancer_name.$cn_genes
    done
    sort -u $cntype > x
    mv x $cntype
done
#sort -u $cn_genes > x
#mv x $cn_genes

mut_table=sig_genes.txt
mut_genes=mut_genes
for mut_res in `ls analyses__2013_09_23/*/20130923/$mut_prefix*$mut_suffix/*.$mut_table`; do
    #awk -F "\t" '{if($20 < .25){print $2}}' $mut_res  | sort >> $mut_table
    #cut -f 2 $mut_res | sed '1d' >> $mut_genes
    cancer_name=`echo $mut_res | cut -d "/" -f 2`
    ln -s `pwd`/$mut_res ~/mendelian/cancer_alterations/$cancer_name.mut
done
sort -u $mut_table > x
mv x $mut_table
sort -u $mut_genes > x
mv x $mut_genes

ln -s `pwd`/sig_genes.txt ~/mendelian/mut_info/mutation_all_cancer
ln -s `pwd`/table_amp.conf_99.txt  ~/mendelian/mut_info/amp_all_cancer
ln -s `pwd`/table_del.conf_99.txt  ~/mendelian/mut_info/del_all_cancer
ln -s `pwd`/$cn_genes ~/mendelian/mut_info/cn_genes_all
ln -s `pwd`/$mut_genes ~/mendelian/mut_info/mut_genes_all


#################
cn_suffix=CopyNumber_Gistic2.Level_4.2013092300.0.0
mut_suffix="-T*.MutSigNozzleReportMerged.Level_4.2013092300.0.0"
mut_prefix="gdac.broadinstitute.org_"

del_table=table_del.conf_99.txt
amp_table=table_amp.conf_99.txt
cn_genes=all_data_by_genes.txt
cn_thres=all_thresholded.by_genes.txt
for cntype in {$del_table,$amp_table}; do
    for cn_res in `ls analyses__2013_09_23/*/20130923/$mut_prefix*$cn_suffix/$cntype`; do
        cancer_name=`echo $cn_res | cut -d "/" -f 2`
        f=~/mendelian/cancer_alterations/$cancer_name.$cntype
        if [ ! -e $f ]; then ln -s `pwd`/$cn_res $f; fi
        f=~/mendelian/cancer_alterations/$cancer_name.$cn_genes
        if [ ! -e $f ]; then ln -s `pwd`/`dirname $cn_res`/$cn_genes $f; fi
        f=~/mendelian/cancer_alterations/$cancer_name.$cn_thres
        if [ ! -e $f ]; then ln -s `pwd`/`dirname $cn_res`/$cn_thres $f; fi
    done
done

mut_table=sig_genes.txt
for mut_res in `ls analyses__2013_09_23/*/20130923/$mut_prefix*$mut_suffix/*.$mut_table`; do
    cancer_name=`echo $mut_res | cut -d "/" -f 2`
    ln -s `pwd`/$mut_res ~/mendelian/cancer_alterations/$cancer_name.mut
done


cat cancer_gene_census.tsv | tr "\r" "\n" | sed 's/"//g' | sed 's/\.//g' | awk -F '\t' '{
                                                  gsub(/ /,",",$13);
                                                  split($13, alt, ",");
                                                  for (m in alt){print alt[m]; if(alt[m]=="Mis" || alt[m]=="S" || alt[m]=="F" || alt[m]=="N"){print $1 > "census_mutation"}; if(alt[m]=="D"){print $1 > "census_deletion"}; if(alt[m]=="A"){print $1 > "census_amplification"}; };}'

cat cancer_gene_census.tsv | tr "\r" "\n" | sed 's/"//g' | sed 's/\.//g' | awk -F '\t' '{
                                                  gsub(/ /,",",$13);
                                                  split($13, alt, ",");
                                                  for (m in alt){print alt[m]; if(alt[m]=="Mis" || alt[m]=="F" || alt[m]=="N"){print $1 > "census_mutation2"}; if(alt[m]=="D"){print $1 > "census_deletion"}; if(alt[m]=="A"){print $1 > "census_amplification"}; };}'



