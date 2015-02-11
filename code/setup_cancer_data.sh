#!/bin/bash
firehose=analyses__2013_09_23/*/20130923 ## CHANGE THIS LINE FOR DIFFERENT FIREHOSE DOWNLOAD
cancer_links=data_source/cancer_alterations
mkdir $cancer_links
cn_suffix=CopyNumber_Gistic2.Level_4.*
del_table=table_del.conf_99.txt
amp_table=table_amp.conf_99.txt
cn_genes=all_data_by_genes.txt
cn_thres=all_thresholded.by_genes.txt
for cntype in {$del_table,$amp_table}; do
    for cn_res in `ls $firehose/$mut_prefix*$cn_suffix/$cntype`; do
        cancer_name=`echo $cn_res | cut -d "/" -f 2`
        f=$cancer_links/$cancer_name.$cntype
        #echo $cn_res
        if [ ! -e $f ]; then ln -s `pwd`/$cn_res $f; fi
        f=$cancer_links/$cancer_name.$cn_genes
        if [ ! -e $f ]; then ln -s `pwd`/`dirname $cn_res`/$cn_genes $f; fi
        f=$cancer_links/$cancer_name.$cn_thres
        if [ ! -e $f ]; then ln -s `pwd`/`dirname $cn_res`/$cn_thres $f; fi
    done
done

mut_suffix="-T*.MutSigNozzleReportMerged.Level_4.*"
mut_prefix="gdac.broadinstitute.org_"
mut_table=sig_genes.txt
for mut_res in `ls -1 $firehose/$mut_prefix*$mut_suffix/*.$mut_table`; do
    cancer_name=`echo $mut_res | cut -d "/" -f 2`
    #echo $mut_res
    ln -s `pwd`/$mut_res $cancer_links/$cancer_name.mut
done
echo DONE
