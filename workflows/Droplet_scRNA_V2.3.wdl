workflow main{
	String Outdir
	String Root
	String Env
	String Species
	String cDNA_Fastq1
	String cDNA_Fastq2
	String Oligo_Fastq1
	String Oligo_Fastq2
	String BeadsBarcode
	String OligoBarcode
	String Refdir
	String Gtf
	String SampleName
	String Oligo_type8
	call makedir{
		input:
		outdir=Outdir
	}
    call parseFastq as cDNAparseFastq{
		input:
	        outdir=makedir.dir,
        	root=Root,
	        fastq1=cDNA_Fastq1,
        	fastq2=cDNA_Fastq2,
	        barcode=BeadsBarcode,
		barcode_counts_raw="cDNA_barcode_counts_raw.txt",
		readFq="cDNA_reads.fq",
		report="cDNA_sequencing_report.csv"
	}
	call parseFastq as M280parseFastq{
		input:
		outdir=makedir.dir,
		root=Root,
		fastq1=Oligo_Fastq1,
		fastq2=Oligo_Fastq2,
		barcode=OligoBarcode,
		barcode_counts_raw="Index_barcode_counts_raw_T1.txt",
		readFq="Index_reads_T1.fq",
		report="Index_sequencing_report_T1.csv"
	}
	call cDNAanno{
		input:
		outdir=makedir.dir,
		root=Root,
		refdir=Refdir,
		gtf=Gtf,
		fastq=cDNAparseFastq.fastqFile,
		cellBarcode=cDNAparseFastq.barcodeFile
	}
	call getM280UMI{
		input:
		outdir=makedir.dir,
		root=Root,
		env=Env,
		Oligo_fastq=M280parseFastq.fastqFile,
		sampleName=SampleName,
		oligo_type8=Oligo_type8,
		finalbam=cDNAanno.finalbam,
		beads_stat=cDNAanno.beads_stat
	}
	call count{
		input:
		outdir=makedir.dir,
		root=Root,
		anno_decon=getM280UMI.anno_decon
	}
	call QC{
		input:
		outdir=makedir.dir,
		root=Root,
		matrixdir=count.matrixdir,
		sampleName=SampleName
	}
	call Cluster{
		input:
		outdir=makedir.dir,
		root=Root,
		RDSlist=QC.RDSlist
	}
	call report{
		input:
		outdir=makedir.dir,
		root=Root,
		cell_report_2=Cluster.cell_report_2,
		sampleName=SampleName,
		species=Species,
		cutoff=getM280UMI.cut,
		DNAsequencing=cDNAparseFastq.sequencing,
		M280sequencing=M280parseFastq.sequencing,
		mapping=cDNAanno.mapping,
		annostat=cDNAanno.annostat,
		M280merge=getM280UMI.M280merge,
		cluster=Cluster.cluster,
		marker=Cluster.marker,
		mergestat=getM280UMI.mergestat,
		matrixdir=count.matrixdir,
		QCfig=QC.QCfig,		
		rm_alnSAM=cDNAanno.rm_alnSAM,
		rm_finalbam=cDNAanno.finalbam,
		rm_alnBAM=cDNAanno.rm_alnBAM,
		rm_annoBAM=cDNAanno.rm_annoBAM,
		rm_anno_decon=getM280UMI.anno_decon,
		rm_cDNAFq=cDNAparseFastq.fastqFile,
		rm_M280Fq=M280parseFastq.fastqFile
	}	
}

## task.14 report ##
task report{
	String outdir
	String root
	String cell_report_2
	String sampleName
	String species
	String cutoff
	String DNAsequencing
	String M280sequencing
	String mapping
	String annostat
	String M280merge
	String cluster
	String marker
	String mergestat
	String matrixdir
	String rm_alnSAM
	String rm_alnBAM
	String rm_annoBAM
	String rm_finalbam
	String rm_cDNAFq
	String rm_M280Fq
	String rm_anno_decon
	String QCfig
	command<<<
	    if [ -f ${outdir}/symbol/07.report_sigh.txt ];then
		echo "report node success"
	    else
		echo "SampleName,${sampleName}" >${outdir}/07.report/1.cell_report.csv &&
		echo "Species,${species}" >>${outdir}/07.report/1.cell_report.csv &&
		`Rscript ${root}/scripts/cell_report_1.R -I ${mergestat} -M ${matrixdir} -O ${outdir}/07.report/` &&
		`cat ${outdir}/07.report/cell_report_1.csv ${cell_report_2} >>${outdir}/07.report/1.cell_report.csv` &&
		`rm ${outdir}/07.report/cell_report_1.csv` &&
		`cp ${cutoff} ${outdir}/07.report/2.cut.off.csv` &&
		`cp ${DNAsequencing} ${outdir}/07.report/3.cDNA_sequencing_report.csv` &&
		`cp ${M280sequencing} ${outdir}/07.report/3.Index_sequencing_report_T1.csv` &&
		`cp ${mapping} ${outdir}/07.report/4.alignment_report.csv` &&
		`cp ${annostat} ${outdir}/07.report/5.anno_report.csv` &&
		`cp ${M280merge} ${outdir}/07.report/6.${sampleName}_CellNumber_merge.png` &&
		`cp ${cluster} ${outdir}/07.report/9.cluster.csv` &&
		`cp ${marker} ${outdir}/07.report/10.marker.csv` &&
		`cp ${QCfig} ${outdir}/07.report/7.${sampleName}_raw_QCplot.png` &&
		`python ${root}/scripts/pre_process-v1.0.1.py --outPath ${outdir} >/dev/null 2>&1 ` &&
		`python ${root}/scripts/generate-v1.0.1.py --outPath ${outdir} --ID ${sampleName} --htmlTemplate ${root}/scripts/template-v1.0.1.html` &&
#		`rm ${rm_alnSAM} ${rm_alnBAM} ${rm_annoBAM} ${rm_finalbam} ${rm_anno_decon}` &&
#		`rm ${rm_cDNAFq} ${rm_M280Fq}` &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/07.report_sigh.txt
	    fi
	>>>
}

## task.1 Create directory ##
task makedir{
	String outdir
	command{
		mkdir -p ${outdir}
		mkdir -p ${outdir}/00.log
		mkdir -p ${outdir}/01.DataStat
		mkdir -p ${outdir}/02.cDNAAnno
		mkdir -p ${outdir}/03.M280UMI_stat
		mkdir -p ${outdir}/04.Matrix
		mkdir -p ${outdir}/05.QCCluster
		mkdir -p ${outdir}/07.report
		mkdir -p ${outdir}/symbol
	}
	output{
		String dir="${outdir}"
	}
}

## task.2 parseFastq ##
task parseFastq{
	String outdir
	String root
	String fastq1
	String fastq2
	String barcode
	String barcode_counts_raw
	String readFq
	String report
	command<<<
	    if [ -f ${outdir}/symbol/01.DataStat_sigh.txt ];then
		echo "01.DataStat node success"
	    else
		${root}/software/scRNA_parse -t 4 -f -q 4 -dropN -config ${barcode} -cbdis ${outdir}/01.DataStat/${barcode_counts_raw} -1 ${outdir}/01.DataStat/${readFq} -report ${outdir}/01.DataStat/${report} ${fastq1} ${fastq2} &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/01.DataStat_sigh.txt
	    fi
	>>>
	output{
    	File fastqFile="${outdir}/01.DataStat/${readFq}"
		File barcodeFile="${outdir}/01.DataStat/${barcode_counts_raw}"
		File sequencing="${outdir}/01.DataStat/${report}"
	}
}


## task.3 1)Align, 2)add Tag in SAM, 3)MAPQ adjustment 4)sort, 5)Annotate gene and transcript with GTF, 6)UMI or Barcode correction ##
task cDNAanno{
	String fastq
	String outdir
	String refdir
	String root
	String gtf
	String cellBarcode
	command<<<
	    if [ -f ${outdir}/symbol/02.cDNAAnno_sigh.txt ];then
		echo "02.cDNAAnno node success"
	    else
		${root}/software/STAR --limitOutSJcollapsed 2000000 --outStd SAM --outSAMunmapped Within --runThreadN 10 --genomeDir ${refdir} --readFilesIn ${fastq} --outFileNamePrefix  ${outdir}/02.cDNAAnno/ 1> ${outdir}/02.cDNAAnno/aln.sam &&
		${root}/software/PISA sam2bam -adjust-mapq -gtf ${gtf} -o ${outdir}/02.cDNAAnno/aln.bam -report ${outdir}/02.cDNAAnno/alignment_report.csv ${outdir}/02.cDNAAnno/aln.sam &&
		${root}/software/PISA anno -gtf ${gtf} -o ${outdir}/02.cDNAAnno/anno.bam -report ${outdir}/02.cDNAAnno/anno_report.csv ${outdir}/02.cDNAAnno/aln.bam &&
		${root}/software/PISA corr -tag UR -new-tag UB -tags-block CB,GN -@ 10 -o ${outdir}/02.cDNAAnno/final.bam ${outdir}/02.cDNAAnno/anno.bam &&
		${root}/software/PISA attrcnt -cb CB -tags UB,GN -@ 5 -dedup -o ${outdir}/02.cDNAAnno/beads_stat.txt -list ${cellBarcode} ${outdir}/02.cDNAAnno/final.bam &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/02.cDNAAnno_sigh.txt
	    fi
	>>>
	output{
		String beads_stat="${outdir}/02.cDNAAnno/beads_stat.txt"
		File mapping="${outdir}/02.cDNAAnno/alignment_report.csv"
		File annostat="${outdir}/02.cDNAAnno/anno_report.csv"
		File finalbam="${outdir}/02.cDNAAnno/final.bam"
		File rm_alnSAM="${outdir}/02.cDNAAnno/aln.sam"
		File rm_alnBAM="${outdir}/02.cDNAAnno/aln.bam"
		File rm_annoBAM="${outdir}/02.cDNAAnno/anno.bam"
		
	}
}

## task.4 getM280 ## 
task getM280UMI{
	String outdir
	String root
	String env
	String Oligo_fastq
	String sampleName
	String oligo_type8
	String finalbam
	String beads_stat
	command<<<
	    if [ -f ${outdir}/symbol/03.M280UMI_stat_sigh.txt ];then
		echo "03.M280UMI_stat node success"
	    else
		Rscript ${root}/scripts/scRNA_cell_calling_v1.2.2.R -i ${beads_stat} -o ${outdir}/03.M280UMI_stat/ -e 0 -f 0 &&
		`cat  ${outdir}/01.DataStat/cDNA_barcode_counts_raw.txt |awk '{print $1}' >${outdir}/03.M280UMI_stat/beads_barcode_all.txt` &&
		export LD_LIBRARY_PATH="${env}/lib/:$LD_LIBRARY_PATH" && ${root}/software/mergeBarcodes -b ${outdir}/03.M280UMI_stat/beads_barcode_all.txt -f ${Oligo_fastq} -n ${sampleName} -o ${outdir}/03.M280UMI_stat/ && unset LD_LIBRARY_PATH &&
		${root}/software/s1.get.similarityOfBeads -n 4 ${sampleName} ${outdir}/03.M280UMI_stat/${sampleName}_CB_UB_count.txt ${outdir}/03.M280UMI_stat/beads_barcodes.txt ${oligo_type8} ${outdir}/03.M280UMI_stat/Similarity.all.csv ${outdir}/03.M280UMI_stat/Similarity.droplet.csv ${outdir}/03.M280UMI_stat/Similarity.droplet.filtered.csv &&
		python ${root}/scripts/s2.get.combinedListOfBeads.py ${outdir}/03.M280UMI_stat/Similarity.droplet.filtered.csv ${outdir}/03.M280UMI_stat/${sampleName}_combined_list.txt &&
		perl ${root}/scripts/fishInWinter.fq.gz.pl -bf table -ff table ${outdir}/03.M280UMI_stat/${sampleName}_combined_list.txt ${outdir}/03.M280UMI_stat/beads_barcodes.txt  --except >${outdir}/03.M280UMI_stat/${sampleName}_N1.beads  &&
		`awk -F '_N'  '{if ($2<=9) print $1"_N"$2}' ${outdir}/03.M280UMI_stat/${sampleName}_combined_list.txt >${outdir}/03.M280UMI_stat/${sampleName}_combined_list.filter.txt` &&
		`num=$(tail -n 1 ${outdir}/03.M280UMI_stat/${sampleName}_combined_list.filter.txt | cut -f 2 | awk -F"_" '{print $1}' | sed 's/CELL//') && awk -v N=$num '{print $1"\tCELL"N+NR"_N1"}' ${outdir}/03.M280UMI_stat/${sampleName}_N1.beads > ${outdir}/03.M280UMI_stat/${sampleName}_N1.cell &&  cat ${outdir}/03.M280UMI_stat/${sampleName}_N1.cell  ${outdir}/03.M280UMI_stat/${sampleName}_combined_list.filter.txt >${outdir}/03.M280UMI_stat/${sampleName}_barcodeTranslate.txt` &&
		${root}/software/tagAdd -n 4 -bam ${finalbam} -file ${outdir}/03.M280UMI_stat/${sampleName}_barcodeTranslate.txt -out ${outdir}/03.M280UMI_stat/anno_decon.bam -tag_check CB:Z: -tag_add DB:Z: -root ${root} &&
		perl ${root}/scripts/cell_stat.pl -c ${beads_stat} -m ${outdir}/03.M280UMI_stat/${sampleName}_barcodeTranslate.txt -o ${outdir}/03.M280UMI_stat/ &&
		Rscript ${root}/scripts/CellMergeStat.R -I ${outdir}/03.M280UMI_stat/${sampleName}_barcodeTranslate.txt -O ${outdir}/03.M280UMI_stat/ -n ${sampleName} &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/03.M280UMI_stat_sigh.txt
	    fi
	>>>
	output{
		File anno_decon="${outdir}/03.M280UMI_stat/anno_decon.bam"
		File M280merge="${outdir}/03.M280UMI_stat/${sampleName}_CellNumber_merge.png"
		File cut="${outdir}/03.M280UMI_stat/cutoff.csv"
		File mergestat="${outdir}/03.M280UMI_stat/merge_cell.stat"
	}

}

## task.5 get Raw matrix ##
task count{
		String root
		String outdir
		String anno_decon
		command<<<
		    if [ -f ${outdir}/symbol/04.Matrix_sigh.txt ];then
			echo "04.Matrix node success"
		    else
			${root}/software/PISA count -@ 10 -tag DB -anno-tag GN -umi UB -outdir ${outdir}/04.Matrix ${anno_decon} &&
			echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/04.Matrix_sigh.txt
		    fi
		>>>
		output{
			String matrixdir = "${outdir}/04.Matrix"
		}
}


## task.6 QC ##
task QC{
	String outdir
	String root
	String matrixdir
	String ?clusterdim
	Float ?doublepercentage
	String ?mitpercentage
	String ?minfeatures
	String sampleName
	command<<<
	    if [ -f ${outdir}/symbol/05.QC_sigh.txt ];then
		echo "05.QC node success"
	    else
		mkdir -p ${outdir}/05.QCCluster/QC &&
		Rscript ${root}/scripts/step1.QC.dir.ggplot.R -I ${matrixdir} -D ${default=20 clusterdim} -MP ${default=1 mitpercentage} -F ${default=200 minfeatures} -B ${sampleName} -O ${outdir}/05.QCCluster/ &&
		echo "${outdir}/05.QCCluster/QC/${sampleName}_QCobject.RDS" >${outdir}/05.QCCluster/QC/RDS.list &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/05.QC_sigh.txt
	    fi
	>>>
	output{
		File RDSlist="${outdir}/05.QCCluster/QC/RDS.list"
		File QCfig="${outdir}/05.QCCluster/QC/${sampleName}_raw_QCplot.png"
	}

}

## task.8 Cluster ##
task Cluster{
	String outdir
	String root
	String RDSlist
	String ?clusterdim
	String ?PCusage
	Float ?resolution=0.5
	command<<<
	    if [ -f ${outdir}/symbol/05.Cluster_sigh.txt ];then
		echo "05.Cluster node success"
	    else
		mkdir -p ${outdir}/05.QCCluster/Clustering &&
		Rscript ${root}/scripts/step2.1.BC.R -I ${RDSlist} -D ${default=20 clusterdim} -PC ${default=50 PCusage} -RES ${resolution}  -O ${outdir}/05.QCCluster/ &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/05.Cluster_sigh.txt
	    fi
	>>>
		output{
			File clusterRDS="${outdir}/05.QCCluster/Clustering/clustering_object.RDS"
			File cell_report_2="${outdir}/05.QCCluster/Clustering/cell_report_2.csv"
			File cluster="${outdir}/05.QCCluster/Clustering/cluster.csv"
			File marker="${outdir}/05.QCCluster/Clustering/marker.csv"
		}
}

## task.9 marker gene GO KEGG enrichment ##
task GOKEGG{
        String outdir
        String root
        String RDSlist
        String species
        command<<<
	    if [ -f ${outdir}/symbol/06.markerGeneGOKEGG_sigh.txt ];then
		echo "06.markerGeneGOKEGG node success"
	    else
                Rscript ${root}/scripts/GOKEGG_enrichment.R -d ${RDSlist} -o ${outdir}/06.markerGeneGOKEGG -s ${species} &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/06.markerGeneGOKEGG_sigh.txt
	    fi
        >>>
		output{
			File marker="${outdir}/06.markerGeneGOKEGG/marker.csv"
		}
}

## task.10 Record the start time of step ##
task stepstartlog{
	String outdir
	String sampleid
	String step
	command{
		echo "[INFO ] `date +%F` `date +%T` | ${sampleid} - ${step} start" >>${outdir}/000.Log/run.log
	}
	output{
		String success="done"
	}
}

## task.11 Record the finish time of step ##
task stepfinishlog{
	String outdir
	String step
	String sampleid
	command{
		echo "[INFO ] `date +%F` `date +%T` | ${sampleid} - ${step} completed" >>${outdir}/000.Log/run.log
	}
	output{
		String success="finish"
	}
}

## task.12 Return run code  ##
task sum{
	Array[Int] tar
	command<<<
		perl -we 'my $aa="${sep='' tar}";$aa=~s/0//g;if($aa){print "1\n";}else{print "0\n";}'
	>>>
	output{
		Int rc=read_lines(stdout())[0]
	}
}

## task.13 Step control
task DealStep{
	Map[String,Int] map
	command<<<
		perl -we 'my %ex;while(<>){my @aa=split;$ex{$aa[0]}=$aa[1]} foreach my $stp(("parseFastq", "getFinalBam", "countMatrix")){if(exists $ex{$stp}){print "$stp\t$ex{$stp}\n";}else{print "$stp\t0\n";}}' ${write_map(map)}
	>>>
	output{
		Map[String,Int] rc=read_map(stdout())
	}
}

