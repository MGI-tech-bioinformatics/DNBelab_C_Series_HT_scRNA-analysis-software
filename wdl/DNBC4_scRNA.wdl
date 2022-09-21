workflow main{
	String Outdir
	String Root
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
	Int? expectCellNum = 3000
	String? calling_method = 'emptydrops'
	Int? forceCellNum = 0
	Int? UMI_min = 1000
	Boolean? Intron = true
	Int? clusterdim = 20
	Float? doublepercentage = 0.05
	Int? mitpercentage = 15
	Int? minfeatures = 200
	String? mtgenes = 'auto'
	Int? PCusage = 50
	Float? resolution = 0.5

	call makedir{
		input:
		outdir=Outdir
	}
	call process_para{
		input:
		outdir=makedir.dir,
		root=Root,
		cDNA_fastq1=cDNA_Fastq1,
		cDNA_fastq2=cDNA_Fastq2,
		cDNA_barcode=BeadsBarcode,
		oligo_fastq1=Oligo_Fastq1,
		oligo_fastq2=Oligo_Fastq2,
		oligo_barcode=OligoBarcode,
		refdir=Refdir,
		gtf=Gtf,
		oligo_type8=Oligo_type8,
	}
	call oligo_parse{
		input:
		oligo_para=process_para.oligo_para,
		root=Root,
		outdir=makedir.dir,
	}
	call cDNAanno{
		input:
		outdir=makedir.dir,
		root=Root,
		cDNA_barcode=BeadsBarcode,
		cDNA_para=process_para.cDNA_para,
		refdir=Refdir,
		gtf=Gtf,
		intron=Intron,
		sampleName=SampleName,
	}
	call getM280UMI{
		input:
		outdir=makedir.dir,
		root=Root,
		expectCell=expectCellNum,
		forceCell=forceCellNum,
		method=calling_method,
		min_umi=UMI_min,
		Oligo_fastq=oligo_parse.oligoFastq,
		sampleName=SampleName,
		oligo_type8=Oligo_type8,
		finalbam=cDNAanno.finalbam,
		rawMatrix=cDNAanno.rawMatrix
	}
	call saturation{
		input:
		outdir=makedir.dir,
		root=Root,
		anno_decon=getM280UMI.anno_decon
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
		sampleName=SampleName,
		clusterdim=clusterdim,
		doublepercentage=doublepercentage,
		mitpercentage=mitpercentage,
		minfeatures=minfeatures,
		mtgenes=mtgenes
	}
	call Cluster{
		input:
		outdir=makedir.dir,
		QCdir=QC.QCdir,
		root=Root,
		clusterdim=clusterdim,
		species=Species,
		PCusage=PCusage,
		resolution=resolution
	}
	call report{
		input:
		outdir=makedir.dir,
		root=Root,
		cell_report=Cluster.cell_report,
		cellcount=count.cellcount,
		sampleName=SampleName,
		species=Species,
		cutoff=getM280UMI.cut,
		cDNAReport=cDNAanno.cDNAReport,
		oligoReport=oligo_parse.oligoReport,
		mapping=cDNAanno.mapping,
		annostat=cDNAanno.annostat,
		M280merge=getM280UMI.M280merge,
		cluster=Cluster.cluster,
		marker=Cluster.marker,
		QCfig=QC.QCfig,
		intron=Intron,
		# cell_type=Cluster.cell_type,
		saturation=saturation.satfile
	}
	if (Intron==true) {
        call splice_matrix{
			input: 
			outdir=makedir.dir,
			root=Root,
			#anno_decon=getM280UMI.anno_decon,
			attachment_dir=report.attachment_dir
		}
    }
}

## task.14 report ##
task report{
	String outdir
	String root
	String cell_report
	String cellcount
	String sampleName
	String species
	String cutoff
	String cDNAReport
	String oligoReport
	String mapping
	String annostat
	String M280merge
	String cluster
	String marker
	String saturation
	String QCfig
	Boolean? intron
	command<<<
        if [ -f ${outdir}/symbol/07.report_sigh.txt ];then
            echo "report node success"
	else
		python ${root}/DNBC4tools/rna/pre_process.py --outPath ${outdir} --sample ${sampleName} &&
		python ${root}/DNBC4tools/rna/generate_report.py --outPath ${outdir} --name ${sampleName} --htmlTemplate ${root}/DNBC4tools/template/template.html --species ${species} --intron ${true="true" false="false" intron } &&
		python ${root}/DNBC4tools/rna/report_output.py --indir ${outdir} &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/07.report_sigh.txt
	fi
	>>>
	output{
		String attachment_dir="${outdir}/output/attachment"
	}
}

## task.1 Create directory ##
task makedir{
	String outdir
	command{
		mkdir -p ${outdir}
		mkdir -p ${outdir}/01.data/raw_matrix
		mkdir -p ${outdir}/02.count/filter_matrix
		mkdir -p ${outdir}/03.analysis/QC
		mkdir -p ${outdir}/03.analysis/Clustering
		mkdir -p ${outdir}/04.report
        mkdir -p ${outdir}/symbol
	}
	output{
		String dir="${outdir}"
	}
}
## task.0 bc_para ##
task process_para{
	String cDNA_fastq1
	String cDNA_fastq2
	String cDNA_barcode
	String outdir
	String root
	String oligo_fastq1
	String oligo_fastq2
	String oligo_barcode
	String refdir
	String gtf
	String oligo_type8
	command<<<
		python ${root}/DNBC4tools/rna/bc_para.py --outdir ${outdir} --cDNAfastq1 ${cDNA_fastq1} --cDNAfastq2 ${cDNA_fastq2} --cDNAconfig ${cDNA_barcode} --oligofastq1 ${oligo_fastq1} --oligofastq2 ${oligo_fastq2} --oligoconfig ${oligo_barcode} --thread 10 --star_index ${refdir} --gtf ${gtf} --oligotype ${oligo_type8}
	>>>
	output{
		File cDNA_para="${outdir}/01.data/cDNA_para"
		File oligo_para="${outdir}/01.data/oligo_para"		
	}
}

## task.1 oligo parse
task oligo_parse{
	String oligo_para
	String root
	String outdir
	command<<<
        if [ -f ${outdir}/symbol/01.oligoparse_sigh.txt ];then
            	echo "01.oligoparse node success"
        else
		${root}/DNBC4tools/soft/parseFq ${oligo_para} &&
        	echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/01.oligoparse_sigh.txt
        fi
	>>>
	output{
		File oligoFastq="${outdir}/01.data/Index_reads.fq.gz"
		File oligoBarcode="${outdir}/01.data/Index_barcode_counts_raw.txt"
		File oligoReport="${outdir}/01.data/Index_sequencing_report.csv"		
	}	
}

## task.2 1)Align, add Tag in SAM, MAPQ adjustment, sort, 2)Annotate gene and transcript with GTF, 3)UMI or Barcode correction ##
task cDNAanno{
	String cDNA_para
	String outdir
	String refdir
	String root
	String gtf
	Boolean? intron
	String sampleName
	String cDNA_barcode
	command<<<
        if [ -f ${outdir}/symbol/02.cDNAAnno_sigh.txt ];then
            	echo "02.cDNAAnno node success"
        else
		${root}/DNBC4tools/soft/scStar --outSAMattributes singleCell --outSAMtype BAM Unsorted --genomeDir ${refdir} --outFileNamePrefix ${outdir}/01.data/ --stParaFile ${cDNA_para} --outSAMmode NoQS --runThreadN 10 --limitOutSJcollapsed 10000000 --limitIObufferSize 350000000 &&
		${root}/DNBC4tools/soft/Anno -I ${outdir}/01.data/Aligned.out.bam -a ${gtf} -L ${outdir}/01.data/cDNA_barcode_counts_raw.txt -o ${outdir}/01.data -c 10 -m chrM -B ${cDNA_barcode} ${true="--intron" false="" intron } --anno 1 &&
		rm -rf ${outdir}/01.data/Aligned.out.bam &&
		samtools sort -@ 10 ${outdir}/01.data/final.bam -o ${outdir}/01.data/final_sorted.bam &&
		rm -rf ${outdir}/01.data/final.bam &&
		${root}/DNBC4tools/soft/PISA count -@ 10 -cb CB -anno-tag GN -umi UB -outdir ${outdir}/01.data/raw_matrix ${outdir}/01.data/final_sorted.bam &&
        echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/02.cDNAAnno_sigh.txt
    fi
	>>>
	output{
		File beads_stat="${outdir}/01.data/beads_stat.txt"
		File mapping="${outdir}/01.data/alignment_report.csv"
		File annostat="${outdir}/01.data/anno_report.csv"
		File finalbam="${outdir}/01.data/final_sorted.bam"
		File rawMatrix="${outdir}/01.data/raw_matrix"
		File cDNABarcode="${outdir}/01.data/cDNA_barcode_counts_raw.txt"
		File cDNAReport="${outdir}/01.data/cDNA_sequencing_report.csv"	
	}
}

## task.3 cell calling and getM280 ## 
task getM280UMI{
	String outdir
	String root
	String Oligo_fastq
	String sampleName
	String oligo_type8
	String finalbam
	String rawMatrix
	String? method
	Int? expectCell
	Int? forceCell
	Int? min_umi
	command<<<
	if [ -f ${outdir}/symbol/03.M280UMI_stat_sigh.txt ];then
		echo "03.M280UMI_stat node success"
	else
		Rscript ${root}/DNBC4tools/rna/cell_calling.R --matrix ${outdir}/01.data/raw_matrix --outdir ${outdir}/02.count/ --method ${method} --expectcells ${expectCell} --forcecells ${forceCell} --minumi ${min_umi} &&
		python ${root}/DNBC4tools/rna/get_barcode.py --raw ${outdir}/01.data/cDNA_barcode_counts_raw.txt --select ${outdir}/02.count/beads_barcodes_hex.txt -o ${outdir}/02.count &&
		${root}/DNBC4tools/soft/mergeBarcodes -b ${outdir}/02.count/beads_barcode_all.txt -f ${Oligo_fastq} -n ${sampleName} -o ${outdir}/02.count/ &&
		${root}/DNBC4tools/soft/s1.get.similarityOfBeads -n 4 ${sampleName} ${outdir}/02.count/${sampleName}_CB_UB_count.txt ${outdir}/02.count/beads_barcodes.txt ${oligo_type8} ${outdir}/02.count/Similarity.all.csv ${outdir}/02.count/Similarity.droplet.csv ${outdir}/02.count/Similarity.droplet.filtered.csv &&
		python ${root}/DNBC4tools/rna/combinedListOfBeads.py --similarity_droplet ${outdir}/02.count/Similarity.droplet.csv --beads_list ${outdir}/02.count/beads_barcodes.txt --combined_list ${outdir}/02.count/${sampleName}_combined_list.txt &&
		python ${root}/DNBC4tools/rna/cellMerge.py --indir ${outdir}/02.count --name ${sampleName} &&
		${root}/DNBC4tools/soft/tagAdd -n 4  -bam ${finalbam} -file ${outdir}/02.count/${sampleName}_barcodeTranslate_hex.txt -out ${outdir}/02.count/anno_decon_sorted.bam -tag_check CB:Z: -tag_add DB:Z: &&
        echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/03.M280UMI_stat_sigh.txt
    fi
	>>>
	output{
		File M280merge="${outdir}/02.count/cellNumber_merge.png"
		File cut="${outdir}/02.count/cutoff.csv"
		File anno_decon = "${outdir}/02.count/anno_decon_sorted.bam"
    }
}

## task.4-1 count_matrix ##
task count{
	String outdir
	String root
	String anno_decon
	command<<<
	if [ -f ${outdir}/symbol/04.count_matrix_sigh.txt ];then
		echo "count_matrix node success"
	else
		${root}/DNBC4tools/soft/PISA count -one-hit -@ 10 -cb DB -anno-tag GN -umi UB -list ${outdir}/02.count/cell.id -outdir ${outdir}/02.count/filter_matrix ${outdir}/02.count/anno_decon_sorted.bam  &&
		Rscript ${root}/DNBC4tools/rna/cell_report.R -M ${outdir}/02.count/filter_matrix -O ${outdir}/02.count/ &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/04.count_matrix_sigh.txt
    fi
	>>>
	output{
		String matrixdir="${outdir}/02.count/filter_matrix"
		File cellcount='${outdir}/02.count/cellCount_report.csv'
    }
}

## task.intron ##
task splice_matrix{
	String outdir
	String root
	#String anno_decon
	String attachment_dir
	command<<<
	if [ -f ${outdir}/symbol/06.splice_matrix_sigh.txt ];then
		echo "splice_matrix node success"
	else
		mkdir -p ${outdir}/output/attachment/splice_matrix
		mkdir -p ${outdir}/output/attachment/RNAvelocity_matrix
		${root}/DNBC4tools/soft/PISA count -one-hit -@ 10 -cb DB -ttype E,S -anno-tag GN -umi UB -list ${outdir}/02.count/cell.id -outdir ${outdir}/output/attachment/splice_matrix ${outdir}/output/anno_decon_sorted.bam  &&
		${root}/DNBC4tools/soft/PISA count -one-hit -@ 10 -cb DB -velo -anno-tag GN -umi UB -list ${outdir}/02.count/cell.id -outdir ${outdir}/output/attachment/RNAvelocity_matrix ${outdir}/output/anno_decon_sorted.bam  &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/06.splice_matrix_sigh.txt
	fi
	>>>
}

## task.4-2 saturation ##
task saturation{
	String outdir
	String root
	String anno_decon
	command<<<
	if [ -f ${outdir}/symbol/04.saturation_sigh.txt ];then
		echo "saturation node success"
	else
		#samtools sort -@ 10 ${outdir}/02.count/anno_decon.bam -o ${outdir}/02.count/anno_decon_sorted.bam &&
		samtools index -@ 10 ${outdir}/02.count/anno_decon_sorted.bam &&
		python ${root}/DNBC4tools/rna/saturation.py  -i ${outdir}/02.count/anno_decon_sorted.bam -o ${outdir}/02.count -f ${outdir}/02.count/cellCount_report.csv --quality 20 --threads 10 &&
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/04.saturation_sigh.txt
	fi
	>>>
	output{
		File satfile="${outdir}/02.count/saturation.xls"
		### File anno_sort='${outdir}/02.count/anno_decon_sorted.bam'
	}
}	

## task.5_1 QC ##
task QC{
	String outdir
	String root
	String matrixdir
	Int? clusterdim
	Float? doublepercentage
	Int? mitpercentage
	Int? minfeatures
	String? mtgenes
	String sampleName
	command<<<
	if [ -f ${outdir}/symbol/05.QC_sigh.txt ];then
		echo "05.QC node success"
	else
		Rscript ${root}/DNBC4tools/rna/QC_analysis.R -I ${matrixdir} -D ${clusterdim} -P ${doublepercentage} -M ${mtgenes} -MP ${mitpercentage} -F ${minfeatures} -B ${sampleName} -O ${outdir}/03.analysis/ &&
        	echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/05.QC_sigh.txt
        fi
	>>>
	output{
		File QCfig="${outdir}/03.analysis/QC/raw_QCplot.png"
		String QCdir = "${outdir}/03.analysis/QC"
	}

}

## task.5_2 Cluster ##
task Cluster{
	String outdir
	String root
	String QCdir
	Int? clusterdim
	Int? PCusage
	String species
	Float? resolution
	command<<<
	if [ -f ${outdir}/symbol/05.Cluster_sigh.txt ];then
		echo "05.Cluster node success"
	else
		Rscript ${root}/DNBC4tools/rna/Cluster_analysis.R -I ${QCdir} -D ${clusterdim} -PC ${PCusage} -RES ${resolution} -O ${outdir}/03.analysis/ -SP ${species} &&
        	echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/05.Cluster_sigh.txt
        fi
	>>>
	output{
		# File clusterRDS="${outdir}/03.analysis/Clustering/clustering_annotation_object.RDS"
		File cell_report="${outdir}/03.analysis/Clustering/cell_report.csv"
		File cluster="${outdir}/03.analysis/Clustering/cluster.csv"
		File marker="${outdir}/03.analysis/Clustering/marker.csv"
		# File cell_type="${outdir}/03.analysis/Clustering/cluster_annotation.png"
	}
}

task bam2cram{
	String root
    	String outdir
    	String finalbam
		String matrixdir
    	String fai
	command<<<
	if [ -f ${outdir}/symbol/06.bam2cram.txt ];then
		echo "06.bam2cram node success"
	else
		samtools view -@ 10 -C -T ${fai} ${finalbam} > ${outdir}/01.data/raw_feature.cram &&
        samtools view -@ 10 -C -T ${fai} ${outdir}/02.count/anno_decon.bam > ${outdir}/02.count/filtered_feature.cram &&
		rm -rf ${outdir}/02.count/anno_decon.bam
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/06.bam2cram.txt
        fi
    	>>>
	output{
		String finalcram = "${outdir}/01.data/raw_feature.cram"
		String annocram = "${outdir}/02.count/filtered_feature.cram"
        }
}
