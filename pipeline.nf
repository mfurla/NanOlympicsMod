#!/usr/bin/env nextflow
/*
========================================================================================
                         mfurla/NanOlympicsMod
========================================================================================
 mfurla/NanOlympicsMod analysis pipeline.
 #### Homepage / Documentation
 https://github.com/mfurla/NanOlympicsMod
----------------------------------------------------------------------------------------
*/

def helpMessage() {
        log.info"""
    Usage:
    nextflow -c pipeline.conf run pipeline.nf --samples="/path/to/samples.txt" --resultsDir="/path/to/resultsDir"
    Mandatory arguments which may be specified in the pipeline.conf file

        --samples                                                Path to the tab-separated sample file including sample name, condition and path to base-called fast5 folder
        --test_condition                                         Condition that we are interested to profile (e. g. 'WT')
        --resultsDir                                             Path to a folder where to store results
        --fast5_slot                                             FAST5 slot containing the basecalled bases
        --fast5_slot_id                                          FAST5 slot containing the basecalled bases (redundant)
        --tombo_slot                                             FAST5 slot containing the resquiggled data
        --tombo_subslot                                          FAST5 slot containing the resquiggled data
        --transcriptome_fasta                                    Path to the transcriptome fasta file
        --transcriptome_fai                                      Path to the transcriptome fasta index file
        --genome_fasta                                           Path to the genome fasta file
        --genome_fai                                             Path to the genome fasta index file
        --genes2transcripts                                      Path to gene-to-transcripts file for Nanom6A
        --transcriptomebed                                       Path to transcripts bed12 file
        --genesbed                                               Path to genes bed file
        --gtf                                                    Path to genome annotation gtf file
        --nanom6AP                                               nanom6A probability thresholds for PR curve plotting
        --yanocompFDR                                            yanocomp FDR threshold
        --differrFDR                                             differr FDR threshold
        --drummerPval                                            drummer Pvalue threshold
        --epinanoErrorSumErr                                     epinanoError threshold sum of errors
        --epinanoErrorResiduals                                  epinanoError threshold residuals
        --postprocessingScript                                   Path to postprocessing R script
        --statisticalAnalysis                                    Path to statistical_analysis R script
        --binLength                                              Size of windows for genome binning
        --threshold                                              Set of thresholds to use for the filtering of m6A sites (choose between 'default' and 'relaxed')
        --peaksfile                                              Path to bed file with set of m6A gold-standard peaks
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Input of sample names, conditions, and FAST5s path.
Channel
	.fromPath( params.samples )
    .splitCsv(header: true, sep:'\t')
    .map{ row-> tuple(row.SampleName, row.Condition, file(row.DataPath)) }
    .set{multi2single_annot}

Channel
	.fromPath( params.samples )
    .splitCsv(header: true, sep:'\t')
    .map{ row-> tuple(row.SampleName, row.Condition, file(row.DataPath)) }
    .set{tombo_annot}

// Input of transcriptome fasta.
Channel
	.fromPath(params.transcriptome_fasta, checkIfExists:true)
	.into{transcriptome_fasta_minimap2;transcriptome_fasta_mines;transcriptome_fasta_nanom6a;transcriptome_fasta_differr;transcriptome_fasta_tombo1;transcriptome_fasta_tombo2;transcriptome_fasta_tombo3;transcriptome_fasta_nanopolish;transcriptome_fasta_dena;transcriptome_fasta_epinanoSVM;transcriptome_fasta_epinanoError;transcriptome_fasta_nanopolish1;transcriptome_fasta_xpore;transcriptome_fasta_nanocompore2}

// Input of transcriptome fasta index.
Channel
	.fromPath(params.transcriptome_fai, checkIfExists:true)
	.into{transcriptome_fai_minimap2;transcriptome_fai_mines;transcriptome_fai_nanom6a;transcriptome_fai_differr;transcriptome_fai_tombo1;transcriptome_fai_tombo2;transcriptome_fai_tombo3;transcriptome_fai_nanopolish;transcriptome_fai_dena;transcriptome_fai_epinanoSVM;transcriptome_fai_epinanoError;transcriptome_fai_nanopolish1;transcriptome_fai_xpore;transcriptome_fai_nanocompore2}

// Input of genome fasta.
Channel
	.fromPath(params.genome_fasta, checkIfExists:true)
	.into{genome_fasta_minimap2;genome_fasta_nanom6a;genome_fasta_differr;genome_fasta_eligos;genome_fasta_nanopolish2;genome_fasta_epinanoSVM;genome_fasta_epinanoError;genome_fasta_nanodoc;genome_fasta_drummer;genome_fasta_nanopolish1;genome_fasta_mines}

// Input of genome fasta index.
Channel
	.fromPath(params.genome_fai, checkIfExists:true)
	.into{genome_fai_minimap2;genome_fai_nanom6a;genome_fai_differr;genome_fai_eligos;genome_fai_nanopolish2;genome_fai_epinanoSVM;genome_fai_epinanoError;genome_fai_nanodoc;genome_fai_drummer;genome_fai_nanopolish1;genome_fai_mines}

// Input of genome bed.
Channel
	.fromPath(params.genomebed, checkIfExists:true)
	.into{bed_eligos;bed_nanodoc}

// Input of genome bed.
Channel
	.fromPath(params.transcriptomebed, checkIfExists:true)
	.into{bed_nanocompore}

// Input of genome gtf.
Channel
	.fromPath(params.gtf, checkIfExists:true)
	.into{gtf_xpore;gtf_yanocomp}

// From multiple read FAST5s to single read FAST5s.
process multi2single {
    input:
	tuple val(sample),val(condition),val(fast5) from multi2single_annot

    output:
    	tuple val(sample), val(condition) into singleReadFAST5_fastq
    	tuple val(condition), val(sample) into singleReadFAST5_tombo1
    	tuple val(condition), val(sample) into singleReadFAST5_nanodoc

    script:
    if(params.multi2single)
    """
    	mkdir -p ${params.resultsDir}
    	mkdir -p ${params.resultsDir}/${condition}
    	mkdir -p ${params.resultsDir}/${condition}/${sample}
    	mkdir -p ${params.resultsDir}/${condition}/${sample}/FAST5/

    	multi_to_single_fast5 --recursive --threads ${task.cpus} --input_path ${fast5} --save_path ${params.resultsDir}/${condition}/${sample}/FAST5/
    """
	else
	"""
	mkdir -p ${params.resultsDir}
    	mkdir -p ${params.resultsDir}/${condition}
    	mkdir -p ${params.resultsDir}/${condition}/${sample}
    	mkdir -p ${params.resultsDir}/${condition}/${sample}/FAST5/
    	mkdir -p ${params.resultsDir}/${condition}/${sample}/FAST5/0/


        #ls -trlh ${fast5}
        f5=\$(find ${fast5} | grep \"\\.fast5\");
        for single_f5 in \$f5; do
          cp \$single_f5 ${params.resultsDir}/${condition}/${sample}/FAST5/0/;
        done

    """
}

// Extract FASTQs from single read FAST5s files.
process fastq {
    input:
		tuple val(sample),val(condition) from singleReadFAST5_fastq

    output:
		tuple val(sample), val(condition) into singleReadFASTQ
    
    script:
    if(params.fastq)
    """
    	mkdir -p ${params.resultsDir}/${condition}/${sample}/FASTQ/
        fast5_dir=\$(find ${params.resultsDir}/${condition}/${sample}/FAST5/ -maxdepth 1 -type d)
        if [ -f ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq ] ; then
          rm ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq
        fi
        for d in \$fast5_dir; do
    	  poretools fastq --group ${params.fast5_slot_id} \$d >> ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq
        done
		mkdir -p ${params.resultsDir}/${condition}/FASTQ/
		mkdir -p ${params.resultsDir}/${condition}/FASTQ/${sample}/
		ln -sf ${params.resultsDir}/${condition}/${sample}/FASTQ/ ${params.resultsDir}/${condition}/FASTQ/${sample}/
    """
	else
	"""
		mkdir -p ${params.resultsDir}/${condition}/FASTQ/
		mkdir -p ${params.resultsDir}/${condition}/FASTQ/${sample}/
		ln -sf ${params.resultsDir}/${condition}/${sample}/FASTQ/ ${params.resultsDir}/${condition}/FASTQ/${sample}/
    """
}

// Genomic and transcriptomic alignment.
process minimap2 {
    input:
		tuple val(sample),val(condition) from singleReadFASTQ

		each file('transcriptome.fa') from transcriptome_fasta_minimap2
		each file('transcriptome.fa.fai') from transcriptome_fai_minimap2

		each file('genome.fa') from genome_fasta_minimap2
		each file('genome.fa.fai') from genome_fai_minimap2

    output:
    	tuple val(condition), val(sample) into minimap2_nanopolish1
		tuple val(condition), file('minimap.filt.sortT.bam'), file('minimap.sortG.bam') into minimap2_minimap2Merge

    	tuple val(condition), file('minimap.filt.sortT.bam'), file('minimap.filt.sortT.bam.bai'), file('minimap.sortG.bam'), file('minimap.sortG.bam.bai') into minimap2_differr
    	tuple val(condition), file('minimap.sortG.bam'), file('minimap.sortG.bam.bai') into minimap2_eligos
    	tuple val(condition), file('minimap.filt.sortT.bam'), file('minimap.filt.sortT.bam.bai'), file('minimap.sortG.bam'), file('minimap.sortG.bam.bai') into minimap2_drummer

    script:
    if(params.minimap2)
    """
		/bin/minimap2/minimap2 -x map-ont -k14 -t ${task.cpus} -a transcriptome.fa ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq | samtools view -hSb | samtools sort -@ ${task.cpus} -o minimapT.bam
		samtools view minimapT.bam -bh -t transcriptome.fa.fai -F 2324 | samtools sort -@ ${task.cpus} -o minimap.filt.sortT.bam
		samtools index -@ ${task.cpus} minimap.filt.sortT.bam

		mkdir -p ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/
		cp minimapT.bam ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.bam
		cp minimap.filt.sortT.bam ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam
		cp minimap.filt.sortT.bam.bai ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam.bai

		/bin/minimap2/minimap2 -ax splice -k14 -t ${task.cpus} genome.fa ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq | samtools view -hSb | samtools sort -@ ${task.cpus} -o minimapG.bam
		samtools view minimapG.bam -bh -t genome.fa.fai -F 2308 | samtools sort -@ ${task.cpus} -o minimap.sortG.bam
		samtools index -@ ${task.cpus} minimap.sortG.bam

		mkdir -p ${params.resultsDir}/${condition}/${sample}/genomeAlignment/
		cp minimapG.bam ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.bam
		cp minimap.sortG.bam ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.sort.bam
		cp minimap.sortG.bam.bai ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.sort.bam.bai
    """
	else
	"""
		ln -s ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam ./minimap.filt.sortT.bam
		ln -s ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam.bai ./minimap.filt.sortT.bam.bai
		
		ln -s ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.sort.bam ./minimap.sortG.bam
		ln -s ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.sort.bam.bai ./minimap.sortG.bam.bai
    """
}

// From a single channel for all the alignments to one channel for each condition.
ni_test_minimap2Merge=Channel.create()
ni_other_minimap2Merge=Channel.create()
minimap2_minimap2Merge.groupTuple(by:0)
	.choice( ni_test_minimap2Merge, ni_other_minimap2Merge ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// From multiple to single bam files.
process minimap2Merge {
    input:
	    tuple val('condition1'), file('minimap.filt.sortT.1.*.bam'), file('minimap.sortG.1.*.bam') from ni_test_minimap2Merge
	    tuple val('condition2'), file('minimap.filt.sortT.2.*.bam'), file('minimap.sortG.2.*.bam') from ni_other_minimap2Merge
    
    output:
		tuple file('minimap.filt.sort.1.bam'), file('minimap.filt.sort.1.bam.bai'), file('minimap.filt.sort.2.bam'), file('minimap.filt.sort.2.bam.bai') into minimap2Merge_dena
		tuple val(condition1), file('minimap.filt.sort.1.bam'), file('minimap.filt.sort.1.bam.bai'), file('minimap.sort.1.bam'), file('minimap.sort.1.bam.bai') into minimap2Merge_1_epinanoSVM
		tuple val(condition2), file('minimap.filt.sort.2.bam'), file('minimap.filt.sort.2.bam.bai'), file('minimap.sort.2.bam'), file('minimap.sort.2.bam.bai') into minimap2Merge_2_epinanoSVM
		tuple val(condition1), file('minimap.filt.sort.1.bam'), file('minimap.filt.sort.1.bam.bai'), file('minimap.sort.1.bam'), file('minimap.sort.1.bam.bai') into minimap2Merge_1_epinanoError
		tuple val(condition2), file('minimap.filt.sort.2.bam'), file('minimap.filt.sort.2.bam.bai'), file('minimap.sort.2.bam'), file('minimap.sort.2.bam.bai') into minimap2Merge_2_epinanoError

    script:
    if(params.minimap2Merge)
    """
		cat ${params.resultsDir}/${condition1}/FASTQ/*/FASTQ/singleReadsFASTQ.fastq > ${params.resultsDir}/${condition1}/FASTQ/singleReadsFASTQ.fastq
		cat ${params.resultsDir}/${condition2}/FASTQ/*/FASTQ/singleReadsFASTQ.fastq > ${params.resultsDir}/${condition2}/FASTQ/singleReadsFASTQ.fastq

    	mkdir -p ${params.resultsDir}/${condition1}/
		mkdir -p ${params.resultsDir}/${condition1}/transcriptomeAlignment/
		samtools merge -f ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.bam minimap.filt.sortT.1.*.bam
		samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.sort.bam ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.bam
		samtools index -@ ${task.cpus} ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.sort.bam

		mkdir -p ${params.resultsDir}/${condition1}/genomeAlignment/
		samtools merge -f ${params.resultsDir}/${condition1}/genomeAlignment/minimap.bam minimap.sortG.1.*.bam
		samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition1}/genomeAlignment/minimap.sort.bam ${params.resultsDir}/${condition1}/genomeAlignment/minimap.bam
		samtools index -@ ${task.cpus} ${params.resultsDir}/${condition1}/genomeAlignment/minimap.sort.bam

    	mkdir -p ${params.resultsDir}/${condition2}/
		mkdir -p ${params.resultsDir}/${condition2}/transcriptomeAlignment/
		samtools merge -f ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.bam minimap.filt.sortT.2.*.bam
		samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.sort.bam ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.bam
		samtools index -@ ${task.cpus} ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.sort.bam

		mkdir -p ${params.resultsDir}/${condition2}/genomeAlignment/
		samtools merge -f ${params.resultsDir}/${condition2}/genomeAlignment/minimap.bam minimap.sortG.2.*.bam
		samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition2}/genomeAlignment/minimap.sort.bam ${params.resultsDir}/${condition2}/genomeAlignment/minimap.bam
		samtools index -@ ${task.cpus} ${params.resultsDir}/${condition2}/genomeAlignment/minimap.sort.bam

		ln -s ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.sort.bam minimap.filt.sort.1.bam
		ln -s ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.sort.bam.bai minimap.filt.sort.1.bam.bai

		ln -s ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.sort.bam minimap.filt.sort.2.bam
		ln -s ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.sort.bam.bai minimap.filt.sort.2.bam.bai

		ln -s ${params.resultsDir}/${condition1}/genomeAlignment/minimap.sort.bam minimap.sort.1.bam
		ln -s ${params.resultsDir}/${condition1}/genomeAlignment/minimap.sort.bam.bai minimap.sort.1.bam.bai

		ln -s ${params.resultsDir}/${condition2}/genomeAlignment/minimap.sort.bam minimap.sort.2.bam
		ln -s ${params.resultsDir}/${condition2}/genomeAlignment/minimap.sort.bam.bai minimap.sort.2.bam.bai
    """
	else
	"""
		ln -s ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.sort.bam minimap.filt.sort.1.bam
		ln -s ${params.resultsDir}/${condition1}/transcriptomeAlignment/minimap.filt.sort.bam.bai minimap.filt.sort.1.bam.bai

		ln -s ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.sort.bam minimap.filt.sort.2.bam
		ln -s ${params.resultsDir}/${condition2}/transcriptomeAlignment/minimap.filt.sort.bam.bai minimap.filt.sort.2.bam.bai

		ln -s ${params.resultsDir}/${condition1}/genomeAlignment/minimap.sort.bam minimap.sort.1.bam
		ln -s ${params.resultsDir}/${condition1}/genomeAlignment/minimap.sort.bam.bai minimap.sort.1.bam.bai

		ln -s ${params.resultsDir}/${condition2}/genomeAlignment/minimap.sort.bam minimap.sort.2.bam
		ln -s ${params.resultsDir}/${condition2}/genomeAlignment/minimap.sort.bam.bai minimap.sort.2.bam.bai
	"""
}

// From a single channel for all the alignments to one channel for each condition.
ni_test_tombo1=Channel.create()
ni_other_tombo1=Channel.create()
singleReadFAST5_tombo1.groupTuple(by:0)
	.choice( ni_test_tombo1, ni_other_tombo1 ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// Resquiggle for each condition.
process tombo1 {
    input:
		tuple val('condition1'), val('samples') from ni_test_tombo1
		tuple val('condition2'), val('samples') from ni_other_tombo1

		each file('transcriptome.fa') from transcriptome_fasta_tombo1
		each file('transcriptome.fa.fai') from transcriptome_fai_tombo1

    output:
		tuple val(condition1) into condition1_tombo1_tombo2
		tuple val(condition2) into condition2_tombo1_tombo2

		tuple val(condition1) into condition1_tombo1_tombo3
		tuple val(condition2) into condition2_tombo1_tombo3

		tuple val(condition1) into condition1_tombo1_nanom6a
		tuple val(condition2) into condition2_tombo1_nanom6a

		tuple val(condition1) into condition1_tombo1_dena
		tuple val(condition2) into condition2_tombo1_dena

    script:
    if(params.tombo1)
    """
		/bin/miniconda3/bin/tombo resquiggle ${params.resultsDir}/${condition1}/ --ignore-read-locks --overwrite --basecall-group ${params.fast5_slot} transcriptome.fa --processes ${task.cpus} --fit-global-scale --include-event-stdev --failed-reads-filename ${params.resultsDir}/${condition1}/failedReads.txt

		/bin/miniconda3/bin/tombo resquiggle ${params.resultsDir}/${condition2}/ --ignore-read-locks --overwrite --basecall-group ${params.fast5_slot} transcriptome.fa --processes ${task.cpus} --fit-global-scale --include-event-stdev --failed-reads-filename ${params.resultsDir}/${condition2}/failedReads.txt
    """
	else
	"""
		echo "Skipped"
    """
}

// RNA modifications detection with Tombo denovo
process tombo2 {
    input:
	    tuple val('condition1') from condition1_tombo1_tombo2
	    tuple val('condition2') from condition2_tombo1_tombo2

		each file('transcriptome.fa') from transcriptome_fasta_tombo2
		each file('transcriptome.fa.fai') from transcriptome_fai_tombo2
    
    output:
	    tuple val(condition1) into condition1_tombo2_mines
	    tuple val(condition2) into condition2_tombo2_mines
	    val('flagtombo2') into tombo2_postprocessing

    script:
    if(params.tombo2)
    """
    	mkdir -p ${params.resultsDir}/${condition1}/tomboDenovo

		/bin/miniconda3/bin/tombo detect_modifications de_novo --per-read-statistics-basename ${params.resultsDir}/${condition1}/tomboDenovo/Per_read_Stats_Filename --fast5-basedirs ${params.resultsDir}/${condition1}/ --statistics-file-basename ${params.resultsDir}/${condition1}/tomboDenovo/Stats_Filename

		/bin/miniconda3/bin/tombo text_output browser_files --fast5-basedirs ${params.resultsDir}/${condition1}/ --statistics-filename ${params.resultsDir}/${condition1}/tomboDenovo/Stats_Filename.tombo.stats --browser-file-basename ${params.resultsDir}/${condition1}/tomboDenovo/output_filename --file-types statistic fraction dampened_fraction coverage valid_coverage signal 

#    	mkdir -p ${params.resultsDir}/${condition2}/tomboDenovo

#		/bin/miniconda3/bin/tombo detect_modifications de_novo --per-read-statistics-basename ${params.resultsDir}/${condition1}/tomboDenovo/Per_read_Stats_Filename --fast5-basedirs ${params.resultsDir}/${condition2}/ --statistics-file-basename ${params.resultsDir}/${condition2}/tomboDenovo/Stats_Filename

#		/bin/miniconda3/bin/tombo text_output browser_files --fast5-basedirs ${params.resultsDir}/${condition2}/ --statistics-filename ${params.resultsDir}/${condition2}/tomboDenovo/Stats_Filename.tombo.stats --browser-file-basename ${params.resultsDir}/${condition2}/tomboDenovo/output_filename --file-types statistic fraction dampened_fraction coverage valid_coverage signal 
    """
	else
	"""
		echo "Skipped"
    """
}

// RNA modifications detection with Tombo comparing samples
process tombo3 {
    input:
	    tuple val('condition1') from condition1_tombo1_tombo3
	    tuple val('condition2') from condition2_tombo1_tombo3

    output:
    	val('flagtombo3') into tombo3_postprocessing
    script:
    if(params.tombo3)
    """
    	mkdir -p ${params.resultsDir}/tomboComparison/

    	/bin/miniconda3/bin/tombo detect_modifications level_sample_compare \
    	--fast5-basedirs ${params.resultsDir}/${condition1}/ \
		--alternate-fast5-basedirs ${params.resultsDir}/${condition2}/ \
		--minimum-test-reads 50 \
		--processes ${task.cpus} --statistics-file-basename ${params.resultsDir}/tomboComparison/sample.level_samp_comp_detect \
		--store-p-value

		/bin/miniconda3/bin/tombo text_output browser_files --statistics-filename ${params.resultsDir}/tomboComparison/sample.level_samp_comp_detect.tombo.stats \
    	--browser-file-basename ${params.resultsDir}/tomboComparison/sample.level_samp_comp_detect --file-types statistic
    """
	else
	"""
        echo "Skipped"
    """
}

// RNA modifications detection with nanom6A
process nanom6a {
    input:
		tuple val('condition1') from condition1_tombo1_nanom6a
		tuple val('condition2') from condition2_tombo1_nanom6a

		each file('genome.fa') from genome_fasta_nanom6a
		each file('genome.fa.fai') from genome_fai_nanom6a

		each file('transcriptome.fa') from transcriptome_fasta_nanom6a
		each file('transcriptome.fa.fai') from transcriptome_fai_nanom6a

    output:
		tuple file('genome.dict'), file('transcriptome.dict') into picard_epinanoSVM
		tuple file('genome.dict'), file('transcriptome.dict') into picard_epinanoError
		val('flagnanom6a') into nanom6a_postprocessing
    
    script:
    if(params.nanom6a)
    """
    	java -jar /picard/build/libs/picard.jar CreateSequenceDictionary -R genome.fa -O genome.dict
		java -jar /picard/build/libs/picard.jar CreateSequenceDictionary -R transcriptome.fa -O transcriptome.dict

		mkdir -p ${params.resultsDir}/${condition1}/nanom6a/

		find ${params.resultsDir}/${condition1}/ -name "*.fast5" > ${params.resultsDir}/${condition1}/nanom6a/files.txt
		/nanom6A_2021_10_22/bin/extract_raw_and_feature_fast --cpu=${task.cpus} --fl=${params.resultsDir}/${condition1}/nanom6a/files.txt -o ${params.resultsDir}/${condition1}/nanom6a/result --clip=10 --basecall_group ${params.tombo_slot} --basecall_subgroup ${params.tombo_subslot}

		for prob in ${params.nanom6AP};do /nanom6A_2021_10_22/bin/predict_sites --cpu ${task.cpus} -i ${params.resultsDir}/${condition1}/nanom6a/result -o ${params.resultsDir}/${condition1}/nanom6a/result_final -r transcriptome.fa -g genome.fa -b ${params.genes2transcripts} --model /nanom6A_2021_10_22/bin/model/ --proba \$prob; done

#		mkdir -p ${params.resultsDir}/${condition2}/nanom6a/

#		find ${params.resultsDir}/${condition2}/ -name "*.fast5" > ${params.resultsDir}/${condition2}/nanom6a/files.txt
#		/nanom6A_2021_10_22/bin/extract_raw_and_feature_fast --cpu=${task.cpus} --fl=${params.resultsDir}/${condition2}/nanom6a/files.txt -o ${params.resultsDir}/${condition2}/nanom6a/result --clip=10 --basecall_group ${params.tombo_slot} --basecall_subgroup ${params.tombo_subslot}

#		for prob in ${params.nanom6AP};do /nanom6A_2021_10_22/bin/predict_sites --cpu ${task.cpus} -i ${params.resultsDir}/${condition2}/nanom6a/result -o ${params.resultsDir}/${condition2}/nanom6a/result_final -r transcriptome.fa -g genome.fa -b ${params.genes2transcripts} --model /nanom6A_2021_10_22/bin/model/ --proba \$prob; done
    """
	else
	"""
    	java -jar /picard/build/libs/picard.jar CreateSequenceDictionary R=genome.fa O=genome.dict
		java -jar /picard/build/libs/picard.jar CreateSequenceDictionary R=transcriptome.fa O=transcriptome.dict
    """
}

// From a single channel for all the alignments to one channel for each condition.
ni_test=Channel.create()
ni_other=Channel.create()
minimap2_differr.groupTuple(by:0)
	.choice( ni_test, ni_other ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with differr
process differr {
    input:
	    tuple val('condition1'), file('minimap.filt.sortT.1.*.bam'), file('minimap.filt.sortT.1.*.bam.bai'), file('minimap.sortG.1.*.bam'), file('minimap.sortG.1.*.bam.bai') from ni_test
	    tuple val('condition2'), file('minimap.filt.sortT.2.*.bam'), file('minimap.filt.sortT.2.*.bam.bai'), file('minimap.sortG.2.*.bam'), file('minimap.sortG.2.*.bam.bai') from ni_other

		each file('genome.fa') from genome_fasta_differr
		each file('genome.fa.fai') from genome_fai_differr

		each file('transcriptome.fa') from transcriptome_fasta_differr
		each file('transcriptome.fa.fai') from transcriptome_fai_differr

    output:
    	val('flagdifferr') into differr_postprocessing
    script:
    if(params.differr)
    """
    	mkdir -p ${params.resultsDir}/differr/

		differr -p ${task.cpus} \
		\$(for file in minimap.sortG.2*.bam; do echo -a \$file; done) \
		\$(for file in minimap.sortG.1*.bam; do echo -b \$file; done) \
		-r genome.fa -o ${params.resultsDir}/differr/differrOut.bed \
		-f ${params.differrFDR}

    """
	else
	"""
        echo "Skipped"
    """
}

// From a single channel for all the alignments to one channel for each condition.
ni_test_eligos=Channel.create()
ni_other_eligos=Channel.create()
minimap2_eligos.groupTuple(by:0)
	.choice( ni_test_eligos, ni_other_eligos ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with eligos
process eligos {
    input:
	    tuple val('condition1'), file('minimap.sortG.1.*.bam'), file('minimap.sortG.1.*.bam.bai') from ni_test_eligos
	    tuple val('condition2'), file('minimap.sortG.2.*.bam'), file('minimap.sortG.2.*.bam.bai') from ni_other_eligos

		each file('genome.fa') from genome_fasta_eligos
		each file('genome.fa.fai') from genome_fai_eligos

		each file('genome.bed') from bed_eligos

    output:
    	val('flageligos') into eligos_postprocessing
    script:
    if(params.eligos)
    """
    	# Merged replicates comparison
    	mkdir -p ${params.resultsDir}/eligos/merged/

    	samtools merge minimap.sortG.1.bam minimap.sortG.1.*.bam
    	samtools merge minimap.sortG.2.bam minimap.sortG.2.*.bam

    	samtools index -@ ${task.cpus} minimap.sortG.1.bam
    	samtools index -@ ${task.cpus} minimap.sortG.2.bam

		/eligos2-v2.1.0/eligos2 pair_diff_mod -t ${task.cpus} \
		-tbam minimap.sortG.1.bam \
		-cbam minimap.sortG.2.bam \
		-reg genome.bed \
		-ref genome.fa \
		-o ${params.resultsDir}/eligos/merged/

    	# Splitted replicates comparison
    	mkdir -p ${params.resultsDir}/eligos/

		/eligos2-v2.1.0/eligos2 pair_diff_mod -t ${task.cpus} \
		-tbam minimap.sortG.1.*.bam \
		-cbam minimap.sortG.2.*.bam \
		-reg genome.bed \
		-ref genome.fa \
		-o ${params.resultsDir}/eligos/
    """
	else
	"""
        echo "Skipped"
    """
}

// RNA modifications detection with mines for each condition
process mines {
    input:
		tuple val('condition1') from condition1_tombo2_mines
		tuple val('condition2') from condition2_tombo2_mines

		each file('transcriptome.fa') from transcriptome_fasta_mines
		each file('transcriptome.fa.fai') from transcriptome_fai_mines

    output:
    	val('flagmines') into mines_postprocessing
    script:
    if(params.mines)
    """
		mkdir -p ${params.resultsDir}/${condition1}/mines/
#		mkdir -p ${params.resultsDir}/${condition2}/mines/

		wig2bed < ${params.resultsDir}/${condition1}/tomboDenovo/output_filename.fraction_modified_reads.plus.wig > ${params.resultsDir}/${condition1}/mines/output_filename.fraction_modified_reads.plus.wig.bed
		
#		wig2bed < ${params.resultsDir}/${condition2}/tomboDenovo/output_filename.fraction_modified_reads.plus.wig > ${params.resultsDir}/${condition2}/mines/output_filename.fraction_modified_reads.plus.wig.bed

		python3 /MINES/cDNA_MINES.py --fraction_modified ${params.resultsDir}/${condition1}/mines/output_filename.fraction_modified_reads.plus.wig.bed --coverage ${params.resultsDir}/${condition1}/tomboDenovo/output_filename.coverage.plus.bedgraph --output ${params.resultsDir}/${condition1}/mines/m6A_output_filename.bed --ref transcriptome.fa --kmer_models /MINES/Final_Models/names.txt
		
#		python3 /MINES/cDNA_MINES.py --fraction_modified ${params.resultsDir}/${condition2}/mines/output_filename.fraction_modified_reads.plus.wig.bed --coverage ${params.resultsDir}/${condition2}/tomboDenovo/output_filename.coverage.plus.bedgraph --output ${params.resultsDir}/${condition2}/mines/m6A_output_filename.bed --ref transcriptome.fa --kmer_models /MINES/Final_Models/names.txt

   """
	else
	"""
        echo "Skipped"
    """
}

// RNA modifications detection with dena for each condition
process dena {
    input:
    	tuple file('minimap.filt.sort.1.bam'), file('minimap.filt.sort.1.bam.bai'), file('minimap.filt.sort.2.bam'), file('minimap.filt.sort.2.bam.bai') from minimap2Merge_dena

		tuple val('condition1') from condition1_tombo1_dena
		tuple val('condition2') from condition2_tombo1_dena

		each file('transcriptome.fa') from transcriptome_fasta_dena
		each file('transcriptome.fa.fai') from transcriptome_fai_dena

    output:
    	val('flagdena') into dena_postprocessing
    script:
    if(params.dena)
    """
		mkdir -p ${params.resultsDir}/${condition1}/dena/

		python3 /DENA/step4_predict/LSTM_extract.py get_pos --fasta transcriptome.fa --motif 'RRACH' --output ${params.resultsDir}/${condition1}/dena/candidate_predict_pos.txt

		python3 /DENA/step4_predict/LSTM_extract.py predict --fast5 ${params.resultsDir}/${condition1} --corr_grp ${params.tombo_slot} --bam minimap.filt.sort.1.bam --sites ${params.resultsDir}/${condition1}/dena/candidate_predict_pos.txt --label "dena_label" --windows 2 2 --processes ${task.cpus}

		mv *_tmp ${params.resultsDir}/${condition1}/dena/

		python3 /DENA/step4_predict/LSTM_predict.py -i ${params.resultsDir}/${condition1}/dena/ -m /DENA/denaModels/ -o ${params.resultsDir}/${condition1}/dena/ -p "dena_label"

		# mkdir -p ${params.resultsDir}/${condition2}/dena/

		# python3 /DENA/step4_predict/LSTM_extract.py get_pos --fasta transcriptome.fa --motif 'RRACH' --output ${params.resultsDir}/${condition2}/dena/candidate_predict_pos.txt

		# python3 /DENA/step4_predict/LSTM_extract.py predict --fast5 ${params.resultsDir}/${condition2} --corr_grp ${params.tombo_slot} --bam minimap.filt.sort.2.bam --sites ${params.resultsDir}/${condition2}/dena/candidate_predict_pos.txt --label "dena_label" --windows 2 2 --processes ${task.cpus}

		# mv *_tmp ${params.resultsDir}/${condition2}/dena/

		# python3 /DENA/step4_predict/LSTM_predict.py -i ${params.resultsDir}/${condition2}/dena/ -m /DENA/denaModels/ -o ${params.resultsDir}/${condition2}/dena/ -p "dena_label"
    """
	else
	"""
        echo "Skipped"
    """
}

// RNA modifications detection with epinano in SVM mode for each condition
process epinanoSVM {
    input:
		tuple val('condition1'), file('minimap.filt.sort.1.bam'), file('minimap.filt.sort.1.bam.bai'), file('minimap.sort.1.bam'), file('minimap.sort.1.bam.bai') from minimap2Merge_1_epinanoSVM
		tuple val('condition2'), file('minimap.filt.sort.2.bam'), file('minimap.filt.sort.2.bam.bai'), file('minimap.sort.2.bam'), file('minimap.sort.2.bam.bai') from minimap2Merge_2_epinanoSVM

		each file('genome.fa') from genome_fasta_epinanoSVM
		each file('genome.fa.fai') from genome_fai_epinanoSVM

		each file('transcriptome.fa') from transcriptome_fasta_epinanoSVM
		each file('transcriptome.fa.fai') from transcriptome_fai_epinanoSVM

		tuple file('genome.fa.dict'), file('transcriptome.fa.dict') from picard_epinanoSVM
    output:
    	val('flagepinanoSVM') into epinanoSVM_postprocessing
    script:
    if(params.epinanoSVM)
    """
    	mkdir -p ${params.resultsDir}/${condition1}/epinanoSVM/
        mkdir -p ${params.resultsDir}/${condition2}/epinanoSVM/

    	samtools view -F16 minimap.sort.1.bam > minimap.sort.1.plus.sam
        samtools view -f16 minimap.sort.1.bam > minimap.sort.1.minus.sam
        samtools view -F16 minimap.sort.2.bam > minimap.sort.2.plus.sam
        samtools view -f16 minimap.sort.2.bam > minimap.sort.2.minus.sam

        if [[ -s minimap.sort.1.plus.sam ]]; then
           if [[ -s minimap.sort.1.minus.sam ]]; then
               /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -n ${task.cpus} -T g -R genome.fa -b minimap.sort.1.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar
               /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.1.plus_strand.per.site.csv 5
	       /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.1.minus_strand.per.site.csv 5
	       /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.1.plus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix plus_mod_prediction
	       /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.1.minus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix minus_mod_prediction

	       mv minimap.sort.1.plus_strand.per.site.csv minimap.sort.1.minus_strand.per.site.csv minimap.sort.1.plus_strand.per.site.5mer.csv minimap.sort.1.minus_strand.per.site.5mer.csv plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv minus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv ${params.resultsDir}/${condition1}/epinanoSVM/

           else
               /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -n ${task.cpus} -R genome.fa -b minimap.sort.1.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar
               /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.1.plus_strand.per.site.csv 5
	       /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.1.plus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix plus_mod_prediction

	       mv minimap.sort.1.plus_strand.per.site.csv minimap.sort.1.plus_strand.per.site.5mer.csv plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv  ${params.resultsDir}/${condition1}/epinanoSVM/
           fi
        else
          echo "No reads mapped for minimap.1.sort.bam"
        fi

#        if [[ -s minimap.sort.2.plus.sam ]]; then
#           if [[ -s minimap.sort.2.minus.sam ]]; then
#               /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -n ${task.cpus} -T g -R genome.fa -b minimap.sort.2.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar
#               /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.2.plus_strand.per.site.csv 5
#	       /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.2.minus_strand.per.site.csv 5
#	       /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.2.plus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix plus_mod_prediction
#               /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.2.minus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix minus_mod_prediction

#	       mv minimap.sort.2.plus_strand.per.site.csv minimap.sort.2.minus_strand.per.site.csv minimap.sort.2.plus_strand.per.site.5mer.csv minimap.sort.2.minus_strand.per.site.5mer.csv plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv minus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv ${params.resultsDir}/${condition2}/epinanoSVM/

        #    else
        #        /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -n ${task.cpus} -R genome.fa -b minimap.sort.2.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar
        #        /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.2.plus_strand.per.site.csv 5
	       # /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.2.plus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix plus_mod_prediction

	       # mv minimap.sort.2.plus_strand.per.site.csv minimap.sort.2.plus_strand.per.site.5mer.csv plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv  ${params.resultsDir}/${condition2}/epinanoSVM/
        #    fi
        # else
        #   echo "No reads mapped for minimap.2.sort.bam"
        # fi

    """
	else
	"""
        echo "Skipped"
    """
}

// RNA modifications detection with epinano in Error mode
process epinanoError {
    input:
		tuple val('condition1'), file('minimap.filt.sort.1.bam'), file('minimap.filt.sort.1.bam.bai'), file('minimap.sort.1.bam'), file('minimap.sort.1.bam.bai') from minimap2Merge_1_epinanoError
		tuple val('condition2'), file('minimap.filt.sort.2.bam'), file('minimap.filt.sort.2.bam.bai'), file('minimap.sort.2.bam'), file('minimap.sort.2.bam.bai') from minimap2Merge_2_epinanoError

		each file('genome.fa') from genome_fasta_epinanoError
		each file('genome.fa.fai') from genome_fai_epinanoError

		each file('transcriptome.fa') from transcriptome_fasta_epinanoError
		each file('transcriptome.fa.fai') from transcriptome_fai_epinanoError

		tuple file('genome.fa.dict'), file('transcriptome.fa.dict') from picard_epinanoError
    output:
    	val('flagepinanoError') into epinanoError_postprocessing
    script:
    if(params.epinanoError)
    """
    	mkdir -p ${params.resultsDir}/epinanoError/
		mkdir -p ${params.resultsDir}/epinanoError/minus/
		mkdir -p ${params.resultsDir}/epinanoError/plus/

		samtools view -F16 minimap.sort.1.bam > minimap.sort.1.plus.sam
                samtools view -f16 minimap.sort.1.bam > minimap.sort.1.minus.sam
                samtools view -F16 minimap.sort.2.bam > minimap.sort.2.plus.sam
                samtools view -f16 minimap.sort.2.bam > minimap.sort.2.minus.sam

        if [[ -s minimap.sort.1.plus.sam ]]; then
           if [[ -s minimap.sort.1.minus.sam ]]; then
               /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -R genome.fa -b minimap.sort.1.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar --type g -n ${task.cpus}
               /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.1.plus*site.csv --out minimap.sort.1.plus.sumErrOut.csv
               /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.1.minus*site.csv --out minimap.sort.1.minus.sumErrOut.csv
           else
               /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -R genome.fa -b minimap.sort.1.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar -n ${task.cpus}
               /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.1.plus*site.csv --out minimap.sort.1.plus.sumErrOut.csv
           fi
        else
           "No reads mapped for minimap.sort.1.bam"
        fi

        if [[ -s minimap.sort.2.plus.sam ]]; then
           if [[ -s minimap.sort.2.minus.sam ]]; then
             /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -R genome.fa -b minimap.sort.2.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar --type g -n ${task.cpus}
             /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.2.plus*site.csv --out minimap.sort.2.plus.sumErrOut.csv
             /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.2.minus*site.csv --out minimap.sort.2.minus.sumErrOut.csv
           else
             /usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -R genome.fa -b minimap.sort.2.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar -n ${task.cpus}
             /usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.2.plus*site.csv --out minimap.sort.2.plus.sumErrOut.csv
           fi
        else
           "No reads mapped for minimap.sort.2.bam"
        fi


        if [[ -f minimap.sort.1.plus.sumErrOut.csv && -f minimap.sort.2.plus.sumErrOut.csv ]]; then
            /bin/miniconda3/bin/Rscript /EpiNano-Epinano1.2.1/Epinano_DiffErr.R -k minimap.sort.2.plus.sumErrOut.csv -w minimap.sort.1.plus.sumErrOut.csv -d ${params.epinanoErrorSumErr} -t 3 -p -o diffErr -f sum_err

            if [[ -f diffErr.delta-sum_err.prediction.csv ]]; then
				mv *diffErr* ${params.resultsDir}/epinanoError/plus/
			fi

			mv *plus_strand.per.site.csv ${params.resultsDir}/epinanoError/plus/
			mv *plus.sumErrOut.csv ${params.resultsDir}/epinanoError/plus/
		fi

		if [[ -f minimap.sort.1.minus.sumErrOut.csv && -f minimap.sort.2.minus.sumErrOut.csv ]]; then
			/bin/miniconda3/bin/Rscript /EpiNano-Epinano1.2.1/Epinano_DiffErr.R -k minimap.sort.2.minus.sumErrOut.csv -w minimap.sort.1.minus.sumErrOut.csv -d ${params.epinanoErrorSumErr} -t 3 -p -o diffErr -f sum_err

			if [[ -f diffErr.delta-sum_err.prediction.csv ]]; then
				mv *diffErr* ${params.resultsDir}/epinanoError/minus/
			fi

			mv *minus_strand.per.site.csv ${params.resultsDir}/epinanoError/minus/
			mv *minus.sumErrOut.csv ${params.resultsDir}/epinanoError/minus/
		fi


    """
	else
	"""
        echo "Skipped"
    """
}

// From a single channel for all the alignments to one channel for each condition.
ni_test_nanodoc=Channel.create()
ni_other_nanodoc=Channel.create()
singleReadFAST5_nanodoc.groupTuple(by:0)
	.choice( ni_test_nanodoc, ni_other_nanodoc ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with nanodoc for each condition
process nanodoc {
    input:
	    tuple val('condition1'), val('samples') from ni_test_nanodoc
	    tuple val('condition2'), val('samples') from ni_other_nanodoc

		each file('genome.fa') from genome_fasta_nanodoc
		each file('genome.fa.fai') from genome_fai_nanodoc
		each file('genome.bed') from bed_nanodoc
    output:
    	val('flagnanodoc') into nanodoc_postprocessing
    script:
    if(params.nanodoc)
    """
    	mkdir -p ${params.resultsDir}/nanodoc/
    	mkdir -p ${params.resultsDir}/nanodoc/${condition1}_output/
    	mkdir -p ${params.resultsDir}/nanodoc/${condition2}_output/

		/bin/miniconda3/bin/python /nanoDoc/src/nanoDoc.py formatfile -i ${params.resultsDir}/${condition1}/ -o . -r genome.fa -t ${task.cpus}
                mv index.txt ${params.resultsDir}/nanodoc/${condition1}_output/
                pq_c1=\$(find . -maxdepth 1| grep "\\.pq")
		mv \$pq_c1 ${params.resultsDir}/nanodoc/${condition1}_output/
                /bin/miniconda3/bin/python /nanoDoc/src/nanoDoc.py formatfile -i ${params.resultsDir}/${condition2}/ -o . -r genome.fa -t ${task.cpus}
                mv index.txt ${params.resultsDir}/nanodoc/${condition2}_output/
                pq_c2=\$(find . -maxdepth 1| grep "\\.pq")
                mv \$pq_c2 ${params.resultsDir}/nanodoc/${condition2}_output/
		cat genome.bed | while read line; do chr=\$(echo \$line | cut -d' ' -f1); start=\$(echo \$line | cut -d' ' -f2); end=\$(echo \$line | cut -d' ' -f3); /bin/miniconda3/bin/python /nanoDoc/src/nanoDoc.py analysis -w /nanoDoc/weight5mer/ -p /nanoDoc/param20.txt -r genome.fa -rraw ${params.resultsDir}/nanodoc/${condition2}_output/ -traw ${params.resultsDir}/nanodoc/${condition1}_output/ -chrom \$chr --start \$start --end \$end -o "nanoDoc_results_"\$chr"_"\$start"_"\$end".txt"; done

		mv nanoDoc_results_*.txt ${params.resultsDir}/nanodoc/
    """
	else
	"""
        echo "Skipped"
    """
}

// From a single channel for all the alignments to one channel for each condition.
ni_test_drummer=Channel.create()
ni_other_drummer=Channel.create()
minimap2_drummer.groupTuple(by:0)
	.choice( ni_test_drummer, ni_other_drummer ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with drummer
process drummer {
    input:
	    tuple val('condition1'), file('minimap.filt.sortT.1.*.bam'), file('minimap.filt.sortT.1.*.bam.bai'), file('minimap.sortG.1.*.bam'), file('minimap.sortG.1.*.bam.bai') from ni_test_drummer
	    tuple val('condition2'), file('minimap.filt.sortT.2.*.bam'), file('minimap.filt.sortT.2.*.bam.bai'), file('minimap.sortG.2.*.bam'), file('minimap.sortG.2.*.bam.bai') from ni_other_drummer

    output:
    	val('flagdrummer') into drummer_postprocessing
    script:
    if(params.drummer)
    """
		mkdir -p ${params.resultsDir}/drummer/

		mkdir -p ${params.resultsDir}/drummer/DRUMMER
		rm -rf ${params.resultsDir}/drummer/DRUMMER
		
		cut ${params.genome_fai} -f 1 > ${params.resultsDir}/drummer/chromosomes.txt

		mv 	minimap.sortG*.bam ${params.resultsDir}/drummer/
		mv 	minimap.sortG*.bai ${params.resultsDir}/drummer/

 		cd ${params.resultsDir}/drummer/
		cp -r /DRUMMER .
        cd ${params.resultsDir}/drummer/DRUMMER/

		while read -r line; 
		do
		  python3 DRUMMER.py -r ${params.genome_fasta} \
		  -t ${params.resultsDir}/drummer/minimap.sortG.2.*.bam \
		  -c ${params.resultsDir}/drummer/minimap.sortG.1.*.bam \
		  -o ${params.resultsDir}/drummer/DRUMMER/\$line/ \
		  -a exome \
		  -p ${params.drummerPval} \
		  -n \$line \
                  -m true ;
		done < ${params.resultsDir}/drummer/chromosomes.txt || true
    """
	else
	"""
        echo "Skipped"
    """
}

// Resquigling with nanopolish for each condition
process nanopolish1 {
    input:
		tuple val(condition), val(sample) from minimap2_nanopolish1

		each file('transcriptome.fa') from transcriptome_fasta_nanopolish1
		each file('transcriptome.fa.fai') from transcriptome_fai_nanopolish1

		each file('genome.fa') from genome_fasta_nanopolish1
		each file('genome.fa.fai') from genome_fai_nanopolish1
    output:
    	tuple val(condition), val(sample) into nanopolish1_xpore
    	tuple val(condition), val(sample) into nanopolish1_nanocompore1
    	tuple val(condition), val(sample) into nanopolish1_yanocomp1
    	tuple val(condition), val(sample) into nanopolish1_m6anet1

    script:
    if(params.nanopolish1)
    """
		mkdir -p ${params.resultsDir}/${condition}/${sample}/nanopolish/
		#mkdir -p ${params.resultsDir}/${condition}/${sample}/nanopolish/genome/
		mkdir -p ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/

		/usr/bin/nanopolish/nanopolish index -d ${params.resultsDir}/${condition}/${sample}/FAST5/ ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq

		#/usr/bin/nanopolish/nanopolish eventalign --reads ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq --bam ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.sort.bam --genome genome.fa --samples --signal-index --scale-events -n --summary ${params.resultsDir}/${condition}/${sample}/nanopolish/genome/summary.txt --threads ${task.cpus} > ${params.resultsDir}/${condition}/${sample}/nanopolish/genome/eventalign_readName.txt

		#/usr/bin/nanopolish/nanopolish eventalign --reads ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq --bam ${params.resultsDir}/${condition}/${sample}/genomeAlignment/minimap.sort.bam --genome genome.fa --signal-index --scale-events --summary ${params.resultsDir}/${condition}/${sample}/nanopolish/genome/summary.txt --threads ${task.cpus} > ${params.resultsDir}/${condition}/${sample}/nanopolish/genome/eventalign_readIndex.txt

		/usr/bin/nanopolish/nanopolish eventalign --reads ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq --bam ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam --genome transcriptome.fa --samples --signal-index --scale-events -n --summary ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/summary.txt --threads ${task.cpus} > ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readName.txt

		/usr/bin/nanopolish/nanopolish eventalign --reads ${params.resultsDir}/${condition}/${sample}/FASTQ/singleReadsFASTQ.fastq --bam ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam --genome transcriptome.fa --signal-index --scale-events --summary ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/summary.txt --threads ${task.cpus} > ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readIndex.txt

    """
	else
	"""
        echo "Skipped"
    """
}

// Data formatting for xpore for each sample
process xpore1 {
    input:
    	tuple val(condition), val(sample) from nanopolish1_xpore
    	each file('genome.gtf') from gtf_xpore
		each file('transcriptome.fa') from transcriptome_fasta_xpore
		each file('transcriptome.fa.fai') from transcriptome_fai_xpore

    output:
    	tuple val(condition), val(sample) into xpore1_xpore2
    script:
    if(params.xpore)
    """

        mkdir -p ${params.resultsDir}/${condition}/${sample}/xpore/
        xpore dataprep --eventalign ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readIndex.txt --out_dir ${params.resultsDir}/${condition}/${sample}/xpore --gtf_or_gff genome.gtf --transcript_fasta transcriptome.fa --genome
    """
	else
	"""
        echo "Skipped"
    """
}

// From a single channel for all the alignments to one channel for each condition.
ni_test_xpore2=Channel.create()
ni_other_xpore2=Channel.create()
xpore1_xpore2.groupTuple(by:0)
	.choice( ni_test_xpore2, ni_other_xpore2 ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with xpore
process xpore2 {
    input:
	    tuple val('condition1'), val('sample1') from ni_test_xpore2
	    tuple val('condition2'), val('sample2') from ni_other_xpore2
    output:
	    val('flagxpore') into xpore_postprocessing
    script:
    if(params.xpore)
    """
    	mkdir -p ${params.resultsDir}/xpore/

    	echo "data:" > ${params.resultsDir}/xpore/xpore.yaml
    	echo "    "${condition1}":" >> ${params.resultsDir}/xpore/xpore.yaml
        for file in \$(find ${params.resultsDir}/${condition1} -maxdepth 2 | grep "xpore"); do sn=\$(basename \$(dirname \$file)); sd=\$(dirname \$file); echo "      rep"\$sn": "\$sd"/xpore"; done >> ${params.resultsDir}/xpore/xpore.yaml
        echo "    "${condition2}":" >> ${params.resultsDir}/xpore/xpore.yaml
        for file in \$(find ${params.resultsDir}/${condition2} -maxdepth 2 | grep "xpore"); do sn=\$(basename \$(dirname \$file));  sd=\$(dirname \$file); echo "      rep"\$sn": "\$sd"/xpore"; done >> ${params.resultsDir}/xpore/xpore.yaml
        echo "" >> ${params.resultsDir}/xpore/xpore.yaml
	echo "out: "${params.resultsDir}"/xpore" >> ${params.resultsDir}/xpore/xpore.yaml

	xpore diffmod --config ${params.resultsDir}/xpore/xpore.yaml --n_processes ${task.cpus}

	xpore postprocessing --diffmod_dir ${params.resultsDir}/xpore/
    """
	else
	"""
        echo "Skipped"
    """
}

// Data formatting for nanocompore for each sample
process nanocompore1 {
    input:
	    tuple val('condition'), val('sample') from nanopolish1_nanocompore1

    output:
    	tuple val(condition), val(sample), file("out_eventalign_collapse.tsv"), file("out_eventalign_collapse.tsv.idx") into nanocompore1_nanocompore2

    script:
    if(params.nanocompore1)
    """
        mkdir mkdir -p ${params.resultsDir}/${condition}/${sample}/nanocompore/
		nanocompore eventalign_collapse --eventalign ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readName.txt --nthreads ${task.cpus} --outpath ${params.resultsDir}/${condition}/${sample}/nanocompore/ --overwrite

		ln -sf ${params.resultsDir}/${condition}/${sample}/nanocompore/out_eventalign_collapse.tsv out_eventalign_collapse.tsv
		ln -sf ${params.resultsDir}/${condition}/${sample}/nanocompore/out_eventalign_collapse.tsv.idx out_eventalign_collapse.tsv.idx
    """
	else
	"""
		ln -sf ${params.resultsDir}/${condition}/${sample}/nanocompore/out_eventalign_collapse.tsv out_eventalign_collapse.tsv
		ln -sf ${params.resultsDir}/${condition}/${sample}/nanocompore/out_eventalign_collapse.tsv.idx out_eventalign_collapse.tsv.idx
    """
}


// From a single channel for all the alignments to one channel for each condition
ni_test_nanocompore2=Channel.create()
ni_other_nanocompore2=Channel.create()
nanocompore1_nanocompore2.groupTuple(by:0)
	.choice( ni_test_nanocompore2, ni_other_nanocompore2 ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with nanocompore
process nanocompore2 {
    input:
	    tuple val('condition1'), val('sample1'), file('out_eventalign_collapse.1.*.tsv'), file('out_eventalign_collapse.1.*.tsv.idx') from ni_test_nanocompore2
	    tuple val('condition2'), val('sample2'), file('out_eventalign_collapse.2.*.tsv'), file('out_eventalign_collapse.2.*.tsv.idx') from ni_other_nanocompore2

		each file('transcriptome.fa') from transcriptome_fasta_nanocompore2
		each file('transcriptome.fa.fai') from transcriptome_fai_nanocompore2
		each file('transcriptome.bed') from bed_nanocompore

    output:
    	val('flagnanocompore') into nanocompore_postprocessing
    script:
    if(params.nanocompore2)
    """
		IFS=','
		f1=(out_eventalign_collapse.1.*.tsv)
		f2=(out_eventalign_collapse.2.*.tsv)

    	mkdir -p ${params.resultsDir}/nanocompore/

    	nanocompore sampcomp --file_list1 "\${f1[*]}" --file_list2 "\${f2[*]}" \
		--label1 ${condition1} \
		--label2 ${condition2} \
	    --fasta transcriptome.fa \
	    --bed transcriptome.bed \
    	--outpath ${params.resultsDir}/nanocompore/ \
    	--allow_warnings \
    	--logit \
    	--nthreads ${task.cpus} \
    	--overwrite
    """
	else
	"""
        echo "Skipped"
    """
}

// Data formatting for m6anet for each sample
process m6anet1 {
    input:
	    tuple val('condition'), val('sample') from nanopolish1_m6anet1

    output:
    	tuple val(condition), val(sample), val() into m6anet1_m6anet2


    script:
    if(params.m6anet1)
    """
        mkdir -p ${params.resultsDir}/${condition}/${sample}/m6anet/

        m6anet-dataprep --eventalign  ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readIndex.txt \
                --out_dir ${params.resultsDir}/${condition}/${sample}/m6anet --n_processes ${task.cpus}
    """
	else
	"""
		ln -sf ${params.resultsDir}/${condition}/${sample}/m6anet m6anet
    """
}


// From a single channel for all the alignments to one channel for each condition
ni_test_m6anet2=Channel.create()
ni_other_m6anet2=Channel.create()
m6anet1_m6anet2.groupTuple(by:0)
	.choice( ni_test_m6anet2, ni_other_m6anet2 ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with m6anet
process m6anet2 {
    input:
	    tuple val('condition1'), val('sample1') from ni_test_m6anet2

    output:
    	val('flagm6anet') into m6anet_postprocessing
    script:
    if(params.m6anet2)
    """
        mkdir -p ${params.resultsDir}/m6anet
	preprocessing_dirs=\$(find ${params.resultsDir}/${condition1} -maxdepth 2 -type d | grep "m6anet\$")
        m6anet-run_inference --input_dir \$preprocessing_dirs --out_dir ${params.resultsDir}/m6anet --infer_mod_rate --n_processes ${task.cpus}
	
	    zcat ${params.resultsDir}/m6anet/data.result.csv.gz > ${params.resultsDir}/m6anet/data.result.csv
    """
	else
	"""
        echo "Skipped"
    """
}

// Data formatting for yanocomp for each sample
process yanocomp1 {
    input:
	    tuple val('condition'), val('sample') from nanopolish1_yanocomp1
	    each file('genome.gtf') from gtf_yanocomp
    output:
    	tuple val(condition), file('outputT.hdf5'), file('outputG.hdf5') into yanocomp1_yanocomp2

    script:
    if(params.yanocomp1)
    """
    	mkdir -p ${params.resultsDir}/${condition}/${sample}/yanocomp/
		mkdir -p ${params.resultsDir}/${condition}/${sample}/yanocomp/transcriptome/
		mkdir -p ${params.resultsDir}/${condition}/${sample}/yanocomp/genome/

		/bin/miniconda3/envs/yanocomp/bin/yanocomp prep -p ${task.cpus} -e ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readName.txt -h ${params.resultsDir}/${condition}/${sample}/yanocomp/transcriptome/output.hdf5
 
		/bin/miniconda3/envs/yanocomp/bin/yanocomp prep -p ${task.cpus} -e ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readName.txt -h ${params.resultsDir}/${condition}/${sample}/yanocomp/genome/output.hdf5 -g genome.gtf

		ln -s ${params.resultsDir}/${condition}/${sample}/yanocomp/transcriptome/output.hdf5 outputT.hdf5
		ln -s ${params.resultsDir}/${condition}/${sample}/yanocomp/genome/output.hdf5 outputG.hdf5
    """
	else
	"""
		ln -s ${params.resultsDir}/${condition}/${sample}/yanocomp/transcriptome/output.hdf5 outputT.hdf5
		ln -s ${params.resultsDir}/${condition}/${sample}/yanocomp/genome/output.hdf5 outputG.hdf5
    """
}

// From a single channel for all the alignments to one channel for each condition
ni_test_yanocomp2=Channel.create()
ni_other_yanocomp2=Channel.create()
yanocomp1_yanocomp2.groupTuple(by:0)
	.choice( ni_test_yanocomp2, ni_other_yanocomp2 ) { a -> a[0] == params.test_condition ? 0 : 1 } 

// RNA modifications detection with yanocomp
process yanocomp2 {
    input:
	    tuple val('condition1'), file('outputT.1.*.hdf5'), file('outputG.1.*.hdf5') from ni_test_yanocomp2
	    tuple val('condition2'), file('outputT.2.*.hdf5'), file('outputG.2.*.hdf5') from ni_other_yanocomp2

    output:
    	val('flagyanocomp2') into yanocomp2_postprocessing
    script:
    if(params.yanocomp2)
    """
		mkdir -p ${params.resultsDir}/yanocomp/
		
		/bin/miniconda3/envs/yanocomp/bin/yanocomp gmmtest \
		\$(for file in outputG.1.*.hdf5; do echo -c \$file; done) \
		\$(for file in outputG.2.*.hdf5; do echo -t \$file; done) \
		-p ${task.cpus} \
		-o ${params.resultsDir}/yanocomp/yanocomp_output.bed \
		-s ${params.resultsDir}/yanocomp/yanocomp_output.json.gzip \
		-f ${params.yanocompFDR}
    """
	else
	"""
        echo "Skipped"
    """
}

// Processing of each output to obtain bed files
process postprocessing {
    input:
		val('flagyanocomp2') from yanocomp2_postprocessing
		val('flagdena') from dena_postprocessing
		val('flagdrummer') from drummer_postprocessing
		val('flagdifferr') from differr_postprocessing
		val('flagnanom6a') from nanom6a_postprocessing
		val('flagnanocompore') from nanocompore_postprocessing
		val('flageligos') from eligos_postprocessing
		val('flagmines') from mines_postprocessing
		val('flagepinanoSVM') from epinanoSVM_postprocessing
		val('flagepinanoError') from epinanoError_postprocessing
		val('flagxpore') from xpore_postprocessing
		val('flagnanodoc') from nanodoc_postprocessing
		val('flagtombo2') from tombo2_postprocessing
		val('flagtombo3') from tombo3_postprocessing
		val('flagm6anet') from m6anet_postprocessing

    output:

    script:
    if(params.postprocessing)
    """
	mkdir -p ${params.resultsDir}/output_bed_files/
	mkdir -p ${params.resultsDir}/output_statistical/

	Rscript ${params.postprocessingScript} path=${params.resultsDir} genomebed=${params.genesbed} genomegtf=${params.gtf} resultsFolder=${params.resultsDir}/output_bed_files/ mccores=${task.cpus} threshold=${params.threshold} pathdena=${params.test_condition}/dena pathdrummer=drummer pathdifferr=differr pathyanocomp=yanocomp pathmines=${params.test_condition}/mines pathnanocompore=nanocompore patheligos=eligos/merged pathepinanoError=epinanoError pathepinanoSVM=${params.test_condition}/epinanoSVM pathxpore=xpore pathnanodoc=nanodoc pathnanom6a=${params.test_condition}/nanom6a/result_final pathtomboComparison=tomboComparison pathm6anet=m6anet
        Rscript ${params.statisticalAnalysis} bed_folder=${params.resultsDir}/output_bed_files genomegtf=${params.gtf} genesbed=${params.genesbed} resultsFolder=${params.resultsDir}/output_statistical/ mccores=${task.cpus} peaks=${params.peaksfile} binLength=${params.binLength} genomefile=${params.genome_fasta}

    """
	else
	"""
        echo "Skipped"
    """
}
