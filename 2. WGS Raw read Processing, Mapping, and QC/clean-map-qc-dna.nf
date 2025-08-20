#!/usr/bin/env nextflow

// Create a channel from input FASTQ R1 files by inferring R2 and sample name
Channel
  .fromPath(params.input)
  // ensure we only keep _L<lane>_1.fq.gz
  .filter { it.name ==~ /.*_L\d+_1\.fq\.gz$/ }
  .map { read1 ->
      // infer R2 by reusing the same lane captured in R1
      def read2 = file(read1.toString().replaceFirst(/_L(\d+)_1\.fq\.gz$/, '_L$1_2.fq.gz'))

      // sample is everything before "_CKDL"
      def basename = read1.getBaseName()                         // e.g., N18_CKDL250016671-1A_22V5GKLT4_L4_1
      def m = (basename =~ /^(.*?)_CKDL/)
      def sample = m ? m[0][1] : basename.replaceFirst(/_L\d+_[12]$/, '')  // fallback if no CKDL block

      tuple(sample, read1, read2)
  }
  .set { read_pairs }

// FASTP - Trims, filters, merges reads, and generates quality reports
process Fastp {
    tag "$sample"
    cpus 6

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_fastp_R1.fastq.gz"),
          path("${sample}_fastp_R2.fastq.gz"),
          path("${sample}_fastp_ME.fastq.gz"),
          path("${sample}_fastp_SE.fastq.gz"),
          path("${sample}_fastp.json"),
          path("${sample}_fastp.html"),
          path("${sample}_fastp.log")

    // Save outputs to appropriate result directories
    publishDir 'results/fastp_cleaned_reads', mode: 'copy', pattern: "*fastp_{R1,R2,ME,SE}.fastq.gz"
    publishDir 'results/logs/fastp', mode: 'copy', pattern: "*fastp.{json,html,log}"

    script:
    """
    fastp \
        -i $read1 -I $read2 \
        -m \
        --merged_out ${sample}_fastp_ME.fastq.gz \
        --out1 ${sample}_fastp_R1.fastq.gz \
        --out2 ${sample}_fastp_R2.fastq.gz \
        --unpaired1 ${sample}_fastp_SE.fastq.gz \
        --unpaired2 ${sample}_fastp_SE.fastq.gz \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_front_window_size 5 \
        --cut_front_mean_quality 20 \
        --correction \
        -l 30 \
        --thread ${task.cpus} \
        -j ${sample}_fastp.json \
        -h ${sample}_fastp.html \
        2> ${sample}_fastp.log
    """
}

// MultiQC_Fastp - Aggregates Fastp reports using MultiQC
process MultiQC_Fastp {
    tag "multiqc_fastp"
    cpus 1

    input:
    path fastp_logs_files

    output:
    path "multiqc_fastp_report.html"

    publishDir "results/multiqc/fastp", mode: 'copy'

    script:
    """
    multiqc $fastp_logs_files --filename multiqc_fastp_report.html --outdir .
    """
}

// AlignReads - Aligns reads using BWA, merges BAMs, adds read groups, and removes duplicates
process AlignReads {
    tag "$sample"
    cpus 8

    input:
    tuple val(sample),
          path(r1), path(r2), path(me), path(se),
          path(json), path(html), path(log)

    output:
    tuple val(sample),
          path("${sample}_deduped.bam"),
          path("${sample}_deduped.bam.bai"),
          path("${sample}_dedupe_metrics.txt")

    publishDir 'results/bams', mode: 'copy', pattern: "*.{bam,bai}"
    publishDir 'results/logs/picard', mode: 'copy', pattern: "*_dedupe_metrics.txt"

    script:
    """
    RG_LINE=\$(zcat $r1 | head -n 1)
    RGID=\$(echo \$RG_LINE | awk -F ':' '{print \$3}')
    PU=\$(echo \$RG_LINE | awk -F ':' '{print \$4}')

    echo "Extracted RGID: \$RGID"
    echo "Extracted PU: \$PU"
    
    bwa mem -M -t ${task.cpus} ${params.ref} $r1 $r2 | samtools view -Sb - > ${sample}_PE.bam
    bwa mem -M -t ${task.cpus} ${params.ref} $me | samtools view -Sb - > ${sample}_ME.bam
    bwa mem -M -t ${task.cpus} ${params.ref} $se | samtools view -Sb - > ${sample}_SE.bam

    samtools merge ${sample}_merged.bam ${sample}_PE.bam ${sample}_ME.bam ${sample}_SE.bam

    picard -Xmx50g AddOrReplaceReadGroups \
        I=${sample}_merged.bam \
        O=${sample}_merge_addRG.bam \
        SO=coordinate \
        RGID=\$RGID \
        LB=tamias_wgs \
        PL=illumina \
        PU=\$PU \
        SM=${sample} \
        TMP_DIR=${params.tmp}

    picard -Xmx50g MarkDuplicates \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        I=${sample}_merge_addRG.bam \
        O=${sample}_deduped.bam \
        METRICS_FILE=${sample}_dedupe_metrics.txt \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=${params.tmp}

    samtools index ${sample}_deduped.bam

    rm ${sample}_PE.bam ${sample}_ME.bam ${sample}_SE.bam
    rm ${sample}_merged.bam ${sample}_merge_addRG.bam
    """
}

// QualimapQC - Performs quality control on BAM files using Qualimap
process QualimapQC {
    tag "$sample"
    cpus 8
    memory '25 GB'

    input:
    tuple val(sample), path(bam), path(bai), path(metrics)

    output:
    path("${sample}_bamqc") // QC output directory

    publishDir 'results/qualimap', mode: 'copy'

    script:
    """
    mkdir ${sample}_bamqc
    qualimap bamqc --java-mem-size=23G -bam $bam -outdir ${sample}_bamqc -outformat HTML -nt ${task.cpus}
    """
}

// MultiQC_Qualimap - Aggregates Qualimap results using MultiQC
process MultiQC_Qualimap {
    tag "multiqc_qualimap"
    cpus 1

    input:
    path qualimap_dirs

    output:
    path "multiqc_qualimap_report.html" // Combined report

    publishDir "results/multiqc/qualimap", mode: 'copy'

    script:
    """
    multiqc $qualimap_dirs --filename multiqc_qualimap_report.html --outdir .
    """
}

// WORKFLOW block - Orchestrates the execution of processes
workflow {
    fastp_results = read_pairs | Fastp          // Run Fastp
    aligned = fastp_results | AlignReads        // Align reads and remove duplicates
    qualimap_results = aligned | QualimapQC     // Run QC

    // Prepare logs and run MultiQC for Fastp
    fastp_logs = fastp_results.flatMap { sample, r1, r2, me, se, json, html, log -> [json, html, log] }
    fastp_logs_files = fastp_logs.collect()
    MultiQC_Fastp(fastp_logs_files)

    // Run MultiQC for Qualimap results
    qualimap_dirs = qualimap_results.collect()
    MultiQC_Qualimap(qualimap_dirs)
}

// Cleanup work directory on successful completion
workflow.onComplete {
    if (workflow.success) {
        println "Workflow succeeded. Cleaning up work directory..."
        def workDir = file(workDir ?: './work')
        if (workDir.exists()) {
            workDir.listFiles().each { it.deleteDir() }
            println "Deleted contents of work directory: ${workDir}"
        } else {
            println "Work directory not found: ${workDir}"
        }
    } else {
        println "Workflow failed. Work directory preserved for debugging."
    }
}


