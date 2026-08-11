#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { loadAlignmentInputs } from './lib/inputs'

/*
 * Nextflow port of alignment.smk — same commands, relative paths, and naming
 * as the Snakemake workflow. Run from the project directory that contains
 * fastq.yml, sample.list, and bam_input/.
 */

workflow {
    def parsed = loadAlignmentInputs(
        params.fastq_config,
        params.sample_list,
        params.reference_key
    )
    def samples = parsed.samples
    def lanes = parsed.lanes

    new File('logs/cluster/aln').mkdirs()
    new File('logs/cluster/metrics').mkdirs()
    samples.each { new File("logs/cluster/${it}").mkdirs() }

    ch_lanes = Channel.fromList(lanes)

    BWA_MEM(ch_lanes)
    SAMTOOLS_READGROUP(BWA_MEM.out.lane)
    SAMTOOLS_MARKDUP(SAMTOOLS_READGROUP.out.lane)

    SAMTOOLS_MARKDUP.out.lane
        .map { sample, reference, run, lane, index, markdup ->
            tuple(sample, reference, markdup.toString())
        }
        .groupTuple(by: [0, 1])
        .map { sample, reference, markdups ->
            tuple(sample, reference, markdups.sort())
        }
        .set { ch_merge }

    INPUT_READY(ch_merge)
    VALIDATE_SAMFILE(INPUT_READY.out.ready)
    READY_BAM(VALIDATE_SAMFILE.out.validated)
    WRITE_BAM_TABLE(
        Channel.value(samples),
        Channel.value(params.reference_key),
        READY_BAM.out.collect()
    )
}

process BWA_MEM {
    tag "${sample}:${reference}:${run}.${lane}.${index}"
    cpus 4
    memory '32 GB'

    input:
    tuple val(sample), val(reference), val(run), val(lane), val(index), val(r1), val(r2)

    output:
    tuple val(sample), val(reference), val(run), val(lane), val(index),
          path("bam_input/work/${sample}/${reference}/${run}/${lane}/${index}/1.mapped.bam"), emit: lane

    script:
    def outdir = "bam_input/work/${sample}/${reference}/${run}/${lane}/${index}"
    """
    mkdir -p ${outdir}
    bwa mem -M -t ${task.cpus} ${params.reference_fasta} ${r1} ${r2} | samtools view -bS -o bam_input/work/${sample}/${reference}/${run}/${lane}/${index}/1.mapped.bam
    """
}

process SAMTOOLS_READGROUP {
    tag "${sample}:${reference}:${run}.${lane}.${index}"
    cpus 4

    input:
    tuple val(sample), val(reference), val(run), val(lane), val(index), path(mapped)

    output:
    tuple val(sample), val(reference), val(run), val(lane), val(index),
          path("bam_input/work/${sample}/${reference}/${run}/${lane}/${index}/4.qsort.bam"), emit: lane

    script:
    def base = "bam_input/work/${sample}/${reference}/${run}/${lane}/${index}"
    """
    samtools addreplacerg -@ ${task.cpus} \\
        -r 'ID:${run}.${lane}' \\
        -r 'PU:${run}.${lane}.${index}' \\
        -r 'PL:illumina' \\
        -r 'LB:${params.library_key}' \\
        -r 'SM:${sample}' \\
        -o ${base}/2.readgroup.bam ${mapped}
    samtools fixmate -m -@ ${task.cpus} ${base}/2.readgroup.bam ${base}/3.fixmate.bam
    samtools sort -@ ${task.cpus} -o ${base}/4.qsort.bam ${base}/3.fixmate.bam
    """
}

process SAMTOOLS_MARKDUP {
    tag "${sample}:${reference}:${run}.${lane}.${index}"
    cpus 4

    input:
    tuple val(sample), val(reference), val(run), val(lane), val(index), path(qsort)

    output:
    tuple val(sample), val(reference), val(run), val(lane), val(index),
          path("bam_input/work/${sample}/${reference}/${run}/${lane}/${index}/5.markdup.bam"), emit: lane

    script:
    def base = "bam_input/work/${sample}/${reference}/${run}/${lane}/${index}"
    """
    samtools markdup -s -f ${base}/5.stats.txt -@ ${task.cpus} ${qsort} ${base}/5.markdup.bam
    """
}

process INPUT_READY {
    tag "${sample}:${reference}"
    cpus 4

    input:
    tuple val(sample), val(reference), path(markdup_bams)

    output:
    tuple val(sample), val(reference),
          path("bam_input/work/${sample}/${reference}/input.bam"), emit: ready

    script:
    def bamList = markdup_bams.join(' ')
    """
    samtools merge -f -@ ${task.cpus} bam_input/work/${sample}/${reference}/input.bam ${bamList}
    samtools index bam_input/work/${sample}/${reference}/input.bam
    """
}

process VALIDATE_SAMFILE {
    tag "${sample}:${reference}"

    input:
    tuple val(sample), val(reference), path(input_bam)

    output:
    tuple val(sample), val(reference), path(input_bam),
          path("metrics/${reference}/${sample}/validation_data.table"), emit: validated

    script:
    """
    set +e
    mkdir -p metrics/${reference}/${sample}
    java -Xmx10240m -jar ${params.picard_jar} ValidateSamFile I=${input_bam} O=metrics/${reference}/${sample}/validation_data.table MODE=SUMMARY
    exitcode=\$?
    if [ \$exitcode -ne 0 ]; then exit 1; fi
    exit 0
    """
}

process READY_BAM {
    tag "${sample}:${reference}"

    input:
    tuple val(sample), val(reference), path(input_bam), path(validation_table)

    output:
    path("bam_input/final/${sample}/${sample}.${reference}.bam")
    path("bam_input/final/${sample}/${sample}.${reference}.bam.bai")

    script:
    """
    rsync -v ${input_bam} bam_input/final/${sample}/${sample}.${reference}.bam
    samtools index bam_input/final/${sample}/${sample}.${reference}.bam
    """
}

process WRITE_BAM_TABLE {
    input:
    val samples
    val ref_key
    val _done

    output:
    path('bam.table')

    script:
    def lines = samples.collect { s -> "${s}\tbam_input/final/${s}/${s}.${ref_key}.bam" }.join('\n')
    """
    cat <<'EOF' > bam.table
${lines}
EOF
    """
}
