nextflow.enable.dsl=2

params.ref_fasta = "../data/GCF_000020105.1_ref/GCF_000020105.1_ASM2010v1_genomic.fna"
params.ref_mmi = "../data/GCF_000020105.1_ref/GCF_000020105.1.mmi"
params.sra_path = "../data/sra/*.sra"
params.outdir = "../aligner_results"

log.info """\
    SRA 2 FASTA ALIGNMENT PIPELINE
    ===================================
    reference   : ${params.ref_fasta}
    sra path    : ${params.sra_path}
    outdir      : ${params.outdir}
    """
    .stripIndent()

workflow {

    ref_fasta = Channel.value(file(params.ref_fasta))
    ref_mmi = Channel.value(file(params.ref_mmi))

    sra_files = Channel.fromPath(params.sra_path).map {
        file -> def sra_id = file.baseName
        tuple(sra_id, file)
        }

    read_pairs = dump_reads(sra_files)
    
    sam_files = align(read_pairs, ref_mmi)
    sorted_bams = sort_bam(sam_files)
    indexed_bams = index_bam(sorted_bams)

    generate_consensus(indexed_bams, ref_fasta).view()

    //aligned_bams = extract_aligned(indexed_bams)
    //fastqs = bam_to_fastq(aligned_bams)
    //fastq_to_fasta(fastqs).view()

    multiqc(fastqc(read_pairs).collect())
}
    

process build_index {
    tag "index"
    publishDir params.outdir, mode: 'copy'

    //container 'quay.io/biocontainers/minimap2:2.26--h5bf99c6_0' 

    input:
    path ref_fasta

    output:
    path "GCF_000020105.1.mmi"

    """
    minimap2 -d GCF_000020105.1.mmi $ref_fasta
    """
}

process dump_reads {
    tag { sra_id }

    //container 'ncbi/sra-tools:2.11.2' 

    input:
    tuple val(sra_id), path(file)

    output:
    tuple val(sra_id), path("*_1.fastq"), path("*_2.fastq")

    """
    fasterq-dump --split-files $file
    """
}

process fastqc {
    tag { sra_id }
    publishDir "${params.outdir}/reports", mode: 'copy' 
    
    //container 'staphb/fastqc:0.11.9' 

    input:
    tuple val(sra_id), path(r1), path(r2)

    output:
    //tuple val(sra_id), path("*_fastqc.*")
    path("*_fastqc.*")

    """
    fastqc $r1 $r2
    """
}

process multiqc {
    publishDir "${params.outdir}/reports", mode: 'copy'

    //container 'quay.io/biocontainers/multiqc:1.9--py_1'

    input:
    path fastqc_files 

    output:
    path "multiqc_report.html"

    """
    multiqc $fastqc_files -o .
    """
}

process align {
    tag { sra_id }

    //container 'quay.io/biocontainers/minimap2:2.26--h5bf99c6_0' 

    input:
    tuple val(sra_id), path(r1), path(r2)
    path ref_mmi

    output:
    tuple val(sra_id), path("*.sam")
    
    """
    minimap2 -ax sr $ref_mmi $r1 $r2 > ${sra_id}.sam
    """
}

process sort_bam {
    tag { sra_id }

    input:
    tuple val(sra_id), path(sam_file)

    output:
    tuple val(sra_id), path("*.sorted.bam")

    """
    samtools view -bS $sam_file | samtools sort -o ${sam_file.baseName}.sorted.bam
    """
}

process index_bam {
    tag { sra_id }

    input:
    tuple val(sra_id), path(bam_file)

    output:
    tuple val(sra_id), path(bam_file), path("*.bai")

    """
    samtools index $bam_file
    """
}

process extract_aligned {
    tag { sra_id }

    input:
    tuple val(sra_id), path(bam_file), path(bai_file)

    output:
    tuple val(sra_id), path("*.aligned.bam")

    """
    samtools view -b -F 4 $bam_file > ${bam_file.baseName}.aligned.bam
    """
}

process bam_to_fastq {
    tag { sra_id }

    input:
    tuple val(sra_id), path(bam_file)

    output:
    tuple val(sra_id), path("*.fastq")

    """
    samtools fastq $bam_file > ${bam_file.baseName}.fastq
    """
}

process fastq_to_fasta {
    tag { sra_id }
    publishDir "${params.outdir}/fasta", mode: 'copy'

    input:
    tuple val(sra_id), path(fastq_file)
    
    output:
    path "*.fasta"

    """
    seqkit fq2fa $fastq_file > ${sra_id}.fasta
    """
}

process generate_consensus {
    tag { sra_id }
    publishDir "${params.outdir}/fasta", mode: 'copy'

    input:
    // Takes the tuple from index_bam
    tuple val(sra_id), path(bam_file), path(bai_file)
    path(ref_fasta)

    output:
    path "*.consensus.fasta"

    script:
    """
    # 1. Generate VCF file for variants/consensus calls
    bcftools mpileup -Ou -f $ref_fasta $bam_file | bcftools call -c -Oz -o calls.vcf.gz
    bcftools index calls.vcf.gz

    # 2. Generate the consensus sequence (replace N for low coverage regions)
    cat $ref_fasta | bcftools consensus calls.vcf.gz > ${sra_id}.consensus.fasta

    """
}


