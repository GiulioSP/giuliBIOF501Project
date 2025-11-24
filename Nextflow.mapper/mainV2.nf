nextflow.enable.dsl=2

params.ref_fasta = "inputs/GCF_000020105.1_ref/GCF_000020105.1_ASM2010v1_genomic.fna"
params.ref_mmi = "inputs/GCF_000020105.1_ref/GCF_000020105.1.mmi"
params.sra_path = "inputs/sra_list.txt"
params.outdir = "outputs/mapper_results"

log.info """\
    SRA 2 sraFASTA Mapping Pipeline
    ===================================
    reference   : ${params.ref_fasta}
    sra path    : ${params.sra_path}
    outdir      : ${params.outdir}
    """
    .stripIndent()

workflow {

    ref_fasta = Channel.value(file(params.ref_fasta))
    ref_mmi = Channel.value(file(params.ref_mmi))

//    sra_files = Channel.fromPath(params.sra_path).map {
//        file -> def sra_id = file.baseName
//        tuple(sra_id, file)
//        }

  //  sra_id = Channel.fromPath(params.sra_path) //.splitText().map{it.trim()}.filter{it.length()>0}

    sra_id = Channel.fromPath(file(params.sra_path)).splitText().map{ it.trim() }.filter{ it.length()>0 }

    sra_files = prefetch_sra(sra_id)

    read_pairs = dump_reads(sra_files)
    
    sam_files = align(read_pairs, ref_mmi)
    
    sorted_bams = sort_bam(sam_files)
    
    indexed_bams = index_bam(sorted_bams)

    generate_consensus(indexed_bams, ref_fasta).view()

    multiqc(fastqc(read_pairs).collect())
}

 process prefetch_sra {
    tag { sra_id }
    //publishDir params.outdir, mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("*/*.sra")

    script:
    """
    prefetch ${sra_id}
    """
}

process dump_reads {
    tag { sra_id }
    //publishDir params.outdir, mode: 'copy'

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

    input:
    tuple val(sra_id), path(r1), path(r2)

    output:
    path("*_fastqc.*")

    """
    fastqc $r1 $r2
    """
}

process multiqc {
    publishDir "${params.outdir}/reports", mode: 'copy'

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
    //publishDir params.outdir, mode: 'copy'

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
    //publishDir params.outdir, mode: 'copy'

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
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sra_id), path(bam_file)

    output:
    tuple val(sra_id), path(bam_file), path("*.bai")

    """
    samtools index $bam_file
    """
}

process generate_consensus {
    tag { sra_id }
    publishDir "${params.outdir}/fasta", mode: 'copy'

    input:
    tuple val(sra_id), path(bam_file), path(bai_file)
    path(ref_fasta)

    output:
    path "*.consensus.fasta"

    script:
    """
    bcftools mpileup -Ou -f $ref_fasta $bam_file | bcftools call -c -Oz -o calls.vcf.gz
    bcftools index calls.vcf.gz
    cat $ref_fasta | bcftools consensus calls.vcf.gz > ${sra_id}.consensus.fasta
    """
}
