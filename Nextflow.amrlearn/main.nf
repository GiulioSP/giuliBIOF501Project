nextflow.enable.dsl=2

params.ref_gbff_path = "inputs/GCF_000020105.1_ref/GCF_000020105.1.gbff"
//params.fasta_dir = "outputs/mapper_results/fasta"
params.vcf = "outputs/mapper_results/parsnp/parsnp.vcf"
params.antibiotics = "inputs/PRJNA776899.antibiotics.txt"
params.project = "PRJNA776899"
params.outdir = "outputs/amrlearn_results"


log.info """\
    SRA 2 VCF ALIGNMENT PIPELINE
    ===================================
    reference                   : ${params.ref_gbff_path}
    vcf path                    : ${params.vcf}
    project name                : ${params.project}
    antibiotic susceptibility   : ${params.antibiotics}
    output directory            : ${params.outdir}
    """
    .stripIndent()

workflow {

    ref_gbff_ch = Channel.value(file(params.ref_gbff_path))
    antibiotics_ch = Channel.value(file(params.antibiotics))
    vcf_ch = Channel.value(file(params.vcf))

    ref_tab = gbff2tab(ref_gbff_ch)

    snp_count = vcf2snp(ref_tab, vcf_ch)

    features = feature2target(snp_count, antibiotics_ch)

    ml_learn(features, ref_tab, antibiotics_ch)
}

process gbff2tab {

    tag "${params.project}"
    publishDir params.outdir, mode: 'copy'

    input:
    path ref_gbff

    output:
    path "*.tab"

    """
    1.gbff2tab.py $ref_gbff ${ref_gbff.baseName}.tab 
    """
}

process vcf2snp {
    tag "$params.project"
    publishDir params.outdir, mode: 'copy'

    input:
    path ref_tab
    path vcf

    output:
    path "*.snp_count.txt"

    """
    2.vcf2snp.py $ref_tab $vcf ${params.project}.snp_count.txt
    """
}

process feature2target {
    tag "$params.project"
    publishDir params.outdir, mode: 'copy'

    input:
    path snp_count
    path antibiotics

    output:
    path "*.feature2target.txt"

    """
    3.feature2target.py $snp_count $antibiotics ${params.project}.feature2target.txt
    """
}

process ml_learn {
    tag "$params.project"
    publishDir params.outdir, mode: 'copy'

    input:
    path features
    path ref_tab
    path antibiotics

    output:
    path "${params.project}.learn/"    

    """
    4.AMR_Learn_linearV2.py $features $ref_tab $antibiotics ${params.project}.learn
    """
}

