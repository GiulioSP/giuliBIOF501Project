nextflow.enable.dsl=2

params.ref_gbff_path = "inputs/GCF_000020105.1_ref/GCF_000020105.1.gbff"
params.fasta_dir = "outputs/mapper_results/fasta"
params.antibiotics = "inputs/PRJNA776899.antibiotics.txt"
params.threshold = 0.08 //threshold for filtering absolute coefficient
params.project = "PRJNA776899"
params.outdir = "outputs/amrlearn_results"


log.info """\
    SRA 2 FASTA ALIGNMENT PIPELINE
    ===================================
    reference                   : ${params.ref_gbff_path}
    fasta path                  : ${params.fasta_dir}
    project name                : ${params.project}
    antibiotic susceptibility   : ${params.antibiotics}
    threshold for coefficients  : ${params.threshold}
    output directory            : ${params.outdir}
    """
    .stripIndent()

workflow {

    ref_gbff_ch = Channel.value(file(params.ref_gbff_path))
    antibiotics_ch = Channel.value(file(params.antibiotics))
    fasta_ch = Channel.value(file(params.fasta_dir))

    ref_tab = gbff2tab(ref_gbff_ch)

    vcf = parsnp_vcf(ref_gbff_ch, fasta_ch)

    snp_count = vcf2snp(ref_tab, vcf)

    features = feature2target(snp_count, antibiotics_ch)

    ml_learn(features, ref_tab, antibiotics_ch, params.threshold)
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

process parsnp_vcf {
    tag "$params.project"
    publishDir params.outdir, mode: 'copy'

    input:
    path ref_gbff_ch
    path fasta_ch
    
    output:
    path "*.parsnp/*.vcf"

    //-p 6
    //-o ${params.outdir}/parsnp 
    """
    parsnp -g $ref_gbff_ch -d $fasta_ch/ -o ${params.project}.parsnp --vcf
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
    val threshold

    output:
    path "${params.project}.learn/"    

    """
    4.AMR_Learn_linearV2.py $features $ref_tab $antibiotics ${params.project}.learn $threshold
    """
}

