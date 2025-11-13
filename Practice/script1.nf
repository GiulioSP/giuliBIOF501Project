params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

//println "reads: $params.reads"

log.info """
reads: $params.reads
transcriptone file: $params.transcriptome_file
multiqc: $params.multiqc
outdir: $params.outdir 
"""

/*
 * define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}