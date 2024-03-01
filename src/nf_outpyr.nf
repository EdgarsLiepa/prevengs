#!/usr/bin/env nextflow

params.input_file = '/home/ubuntu/projects/bioinformatics/prevengs/data/feature_table_BKUS.csv'
params.output_dir = '/home/ubuntu/projects/bioinformatics/prevengs/outpyr_out'
params.jobs = 10
params.preprocessing = 'normalize_sf'


process run_outpyr {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path(input_file)

    """
    # install this while container not fixed
    apt-get install python3-tk -y
    
    echo "Running outpyr"
    
    outpyr -j ${params.jobs} -p ${params.preprocessing} -n ${input_file}
    """
}

workflow {
    input = Channel.fromPath(params.input_file)
    run_outpyr(input)
}
