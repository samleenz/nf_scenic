#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* This pipeline executes the three pySCENIC steps:
*  1. GRN Inference (pySCENIC grn)
*  2. Motif Enrichment / Context (pySCENIC ctx)
*  3. Regulon Activity Scoring (pySCENIC aucell)
*
* It uses an Apptainer container specified by the 'params.container' parameter.
* Modify the input files and command options to suit your dataset and pySCENIC version.
*/


// Process for GRN Inference using pySCENIC grn
process GRN {
    container params.container
    cpus 16
    memory "100.G"
    time "3h"

    input:
    path(expr)
    path(genes)
    val(seed)

    output:
    tuple val(seed), path("grn_output.tsv")

    script:
    """
    pyscenic grn \
        -o grn_output.tsv \
        --num_workers ${task.cpus} \
        --method grnboost2 \
        --seed ${seed} \
        ${expr} \
        ${genes} 
    """
}

// Process for motif enrichment / context using pySCENIC ctx
process CTX {
    container params.container
    cpus 16
    memory "100.G"
    time "3h"


    input:
    tuple val(seed), path(grn_file)
    path(expr)
    path(gene_motifs)
    path(motifs)

    output:
    tuple val(seed), path("ctx_output.tsv")

    script:
    """
    pyscenic ctx \
        ${grn_file} \
        ${gene_motifs} \
        --annotations_fname ${motifs} \
        --expression_mtx_fname ${expr}\
        -o ctx_output.tsv \
        --mask_dropouts \
        --num_workers ${task.cpus}
    """
}

// Process for scoring regulon activity using pySCENIC aucell
process AUCell {
    container params.container
    cpus 12
    memory "100.G"
    time "3h"
    publishDir "results/${params.project}/${seed}", mode: 'copy'

    input:
    tuple val(seed), path(ctx_file)
    path(expr)
    
    output:
    path(output)

    script:
    output = "auc_output.tsv"
    """
    pyscenic aucell \
        ${expr} \
        ${ctx_file} \
        -o ${output} \
        --num_workers ${task.cpus}
    """
}

// // Generate 50 random seeds
seeds_ch = Channel.of( 1000..9999 )
    .randomSample( ${params.nruns} )

// Define the workflow by chaining the processes
workflow {
    grn_ch = GRN(
        file(params.expr),
        file(params.genes),
        seeds_ch
    )   // Run GRN inference
    
    ctx_ch = CTX(
        grn_ch, 
        file(params.expr),
        file(params.genes_motifs), 
        file(params.motifs)
    )   // Run motif enrichment using the GRN output
    
    auc_ch = AUCell(
        ctx_ch,
        file(params.expr)
    ) // Run regulon activity scoring using the expression matrix and CTX output
}
