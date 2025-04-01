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
    time "4h"

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
    time "30m"


    input:
    tuple val(seed), path(grn_file)
    path(expr)
    path(gene_motifs)
    path(motifs)

    output:
    path(output)

    script:
    output = "ctx_output_${seed}.tsv"
    """
    pyscenic ctx \
        ${grn_file} \
        ${gene_motifs} \
        --annotations_fname ${motifs} \
        --expression_mtx_fname ${expr}\
        -o ${output} \
        --mask_dropouts \
        --num_workers ${task.cpus}
    """
}

// Process to generate HC regulons from the ctx output
process HCRegulons {
    container params.container_r
    cpus 1
    memory "32G"
    time "30m"

    input:
    path(ctx_files)

    output:
    path(output)

    script:
    output = "hc_regulons.gmt"
    """
    ./combine_ctx.R ${ctx_files}}
    """
}

// Process for scoring regulon activity using pySCENIC aucell
process AUCell {
    container params.container
    cpus 12
    memory "100.G"
    time "30m"
    publishDir "results/${params.project}/${seed}", mode: 'copy'

    input:
    path(ctx_file)
    path(expr)
    
    output:
    path(output)

    script:
    output = "auc_output.loom"
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
    .randomSample( params.nruns, params.myseed )

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
    
    hc_ch = HCRegulons(
        ctx_ch.collect()
    ) // Generate HC regulons from the CTX output

    auc_ch = AUCell(
        hc_ch,
        file(params.expr)
    ) // Run regulon activity scoring using the expression matrix and HC regulons
}
