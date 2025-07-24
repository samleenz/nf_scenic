# pySCENIC Nextflow Pipeline

This repository contains a Nextflow DSL2 pipeline to run a complete pySCENIC workflow inside an Apptainer container. The pipeline covers the following steps:

- **GRN Inference:** Identifies gene regulatory networks using `pyscenic grn`.
- **Motif Enrichment/Context:** Performs motif enrichment with `pyscenic ctx`.
- **Regulon Activity Scoring:** Scores regulon activity using `pyscenic aucell`.

> **Note:** This repository is private.

## Prerequisites

- [Nextflow](https://www.nextflow.io/) installed.
- [Apptainer](https://apptainer.org/) installed.
- An Apptainer container image with pySCENIC and all necessary dependencies.
- Input data files as defined in the pipeline (`data/expression_matrix.csv`, transcription factor list, and motif databases).

## Repository Structure

```plaintext
.
├── main.nf          # Nextflow pipeline script
├── README.md        # This README file
└── databases        # motif, and TF files for SCENIC
```

## Getting Started

1. **Clone the Repository**  
   Ensure you have access to the private repo:
   ```bash
   git clone git@github.com:samleenz/nf_scenic.git
   cd nf_scenic
   ```

2. **Prepare the Input Data**  
   Place your required files:
   - Expression matrix: `data/expression_data.loom`
   - Transcription factors: `databases/hs_hgnc_tfs.txt`
   - Gene motif rankings: `databases/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
   - Motif annotations: `/databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`

3. **Run the Pipeline**  
   Launch the pipeline with:
   ```bash
   nextflow run main.nf
   ```
   Adjust parameters either by editing `main.nf` or by passing options via the command line.

## Configuration

- **Container:**  
  The pipeline uses an Apptainer container defined in `params.container`. Update the path to your container image as needed.

- **Resource Allocation:**  
  CPU and runtime specifications are set in each process (GRN, CTX, AUCell). Modify these values according to your environment.

- **Reproducibility:**  
  A fixed seed (`myseed = 12345`) is used for reproducibility. Change the seed if needed.

## Additional Notes

- Verify that the Apptainer container has pySCENIC properly installed and configured.
- Adapt input file paths and parameters to suit your dataset and specific analysis requirements.


