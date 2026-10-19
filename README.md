# Duck circovirus recombination analysis

Code and results for the recombination analysis of 444 complete duck circovirus genomes (GenBank, retrieved 18 October 2026).

## Contents
- `step1_fetch.py` ... `step28_remap.py`: analysis scripts, run in numeric order from the repository root (each script's header or code shows its inputs and outputs). The scripts use relative paths.
- `gardrun.sh`, `gard_mpi*.sh`, `step16_hyphy.sh`, `seq3/run_3seq.sh`: launchers for GARD, HyPhy and 3SEQ.
- `tree_aln.fasta`: the 444-genome alignment (1,980 columns). `setA_aln_raw.fasta`: the untrimmed MAFFT alignment.
- `events_table_HM.tsv`: the 14 extreme genomes with parents, breakpoints (HM162346.1 numbering) and mosaic fit.
- `recomb_candidates.tsv`, `gard_breakpoints.tsv`, `triplets.txt`, `phi_test.txt`, `phi_validation.txt`, `seq3/`: screening and confirmation outputs.
- `fig_breakpoint_map_v3.png/.pdf`: breakpoint map.

## Software
MAFFT 7.526, IQ-TREE 3.1.4, HyPhy 2.5.103, 3SEQ 1.8.0, FastTree 2.2.0, Biopython 1.88.

## Data provenance
Sequences are public GenBank records; accession lists are in `setA_ids.txt` and `uids.json`.

## Citation
See `CITATION.cff`. Archived release DOI: to be added after the Zenodo release.

## Archive and citation

Archived on Zenodo: [10.5281/zenodo.23463655](https://doi.org/10.5281/zenodo.23463655)

If you use these data or scripts, please cite this archive (see `CITATION.cff`).

## Licence

Code: MIT (`LICENSE`). Data, results and figures: CC BY 4.0 (`LICENSE-DATA.md`).
