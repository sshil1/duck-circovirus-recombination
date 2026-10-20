# Duck circovirus recombination analysis

Code and results for the recombination analysis of 444 complete duck circovirus genomes (GenBank, retrieved 18 October 2026).

## Contents
- `step1_fetch.py` ... `step28_remap.py`: analysis scripts, run in numeric order from the repository root (each script's header or code shows its inputs and outputs). The scripts use relative paths.
- `gardrun.sh`, `gard_mpi*.sh`, `step16_hyphy.sh`, `seq3/run_3seq.sh`: launchers for GARD, HyPhy and 3SEQ.
- `tree_aln.fasta`: the 444-genome alignment (1,980 columns). `setA_aln_raw.fasta`: the untrimmed MAFFT alignment.
- `events_table_HM.tsv`: the 14 extreme genomes with parents, breakpoints (HM162346.1 numbering) and mosaic fit.
- `recomb_candidates.tsv`, `gard_breakpoints.tsv`, `triplets.txt`, `phi_test.txt`, `phi_validation.txt`, `seq3/`: screening and confirmation outputs.
- `fig_breakpoint_map_v3.png/.pdf`: breakpoint map.

- `additional_analyses/`: scripts, logs and summary tables for the group 1 monophyly and subgroup tests, replication-origin motif, recurrence permutation test, homoplasy-index negative control and selection sensitivity analysis (v1.0.2; see below).

## Software
MAFFT 7.526, IQ-TREE 3.1.4, HyPhy 2.5.103, 3SEQ 1.8.0, FastTree 2.2.0, Biopython 1.88.

## Additional analyses (v1.0.2)

These scripts and logs support the monophyly, subgroup, replication-origin, recurrence and selection-sensitivity analyses added in the revised manuscript. Run all server scripts from `~/dcv` inside `conda activate cavphylo` (IQ-TREE 3, HyPhy 2.5.103); inputs are `tree_aln.fasta`, `recomb_candidates.tsv` and the clustering files already in the repository.

| Script | Purpose | Main output |
|---|---|---|
| `additional_analyses/scripts/g1_monophyly.sh` | Builds IQ-TREE trees for arc I (alignment columns 80–570, TIM2+R3) and arc O (columns 730–1940, TIM3+F+R4). Its midpoint-rooted test is superseded by `g1_post.sh`. | `g1_out/` trees |
| `additional_analyses/scripts/g1_post.sh` | Unrooted monophyly test of the ten group 1 genomes, then constrained tree and AU test (10,000 RELL resamplings). | `g1_post.log` |
| `additional_analyses/scripts/g1_subgroups.sh` | Constrained AU tests for Henan 6, Fujian 4 and cluster C in both arcs. | `g1_out/subgroups_summary.tsv`, `g1_subgroups.log` |
| `additional_analyses/scripts/g1_clusters_AB.sh` | Same test for clusters A and B. | `g1_out/clustersAB_summary.tsv`, `g1_clusters_AB.log` |
| `additional_analyses/scripts/origin_motif.py` | NANTATTAC matches (both strands), stem-loops, and breakpoint-to-nonamer distances against 20,000 random placements. | `origin_motif_result.txt` |
| `additional_analyses/scripts/hotspot_test.py` | Breakpoint recurrence permutation test (nine events; uniform, informative-site and rigid-arc-rotation nulls). | stdout |
| `additional_analyses/scripts/hotspot_merged.py` | Same test with the two Henan clusters merged (eight events). | stdout |
| `additional_analyses/scripts/phi_negctl.py` | Recombination-free simulations for the pairwise homoplasy index. | `negctl_summary.txt` |
| `additional_analyses/scripts/sel_noflag.sh` | Selection analysis after removing all 59 flagged genomes (FUBAR, FEL, MEME, BUSTED). | `sel_noflag/summary.txt` |
| `additional_analyses/scripts/step17_summary.py` | Summary helper used by `sel_noflag.sh`, with the path changed to `~/dcv/sel_noflag`. | |

Headline results (see the logs for the exact output): the ten group 1 genomes are not monophyletic in either arc (AU p = 8.5×10⁻³⁴ and 1.5×10⁻⁴⁹); Henan 6 is not rejected in arc I (p = 0.19) but is in arc O (p = 3.4×10⁻⁶⁴); Fujian 4 is rejected in both arcs; clusters A, B and C are each consistent with a single event; the single NANTATTAC match lies at position 600; and five of the seven multiply supported selected sites remain after removing the 59 flagged genomes.

Checksums for every file in `additional_analyses/` are in `additional_analyses/MANIFEST.sha256`. Verify with `cd additional_analyses && sha256sum -c MANIFEST.sha256`. The reference list checked against PubMed is in `additional_analyses/docs/references_verified.tsv`.

## Data provenance
Sequences are public GenBank records; accession lists are in `setA_ids.txt` and `uids.json`.

## Citation
See `CITATION.cff`. Archived release DOI: to be added after the Zenodo release.

## Archive and citation

Archived on Zenodo: [10.5281/zenodo.23463655](https://doi.org/10.5281/zenodo.23463655)

If you use these data or scripts, please cite this archive (see `CITATION.cff`).

## Licence

Code: MIT (`LICENSE`). Data, results and figures: CC BY 4.0 (`LICENSE-DATA.md`).
