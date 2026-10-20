## Additional analyses (v1.0.2)

These scripts and logs support the monophyly, subgroup, replication-origin, recurrence and selection-sensitivity analyses added in the revised manuscript. Run all server scripts from `~/dcv` inside `conda activate cavphylo` (IQ-TREE 3, HyPhy 2.5.103); inputs are `tree_aln.fasta`, `recomb_candidates.tsv` and the clustering files already in the repository.

| Script | Purpose | Main output |
|---|---|---|
| `scripts/g1_monophyly.sh` | Builds IQ-TREE trees for arc I (alignment columns 80–570, TIM2+R3) and arc O (columns 730–1940, TIM3+F+R4). Its midpoint-rooted test is superseded by `g1_post.sh`. | `g1_out/` trees |
| `scripts/g1_post.sh` | Unrooted monophyly test of the ten group 1 genomes, then constrained tree and AU test (10,000 RELL resamplings). | `g1_post.log` |
| `scripts/g1_subgroups.sh` | Constrained AU tests for Henan 6, Fujian 4 and cluster C in both arcs. | `g1_out/subgroups_summary.tsv`, `g1_subgroups.log` |
| `scripts/g1_clusters_AB.sh` | Same test for clusters A and B. | `g1_out/clustersAB_summary.tsv`, `g1_clusters_AB.log` |
| `scripts/origin_motif.py` | NANTATTAC matches (both strands), stem-loops, and breakpoint-to-nonamer distances against 20,000 random placements. | `origin_motif_result.txt` |
| `scripts/hotspot_test.py` | Breakpoint recurrence permutation test (nine events; uniform, informative-site and rigid-arc-rotation nulls). | stdout |
| `scripts/hotspot_merged.py` | Same test with the two Henan clusters merged (eight events). | stdout |
| `scripts/phi_negctl.py` | Recombination-free simulations for the pairwise homoplasy index. | `negctl_summary.txt` |
| `scripts/sel_noflag.sh` | Selection analysis after removing all 59 flagged genomes (FUBAR, FEL, MEME, BUSTED). | `sel_noflag/summary.txt` |
| `scripts/step17_summary.py` | Summary helper used by `sel_noflag.sh`, with the path changed to `~/dcv/sel_noflag`. | |

Headline results (see the logs for the exact output): the ten group 1 genomes are not monophyletic in either arc (AU p = 8.5×10⁻³⁴ and 1.5×10⁻⁴⁹); Henan 6 is not rejected in arc I (p = 0.19) but is in arc O (p = 3.4×10⁻⁶⁴); Fujian 4 is rejected in both arcs; clusters A, B and C are each consistent with a single event; the single NANTATTAC match lies at position 600; and five of the seven multiply supported selected sites remain after removing the 59 flagged genomes.

Checksums for every file in this folder are in `MANIFEST.sha256`. Verify with `sha256sum -c MANIFEST.sha256`.
