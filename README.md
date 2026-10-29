# Hybrid Drosophila genome dynamics

Code used for the ancestry and forward-simulation analyses in the D. simulans x D. sechellia experimental evolution study.

## Files

```
R/
  01_marker_ancestry.R          ancestry-marker time series and effective selection parameters
  02_illumina_ancestry.R        generation-46 Illumina ancestry and repeatability analyses
  03_simulation_sensitivity.R       drift, recombination, population-size and sampling sensitivity
  04_genome_hitchhiking.R       optional sparse genome-wide linked-selection analysis

slim/
  make_parental_state_13MB.slim
  fly_hybrid_3R_fullhistory.slim
  fly_hybrid_genome_hitchhiking.slim

bash/
  make_parental_state.sh
  run_main_sims.sh
  run_architecture_sims.sh
  run_genome_hitchhiking.sh
  run_all.sh

data/marker/
  fly_anc-pos_simellia.csv
  fly_anc-pos_sechelulans.csv
```

Generated files are written to `results/` and are not tracked by git.

## Requirements

R packages:

- data.table
- ggplot2
- reshape2
- RColorBrewer
- zoo
- patchwork

Forward simulations use SLiM. Shell scripts assume `slim` is on the PATH; another executable can be supplied with `SLIM=/path/to/slim`.

## Marker ancestry analysis

The marker-frequency input files are included in `data/marker`.

```bash
Rscript R/01_marker_ancestry.R
```

The effective marker-level selection parameter is calculated as the per-generation change in log allele odds:

`m = [logit(pn) - logit(p0)] / n`

This is an effective ancestry-frequency parameter and is not interpreted as the genotype-level fitness coefficient of a causal mutation.

## Illumina ancestry analysis

The large read-count input is not included in this repository. Run with the path to `readcount_by_bottle_Illumina_data`:

```bash
Rscript R/02_illumina_ancestry.R /path/to/readcount_by_bottle_Illumina_data
```

The script includes the parental-SNP QC, 300-kb and 30-kb ancestry windows, bottle-level ancestry, 3R/chr4 peak summaries, repeatability scans and correlation analyses.

## Forward simulations

The full-history model uses one shared parental state, followed by an explicit F1-female x D. simulans backcross. Generation 0 is the post-backcross bottle and therefore starts at approximately 25% D. sechellia ancestry. Recombination occurs in females only.

Generate the parental state once:

```bash
bash bash/make_parental_state.sh
```

Run the main simulation sensitivity analyses:

```bash
bash bash/run_main_sims.sh
```

This runs:

- 1,000 neutral replicates at N=50
- recombination multipliers 0, 0.25, 0.5, 1 and 2 under the AE20 model
- population-size sensitivity at N=10, 20, 30 and 50

The model used for the recombination and population-size sensitivity has linked deleterious effects at the two flanking loci (`s=-0.20`) and beneficial effects at 5.9 and 12.9 Mb (`s=0.145` and `0.078`). These are simulation parameters used to test the linked-selection/recombination mechanism, not direct causal estimates for the marker sites.

Summarize the simulations:

```bash
Rscript R/03_simulation_sensitivity.R
```

This also reconstructs the historical one-/five-female genotyping scheme from the `obs_*` columns and compares sampled frequencies with the true simulated bottle frequencies.

## Linkage-architecture diagnostics

The architecture grid used while evaluating the reviewer alternative can be rerun with:

```bash
bash bash/run_architecture_sims.sh
```

It includes the original five-locus model, the dominance diagnostic, corrected positive-effect parameters, one- and two-deleterious-locus arrangements, and the AE linked-load sensitivity. These simulations are intended to test mechanism and parameter dependence rather than fit a unique causal architecture.

## Optional genome-wide sparse-selection sensitivity

```bash
bash bash/run_genome_hitchhiking.sh 1000
Rscript R/04_genome_hitchhiking.R
```

This model asks whether a small number of deleterious introgressed loci can drive broad loss of otherwise neutral ancestry through linkage.

## Generation labels

The simulation variable `generation` counts generations after construction of the post-backcross bottle, with the post-backcross state as generation 0. In the historical sampling schedule used for the `obs_*` columns, empirical N2 therefore maps to simulation generation 1, N11 to generation 10, and N45 to generation 44.
