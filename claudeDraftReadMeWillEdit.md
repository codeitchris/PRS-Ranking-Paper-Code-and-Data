# A Literature-Based Benchmarking Database for Polygenic Risk Score Methods

This repository contains the data and code for **"Evaluating Methods for Developing Polygenic Risk Scores with Spectral Ranking Inferences based on Literature Data."** Currently available as a preprint on [medRxiv](https://www.medrxiv.org/content/10.64898/2026.03.01.26347258v1).

We built a literature-derived benchmarking database that pools published head-to-head comparisons of 14 polygenic risk score (PRS) construction methods, drawn from both the original method-development papers and independent applied/benchmarking studies. A spectral ranking algorithm with bootstrap-based uncertainty quantification is then used to turn those scattered pairwise comparisons into an overall ranking (with confidence intervals) of the methods — either across all available data, or restricted to a specific phenotype, a specific subset of methods, or a chosen minimum number of comparisons.

## Repository structure

```
├── R/
│   └── prsRankingCalculation.R   # Spectral ranking algorithm, CI computation, figures/tables
├── python/
│   ├── methodRanking.ipynb       # Builds comparison matrices from the method-development papers
│   └── appliedRankings.ipynb     # Builds comparison matrices from the applied/benchmarking papers and creates combined method and applied dataset
├── data/
│   ├── PRSMethodPapersGitHub.xlsx   # Curated data table: method-development paper comparisons
│   └── PRSPaperAppliedGitHub.xlsx   # Curated data table: applied/benchmarking paper comparisons
└── README.md
```

## The curated data table (supplementary table)

The two Excel workbooks in `data/` **are** the curated literature database — every sheet was hand-extracted from a published paper.

- **`PRSMethodPapersGitHub.xlsx`** — comparisons reported in the original papers that introduced each of the 14 PRS methods.
- **`PRSPaperAppliedGitHub.xlsx`** — comparisons reported in independent applied and benchmarking studies (papers that were *not* introducing a new method, just comparing existing ones). Its `Application and benchmarking` sheet lists every source paper and link.

Each workbook has one reference sheet and then one data sheet per source paper:

| Sheet | Contents |
|---|---|
| `models` (method-papers workbook) / `Method papers` (applied workbook) | The 14 methods, with `Method`, `Other Names`, and a link to the `Method paper` that introduced it. |
| One sheet per paper (e.g. `AUC-DongJunKim`) | First column = phenotype/cohort/dataset label as reported by that paper; remaining columns = one per PRS method, with the performance metric (AUC, R², etc.) that paper reported for that method on that phenotype. A blank cell means that paper didn't compare that method. |

## How the ranking pipeline works

1. **Python (`python/`)** turns each paper's raw performance table into two matrices per comparison:
   - **AA0** — a 0/1 indicator of which two methods were compared in that row.
   - **WW0** — a 0/1 indicator of which of those two methods won (had the better reported metric).

   `appliedRankings.ipynb` additionally normalizes each paper's raw phenotype label onto one of ~85 standardized phenotype names (via the `traitMap` dictionary) and builds a **phenotype-specific** AA0/WW0 pair for each one, stored in a `traits` dictionary.

2. **R (`R/prsRankingCalculation.R`)** takes an AA0/WW0 pair and runs the spectral ranking algorithm (vanilla spectral estimator + a two-stage estimator), with a weighted bootstrap to produce confidence intervals for each method's rank. It also generates the rank plot, head-to-head table, trait-comparison table, and violin plot used in the paper.

## Tutorial: running the ranking analysis

### Requirements

- **Python 3** with `pandas`, `numpy`, `openpyxl`, and Jupyter
- **R** with `readr`, `dplyr`, `ggplot2`, `gt`, `tidyverse`, `ggstatsplot`, `gghalves`, `paletteer`, `colorBlindness`, `scales`

### Step 1 — Build the comparison matrices

Pick the notebook that matches the data source you want to rank:

- `python/methodRanking.ipynb` → reads `data/PRSMethodPapersGitHub.xlsx` → writes `aa0MethodPapers.csv`, `ww0MethodPapers.csv`, `comparisonsMethod.csv`
- `python/appliedRankings.ipynb` → reads `data/PRSPaperAppliedGitHub.xlsx` → writes `aa0AppliedPapers.csv`, `ww0AppliedPapers.csv`, `comparisonsApplied.csv`, `traitCompAppliedDf.csv`, `violin.csv`, and builds the in-memory `traits` dictionary (phenotype-specific matrices)

Open the notebook, replace the `'YOUR PATH TO ...'` placeholder(s) in the `pd.read_excel(...)` cell near the top with the local path to the matching `.xlsx` file, then run all cells.

### Step 2 — Run the spectral ranking in R

Open `R/prsRankingCalculation.R` and set the three `read_csv(file.path("YOUR PATH TO ..."))` calls to point at the CSVs you want to rank — e.g. `aa0MethodPapers.csv` / `ww0MethodPapers.csv` for a ranking over the method-development papers, or `aa0AppliedPapers.csv` / `ww0AppliedPapers.csv` for the applied/benchmarking papers. Then source the script.

The result is a matrix `RR2` with one column per method and these rows:

| Row | Meaning |
|---|---|
| 1–2 | `theta.hat`, rank (vanilla spectral method) |
| 3–4 | Two-sided 95% CI for rank (left, right) |
| 5 | Left-sided CI for rank |
| 6 | Uniform left-sided CI for rank |
| 7–12 | Same five quantities, using the two-stage theta estimator |

The script also produces a ranked forest plot with CIs, a head-to-head `gt` table, a trait-comparison count table, and (for the applied data) a violin plot of normalized per-phenotype ranks.

### Customizing your ranking

#### Restrict to specific methods

There's no command-line flag for this yet — filter the comparison matrices before ranking:

1. Drop the column(s) for the method(s) you want to exclude from `AA0`/`WW0`.
2. Drop any row where that method had a `1` in `AA0`, since that comparison no longer makes sense without it. `appliedRankings.ipynb` already has a helper for exactly this — `cleanUp(AA, WW, col)` (used internally to drop under-threshold methods; you can call it on any method name you want removed).
3. Save the filtered `AA0`/`WW0` to CSV and point the R script at those files.

#### Rank a specific phenotype

Phenotype-specific ranking is currently only available for the **applied/benchmarking** dataset — `methodRanking.ipynb` doesn't split by phenotype, since it draws from the 14 method-development papers rather than independent applications.

In `appliedRankings.ipynb`, `traits[<phenotype name>]` holds `[AA0, WW0]` already restricted to that phenotype (see `traitMap` for the full mapping from each paper's raw phenotype label — e.g. `"BrCa-UKB"` — to its standardized name — e.g. `"Breast Cancer"`).

- **Quick ranking, no CIs** (matches the paper's violin plot): call `rankingCalculation(traits["<Phenotype>"][0], traits["<Phenotype>"][1], "<Phenotype>")`.
- **Full ranking with confidence intervals**: export that phenotype's matrices — `traits["<Phenotype>"][0].to_csv("aa0_<phenotype>.csv", index=False)` and the same for `[1]` as `ww0_<phenotype>.csv` — then point the R script's `aCSV`/`wCSV` reads at those two files.

#### Set the minimum number of comparisons required

This threshold lives inside `rankingCalculation()` in `appliedRankings.ipynb`:

```python
if sum(AA[x]) < 10:
    removes.append(x)
```

Any method with fewer than this many total comparisons (within the phenotype/dataset being ranked) is dropped before ranking, since too few comparisons can't support a reliable rank. To change the minimum, edit the `10` and re-run the notebook to regenerate the CSVs before ranking in R. The method-development ranking (`methodRanking.ipynb` / R script directly) doesn't use a separate threshold, since every method in that pool already comes from its own dedicated method paper.
