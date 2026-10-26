# main.py – rebuilding the GSE29985 and GSE46791 gene lists

## The problem
The `GSE29985.xlsx` and `GSE46791.xlsx` workbooks in the `Datasets` folders held the same table under two different names. Both were the GSE29985 list, so GSE46791 had no analysis of its own. In both, the Upstream and Downstream sheets were unfiltered copies of the full list.

`main.py` rebuilds both workbooks from each dataset's own GEO series-matrix `.txt` files, so each file reflects its own data.

## What it produces
Each `<GSE>.xlsx` has four sheets:

| Sheet | Content |
|---|---|
| Annotate Gene List | `ID`, `Gene.symbol`, `Gene.title`, `log2(fold change)`, `=-LOG10(Pvalue)` for every probe passing the cutoff |
| Upstream | rows of the list with `log2(fold change)` **> 0** |
| Downstream | rows of the list with `log2(fold change)` **< 0** |
| Gene Overlap | `Gene.symbol` column of the annotated list (no header) |

The `=-LOG10(Pvalue)` header is stored as text, not as a formula.

## Inputs
- The series-matrix files in each dataset folder:
  - GSE29985: `GSE29985-GPL4128_series_matrix.txt` and `GSE29985-GPL4129_series_matrix.txt`
  - GSE46791: `GSE46791_series_matrix.txt`
- Platform annotation from GEO: `GPL4128.txt`, `GPL4129.txt`, `GPL6885.txt` (SOFT "full" view, which contains the probe → gene table).
- `Mus_musculus.gene_info.gz` from NCBI, used only to fill missing gene titles.

The series-matrix files hold only probe IDs and values. Gene symbols and titles have to come from the platform annotation.

## Method
1. Read the expression table from each series matrix.
2. Compare two groups per probe (**log2(fold change) = mean(B) − mean(A)**):

   | Dataset | A (reference) | B |
   |---|---|---|
   | GSE46791 | NSC (GSM1138448–52) | Neuron (GSM1138453–55) |
   | GSE29985 | Forebrain (GSM742103–05, GSM742109–11) | ARX-transfected N2a cells (GSM742100–02, GSM742106–08) |

3. P-values come from a limma-style moderated t-test (empirical-Bayes variance shrinkage), implemented in `moderated_t`. This is the statistic GEO2R reports.
4. Join the platform annotation, drop probes without a gene symbol, and keep probes with **p < 0.0009** (−log10 P ≥ 3.05).
5. Fill missing `Gene.title` values from the NCBI gene info (official symbol, otherwise an unambiguous old synonym). Platform gene symbols are kept as published.
6. Write the four sheets.

GSE29985 is analysed per array (GPL4128 and GPL4129 each have 3 + 3 samples) and the two results are concatenated. Probe IDs are therefore not unique across the whole sheet: the same number can appear once per array.

## Validation
On GPL4128 alone, the code reproduces the original GSE29985 GEO2R list. Fold changes agree to within 0.0005 and −log10 P values to three decimals, which is how the cutoff and the moderated t-test were identified. 4,228 of the original 4,264 probes were recovered; the rest have no gene symbol in the platform annotation.

## Running it
```
uv add pandas openpyxl scipy        # or: pip install pandas openpyxl scipy
python main.py <annotation_dir> [datasets_dir]
```
`datasets_dir` defaults to the folder holding `main.py`. It must contain one folder per dataset whose name starts with the GSE number. The rebuilt `<GSE>.xlsx` is written into that folder, replacing the existing file.

## Known limitations
- The comparison direction is a choice. Reversing A and B swaps Upstream and Downstream.
- Some `Gene.title` values are blank (about 420 in GSE29985 and 110 in GSE46791), where no NCBI title could be matched.
- Titles come from current NCBI records, so wording can differ from older GEO2R output. Symbols follow the platform and may be outdated for some genes.
- "Upstream" and "Downstream" here follow the sign of the log2 fold change, as defined for this project.
