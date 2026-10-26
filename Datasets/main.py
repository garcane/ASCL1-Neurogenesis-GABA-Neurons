"""Rebuild the GSE29985 / GSE46791 workbooks from their GEO series-matrix files.

Usage: python main.py <annotation_dir> [datasets_dir]
<annotation_dir> holds GPL4128.txt, GPL4129.txt, GPL6885.txt (GEO SOFT "full" downloads) and
Mus_musculus.gene_info.gz (NCBI). [datasets_dir] defaults to the folder containing this script and must
hold one folder per dataset whose name starts with the GSE number (e.g. "GSE46791 An epigenetic ...");
each folder holds that dataset's *_series_matrix.txt files and receives the rebuilt <GSE>.xlsx.
See main.md for background.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

DATA = Path(__file__).parent  # overridden by the optional datasets_dir argument
P_CUTOFF = 0.0009  # reproduces the old GSE29985 list (min -log10 P = 3.047)

SYMBOL_COLS = ["GENE_SYMBOL", "Symbol", "Gene Symbol", "Gene symbol"]
TITLE_COLS = ["GENE_NAME", "Definition", "Gene Title", "Gene title"]

# (series-matrix files, annotation platforms, group A (reference), group B)
DATASETS = {
    "GSE46791": {
        "matrices": ["GSE46791_series_matrix.txt"],
        "platforms": ["GPL6885"],
        "A": [f"GSM{n}" for n in range(1138448, 1138453)],  # NSC
        "B": [f"GSM{n}" for n in range(1138453, 1138456)],  # Neuron
    },
    "GSE29985": {
        "matrices": [
            "GSE29985-GPL4128_series_matrix.txt",
            "GSE29985-GPL4129_series_matrix.txt",
        ],
        "platforms": ["GPL4128", "GPL4129"],
        "B": ["GSM742100", "GSM742101", "GSM742102", "GSM742106", "GSM742107", "GSM742108"],  # ARX N2a cells
        "A": ["GSM742103", "GSM742104", "GSM742105", "GSM742109", "GSM742110", "GSM742111"],  # forebrain
    },
}


def read_table(path, begin, end):
    """Read the tab-separated table between two marker lines of a GEO text file."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = []
        inside = False
        for line in fh:
            if not inside:
                inside = line.startswith(begin)
            elif line.startswith(end):
                break
            else:
                lines.append(line)
    from io import StringIO

    return pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str, quotechar='"', low_memory=False)


def trigamma_inverse(x):
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = polygamma(1, y)
        step = tri * (1 - tri / x) / polygamma(2, y)
        y += step
        if abs(step) < 1e-10:
            break
    return y


def moderated_t(a, b):
    """limma-style two-group moderated t-test (what GEO2R reports). Returns (log2FC, p-value) Series."""
    n1, n2 = a.shape[1], b.shape[1]
    df = n1 + n2 - 2
    s2 = (a.var(axis=1, ddof=1) * (n1 - 1) + b.var(axis=1, ddof=1) * (n2 - 1)) / df
    beta = b.mean(axis=1) - a.mean(axis=1)
    # empirical-Bayes prior variance (limma squeezeVar / fitFDist)
    e = np.log(s2[s2 > 0]) - digamma(df / 2) + np.log(df / 2)
    evar = ((e - e.mean()) ** 2).sum() / (len(e) - 1) - polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * trigamma_inverse(evar)
        s02 = np.exp(e.mean() + digamma(d0 / 2) - np.log(d0 / 2))
        post = (d0 * s02 + df * s2) / (d0 + df)
    else:
        d0, post = np.inf, pd.Series(np.exp(e.mean()), index=s2.index)
    t = beta / np.sqrt(post * (1 / n1 + 1 / n2))
    return beta, pd.Series(2 * stats.t.sf(np.abs(t), d0 + df), index=s2.index)


def load_expression(matrix):
    expr = read_table(matrix, "!series_matrix_table_begin", "!series_matrix_table_end").set_index("ID_REF")
    return expr.apply(pd.to_numeric, errors="coerce")


def load_annotation(ann_dir, gpl):
    t = read_table(Path(ann_dir) / f"{gpl}.txt", "!platform_table_begin", "!platform_table_end")
    sym = next(c for c in SYMBOL_COLS if c in t.columns)
    title = next(c for c in TITLE_COLS if c in t.columns)
    ann = t[["ID", sym, title]].set_axis(["ID", "Gene.symbol", "Gene.title"], axis=1)
    # Illumina "Definition" reads "Mus musculus <title> (<symbol>), mRNA." -> keep just <title>
    ann["Gene.title"] = ann["Gene.title"].str.replace(r"^Mus musculus\s+|\s*\([^()]*\),\s*[^,]*$", "", regex=True)
    return ann.drop_duplicates("ID").set_index("ID")


def load_gene_titles(ann_dir):
    """NCBI mouse gene_info -> {symbol: title}, covering official symbols and unambiguous old synonyms."""
    gi = pd.read_csv(Path(ann_dir) / "Mus_musculus.gene_info.gz", sep="	", dtype=str, usecols=["Symbol", "Synonyms", "description"])
    official = gi.drop_duplicates("Symbol").set_index("Symbol")["description"]
    syn = gi.assign(syn=gi["Synonyms"].str.split("|")).explode("syn")
    syn = syn[syn["syn"].notna() & (syn["syn"] != "-") & ~syn["syn"].isin(official.index)]
    syn = syn.drop_duplicates(["syn", "Symbol"]).drop_duplicates("syn", keep=False)
    return pd.concat([official, syn.set_index("syn")["description"]])


def fill_titles(res, ann_dir):
    """Fill missing Gene.title from NCBI gene_info (platform symbols are kept as published)."""
    titles = load_gene_titles(ann_dir)
    missing = res["Gene.title"].isna() | (res["Gene.title"] == "")
    res.loc[missing, "Gene.title"] = res.loc[missing, "Gene.symbol"].map(titles)
    return res


def analyse_platform(matrix, gpl, groups, ann_dir):
    expr = load_expression(matrix)
    a = expr[[c for c in groups["A"] if c in expr.columns]]
    b = expr[[c for c in groups["B"] if c in expr.columns]]
    log2fc, pval = moderated_t(a, b)
    res = pd.DataFrame({"log2(fold change)": log2fc, "p": pval}).dropna()
    res = res.join(load_annotation(ann_dir, gpl), how="inner")
    res = res[res["Gene.symbol"].notna() & (res["Gene.symbol"].str.strip() != "")]
    res = res[res["p"] < P_CUTOFF].copy()
    res["=-LOG10(Pvalue)"] = (-np.log10(res["p"].clip(lower=1e-300))).round(3)
    res["log2(fold change)"] = res["log2(fold change)"].round(3)
    res = fill_titles(res.reset_index(names="ID"), ann_dir)
    if res["ID"].str.fullmatch(r"\d+").all():  # GPL4128/4129 use numeric probe IDs; keep them numeric
        res["ID"] = res["ID"].astype(int)
    return res[["ID", "Gene.symbol", "Gene.title", "log2(fold change)", "=-LOG10(Pvalue)"]]


def dataset_dir(name):
    """The folder under DATA whose name starts with the GSE number."""
    return next(p for p in sorted(DATA.iterdir()) if p.is_dir() and p.name.startswith(name))


def analyse(name, cfg, ann_dir):
    folder = dataset_dir(name)
    parts = [
        analyse_platform(folder / m, gpl, cfg, ann_dir) for m, gpl in zip(cfg["matrices"], cfg["platforms"])
    ]
    return pd.concat(parts, ignore_index=True)


def write_workbook(name, table):
    lfc = "log2(fold change)"
    path = dataset_dir(name) / f"{name}.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        table.to_excel(xw, sheet_name="Annotate Gene List", index=False)
        table[table[lfc] > 0].to_excel(xw, sheet_name="Upstream", index=False)
        table[table[lfc] < 0].to_excel(xw, sheet_name="Downstream", index=False)
        table[["Gene.symbol"]].to_excel(xw, sheet_name="Gene Overlap", index=False, header=False)
        for ws in list(xw.book.worksheets)[:3]:
            ws["E1"].data_type = "s"  # keep the "=-LOG10(Pvalue)" header as text, not a formula
    print(f"{name}: {len(table)} rows, up {int((table[lfc] > 0).sum())}, down {int((table[lfc] < 0).sum())} -> {path}")


def main():
    global DATA
    ann_dir = sys.argv[1]
    if len(sys.argv) > 2:
        DATA = Path(sys.argv[2])
    for name, cfg in DATASETS.items():
        write_workbook(name, analyse(name, cfg, ann_dir))


if __name__ == "__main__":
    main()
