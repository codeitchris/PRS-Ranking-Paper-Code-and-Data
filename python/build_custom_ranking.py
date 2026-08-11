"""
build_custom_ranking.py

Builds a customized AA0/WW0 comparison-matrix pair (ready to feed into
R/prsRankingCalculation.R) from the curated literature workbooks, filtered by:

  --methods       a subset of PRS methods to include
  --phenotype     a single standardized phenotype to rank (applied data only)
  --min-comparisons  drop any method with fewer than this many comparisons
                      in the selected data (default: 10, matching the
                      threshold used in appliedRankings.ipynb)

This reproduces the matrix-building logic from python/methodRanking.ipynb and
python/appliedRankings.ipynb (getInfo / getInfoCohort / traitMap / cleanUp),
generalized into a single reusable, parameterized script instead of requiring
hand edits to the notebooks each time.

USAGE
-----
Rank all methods, all applied-paper data, default threshold:
    python build_custom_ranking.py --data-type applied \
        --workbook ../data/PRSPaperAppliedGitHub.xlsx --outdir ./out

Rank only a subset of methods:
    python build_custom_ranking.py --data-type applied \
        --workbook ../data/PRSPaperAppliedGitHub.xlsx \
        --methods "LDpred2,PRS-CS,SBayesR,lassosum2" --outdir ./out

Rank a single phenotype (applied data only):
    python build_custom_ranking.py --data-type applied \
        --workbook ../data/PRSPaperAppliedGitHub.xlsx \
        --phenotype "Breast Cancer" --outdir ./out

Change the minimum-comparisons threshold:
    python build_custom_ranking.py --data-type applied \
        --workbook ../data/PRSPaperAppliedGitHub.xlsx \
        --min-comparisons 5 --outdir ./out

Method-development paper data (no phenotype splitting available):
    python build_custom_ranking.py --data-type method \
        --workbook ../data/PRSMethodPapersGitHub.xlsx --outdir ./out

Each run writes AA0.csv and WW0.csv to --outdir. Point
R/prsRankingCalculation.R's `aCSV`/`wCSV` reads at those two files to get
ranks and confidence intervals for exactly the slice of data you asked for.
"""

import argparse
import sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Phenotype-label normalization map, copied verbatim from
# python/appliedRankings.ipynb (and python/methodRanking.ipynb, which
# defines the identical dictionary). Only used when --phenotype is given.
# ---------------------------------------------------------------------------
TRAIT_MAP = {
    'AD': "Alzheimer's Disease", 'ADHD': "Attention-deficit/hyperactivity disorder",
    'ADNI123-M1': "Alzheimer's Disease", 'ADNI123-M2': "Alzheimer's Disease",
    'ADNI123-M3': "Alzheimer's Disease", 'ADNI123-M4': "Alzheimer's Disease",
    'ADNI123-M5': "Alzheimer's Disease", 'AFRBMI': "Body Mass Index",
    'AFRHDL': "High-Density Lipoprotein", 'AFRHeight': "Height",
    'AFRLDL': "Low-Density Lipoprotein", 'AFRTC': "Cholesterol",
    'AFRlogTG': "Triglyceride", 'AM': "Age at Menarche", 'AN': "Angina",
    'AS': "Asthma Disease", 'ATH': "Asthma Disease",
    'African AmericanAny CVD': "Cardio Vascular Disease",
    'African AmericanDepression': "Depression",
    'African AmericanHeart metabolic disease burden': "Heart Metabolic Disease Burden",
    'African AmericanHeight': "Height",
    'African AmericanMigraine Diagnosis': "Migraine Diagnosis",
    'African AmericanMorning Person': "Morning Person",
    'African AmericanSBMN': "Sing Back Musical Note",
    'Age-related macular degeneration': "Age Related Macular Degeneration",
    'Alzheimer\u2019s disease': "Alzheimer's Disease", 'Asthma': "Asthma Disease",
    'Asthma-BBJ': "Asthma Disease", 'Asthma-Consortium': "Asthma Disease",
    'Asthma-FinnGen': "Asthma Disease", 'Asthma-UKB': "Asthma Disease",
    'Asthma-UKBB': "Asthma Disease", 'Atrial fibrillation': "Atrial Fibrillation",
    'BFP': "Body Fat Percentage", 'BIP': "Bipolar Disorder",
    'BMD': "Bone Mineral Density", 'BMI': "Body Mass Index",
    'BMR': "Basal Metabolic Rate", 'BRCA': "Breast Cancer", 'BW': "Birth Weight",
    'BioFINDER-M1': "Alzheimer's Disease", 'BioFINDER-M2': "Alzheimer's Disease",
    'BioFINDER-M3': "Alzheimer's Disease", 'BioFINDER-M4': "Alzheimer's Disease",
    'BioFINDER-M5': "Alzheimer's Disease", 'Bipolar disorder': "Bipolar Disorder",
    'Bladder Cancer': "Bladder Cancer", 'Bowel cancer': "Bowel Cancer",
    'BrCa': "Breast Cancer", 'BrCa-BBJ': "Breast Cancer",
    'BrCa-Consortium': "Breast Cancer", 'BrCa-FinnGen': "Breast Cancer",
    'BrCa-UKB': "Breast Cancer", 'BrCa-UKBB': "Breast Cancer",
    'Breast Cancer': "Breast Cancer", 'Breast cancer': "Breast Cancer",
    'CAD': "Coronary Artery Disease", 'CAD-BBJ': "Coronary Artery Disease",
    'CAD-Consortium': "Coronary Artery Disease", 'CAD-FinnGen': "Coronary Artery Disease",
    'CAD-UKB': "Coronary Artery Disease", 'CAD-UKBB': "Coronary Artery Disease",
    'CKD': "Chronic Kidney Disease", 'Cardiovascular disease': "Cardio Vascular Disease",
    'Cataract-BBJ': "Cataracts", 'Cataract-UKB': "Cataracts",
    'Celiac Disease': "Celiac Diesease", 'Celiac disease': "Celiac Diesease",
    'Colorectal Cancer': "Colorectal Cancer",
    'Coronary artery disease': "Coronary Artery Disease",
    'Crohns Disease': "Crohn's Disease", 'Crohn\u2019s disease': "Crohn's Disease",
    'DBP': "Diastolic Blood Pressure", 'DBT': "Diabetes", 'DEP': "Depression",
    'DFI': "Dried Fruit Intake", 'Depression': "Depression",
    'Distolic Blood Pressure': "Diastolic Blood Pressure",
    'EASHDL': "High-Density Lipoprotein", 'EASLDL': "Low-Density Lipoprotein",
    'EASTC': "Cholesterol", 'EASlogTG': "Triglyceride", 'EC': "Eosinophil Count",
    'EOS': "Eosinophil Count", 'ES': "Ever Smoked",
    'East AsianAny CVD': "Cardio Vascular Disease", 'East AsianDepression': "Depression",
    'East AsianHeart metabolic disease burden': "Heart Metabolic Disease Burden",
    'East AsianHeight': "Height", 'East AsianMigraine Diagnosis': "Migraine Diagnosis",
    'East AsianMorning Person': "Morning Person", 'East AsianSBMN': "Sing Back Musical Note",
    'Epithelial ovarian cancer': "Epithelia Ovarian Cancer",
    'FEV': "Forced Expiratory Volume", 'FFI': "Fresh Fruit Intake",
    'FFR': "(Forced Expiratory Volume in 1 second)/(Forced Vital Capacity)",
    'FVC': "Forced Vital Capacity", 'GCSE': "General Certificate of Secondary Education",
    'GO': "Gout", 'Gastric Cancer-BBJ': "Gastric Cancer",
    'Gastric Cancer-UKB': "Gastric Cancer", 'Glaucoma-BBJ': "Glaucoma",
    'Glaucoma-UKB': "Glaucoma", 'Gout': "Gout", 'HA': "Headache",
    'HBP': "High Blood Pressure", 'HC': "Hip Circumference",
    'HDL': "High-Density Lipoprotein", 'HDL Cholesterol': "High-Density Lipoprotein",
    'HDL-UKBB': "High-Density Lipoprotein", 'HDL-constoritum': "High-Density Lipoprotein",
    'HP': "Haematocrit Percentage", 'HT': "Height", 'HbA1c': "Hemoglobin A1c",
    'Height': "Height", 'Height-UKBB': "Height", 'Height-consortium': "Height",
    'Hypercholesterolemia': "Hypercholesterolemia", 'Hypertension': "Hypertension",
    'Hyperthyroidism-BBJ': "Hyperthyroidism", 'Hypertriglyceridemia': "Hypertriglyceridemia",
    'Hypothyroidism-BBJ': "Hypothyroidism", 'Hypothyroidism-UKB': "Hypothyroidism",
    'Hyptertension': "Hypertension", 'IBD': "Inflammatory Bowel Disease",
    'Intelligence': "Intelligence", 'Ischemic stroke': "Ischemic Stroke",
    'Kidney Cancer': "Kidney Cancer", 'LC': "Lymphocyte Count",
    'LDL': "Low-Density Lipoprotein", 'LDL Cholesterol': "Low-Density Lipoprotein",
    'LDL-C': "Low-Density Lipoprotein", 'LDL-UKBB': "Low-Density Lipoprotein",
    'LDL-consortium': "Low-Density Lipoprotein",
    'LatinoAny CVD': "Cardio Vascular Disease", 'LatinoDepression': "Depression",
    'LatinoHeart metabolic disease burden': "Heart Metabolic Disease Burden",
    'LatinoHeight': "Height", 'LatinoMigraine Diagnosis': "Migraine Diagnosis",
    'LatinoMorning Person': "Morning Person", 'LatinoSBMN': "Sing Back Musical Note",
    'Lung Cancer': "Lung Cancer", 'MBI-UKBB': "Body Mass Index",
    'MBI-consortium': "Body Mass Index", 'MC': "Monocyte Count",
    'MCH': "Mean Corpuscular Hemoglobin",
    'MCHC': "Mean Corpuscular Hemoglobin Concentration",
    'MCV': "Mean Corpuscular Volume", 'MDD': "Depression", 'MP': "Morning Person",
    'MY': "Myxedema", 'Melanoma': "Melanoma", 'MultiScler': "Multiple Sclerosis",
    'NC': "Neutrophil Count", 'OA': "Osteoarthritis", 'Obesity': "Obesity",
    'Osteoporosis': "Osteoporosis", 'Osteoporosis-BBJ': "Osteoporosis",
    'Osteoporosis-UKB': "Osteoporosis", 'Ovary Cancer': "Ovarian Cancer",
    'PC': "Platelet Count", 'PLT': "Platelet Count", 'PRCA': "Prostate Cancer",
    'Pancreas Cancer': "Pancreatic Cancer", 'Parkinson\u2019s disease': "Parkinson's Disease",
    'PrCa': "Prostate Cancer", 'PrCa-Consortium': "Prostate Cancer",
    'PrCa-FinnGen': "Prostate Cancer", 'PrCa-UKBB': "Prostate Cancer",
    'Primary open angle glaucoma': "Primary Open Angle Glaucoma",
    'Prostate Cancer': "Prostate Cancer", 'Prostate cancer': "Prostate Cancer",
    'Psoriasis': "Psoriasis", 'QU': "Qualification", 'RA': "Rheumatoid Arthirits",
    'RBC': "Red Blood Cell Count", 'RBCC': "Red Blood Cell Count",
    'RDW': "Red Blood Cell Distribution Width", 'RheuArth': "Rheumatoid Arthirits",
    'Rheumatoid arthritis': "Rheumatoid Arthirits", 'SAF': "Salt Added to Food",
    'SASHDL': "High-Density Lipoprotein", 'SASLDL': "Low-Density Lipoprotein",
    'SASTC': "Cholesterol", 'SASlogTG': "Triglyceride", 'SBP': "Systolic Blood Pressure",
    'SCZ': "Schizophrenia", 'SCZ-MGS': "Schizophrenia", 'SCZ-ISC': "Schizophrenia",
    'CHOL': "Cholesterol", 'TRIG': "Triglyceride", 'SH': "Height", 'SN': "Snoring",
    'SS': "Smoking Status", 'SU': "Sodium in Urine", 'Schizophrenia': "Schizophrenia",
    'Sleep Duration': "Sleep Duration",
    'South AsianAny CVD': "Cardio Vascular Disease", 'South AsianDepression': "Depression",
    'South AsianHeart metabolic disease burden': "Heart Metabolic Disease Burden",
    'South AsianHeight': "Height", 'South AsianMigraine Diagnosis': "Migraine Diagnosis",
    'South AsianMorning Person': "Morning Person", 'South AsianSBMN': "Sing Back Musical Note",
    'Stroke': "Stroke", 'Systemic lupus erythematosus': "Systemic Lupus Erythematosus",
    'Systollic Blood Pressure': "Systolic Blood Pressure", 'T1B': "Type One Balding",
    'T1D': "Type One Diabetes", 'T2D': "Type Two Diabetes",
    'T2D-BBJ': "Type Two Diabetes", 'T2D-Consortium': "Type Two Diabetes",
    'T2D-FinnGen': "Type Two Diabetes", 'T2D-UKB': "Type Two Diabetes",
    'T2D-UKBB': "Type Two Diabetes", 'T2d': "Type Two Diabetes",
    'TA': "Tanning Ability", 'TC': "Cholesterol", 'TC-UKBB': "Cholesterol",
    'TC-consortium': "Cholesterol", 'TE': "Tense", 'TFP': "Trunk Fat Percentage",
    'TG': "Triglyceride", 'Total Cholesterol': "Cholesterol",
    'Total cholesterol': "Cholesterol", 'Triglyceride': "Triglyceride",
    'Triglycerides': "Triglyceride", 'Ulcerative colitis': "Ulcerative Colitis",
    'Urate': "Urate", 'VMS': "Supplementary Vitamin and Mineral",
    'Venous thromboembolic disease': "Venous Thromboembolic Disease",
    'WBC': "White Blood Cell Count", 'WBCC': "White Blood Cell Count",
    'WC': "Waist Circumference", 'WHR': "Waist Hip Ratio",
    'eGFR': "Estimated Glomerular Filtration Rate", 'logTG-UKBB': "Triglyceride",
    'logTG-consortium': "Triglyceride", 'Height-H&R Study': "Height",
    'BMI-H&R Study': "Body Mass Index", 'Height-Estonian BioBank': "Height",
    'BMI-H&R Estonian Biobank': "Height", 'T2D AUC': "Type Two Diabetes",
}

# Sheets in the applied workbook that report per-cohort results for a single
# phenotype rather than being named after that phenotype directly. Mirrors
# the special-cased calls in python/appliedRankings.ipynb.
COHORT_SHEET_PHENOTYPES = {
    'AUC-Alzheimer-YuexuanXu': "Alzheimer's Disease",
    'AUC-Alzheimer-GannaLeonenk': "Alzheimer's Disease",
    'AUC-alz-GWAS-EftychiaBellou': "Alzheimer's Disease",
    'AUC-alz-RD-GWAS-EftychiaBellou': "Alzheimer's Disease",
    'AUC-MDD-ByCohort-GuiyanNi': "Depression",
    'AUC-SCZ-ByCohort-GuiyanNi': "Schizophrenia",
}

METADATA_SHEETS = {'models', 'Method papers', 'Application and benchmarking'}


def get_info(df, method_cols, helper, all_aa, all_ww,
             traits=None, trait_hits=None, paper_name=None):
    """
    Reproduces getInfo() from python/methodRanking.ipynb and
    python/appliedRankings.ipynb: for every pair of methods compared in a
    row, add one row to AA0 (comparison indicator) and WW0 (winner
    indicator). If `traits` is given, also route the comparison into the
    matching phenotype-specific matrices.
    """
    n_methods = len(method_cols)
    cols = df.columns
    r, c = df.shape
    aa_rows, ww_rows = [], []
    trait_rows = {}

    for x in range(r):
        for y in range(1, c):
            comparer = cols[y]
            if comparer not in helper:
                continue
            for z in range(y + 1, c):
                comparie = cols[z]
                if comparie not in helper:
                    continue
                v1, v2 = df[comparer].iloc[x], df[comparie].iloc[x]
                if pd.isna(v1) or pd.isna(v2):
                    continue
                diff = v1 - v2
                if diff == 0:
                    continue

                arr = np.zeros(n_methods)
                arr[helper[comparer]] = 1
                arr[helper[comparie]] = 1
                win = np.zeros(n_methods)
                win[helper[comparer] if diff > 0 else helper[comparie]] = 1
                aa_rows.append(arr)
                ww_rows.append(win)

                if traits is not None:
                    phen = df[cols[0]].iloc[x]
                    std_phen = TRAIT_MAP.get(phen)
                    if std_phen is not None and std_phen in traits:
                        if trait_hits is not None and paper_name is not None:
                            trait_hits.setdefault(std_phen, set()).add(paper_name)
                        trait_rows.setdefault(std_phen, ([], []))
                        trait_rows[std_phen][0].append(arr)
                        trait_rows[std_phen][1].append(win)

    if aa_rows:
        all_aa.extend(aa_rows)
        all_ww.extend(ww_rows)
    for phen, (a_list, w_list) in trait_rows.items():
        traits[phen][0].extend(a_list)
        traits[phen][1].extend(w_list)


def load_workbook(path, data_type):
    """Loads the models list and per-paper comparison sheets."""
    models_sheet = 'models' if data_type == 'method' else 'Method papers'
    models = pd.read_excel(path, sheet_name=models_sheet)
    method_cols = list(models['Method'])
    helper = {m: i for i, m in enumerate(method_cols)}

    all_sheets = pd.read_excel(path, sheet_name=None)
    data_sheets = {k: v for k, v in all_sheets.items() if k not in METADATA_SHEETS}
    return method_cols, helper, data_sheets


def build_matrices(method_cols, helper, data_sheets, want_phenotype=None):
    all_aa, all_ww = [], []
    traits = None
    trait_hits = None
    if want_phenotype is not None:
        traits = {want_phenotype: ([], [])}
        trait_hits = {}

    for name, df in data_sheets.items():
        paper_name = name.split('-')[-1]
        get_info(df, method_cols, helper, all_aa, all_ww,
                 traits=traits, trait_hits=trait_hits, paper_name=paper_name)

        if name in COHORT_SHEET_PHENOTYPES and traits is not None:
            phen = COHORT_SHEET_PHENOTYPES[name]
            if phen in traits:
                # cohort-level sheets are keyed by their own phenotype-style
                # first column already, so a second pass under the fixed
                # phenotype label captures them too, matching
                # getInfoCohort()'s behavior in appliedRankings.ipynb
                fake_traits = {phen: traits[phen]}
                get_info(df, method_cols, helper, [], [],
                         traits=fake_traits, trait_hits=trait_hits,
                         paper_name=paper_name)

    aa_df = pd.DataFrame(all_aa if want_phenotype is None else traits[want_phenotype][0],
                          columns=method_cols)
    ww_df = pd.DataFrame(all_ww if want_phenotype is None else traits[want_phenotype][1],
                          columns=method_cols)
    return aa_df, ww_df


def restrict_to_methods(aa_df, ww_df, keep_methods):
    drop_methods = [m for m in aa_df.columns if m not in keep_methods]
    for m in drop_methods:
        rows_with_m = aa_df.index[aa_df[m] == 1]
        aa_df = aa_df.drop(index=rows_with_m)
        ww_df = ww_df.drop(index=rows_with_m)
    aa_df = aa_df.drop(columns=drop_methods).reset_index(drop=True)
    ww_df = ww_df.drop(columns=drop_methods).reset_index(drop=True)
    return aa_df, ww_df


def apply_min_comparisons(aa_df, ww_df, min_comparisons):
    """
    Iteratively drops any method with fewer than `min_comparisons` total
    comparisons, and the rows that reference it -- same loop as
    rankingCalculation() in python/appliedRankings.ipynb.
    """
    while True:
        totals = aa_df.sum(axis=0)
        under = list(totals[totals < min_comparisons].index)
        if not under:
            break
        for m in under:
            rows_with_m = aa_df.index[aa_df[m] == 1]
            aa_df = aa_df.drop(index=rows_with_m)
            ww_df = ww_df.drop(index=rows_with_m)
        aa_df = aa_df.drop(columns=under).reset_index(drop=True)
        ww_df = ww_df.drop(columns=under).reset_index(drop=True)
    return aa_df, ww_df


def combine_aa_ww(aa1, ww1, aa2, ww2, label1, label2):
    """
    Row-stacks two AA0/WW0 pairs into one combined comparison pool -- e.g.
    method-development paper comparisons plus applied/benchmarking paper
    comparisons, as done by hand at the end of appliedRankings.ipynb. Aligns
    columns by method name (union of both, missing entries filled with 0)
    in case the two workbooks ever list methods in different orders.
    """
    all_methods = list(dict.fromkeys(list(aa1.columns) + list(aa2.columns)))
    if list(aa1.columns) != list(aa2.columns):
        print(f"Note: {label1} and {label2} method-column order/contents differ; "
              f"aligning to the union of both ({len(all_methods)} methods).",
              file=sys.stderr)
    aa1 = aa1.reindex(columns=all_methods, fill_value=0)
    ww1 = ww1.reindex(columns=all_methods, fill_value=0)
    aa2 = aa2.reindex(columns=all_methods, fill_value=0)
    ww2 = ww2.reindex(columns=all_methods, fill_value=0)
    aa_df = pd.concat([aa1, aa2], ignore_index=True)
    ww_df = pd.concat([ww1, ww2], ignore_index=True)
    return aa_df, ww_df, all_methods


def main():
    parser = argparse.ArgumentParser(
        description="Build a customized AA0/WW0 pair for R/prsRankingCalculation.R")
    parser.add_argument('--data-type', choices=['method', 'applied', 'combined'], required=True,
                         help="'method' = data/PRSMethodPapersGitHub.xlsx, "
                              "'applied' = data/PRSPaperAppliedGitHub.xlsx, "
                              "'combined' = both, row-stacked into one comparison pool "
                              "(requires --method-workbook and --applied-workbook instead of --workbook)")
    parser.add_argument('--workbook', default=None,
                         help='Path to the .xlsx workbook (for --data-type method or applied)')
    parser.add_argument('--method-workbook', default=None,
                         help='Path to the method-development workbook (for --data-type combined)')
    parser.add_argument('--applied-workbook', default=None,
                         help='Path to the applied/benchmarking workbook (for --data-type combined)')
    parser.add_argument('--methods', default=None,
                         help='Comma-separated list of methods to include (default: all)')
    parser.add_argument('--phenotype', default=None,
                         help='Standardized phenotype name to restrict to '
                              '(applied data only, e.g. "Breast Cancer")')
    parser.add_argument('--min-comparisons', type=int, default=10,
                         help='Drop methods with fewer than this many comparisons '
                              'in the selected data (default: 10)')
    parser.add_argument('--outdir', required=True, help='Directory to write AA0.csv/WW0.csv')
    args = parser.parse_args()

    if args.data_type == 'combined':
        if not args.method_workbook or not args.applied_workbook:
            sys.exit('--data-type combined requires both --method-workbook and --applied-workbook.')
        if args.phenotype is not None:
            sys.exit('--phenotype cannot be combined with --data-type combined -- the '
                      'method-development workbook has no phenotype information to '
                      'restrict to, so a phenotype-specific slice can only be built '
                      'from --data-type applied on its own.')

        m_cols, m_helper, m_sheets = load_workbook(args.method_workbook, 'method')
        a_cols, a_helper, a_sheets = load_workbook(args.applied_workbook, 'applied')
        aa_m, ww_m = build_matrices(m_cols, m_helper, m_sheets)
        aa_a, ww_a = build_matrices(a_cols, a_helper, a_sheets)
        aa_df, ww_df, method_cols = combine_aa_ww(
            aa_m, ww_m, aa_a, ww_a, 'method workbook', 'applied workbook')
    else:
        if not args.workbook:
            sys.exit(f'--workbook is required for --data-type {args.data_type}.')
        if args.phenotype is not None and args.data_type != 'applied':
            sys.exit('--phenotype is only available for --data-type applied '
                      '(the method-development workbook is not split by phenotype).')

        method_cols, helper, data_sheets = load_workbook(args.workbook, args.data_type)
        aa_df, ww_df = build_matrices(method_cols, helper, data_sheets,
                                       want_phenotype=args.phenotype)

    if args.methods:
        keep = [m.strip() for m in args.methods.split(',')]
        unknown = [m for m in keep if m not in method_cols]
        if unknown:
            sys.exit(f"Unknown method name(s): {unknown}. "
                      f"Available methods: {method_cols}")
        aa_df, ww_df = restrict_to_methods(aa_df, ww_df, keep)

    aa_df, ww_df = apply_min_comparisons(aa_df, ww_df, args.min_comparisons)

    if aa_df.shape[1] < 2 or aa_df.shape[0] == 0:
        sys.exit('Not enough data left after filtering to rank anything -- '
                  'try a lower --min-comparisons, a broader --methods list, '
                  'or a different --phenotype.')

    import os
    os.makedirs(args.outdir, exist_ok=True)
    aa_path = os.path.join(args.outdir, 'AA0.csv')
    ww_path = os.path.join(args.outdir, 'WW0.csv')
    aa_df.to_csv(aa_path, index=False)
    ww_df.to_csv(ww_path, index=False)

    print(f"Wrote {aa_path} and {ww_path}")
    print(f"Methods included ({aa_df.shape[1]}): {list(aa_df.columns)}")
    print(f"Comparisons included: {aa_df.shape[0]}")
    print("\nPoint R/prsRankingCalculation.R's aCSV/wCSV reads at these two files "
          "to get ranks and confidence intervals for this selection.")


if __name__ == '__main__':
    main()
