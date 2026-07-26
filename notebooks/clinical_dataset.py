# /// script
# dependencies = ["biopython", "kaleido", "matplotlib", "numpy", "openpyxl", "plotly", "plotly-express", "scikit-posthocs", "scipy", "seaborn", "statannotations", "tabulate", "tqdm"]
# ///

import marimo

__generated_with = "0.23.11"
app = marimo.App()


@app.cell
def _():
    # packages added via marimo's package management: biopython !pip install biopython
    # packages added via marimo's package management: matplotlib !pip install matplotlib
    # packages added via marimo's package management: numpy !pip install numpy
    # packages added via marimo's package management: seaborn !pip install seaborn
    # packages added via marimo's package management: tqdm !pip install tqdm
    # packages added via marimo's package management: plotly !pip install plotly
    # packages added via marimo's package management: scipy !pip install scipy
    # packages added via marimo's package management: plotly.express !pip install plotly.express
    # packages added via marimo's package management: kaleido==0.2.1 !pip install kaleido==0.2.1
    # packages added via marimo's package management: openpyxl !pip install openpyxl
    # packages added via marimo's package management: scikit_posthocs !pip install scikit_posthocs
    return


@app.cell
def _():
    import pandas as pd
    import subprocess
    from Bio.Seq import Seq
    from Bio import SeqIO
    import io

    MAFFT_PATH = "mafft"
    # --- Configuration ---
    # If mafft is not in your system PATH, provide the full path here.
    # e.g., MAFFT_PATH = "/usr/local/bin/mafft" or "C:\\mafft-win\\mafft.bat"

    def get_mafft_consensus(row):
        fwd = (
            str(row["Forward Sequence"]).strip()
            if pd.notna(row["Forward Sequence"])
            else ""
        )
        rev = (
            str(row["Reverse Sequence"]).strip()
            if pd.notna(row["Reverse Sequence"])
            else ""
        )  # 1. Extract and Clean Sequences
        if not fwd and (not rev):
            return ""
        if fwd and (not rev):
            return fwd  # Handle missing data
        if rev and (not fwd):
            try:
                return str(Seq(rev).reverse_complement())
            except:
                return ""
        try:
            rev_rc = str(
                Seq(rev).reverse_complement()
            )  # 2. Prepare Fasta Content for MAFFT
            fasta_content = f">Forward\n{fwd}\n>Reverse_RC\n{rev_rc}\n"
        except Exception:  # Reverse complement the reverse read
            return fwd
        try:
            process = subprocess.Popen(
                [MAFFT_PATH, "--quiet", "--auto", "-"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )  # Create a string that looks like a FASTA file
            (stdout, stderr) = process.communicate(input=fasta_content)
            if process.returncode != 0:
                print(
                    f"MAFFT Error for {_row.get('Patient Code', 'Unknown')}: {stderr}"
                )  # Fallback if RC fails
                return fwd
            aligned_seqs = list(
                SeqIO.parse(io.StringIO(stdout), "fasta")
            )  # 3. Run MAFFT using subprocess
            if len(aligned_seqs) < 2:
                return fwd  # We pass the fasta content directly to stdin using Popen
            seq1 = str(aligned_seqs[0].seq).upper()  # Command: mafft --quiet --auto -
            seq2 = str(
                aligned_seqs[1].seq
            ).upper()  # The '-' tells mafft to read from stdin
            consensus = []
            for b1, b2 in zip(seq1, seq2):
                if b1 == b2:
                    consensus.append(b1)
                elif b1 == "-":
                    consensus.append(b2)
                elif b2 == "-":
                    consensus.append(b1)
                else:  # Send data and get output
                    consensus.append("N")
            return "".join(consensus).replace("-", "")
        except FileNotFoundError:
            print("Error: MAFFT not found. Please install MAFFT or check the path.")
            return fwd  # Fallback
        except Exception as e:
            print(
                f"Error processing {_row.get('Patient Code', 'Unknown')}: {e}"
            )  # 4. Parse Aligned Sequences
            return fwd  # stdout contains the aligned fasta

    df = pd.read_excel("../data/raw/Patient_Metadata.xlsx")
    print("Running MAFFT alignment for each patient...")
    try:
        subprocess.run(
            [MAFFT_PATH, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except FileNotFoundError:
        print(
            "CRITICAL WARNING: 'mafft' command not found. Install it with '!apt-get install mafft' in Colab."
        )
    df["Consensus_Sequence"] = df.apply(get_mafft_consensus, axis=1)
    output_file = "../data/processed/Patient_Consensus_MAFFT.xlsx"
    df.to_excel(output_file, index=False)  # 5. Build Consensus
    # --- Main Script ---
    # 1. Load Data
    # 2. Check if MAFFT is callable
    # 3. Apply Consensus
    # 4. Save
    print(
        f"Done! Saved to {output_file}"
    )  # Gap in seq1 -> Insertion in seq2 (or seq1 missing coverage) -> Keep seq2  # Gap in seq2 -> Keep seq1  # Mismatch -> N  # Join and remove any remaining gap characters (if any artifacts exist)
    return SeqIO, df, pd, subprocess


@app.cell
def _(df):
    fasta_output = "../data/processed/Final_Consensus_Sequences.fasta"
    print(f"Writing FASTA file to {fasta_output}...")
    with open(fasta_output, "w") as f:
        count = 0
        for _i, _row in df.iterrows():
            seq = str(_row["Consensus_Sequence"]).strip()
            if seq and seq.upper() != "NAN":
                p_code = str(_row["Patient Code"]).strip().replace(" ", "_")
                age = str(_row["Age Group"]).strip().replace(" ", "_")
                category = str(_row["Category"]).strip().replace(" ", "_")
                header = f">{p_code}_{age}_{category}"
                f.write(f"{header}\n{seq}\n")
                count = count + 1
    print(f"Done! {count} sequences written to '{fasta_output}'.")
    return


@app.cell
def _(SeqIO, pd):
    import re
    from Bio import SeqUtils

    def calculate_metrics(fasta_file, output_file):
        data = []
        print(f"Reading {fasta_file}...")
        for record in SeqIO.parse(fasta_file, "fasta"):
            header_parts = record.description.split("_")
            if len(header_parts) >= 3:
                patient_code = header_parts[0]
                age_group = header_parts[1]
                category = "_".join(header_parts[2:])
            else:
                patient_code = record.id
                age_group = "Unknown"
                category = "Unknown"
            seq = str(record.seq).upper()
            seq_len = len(seq)
            if seq_len == 0:
                continue
            _metrics = {
                "Patient Code": patient_code,
                "Age Group": age_group,
                "Category": category,
                "seq_len": seq_len,
                "seq": seq,
            }
            _metrics["A"] = seq.count("A")
            _metrics["T"] = seq.count("T")
            _metrics["C"] = seq.count("C")
            _metrics["G"] = seq.count("G")
            dinucl = [
                "AA",
                "AT",
                "AC",
                "AG",
                "TA",
                "TT",
                "TC",
                "TG",
                "CA",
                "CT",
                "CC",
                "CG",
                "GA",
                "GT",
                "GC",
                "GG",
            ]
            for din in dinucl:
                _metrics[f"{din.lower()}_count"] = seq.count(din)
            zap_patterns = {
                "zap4": "CNNNNGNCG",
                "zap5": "CNNNNNGNCG",
                "zap6": "CNNNNNNGNCG",
                "zap7": "CNNNNNNNGNCG",
                "zap8": "CNNNNNNNNGNCG",
            }
            total_zap = 0
            for _name, pattern in zap_patterns.items():
                search_res = SeqUtils.nt_search(seq, pattern)
                count = len(search_res) - 1
                total_zap = total_zap + count
            _metrics["zap_motif"] = total_zap
            drach_count = len(re.findall("[AGT][AG]AC[ACT]", seq))
            _metrics["drach_motif"] = drach_count
            _metrics["cpg_count"] = _metrics["cg_count"]
            _metrics["cpg_per"] = 100 * (_metrics["cpg_count"] / seq_len)
            try:
                _metrics["gc_content"] = SeqUtils.gc_fraction(seq) * 100
            except AttributeError:
                _metrics["gc_content"] = SeqUtils.GC(seq)
            for din in dinucl:
                n1 = din[0]
                n2 = din[1]
                count_n1 = _metrics[n1]
                count_n2 = _metrics[n2]
                din_name = din.lower()
                if count_n1 > 0 and count_n2 > 0:
                    obe = (
                        _metrics[f"{din_name}_count"] * seq_len / (count_n1 * count_n2)
                    )
                else:
                    obe = 0.0
                if din == "CG":
                    _metrics["cpg_obye"] = obe
                else:
                    _metrics[f"{din_name}_obye"] = obe
            data.append(_metrics)
        df_results = pd.DataFrame(data)
        df_results.to_excel(output_file, index=False)
        print(f"Metrics generated for {len(df_results)} sequences.")
        print(f"Saved to {output_file}")
        return df_results

    file_input = "../data/raw/ILBS_ALL.fasta"
    file_output = "../results/tables/ILBS_Sequence_Metrics.xlsx"
    try:
        df_metrics = calculate_metrics(file_input, file_output)
        print(df_metrics.head())
    except FileNotFoundError:
        print(
            f"Error: Could not find file '{file_input}'. Please make sure it is uploaded."
        )
    return


@app.cell
def _():
    # packages added via marimo's package management: statannotations !pip install statannotations
    return


@app.cell
def _(subprocess):
    with open("../data/processed/ILBS_ALL_aligned_linsi.fasta", "w") as out_file:
        subprocess.call(
            [
                "mafft",
                "--localpair",
                "--maxiterate",
                "1000",
                "--thread",
                "-1",
                "../data/raw/ILBS_ALL.fasta",
            ],
            stdout=out_file,
        )
    return


@app.cell
def _(pd):
    from Bio import AlignIO

    def find_mutations_from_alignment(alignment_file, output_file):
        try:
            alignment = AlignIO.read(alignment_file, "fasta")
        except FileNotFoundError:
            print(f"Error: File {alignment_file} not found.")
            return
        reference_record = alignment[0]
        ref_seq = str(reference_record.seq).upper()
        ref_name = reference_record.id
        print(f"Reference Sequence: {ref_name} (Length: {len(ref_seq)})")
        mutation_data = []
        for _i in range(1, len(alignment)):
            query_record = alignment[_i]
            query_seq = str(query_record.seq).upper()
            query_name = query_record.id
            parts = query_name.split("_")
            patient_code = parts[0] if len(parts) > 0 else query_name
            age_group = parts[1] if len(parts) > 1 else "Unknown"
            category = "_".join(parts[2:]) if len(parts) > 2 else "Unknown"
            for pos, (ref_base, query_base) in enumerate(zip(ref_seq, query_seq)):
                display_pos = pos + 1
                if ref_base == query_base:
                    continue
                if ref_base == "-" and query_base == "-":
                    continue
                mut_type = ""
                notation = ""
                if ref_base == "-":
                    mut_type = "Insertion"
                    notation = f"ins_{display_pos}_{query_base}"
                elif query_base == "-":
                    mut_type = "Deletion"
                    notation = f"del_{display_pos}_{ref_base}"
                elif query_base == "N":
                    mut_type = "Ambiguous"
                    notation = f"{ref_base}{display_pos}N"
                else:
                    mut_type = "Substitution"
                    notation = f"{ref_base}{display_pos}{query_base}"
                if mut_type != "":
                    mutation_data.append(
                        {
                            "Patient Code": patient_code,
                            "Age Group": age_group,
                            "Category": category,
                            "Position": display_pos,
                            "Ref_Base": ref_base,
                            "Query_Base": query_base,
                            "Mutation_Type": mut_type,
                            "Notation": notation,
                            "Ref_Context": ref_seq[
                                max(0, pos - 2) : min(len(ref_seq), pos + 3)
                            ],
                        }
                    )
        if mutation_data:
            df_mutations = pd.DataFrame(mutation_data)
            df_mutations.to_excel(output_file, index=False)
            print(f"Successfully found {len(df_mutations)} mutations.")
            print(f"Report saved to: {output_file}")
            print("\nMutation Summary by Type:")
            print(df_mutations["Mutation_Type"].value_counts())
        else:
            print("No mutations found (sequences are identical).")

    alignment_filename = "../data/processed/ILBS_ALL_aligned_linsi.fasta"
    output_filename = "../results/tables/Mutation_Report.xlsx"
    find_mutations_from_alignment(alignment_filename, output_filename)
    return


@app.cell
def _(pd):
    import matplotlib.pyplot as plt
    import seaborn as sns

    df_1 = pd.read_excel("../results/tables/Mutation_Report.xlsx")

    def clean_data(df):
        substitutions = df[df["Mutation_Type"] == "Substitution"].copy()
        indels = df[df["Mutation_Type"].isin(["Insertion", "Deletion"])].copy()
        max_pos = df["Position"].max()
        indels = indels[
            ~((indels["Position"] < 50) | (indels["Position"] > max_pos - 50))
        ]

        def is_homopolymer(row):
            context = str(row["Ref_Context"])
            ref_base = str(row["Ref_Base"])
            if ref_base * 3 in context or any(
                (b * 4 in context for b in ["A", "C", "G", "T"])
            ):
                return True
            return False

        indels = indels[~indels.apply(is_homopolymer, axis=1)]
        indels = indels.sort_values(by=["Patient Code", "Position"])
        indels["prev_pos"] = indels.groupby("Patient Code")["Position"].shift(1)
        indels["gap_group"] = (indels["Position"] - indels["prev_pos"] != 1).cumsum()
        event_lengths = (
            indels.groupby(["Patient Code", "gap_group"])["Position"]
            .count()
            .reset_index(name="Event_Length")
        )
        indels = pd.merge(indels, event_lengths, on=["Patient Code", "gap_group"])
        valid_indels = indels[indels["Event_Length"] % 3 == 0].drop(
            columns=["prev_pos", "gap_group", "Event_Length"]
        )
        return pd.concat([substitutions, valid_indels])

    _df_clean = clean_data(df_1)
    counts_raw = df_1["Patient Code"].value_counts().reset_index()
    counts_raw.columns = ["Patient Code", "Count"]
    counts_raw["Dataset"] = "Raw (With Artifacts)"
    counts_clean = _df_clean["Patient Code"].value_counts().reset_index()
    counts_clean.columns = ["Patient Code", "Count"]
    counts_clean["Dataset"] = "Cleaned (High Confidence)"
    comparison = pd.concat([counts_raw, counts_clean])
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=comparison, x="Dataset", y="Count", palette=["#e74c3c", "#2ecc71"])
    plt.title("Impact of Artifact Removal on Mutation Load")
    plt.ylabel("Mutations per Patient")
    plt.grid(True, alpha=0.3)
    plt.show()
    print(f"Original Rows: {len(df_1)}")
    print(f"Cleaned Rows:  {len(_df_clean)}")
    print(f"Artifacts Removed: {len(df_1) - len(_df_clean)}")
    _df_clean.to_excel("../results/tables/Mutation_Report_Cleaned.xlsx", index=False)
    return plt, sns


@app.cell
def _(pd, plt, sns):
    import sys
    import warnings

    from matplotlib.gridspec import GridSpec
    import numpy as np

    import scikit_posthocs as sp
    from scipy.stats import chi2_contingency, fisher_exact, kruskal

    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.utils import resample

    warnings.filterwarnings("ignore")

    # --- PLOTTING & STYLING PARAMETERS ---
    _LABEL_FS = 18
    _TICK_FS = 15
    _TITLE_FS = 20
    _ANNOT_FS = 14
    _LETTER_FS = 28
    _LEGEND_FS = 14

    plt.rcParams["svg.fonttype"] = "none"
    sns.set_context("paper", font_scale=2.5)
    sns.set_style(
        "ticks",
        {
            "axes.linewidth": 2.0,
            "xtick.major.width": 2.0,
            "ytick.major.width": 2.0,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
        },
    )

    # --- FILE PATHS ---
    _MUTATION_FILE = "../results/tables/Mutation_Report_Cleaned.xlsx"
    _METRICS_FILE = "../results/tables/ILBS_Sequence_Metrics.xlsx"
    _GLOBAL_FILE = "../data/raw/global_dataset.xlsx"
    _OUTPUT_FIG = "../results/figures/Figure6.png"
    _OUTPUT_SVG = "../results/figures/Figure6.svg"
    _OUTPUT_REPORT = "../results/tables/Figure6_stats.txt"

    # --- PALETTES & MAPPINGS ---
    _PALETTE_GROUPS = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD"]
    _RISK_COLORS = ["#4C72B0", "#E8A317", "#D73027"]
    _COLOR_MAP = {
        "ga_obye": "#4C72B0",
        "APOBEC_Load": "#D73027",
        "CpG_Deam_Load": "#E8A317",
    }
    _LABEL_MAP = {
        "ga_obye": "GpA O/E",
        "APOBEC_Load": "G[G>A]T",
        "CpG_Deam_Load": "C[C>T]G",
    }
    _AGE_ABBR = {"Pediatric": "Ped", "Adolescent": "Adol", "Adult": "Adult"}


    # --- HELPER FUNCTIONS ---
    def report_kruskal_dunn(df, group_col, value_col, f):
        clean_df = df.dropna(subset=[group_col, value_col])
        _groups = [
            clean_df[clean_df[group_col] == g][value_col]
            for g in sorted(clean_df[group_col].unique())
        ]
        if len(_groups) > 1:
            (_stat, _p) = kruskal(*_groups)
            f.write("\n==================================================\n")
            f.write(f"--- OMNIBUS KRUSKAL-WALLIS ({value_col}) ---\n")
            f.write(f"Statistic: {_stat:.4f}, p-value: {_p:.4e}\n")
            if _p < 0.05:
                dunn_matrix = sp.posthoc_dunn(
                    clean_df, val_col=value_col, group_col=group_col, p_adjust="fdr_bh"
                )
                f.write("\n--- DUNN'S POST-HOC (Benjamini-Hochberg FDR Adjusted) ---\n")
                f.write(dunn_matrix.to_string())
                f.write("\n")
            return _p
        return 1.0


    def _report_chi2(contingency_table, f):
        (_chi2, _p, _dof, _ex) = chi2_contingency(contingency_table)
        f.write(f"\n--- Chi-Square ---\nStatistic: {_chi2:.4f}, p-value: {_p:.4e}\n")
        return _p


    def get_risk_score(row, cut_ga, cut_apobec, cut_cpg):
        s = 0
        if row["ga_obye"] > cut_ga:
            s += 1
        if row["APOBEC_Load"] > cut_apobec:
            s += 1
        if row["CpG_Deam_Load"] > cut_cpg:
            s += 1
        return "Low" if s == 0 else "Intermediate" if s == 1 else "High"


    def _calculate_bootstrap_stats(data, group_name, n_boot=1000):
        predictors = ["ga_obye", "APOBEC_Load", "CpG_Deam_Load"]
        if len(data) < 8 or len(data["Category"].unique()) < 2:
            return None
        _X = data[predictors].fillna(0)
        _y = data["Category"].apply(lambda x: 1 if x == "ALF" else 0).values
        _scaler = StandardScaler()
        _X_scaled = _scaler.fit_transform(_X)
        coefs = []
        _clf = LogisticRegression(
            penalty="l2", C=1.0, solver="liblinear", class_weight="balanced"
        )
        for _i in range(n_boot):
            (_X_res, _y_res) = resample(_X_scaled, _y, stratify=_y, random_state=_i)
            if len(np.unique(_y_res)) < 2:
                continue
            try:
                _clf.fit(_X_res, _y_res)
                coefs.append(_clf.coef_[0])
            except Exception:
                continue
        _boot_df = pd.DataFrame(coefs, columns=predictors)
        _results = []
        for _col in predictors:
            mean_log = _boot_df[_col].mean()
            lower_log = _boot_df[_col].quantile(0.025)
            upper_log = _boot_df[_col].quantile(0.975)
            _results.append(
                {
                    "Group": group_name,
                    "Marker": _col,
                    "OR": np.exp(mean_log),
                    "Lower": np.exp(lower_log),
                    "Upper": np.exp(upper_log),
                    "Significant": lower_log > 0 or upper_log < 0,
                    "N": len(data),
                }
            )
        return pd.DataFrame(_results)


    # --- DATA LOADING & PROCESSING ---
    print("Loading data...")
    try:
        _df_mut_raw = pd.read_excel(_MUTATION_FILE)
        _df_metrics_raw = pd.read_excel(_METRICS_FILE)
        for df_6 in [_df_mut_raw, _df_metrics_raw]:
            if "Patient Code" in df_6.columns:
                df_6["Patient Code"] = df_6["Patient Code"].astype(str).str.strip()
            if "Age Group" in df_6.columns:
                df_6["Age Group"] = df_6["Age Group"].astype(str).str.strip()

        def _create_short_label(df_in):
            return (
                df_in["Age Group"].map(_AGE_ABBR).fillna(df_in["Age Group"])
                + "-"
                + df_in["Category"]
            )

        _df_mut_p1 = _df_mut_raw.copy()
        _df_mut_p1["ShortGroup"] = _create_short_label(_df_mut_p1)
        _df_metrics_p1 = _df_metrics_raw.copy()
        _df_metrics_p1 = _df_metrics_p1.dropna(subset=["Category", "Age Group"])
        _df_metrics_p1 = _df_metrics_p1[
            _df_metrics_p1["Category"].astype(str).str.lower() != "unknown"
        ]
        _df_metrics_p1["ShortGroup"] = _create_short_label(_df_metrics_p1)
        _df_metrics_p2 = _df_metrics_raw.copy()
        _ctx_col = next(
            (
                _c
                for _c in ["Ref_Context", "Context", "Trinucleotide"]
                if _c in _df_mut_raw.columns
            ),
            None,
        )
        if _ctx_col:

            def _get_sig(row):
                try:
                    ctx = str(row[_ctx_col])
                    if len(ctx) < 3:
                        return "Unknown"
                    return f"{ctx[0]}[{row['Ref_Base']}>{row['Query_Base']}]{ctx[-1]}"
                except Exception:
                    return "Unknown"

            _df_mut_raw["Sig"] = _df_mut_raw.apply(_get_sig, axis=1)
            _apobec_counts = (
                _df_mut_raw[_df_mut_raw["Sig"] == "G[G>A]T"]
                .groupby("Patient Code")
                .size()
            )
            _cpg_counts = (
                _df_mut_raw[_df_mut_raw["Sig"] == "C[C>T]G"]
                .groupby("Patient Code")
                .size()
            )
            _df_metrics_p2["APOBEC_Load"] = (
                _df_metrics_p2["Patient Code"].map(_apobec_counts).fillna(0)
                / _df_metrics_p2["seq_len"]
                * 1000
            )
            _df_metrics_p2["CpG_Deam_Load"] = (
                _df_metrics_p2["Patient Code"].map(_cpg_counts).fillna(0)
                / _df_metrics_p2["seq_len"]
                * 1000
            )
        else:
            for _col in ["APOBEC_Load", "CpG_Deam_Load"]:
                if _col not in _df_metrics_p2.columns:
                    _df_metrics_p2[_col] = 0
        if "ga_obye" not in _df_metrics_p2.columns:
            _df_metrics_p2["ga_obye"] = 0
        _df_model = _df_metrics_p2[
            _df_metrics_p2["Category"].isin(["ALF", "AVH"])
        ].copy()

        print("Running Clinical Models (n=1000)...")
        _stats_global = _calculate_bootstrap_stats(_df_model, "Global", n_boot=1000)
        _stats_ped = _calculate_bootstrap_stats(
            _df_model[
                _df_model["Age Group"].str.contains("Pediatric", case=False, na=False)
            ],
            "Paediatrics",
            n_boot=1000,
        )
        _stats_adol = _calculate_bootstrap_stats(
            _df_model[
                _df_model["Age Group"].str.contains("Adolescent", case=False, na=False)
            ],
            "Adolescents",
            n_boot=1000,
        )
        _stats_adult = _calculate_bootstrap_stats(
            _df_model[
                _df_model["Age Group"].str.contains("Adult", case=False, na=False)
            ],
            "Adults",
            n_boot=1000,
        )
        _valid_stats = [
            s
            for s in [_stats_adult, _stats_adol, _stats_ped, _stats_global]
            if s is not None
        ]
        _all_stats = (
            pd.concat(_valid_stats, ignore_index=True)
            if _valid_stats
            else pd.DataFrame()
        )

        print("Loading Global Dataset...")
        _df_global = pd.read_excel(_GLOBAL_FILE)
        for _col in ["ga_obye", "APOBEC_Load", "CpG_Deam_Load"]:
            if _col not in _df_global.columns:
                if _col == "APOBEC_Load" and "G[G>A]T" in _df_global.columns:
                    _df_global["APOBEC_Load"] = _df_global["G[G>A]T"]
                elif _col == "CpG_Deam_Load" and "C[C>T]G" in _df_global.columns:
                    _df_global["CpG_Deam_Load"] = _df_global["C[C>T]G"]
                else:
                    _df_global[_col] = 0
        if "year" in _df_global.columns:
            _df_global["Era"] = _df_global["year"].apply(
                lambda x: "Contemporary" if x >= 2018 else "Historical"
            )
        else:
            _df_global["Era"] = "Unknown"
    except Exception as e:
        print(f"Data Error: {e}")
        sys.exit(1)


    # --- PLOT GENERATION ---
    print("Generating Plot...")
    _f_rep = open(_OUTPUT_REPORT, "w")
    _fig = plt.figure(figsize=(24, 22))
    _gs = GridSpec(3, 6, figure=_fig, hspace=0.8, wspace=1.2, height_ratios=[1, 1.8, 1])


    def _format_ticks(ax):
        ax.tick_params(axis="both", labelsize=_TICK_FS)


    # Panel A: C>T Enrichment
    _subs_only = _df_mut_p1[_df_mut_p1["Mutation_Type"] == "Substitution"]
    _base_counts = (
        _subs_only.groupby("Patient Code").agg({"ShortGroup": "first"}).reset_index()
    )
    _n_counts_A = _base_counts["ShortGroup"].value_counts()
    _axA = _fig.add_subplot(_gs[0, 0:2])
    _subs_df = _df_mut_p1[_df_mut_p1["Mutation_Type"] == "Substitution"].copy()
    _subs_df["Sig"] = _subs_df["Ref_Base"] + ">" + _subs_df["Query_Base"]
    _data_A = _subs_df.copy()
    _data_A["is_CT"] = (_data_A["Sig"] == "C>T").astype(int)
    _data_A["Label"] = _data_A["ShortGroup"].map(
        {_idx: f"{_idx}\n(n={_c})" for (_idx, _c) in _n_counts_A.items()}
    )
    _means_A = _data_A.groupby("Label")["is_CT"].mean().sort_values(ascending=False)
    _contingency_A = pd.crosstab(_data_A["ShortGroup"], _data_A["is_CT"])
    _p_val_A = _report_chi2(_contingency_A, _f_rep)
    sns.barplot(
        data=_data_A,
        x="Label",
        y="is_CT",
        order=_means_A.index,
        palette="flare_r",
        ax=_axA,
        edgecolor="black",
        linewidth=2.0,
        errorbar="se",
        capsize=0.1,
    )
    _axA.set_title(
        f"C>T Enrichment\n($p$ = {_p_val_A:.2e})", fontweight="bold", fontsize=_TITLE_FS
    )
    _axA.set_ylabel("(C>T / Total Subs) %", fontweight="bold", fontsize=_LABEL_FS)
    _axA.set_xlabel("")
    _axA.tick_params(axis="x", rotation=45)
    _format_ticks(_axA)

    # Panel B: CpG O/E
    _axB = _fig.add_subplot(_gs[0, 2:4])
    if "cpg_obye" in _df_metrics_p1.columns:
        p_b = report_kruskal_dunn(_df_metrics_p1, "ShortGroup", "cpg_obye", _f_rep)
        sns.violinplot(
            data=_df_metrics_p1,
            x="ShortGroup",
            y="cpg_obye",
            order=sorted(_df_metrics_p1["ShortGroup"].unique()),
            palette=_PALETTE_GROUPS,
            ax=_axB,
            inner=None,
            linewidth=2.0,
        )
        sns.swarmplot(
            data=_df_metrics_p1,
            x="ShortGroup",
            y="cpg_obye",
            order=sorted(_df_metrics_p1["ShortGroup"].unique()),
            color=".2",
            alpha=0.6,
            size=5,
            ax=_axB,
        )
        _axB.set_title(
            f"CpG O/E\n($p$ = {p_b:.1e})", fontweight="bold", fontsize=_TITLE_FS
        )
        _axB.set_ylabel("CpG O/E Ratio", fontweight="bold", fontsize=_LABEL_FS)
        _axB.set_xlabel("")
        _axB.tick_params(axis="x", rotation=45)
        _format_ticks(_axB)

    # Panel C: GpA O/E
    _axC = _fig.add_subplot(_gs[0, 4:6])
    if "ga_obye" in _df_metrics_p1.columns:
        p_c = report_kruskal_dunn(_df_metrics_p1, "ShortGroup", "ga_obye", _f_rep)
        sns.violinplot(
            data=_df_metrics_p1,
            x="ShortGroup",
            y="ga_obye",
            order=sorted(_df_metrics_p1["ShortGroup"].unique()),
            palette=_PALETTE_GROUPS,
            ax=_axC,
            inner=None,
            linewidth=2.0,
        )
        sns.swarmplot(
            data=_df_metrics_p1,
            x="ShortGroup",
            y="ga_obye",
            order=sorted(_df_metrics_p1["ShortGroup"].unique()),
            color=".2",
            alpha=0.6,
            size=5,
            ax=_axC,
        )
        _axC.set_title(
            f"GpA O/E\n($p$ = {p_c:.1e})", fontweight="bold", fontsize=_TITLE_FS
        )
        _axC.set_ylabel("GpA O/E Ratio", fontweight="bold", fontsize=_LABEL_FS)
        _axC.set_xlabel("")
        _axC.tick_params(axis="x", rotation=45)
        _format_ticks(_axC)

    # Panel D: Clinical Odds Ratios
    _axD = _fig.add_subplot(_gs[1, 0:3])
    if not _all_stats.empty:
        _plot_df = _all_stats.iloc[::-1].reset_index(drop=True)
        for _i, _row in _plot_df.iterrows():
            _c = _COLOR_MAP.get(_row["Marker"], "black")
            _axD.plot(_row["OR"], _i, "D", color=_c, markersize=10)
            _axD.plot(
                [_row["Lower"], _row["Upper"]], [_i, _i], lw=3.0, color=_c, alpha=0.7
            )
            _txt = f"{_row['OR']:.2f}" + ("*" if _row["Significant"] else "")
            _axD.text(
                _row["Upper"] + 0.2,
                _i,
                _txt,
                va="center",
                fontweight="bold",
                fontsize=_ANNOT_FS,
            )
        _axD.set_yticks(range(len(_plot_df)))
        _axD.set_yticklabels(
            [
                f"{_row['Group']}\n{_LABEL_MAP.get(_row['Marker'], _row['Marker'])}"
                for (_, _row) in _plot_df.iterrows()
            ],
            fontsize=_TICK_FS,
        )
        _axD.axvline(1.0, color="gray", linestyle="--", linewidth=2.0)
        _axD.set_title("Clinical Odds Ratios", fontweight="bold", fontsize=_TITLE_FS)
        _axD.tick_params(axis="x", labelsize=_TICK_FS)

    # Panel E: Feature Importance Heatmap
    _axE = _fig.add_subplot(_gs[1, 3:6])
    if not _all_stats.empty:
        _piv_or = _all_stats[_all_stats["Group"] != "Global"].pivot(
            index="Marker", columns="Group", values="OR"
        )
        _piv_sig = _all_stats[_all_stats["Group"] != "Global"].pivot(
            index="Marker", columns="Group", values="Significant"
        )
        sns.heatmap(
            _piv_or,
            cmap="YlOrRd",
            annot=False,
            ax=_axE,
            vmin=0.5,
            vmax=3.5,
            linewidths=2.0,
            linecolor="black",
        )
        for _y in range(_piv_or.shape[0]):
            for _x in range(_piv_or.shape[1]):
                _val = _piv_or.iloc[_y, _x]
                if pd.isna(_val):
                    continue
                _txt = f"{_val:.2f}" + ("*" if _piv_sig.iloc[_y, _x] else "")
                _axE.text(
                    _x + 0.5,
                    _y + 0.5,
                    _txt,
                    ha="center",
                    va="center",
                    color="black",
                    fontweight="bold",
                    fontsize=_ANNOT_FS,
                )
        _current_labels = [l.get_text() for l in _axE.get_yticklabels()]
        _new_labels = [_LABEL_MAP.get(l, l) for l in _current_labels]
        _axE.set_yticklabels(
            _new_labels, rotation=0, fontsize=_TICK_FS, fontweight="bold"
        )
        _axE.set_xticklabels(
            _axE.get_xticklabels(), fontsize=_TICK_FS, fontweight="bold"
        )
        _axE.set_ylabel("Genome Features", fontweight="bold", fontsize=_LABEL_FS)
        _axE.set_title(
            "Feature Importance Heatmap", fontweight="bold", fontsize=_TITLE_FS
        )

    # Panel F: Clinical Outcome
    _axF = _fig.add_subplot(_gs[2, 0:3])
    _cut_ga = _df_model["ga_obye"].median()
    _cut_apobec = _df_model["APOBEC_Load"].median()
    _cut_cpg = _df_model["CpG_Deam_Load"].median()
    _df_model["Mutational_Signature"] = _df_model.apply(
        lambda r: get_risk_score(r, _cut_ga, _cut_apobec, _cut_cpg), axis=1
    )
    _risk_tab = pd.crosstab(_df_model["Mutational_Signature"], _df_model["Category"])
    _risk_tab = _risk_tab.reindex(["Low", "Intermediate", "High"]).fillna(0)
    try:
        (_odds_risk, _p_risk) = fisher_exact(
            [
                [_risk_tab.loc["High", "ALF"], _risk_tab.loc["High", "AVH"]],
                [_risk_tab.loc["Low", "ALF"], _risk_tab.loc["Low", "AVH"]],
            ]
        )
    except Exception:
        (_odds_risk, _p_risk) = (0, 1.0)
    _risk_tab["ALF_Rate"] = (
        _risk_tab["ALF"] / (_risk_tab["ALF"] + _risk_tab["AVH"]).replace(0, 1) * 100
    )
    _bars = _axF.bar(
        _risk_tab.index,
        _risk_tab["ALF_Rate"],
        color=_RISK_COLORS,
        edgecolor="black",
        linewidth=2.0,
    )
    for _bar in _bars:
        if _bar.get_height() > 0:
            _axF.text(
                _bar.get_x() + _bar.get_width() / 2,
                _bar.get_height() + 2,
                f"{_bar.get_height():.1f}%",
                ha="center",
                fontweight="bold",
                fontsize=_ANNOT_FS,
            )
    _axF.text(
        0.05,
        0.85,
        f"High vs Low:\nOR={_odds_risk:.2f}, $p$={_p_risk:.3f}",
        transform=_axF.transAxes,
        fontsize=_ANNOT_FS,
        bbox=dict(facecolor="white", boxstyle="round", alpha=0.9),
    )
    _axF.set_ylabel("ALF Progression (%)", fontweight="bold", fontsize=_LABEL_FS)
    _axF.set_title("Clinical Outcome", fontweight="bold", fontsize=_TITLE_FS)
    _axF.set_ylim(0, 100)
    _format_ticks(_axF)

    # Panel G: Global Expansion Bootstraps
    _axG = _fig.add_subplot(_gs[2, 3:6])
    print("Running Global Bootstraps (n=10000)...")
    _g_cut_ga = _df_global["ga_obye"].median()
    _g_cut_apobec = _df_global["APOBEC_Load"].median()
    _g_cut_cpg = _df_global["CpG_Deam_Load"].median()
    _df_global["Mutational_Signature"] = _df_global.apply(
        lambda r: get_risk_score(r, _g_cut_ga, _g_cut_apobec, _g_cut_cpg), axis=1
    )
    _plot_global = _df_global.dropna(subset=["year", "Era"]).copy()
    _plot_global = _plot_global[_plot_global["Era"] != "Unknown"]
    _n_boot_g = 10000
    _boot_ors = []
    _data_for_boot = _plot_global[["Era", "Mutational_Signature"]].copy()
    _data_for_boot["is_HighRisk"] = (
        _data_for_boot["Mutational_Signature"] == "High"
    ).astype(int)
    (_chi2, _p_chi, _dof, _ex) = chi2_contingency(
        pd.crosstab(_plot_global["Mutational_Signature"], _plot_global["Era"])
    )
    for _i in range(_n_boot_g):
        _sample = _data_for_boot.sample(frac=1.0, replace=True)
        _ct = pd.crosstab(_sample["is_HighRisk"], _sample["Era"])
        if _ct.shape == (2, 2):
            try:
                _h_post = _ct.loc[1, "Contemporary"]
                _o_post = _ct.loc[0, "Contemporary"]
                _h_pre = _ct.loc[1, "Historical"]
                _o_pre = _ct.loc[0, "Historical"]
                if _o_post * _h_pre > 0:
                    _boot_ors.append(_h_post * _o_pre / (_o_post * _h_pre))
            except Exception:
                continue
    _boot_ors = np.array(_boot_ors)
    if len(_boot_ors) > 0:
        _mean_or = np.mean(_boot_ors)
        _ci_lower = np.percentile(_boot_ors, 2.5)
        _ci_upper = np.percentile(_boot_ors, 97.5)
    else:
        (_mean_or, _ci_lower, _ci_upper) = (0, 0, 0)
    _props = (
        pd.crosstab(
            _plot_global["Era"], _plot_global["Mutational_Signature"], normalize="index"
        )
        * 100
    )
    _era_order = ["Historical", "Contemporary"]
    _era_order = [e for e in _era_order if e in _props.index]
    _props = _props.reindex(_era_order)
    _props = _props.reindex(columns=["Low", "Intermediate", "High"])
    _props.plot(
        kind="bar",
        stacked=True,
        color=_RISK_COLORS,
        edgecolor="black",
        width=0.5,
        ax=_axG,
        linewidth=2.0,
    )
    for _c in _axG.containers:
        _labels = [f"{v.get_height():.1f}%" if v.get_height() > 5 else "" for v in _c]
        _axG.bar_label(
            _c,
            labels=_labels,
            label_type="center",
            fontweight="bold",
            color="white",
            fontsize=_ANNOT_FS,
        )
    _title_text = f"Global Expansion\nOR = {_mean_or:.2f} [{_ci_lower:.2f}-{_ci_upper:.2f}], $p$ = {_p_chi:.1e}"
    _axG.set_title(_title_text, fontweight="bold", fontsize=_TITLE_FS - 2)
    _axG.set_ylabel("Proportion (%)", fontweight="bold", fontsize=_LABEL_FS)
    _axG.set_xlabel("")
    _axG.tick_params(axis="x", rotation=0)
    _format_ticks(_axG)
    (_handles, _labels) = _axG.get_legend_handles_labels()
    _new_labels = ["Basal", "Transitional", "High-Signature"]
    _axG.legend(
        _handles,
        _new_labels,
        bbox_to_anchor=(0.5, -0.25),
        loc="upper center",
        title="Mutational Signatures",
        ncol=3,
        fontsize=_LEGEND_FS,
        title_fontsize=_LEGEND_FS,
        frameon=False,
    )

    # --- FINAL LAYOUT ADJUSTMENTS ---
    all_axes = [_axA, _axB, _axC, _axD, _axE, _axF, _axG]
    panel_letters = ["A", "B", "C", "D", "E", "F", "G"]
    for _ax, _lbl in zip(all_axes, panel_letters):
        _ax.text(
            -0.08,
            1.15,
            _lbl,
            transform=_ax.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="top",
            ha="right",
        )
        sns.despine(ax=_ax)

    plt.tight_layout()
    plt.savefig(_OUTPUT_FIG, dpi=300, bbox_inches="tight")
    plt.savefig(_OUTPUT_SVG, format="svg", bbox_inches="tight")
    _f_rep.close()
    print(f"Done. Saved as {_OUTPUT_FIG} and {_OUTPUT_SVG}")
    plt.show()
    return


if __name__ == "__main__":
    app.run()
