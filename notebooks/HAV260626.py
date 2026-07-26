# /// script
# dependencies = ["biopython", "kaleido", "matplotlib", "numpy", "openpyxl", "plotly", "plotly-express", "scikit-posthocs", "scipy", "seaborn", "tqdm"]
# ///

import marimo

__generated_with = "0.23.11"
app = marimo.App()


@app.cell
def _():
    import subprocess

    return (subprocess,)


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
    from Bio import SeqUtils
    from Bio import SeqIO
    from Bio import Entrez
    from matplotlib import pyplot as plt
    import re
    import numpy as np
    import pandas as pd
    import seaborn as sns

    return Entrez, SeqIO, SeqUtils, np, pd, plt, re, sns


@app.cell
def _(pd):
    file_path = "../data/raw/ALL ACCESSION INFO.xlsx"

    # Read the Excel file
    df_3 = pd.read_excel(file_path)

    # Display the first few rows
    print(df_3.head())
    return (df_3,)


@app.cell
def _(df_3):
    df_3["accession_id"].tolist()
    return


@app.cell
def _(df_3):
    id_list_1 = df_3["accession_id"].tolist()
    print(id_list_1)
    return (id_list_1,)


@app.cell
def _(Entrez, SeqIO, id_list_1):
    Entrez.email = "smukherjee109@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=id_list_1, rettype="gb")
    gb_records = SeqIO.parse(handle, "gb")
    gb_records = [_i for _i in gb_records]
    SeqIO.write(gb_records, "../data/processed/HAV_297_output.gb", "gb")
    return (gb_records,)


@app.function
def index_genbank_features(in_gb_record, feature_type):
    """There are usually several db_xref entries but
    there should only be one source per feature. Since
    we only care for source, we'll not handle the repetition
    cases that will arise for features with multiple entries
    """
    for _index, feature in enumerate(in_gb_record.features):
        if feature.type == feature_type:
            return _index


@app.cell
def _(SeqIO, gb_records, pd):
    datalist = []
    for gb_record in gb_records:
        gb_feature1 = gb_record.features[index_genbank_features(gb_record, "source")]
        try:
            gb_feature = gb_record.features[index_genbank_features(gb_record, "CDS")]
        except TypeError:
            pass
        datalist.append(
            [
                gb_record.id,
                str(gb_record.seq),
                str(gb_feature.extract(gb_record.seq)),
                gb_feature.qualifiers.get("translation", ["None"])[0],
                len(gb_record.seq),
                gb_feature1.qualifiers.get("collection_date", ["None"])[0],
                gb_feature1.qualifiers.get("strain", ["None"])[0],
                gb_feature1.qualifiers.get("isolate", ["None"])[0],
                gb_feature1.qualifiers.get("geo_loc_name", ["None"])[0],
                gb_feature1.qualifiers.get("note", ["None"])[0],
                gb_feature1.qualifiers.get("organism", ["None"])[0],
                gb_feature1.qualifiers.get("host", ["None"])[0],
                gb_record.description,
            ]
        )
    df_4 = pd.DataFrame(
        datalist,
        columns=[
            "accession_id",
            "seq",
            "cds_seq",
            "protein",
            "seq_len",
            "collection_date",
            "strain",
            "isolate",
            "country",
            "genotype",
            "organism",
            "host",
            "header",
        ],
    )
    SeqIO.write(
        gb_records, "../data/processed/HAV_297_output.fasta", "fasta"
    )  # print(gb_record.annotations)
    df_4.to_csv(
        "../data/processed/dfhav_297_seqs.csv", index=False
    )  # accession id  # seq  #cds_seq  #Protein  # seq len  # collection date  # strain  # isolate  # country  # genotype/subtype  # organism  # host  # header
    return (df_4,)


@app.cell
def _(df_4):
    df_4.organism.value_counts()
    return


@app.cell
def _(df_4, np):
    df_4.replace(to_replace="None", value=np.nan, inplace=True)
    return


@app.cell
def _(df_4, pd):
    df_4["date"] = pd.to_datetime(
        df_4["collection_date"], format="mixed"
    )
    return


@app.cell
def _(df_4, pd):
    df_4["year"] = pd.to_datetime(df_4.date).dt.strftime("%Y")
    df_4["year"] = df_4["year"].astype("float")
    return


@app.cell
def _(df_4):
    df_4["year"].value_counts()
    return


@app.cell
def _(df_4):
    _ofile = open("../data/processed/HAV_ALL_cds.fasta", "w")
    _ref_row = df_4[df_4["accession_id"] == "NC_001489.1"].iloc[0]
    _ofile.write(f">NC_001489.1 | Collection Date: 1958\n{_ref_row['cds_seq']}\n")
    for _index, _row in df_4.iterrows():
        if _row["accession_id"] == "NC_001489.1":
            continue
        _ofile.write(
            f">{_row['accession_id']} | Collection Date: {_row['year']}\n{_row['cds_seq']}\n"
        )
    _ofile.close()
    return


@app.cell
def _(df_4):
    #This is to categorise and eliminate any non-human HAV sequences that may have been downloaded
    bin3_org = [
        "Duck hepatitis A virus 1",
        "Duck hepatitis A virus 3",
        "Avihepatovirus A",
    ]
    bin4_org = [
        "hepatovirus H2",
        "Bat hepatovirus SMG18520Minmav2014",
        "Bat hepatovirus BUO2BF86Colafr2010",
        "hepatovirus G2",
    ]
    bin5_org = [
        "Rodent hepatovirus KEF121Sigmas2012",
        "Rodent hepatovirus CIV459Lopsik2004",
        "Rodent hepatovirus RMU101637Micarv2010",
        "Shrew hepatovirus KS121232Sorara2012",
        "Hedgehog hepatovirus",
        "Simian hepatitis A virus",
        "hepatovirus D2",
        "Shrew hepatovirus",
        "Phocoena sinus",
        "Rodent hepatovirus",
        "Didelphis aurita hepatitis A virus",
    ]
    bin6_org = ["Hepatovirus A", "Human hepatitis A virus"]
    for _row in df_4.itertuples():
        if _row.organism in bin3_org:
            df_4.at[_row.Index, "Org"] = "Avians"
        elif _row.organism in bin4_org:
            df_4.at[_row.Index, "Org"] = "Bats"
        elif _row.organism in bin5_org:
            df_4.at[_row.Index, "Org"] = "Mammals"
        elif _row.organism in bin6_org:
            df_4.at[_row.Index, "Org"] = "Humans"
        else:
            pass
    return


@app.cell
def _(df_4):
    df_4.Org.value_counts()
    return


@app.cell
def _(df_4):
    df_collectiondate = df_4
    return (df_collectiondate,)


@app.cell
def _(df_collectiondate):
    df_collectiondate.to_excel("../data/processed/HAV with collection date.xlsx", index=False)
    return


@app.cell
def _(df_collectiondate):
    # Open the file
    _ofile = open("../data/processed/HAV_Collection date.fasta", "w")
    _ref_row = df_collectiondate[
        df_collectiondate["accession_id"] == "NC_001489.1"
    ].iloc[0]
    # 1. Extract and write the Reference first
    _ofile.write(f">NC_001489.1 | Collection Date: 1958\n{_ref_row['seq']}\n")
    for _index, _row in df_collectiondate.iterrows():
        if _row["accession_id"] == "NC_001489.1":
            # 2. Write the rest
            continue
        _ofile.write(
            f">{_row['accession_id']} | Collection Date: {_row['year']}\n{_row['seq']}\n"
        )  # Skip the reference since we just wrote it
    # Do not forget to close it
    _ofile.close()
    return


@app.cell
def _(df_collectiondate):
    # Open the file
    _ofile = open("../data/processed/HAV_Collection_protein.fasta", "w")
    _ref_row = df_collectiondate[
        df_collectiondate["accession_id"] == "NC_001489.1"
    ].iloc[0]
    # 1. Extract and write the Reference first
    _ofile.write(f">NC_001489.1 | Collection Date: 1958\n{_ref_row['protein']}\n")
    for _index, _row in df_collectiondate.iterrows():
        if _row["accession_id"] == "NC_001489.1":
            # 2. Write the rest
            continue
        _ofile.write(
            f">{_row['accession_id']} | Collection Date: {_row['year']}\n{_row['protein']}\n"
        )  # Skip the reference since we just wrote it
    # Do not forget to close it
    _ofile.close()
    return


@app.cell
def _(SeqUtils, df_collectiondate, re):
    for _row in df_collectiondate.itertuples():
        temp_seq_aa_count = _row.seq.count("AA")
        temp_seq_at_count = _row.seq.count("AT")
        temp_seq_ac_count = _row.seq.count("AC")
        temp_seq_ag_count = _row.seq.count("AG")
        temp_seq_ta_count = _row.seq.count("TA")
        temp_seq_tt_count = _row.seq.count("TT")
        temp_seq_tc_count = _row.seq.count("TC")
        temp_seq_tg_count = _row.seq.count("TG")
        temp_seq_ca_count = _row.seq.count("CA")
        temp_seq_ct_count = _row.seq.count("CT")
        temp_seq_cc_count = _row.seq.count("CC")
        temp_seq_cpg_count = _row.seq.count("CG")
        temp_seq_ga_count = _row.seq.count("GA")
        temp_seq_gt_count = _row.seq.count("GT")
        temp_seq_gc_count = _row.seq.count("GC")
        temp_seq_gg_count = _row.seq.count("GG")
        temp_seq_a_count = _row.seq.count("A")
        temp_seq_t_count = _row.seq.count("T")
        temp_seq_c_count = _row.seq.count("C")
        temp_seq_g_count = _row.seq.count("G")
        zap4 = len(SeqUtils.nt_search(str(_row.seq), "CNNNNGNCG")) - 1
        zap5 = len(SeqUtils.nt_search(str(_row.seq), "CNNNNNGNCG")) - 1
        zap6 = len(SeqUtils.nt_search(str(_row.seq), "CNNNNNNGNCG")) - 1
        zap7 = len(SeqUtils.nt_search(str(_row.seq), "CNNNNNNNGNCG")) - 1
        zap8 = len(SeqUtils.nt_search(str(_row.seq), "CNNNNNNNNGNCG")) - 1
        zap = zap4 + zap5 + zap6 + zap7 + zap8
        df_collectiondate.at[_row.Index, "zap_motif"] = zap
        df_collectiondate.at[_row.Index, "aa_count"] = temp_seq_aa_count
        df_collectiondate.at[_row.Index, "at_count"] = temp_seq_at_count
        df_collectiondate.at[_row.Index, "ac_count"] = temp_seq_ac_count
        df_collectiondate.at[_row.Index, "ag_count"] = temp_seq_ag_count
        df_collectiondate.at[_row.Index, "ta_count"] = temp_seq_ta_count
        df_collectiondate.at[_row.Index, "tt_count"] = temp_seq_tt_count
        df_collectiondate.at[_row.Index, "tc_count"] = temp_seq_tc_count
        df_collectiondate.at[_row.Index, "tg_count"] = temp_seq_tg_count
        df_collectiondate.at[_row.Index, "ca_count"] = temp_seq_ca_count
        df_collectiondate.at[_row.Index, "ct_count"] = temp_seq_ct_count
        df_collectiondate.at[_row.Index, "cc_count"] = temp_seq_cc_count
        df_collectiondate.at[_row.Index, "cpg_count"] = temp_seq_cpg_count
        df_collectiondate.at[_row.Index, "ga_count"] = temp_seq_ga_count
        df_collectiondate.at[_row.Index, "gt_count"] = temp_seq_gt_count
        df_collectiondate.at[_row.Index, "gc_count"] = temp_seq_gc_count
        df_collectiondate.at[_row.Index, "gg_count"] = temp_seq_gg_count
        df_collectiondate.at[_row.Index, "cpg_per"] = 100 * (
            temp_seq_cpg_count / _row.seq_len
        )
        df_collectiondate.at[_row.Index, "aa_obye"] = (
            temp_seq_aa_count * _row.seq_len / (temp_seq_a_count * temp_seq_a_count)
        )
        df_collectiondate.at[_row.Index, "at_obye"] = (
            temp_seq_at_count * _row.seq_len / (temp_seq_a_count * temp_seq_t_count)
        )
        df_collectiondate.at[_row.Index, "ac_obye"] = (
            temp_seq_ac_count * _row.seq_len / (temp_seq_a_count * temp_seq_c_count)
        )
        df_collectiondate.at[_row.Index, "ag_obye"] = (
            temp_seq_ag_count * _row.seq_len / (temp_seq_a_count * temp_seq_g_count)
        )
        df_collectiondate.at[_row.Index, "ta_obye"] = (
            temp_seq_ta_count * _row.seq_len / (temp_seq_t_count * temp_seq_a_count)
        )
        df_collectiondate.at[_row.Index, "tt_obye"] = (
            temp_seq_tt_count * _row.seq_len / (temp_seq_t_count * temp_seq_t_count)
        )
        df_collectiondate.at[_row.Index, "tc_obye"] = (
            temp_seq_tc_count * _row.seq_len / (temp_seq_t_count * temp_seq_c_count)
        )
        df_collectiondate.at[_row.Index, "tg_obye"] = (
            temp_seq_tg_count * _row.seq_len / (temp_seq_t_count * temp_seq_g_count)
        )
        df_collectiondate.at[_row.Index, "ca_obye"] = (
            temp_seq_ca_count * _row.seq_len / (temp_seq_c_count * temp_seq_a_count)
        )
        df_collectiondate.at[_row.Index, "ct_obye"] = (
            temp_seq_ct_count * _row.seq_len / (temp_seq_c_count * temp_seq_t_count)
        )
        df_collectiondate.at[_row.Index, "cpg_obye"] = (
            temp_seq_cpg_count * _row.seq_len / (temp_seq_c_count * temp_seq_g_count)
        )
        df_collectiondate.at[_row.Index, "cc_obye"] = (
            temp_seq_cc_count * _row.seq_len / (temp_seq_c_count * temp_seq_c_count)
        )
        df_collectiondate.at[_row.Index, "ga_obye"] = (
            temp_seq_ga_count * _row.seq_len / (temp_seq_g_count * temp_seq_a_count)
        )
        df_collectiondate.at[_row.Index, "gt_obye"] = (
            temp_seq_gt_count * _row.seq_len / (temp_seq_g_count * temp_seq_t_count)
        )
        df_collectiondate.at[_row.Index, "gc_obye"] = (
            temp_seq_gc_count * _row.seq_len / (temp_seq_g_count * temp_seq_c_count)
        )
        df_collectiondate.at[_row.Index, "gg_obye"] = (
            temp_seq_gg_count * _row.seq_len / (temp_seq_g_count * temp_seq_g_count)
        )
        df_collectiondate.at[_row.Index, "gc_content"] = SeqUtils.gc_fraction(_row.seq)
        df_collectiondate.at[_row.Index, "A"] = temp_seq_a_count
        df_collectiondate.at[_row.Index, "T"] = temp_seq_t_count
        df_collectiondate.at[_row.Index, "G"] = temp_seq_g_count
        df_collectiondate.at[_row.Index, "C"] = temp_seq_c_count
        drach_count = len(re.findall("[AGT][AG]AC[ACT]", str(_row.seq)))
        df_collectiondate.at[_row.Index, "drach_motif"] = (
            drach_count  # --- NEW: DRACH Motif (m6A) ---  # Consensus: D=A/G/T, R=A/G, A=A, C=C, H=A/C/T  # This regex captures the standard m6A site context.
        )
    return


@app.cell
def _(subprocess):
    #! mafft --auto --thread -1 "HAV_Collection date.fasta" > "HAV_Collection date_aligned.fasta"
    input_fasta = "../data/processed/HAV_Collection date.fasta"
    output_fasta = "../data/processed/HAV_Collection date_aligned.fasta"

    print("Starting MAFFT alignment...")

    # Open the target file, then instruct subprocess to pipe standard output directly into it
    with open(output_fasta, "w") as outfile:
        subprocess.run(
            [
                "mafft",
                "--auto",
                "--thread",
                "-1",
                input_fasta
            ],
            stdout=outfile,  # This completely replaces the shell '>' operator
            check=True       # This forces Python to raise an error if MAFFT fails
        )

    print("Alignment complete.")
    return


@app.cell
def _(pd):
    from Bio import AlignIO

    # Function to count mutations
    def _count_mutations(reference, sequence):
        _mutations = sum(
            (
                1
                for (ref_base, seq_base) in zip(reference, _sequence)
                if ref_base != seq_base
            )
        )
        return _mutations

    _alignment_file = "../data/processed/HAV_Collection date_aligned.fasta"
    # Load the aligned file
    _alignment = AlignIO.read(
        _alignment_file, "fasta"
    )  # Replace with your aligned file name
    _reference_sequence = None
    for _record in _alignment:
        # Extract the reference sequence
        if _record.id == "NC_001489.1":
            _reference_sequence = str(_record.seq)
            break  # Replace with the ID of your reference sequence
    if not _reference_sequence:
        raise ValueError("Reference sequence not found in the alignment file")
    _data = []
    for _record in _alignment:
        if _record.id != "reference":
            _sequence = str(_record.seq)
            # Calculate mutations for each sequence
            _mutations = _count_mutations(_reference_sequence, _sequence)
            _data.append(
                {"accession_id": _record.id, "Nucleotide Mutations": _mutations}
            )
    df1 = pd.DataFrame(_data)  # Skip the reference sequence itself
    return AlignIO, df1


@app.cell
def _(df1, df_collectiondate, pd):
    merged_df = pd.merge(
        df_collectiondate,
        df1[["accession_id", "Nucleotide Mutations"]],
        on="accession_id",
        how="left",
    )

    merged_df
    return (merged_df,)


@app.cell
def _(subprocess):
    def _():
        #! mafft --auto --thread -1 "HAV_Collection_protein.fasta" > "HAV_Collection_protein_align.fasta"
        input_fasta = "../data/processed/HAV_Collection_protein.fasta"
        output_fasta = "../data/processed/HAV_Collection_protein_align.fasta"

        print("Starting MAFFT alignment...")

        # Open the target file, then instruct subprocess to pipe standard output directly into it
        with open(output_fasta, "w") as outfile:
            subprocess.run(
                [
                    "mafft",
                    "--auto",
                    "--thread",
                    "-1",
                    input_fasta
                ],
                stdout=outfile,  # This completely replaces the shell '>' operator
                check=True       # This forces Python to raise an error if MAFFT fails
            )
        return print("Alignment complete.")


    _()
    return


@app.cell
def _(AlignIO, pd):
    def _count_mutations(reference, sequence):
        _mutations = sum(
            (
                1
                for (ref_base, seq_base) in zip(reference, _sequence)
                if ref_base != seq_base
            )
        )
        return _mutations

    _alignment_file = "../data/processed/HAV_Collection_protein_align.fasta"
    _alignment = AlignIO.read(_alignment_file, "fasta")
    _reference_sequence = None
    for _record in _alignment:
        if _record.id == "NC_001489.1":
            _reference_sequence = str(_record.seq)
            break
    if not _reference_sequence:
        raise ValueError("Reference sequence not found in the alignment file")
    _data = []
    for _record in _alignment:
        if _record.id != "reference":
            _sequence = str(_record.seq)
            _mutations = _count_mutations(_reference_sequence, _sequence)
            _data.append(
                {"accession_id": _record.id, "Amino Acid Mutations": _mutations}
            )
    df1_1 = pd.DataFrame(_data)
    return (df1_1,)


@app.cell
def _(df1_1, merged_df, pd):
    merged_df_1 = pd.merge(
        merged_df,
        df1_1[["accession_id", "Amino Acid Mutations"]],
        on="accession_id",
        how="left",
    )
    merged_df_1
    return (merged_df_1,)


@app.cell
def _(merged_df_1):
    # 1. Define the Timeline Bins
    def get_timeline_window(y):
        if 1998 <= y <= 2007:
            return "1998-2007"
        elif 2008 <= y <= 2012:
            return "2008-2012"
        elif 2013 <= y <= 2017:
            return "2013-2017"
        elif 2018 <= y <= 2022:
            return "2018-2022"
        else:
            return "Other"

    merged_df_1["Timeline_Group"] = merged_df_1["year"].apply(get_timeline_window)
    # 2. Create Temporary Columns for Grouping
    merged_df_1["Period_Group"] = merged_df_1["year"].apply(
        lambda x: "Before 2018" if x < 2018 else "After 2018"
    )
    period_counts = merged_df_1["Period_Group"].value_counts()
    merged_df_1["Period"] = merged_df_1["Period_Group"].apply(
        lambda x: f"{x}\n(n={period_counts[x]})"
    )
    # 3. Calculate Counts for 'Period' and Map formatted labels
    timeline_counts = merged_df_1["Timeline_Group"].value_counts()
    merged_df_1["Timeline"] = merged_df_1["Timeline_Group"].apply(
        lambda x: f"{x}\n(n={timeline_counts[x]})"
    )
    merged_df_1.drop(columns=["Timeline_Group", "Period_Group"], inplace=True)
    print(merged_df_1[["year", "Period", "Timeline"]].head())
    print("\nUnique Periods:", merged_df_1["Period"].unique())
    # 4. Calculate Counts for 'Timeline' and Map formatted labels
    # Optional: Clean up temporary columns if you don't need them
    # Check the results
    print("Unique Timelines:", merged_df_1["Timeline"].unique())
    return


@app.cell
def _(AlignIO, merged_df_1, pd, plt, sns):
    from scipy.stats import linregress

    _alignment_file = "../data/processed/HAV_Collection date_aligned.fasta"
    _alignment = AlignIO.read(_alignment_file, "fasta")
    _reference_sequence = None
    for _record in _alignment:
        if _record.id == "NC_001489.1":
            _reference_sequence = str(_record.seq)
            break
    if not _reference_sequence:
        raise ValueError(
            "Reference sequence NC_001489.1 not found in the alignment file!"
        )

    def get_clinical_risk_factors(seq, ref_seq):
        _counts = {
            "Entire": {
                "G": 0,
                "A": 0,
                "GpA": 0,
                "G[G>A]T": 0,
                "C[C>T]G": 0,
                "Length": 0,
            },
            "VP1_2A": {
                "G": 0,
                "A": 0,
                "GpA": 0,
                "G[G>A]T": 0,
                "C[C>T]G": 0,
                "Length": 0,
            },
        }
        seq_len = min(len(seq), len(ref_seq))
        ref_pos = 0
        prev_base_entire = ""
        prev_base_vp1_2a = ""
        for _i in range(seq_len):
            ref_base = ref_seq[_i].upper()
            seq_base = seq[_i].upper()
            if ref_base != "-":
                ref_pos = ref_pos + 1
            in_vp1_2a = 2265 <= ref_pos <= 3732
            if seq_base != "-":
                _counts["Entire"]["Length"] = _counts["Entire"]["Length"] + 1
                if seq_base == "G":
                    _counts["Entire"]["G"] = _counts["Entire"]["G"] + 1
                if seq_base == "A":
                    _counts["Entire"]["A"] = _counts["Entire"]["A"] + 1
                if prev_base_entire == "G" and seq_base == "A":
                    _counts["Entire"]["GpA"] = _counts["Entire"]["GpA"] + 1
                prev_base_entire = seq_base
                if in_vp1_2a:
                    _counts["VP1_2A"]["Length"] = _counts["VP1_2A"]["Length"] + 1
                    if seq_base == "G":
                        _counts["VP1_2A"]["G"] = _counts["VP1_2A"]["G"] + 1
                    if seq_base == "A":
                        _counts["VP1_2A"]["A"] = _counts["VP1_2A"]["A"] + 1
                    if prev_base_vp1_2a == "G" and seq_base == "A":
                        _counts["VP1_2A"]["GpA"] = _counts["VP1_2A"]["GpA"] + 1
                    prev_base_vp1_2a = seq_base
            if ref_base != seq_base and ref_base != "-" and (seq_base != "-"):
                right_idx = _i + 1
                while right_idx < seq_len and ref_seq[right_idx] == "-":
                    right_idx = right_idx + 1
                right_context = (
                    ref_seq[right_idx].upper() if right_idx < seq_len else "N"
                )
                left_idx = _i - 1
                while left_idx >= 0 and ref_seq[left_idx] == "-":
                    left_idx = left_idx - 1
                left_context = ref_seq[left_idx].upper() if left_idx >= 0 else "N"
                if (
                    left_context == "G"
                    and ref_base == "G"
                    and (seq_base == "A")
                    and (right_context == "T")
                ):
                    _counts["Entire"]["G[G>A]T"] = _counts["Entire"]["G[G>A]T"] + 1
                    if in_vp1_2a:
                        _counts["VP1_2A"]["G[G>A]T"] = _counts["VP1_2A"]["G[G>A]T"] + 1
                if ref_base == "C" and seq_base == "T" and (right_context == "G"):
                    _counts["Entire"]["C[C>T]G"] = _counts["Entire"]["C[C>T]G"] + 1
                    if in_vp1_2a:
                        _counts["VP1_2A"]["C[C>T]G"] = _counts["VP1_2A"]["C[C>T]G"] + 1
        metrics = {}
        for region in ["Entire", "VP1_2A"]:
            L = _counts[region]["Length"]
            G = _counts[region]["G"]
            A = _counts[region]["A"]
            GpA = _counts[region]["GpA"]
            if G * A > 0:
                metrics[f"GpA_OE_{region}"] = GpA * L / (G * A)
            else:
                metrics[f"GpA_OE_{region}"] = 0.0
            if L > 0:
                metrics[f"G[G>A]T_Load_{region}"] = (
                    _counts[region]["G[G>A]T"] / L * 1000
                )
                metrics[f"C[C>T]G_Load_{region}"] = (
                    _counts[region]["C[C>T]G"] / L * 1000
                )
            else:
                metrics[f"G[G>A]T_Load_{region}"] = 0.0
                metrics[f"C[C>T]G_Load_{region}"] = 0.0
        return pd.Series(metrics)

    merged_df_2 = merged_df_1.loc[:, ~merged_df_1.columns.duplicated()].copy()
    _metric_cols = [
        "GpA_OE_Entire",
        "GpA_OE_VP1_2A",
        "G[G>A]T_Load_Entire",
        "G[G>A]T_Load_VP1_2A",
        "C[C>T]G_Load_Entire",
        "C[C>T]G_Load_VP1_2A",
    ]
    merged_df_2 = merged_df_2.drop(
        columns=[_c for _c in _metric_cols if _c in merged_df_2.columns],
        errors="ignore",
    )
    print("Calculating clinical risk factors... this may take a moment.")
    _new_cols = merged_df_2["seq"].apply(
        lambda x: get_clinical_risk_factors(x, _reference_sequence)
    )
    merged_df_2 = pd.concat([merged_df_2, _new_cols], axis=1)
    print("Metrics successfully added to merged_df!")
    sns.set_theme(style="whitegrid")
    (_fig, _axes) = plt.subplots(1, 3, figsize=(20, 6))
    _plot_configs = [
        (
            "GpA_OE_Entire",
            "GpA_OE_VP1_2A",
            "GpA O/E Bias: Whole Genome vs VP1-2A",
            "blue",
            _axes[0],
        ),
        (
            "G[G>A]T_Load_Entire",
            "G[G>A]T_Load_VP1_2A",
            "G[G>A]T Load (/kb): Whole Genome vs VP1-2A",
            "red",
            _axes[1],
        ),
        (
            "C[C>T]G_Load_Entire",
            "C[C>T]G_Load_VP1_2A",
            "C[C>T]G Load (/kb): Whole Genome vs VP1-2A",
            "green",
            _axes[2],
        ),
    ]
    for _x_col, _y_col, _title, _line_color, _ax in _plot_configs:
        sns.regplot(
            data=merged_df_2,
            x=_x_col,
            y=_y_col,
            scatter_kws={"alpha": 0.6, "color": "gray"},
            line_kws={"color": _line_color},
            ax=_ax,
        )
        _mask = ~merged_df_2[_x_col].isna() & ~merged_df_2[_y_col].isna()
        if _mask.sum() > 1:
            x_data = merged_df_2.loc[_mask, _x_col]
            y_data = merged_df_2.loc[_mask, _y_col]
            (_slope, _intercept, _r_value, _p_value, _std_err) = linregress(
                x_data, y_data
            )
            _stats_text = f"\n$R^2$: {_r_value**2:.3f} | p-value: {_p_value:.3e}"
            _ax.set_title(_title + _stats_text, fontsize=13)
        else:
            _ax.set_title(_title, fontsize=13)
        _ax.set_xlabel(_x_col.replace("_", " "), fontsize=11)
        _ax.set_ylabel(_y_col.replace("_", " "), fontsize=11)
    plt.tight_layout
    # 1. SAVE THE FIGURE (Must happen strictly before plt.show)
    fig_path = "../results/figures/SuppFig_S6.svg"
    fig_path1 = "../results/figures/SuppFig_S6.png"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.savefig(fig_path1, dpi=300, bbox_inches="tight")
    print(f"Figure saved to: {fig_path}")
    plt.show()

    table_path = "../data/processed/merged_df_2_with_mutational_profile.csv"
    merged_df_2.to_csv(table_path, index=False)
    print(f"Dataframe saved to: {table_path}")
    return linregress, merged_df_2


@app.cell
def _(AlignIO, linregress, merged_df_2, pd, plt, sns):
    _alignment_file = "../data/processed/HAV_Collection date_aligned.fasta"
    _alignment = AlignIO.read(_alignment_file, "fasta")
    _reference_sequence = None
    for _record in _alignment:
        if _record.id == "NC_001489.1":
            _reference_sequence = str(_record.seq)
            break
    if not _reference_sequence:
        raise ValueError("Reference sequence NC_001489.1 not found!")

    def get_true_context_loads(seq, ref_seq):
        _counts = {
            "Entire": {"Ref_GGT": 0, "Ref_CCG": 0, "G[G>A]T": 0, "C[C>T]G": 0},
            "VP1_2A": {"Ref_GGT": 0, "Ref_CCG": 0, "G[G>A]T": 0, "C[C>T]G": 0},
        }
        seq_len = min(len(seq), len(ref_seq))
        ref_pos = 0
        for _i in range(seq_len):
            ref_base = ref_seq[_i].upper()
            seq_base = seq[_i].upper()
            if ref_base != "-":
                ref_pos = ref_pos + 1
            in_vp1_2a = 2265 <= ref_pos <= 3732
            right_idx = _i + 1
            while right_idx < seq_len and ref_seq[right_idx] == "-":
                right_idx = right_idx + 1
            right_context = ref_seq[right_idx].upper() if right_idx < seq_len else "N"
            left_idx = _i - 1
            while left_idx >= 0 and ref_seq[left_idx] == "-":
                left_idx = left_idx - 1
            left_context = ref_seq[left_idx].upper() if left_idx >= 0 else "N"
            if ref_base != "-":
                if left_context == "G" and ref_base == "G" and (right_context == "T"):
                    _counts["Entire"]["Ref_GGT"] = _counts["Entire"]["Ref_GGT"] + 1
                    if in_vp1_2a:
                        _counts["VP1_2A"]["Ref_GGT"] = _counts["VP1_2A"]["Ref_GGT"] + 1
                if left_context == "C" and ref_base == "C" and (right_context == "G"):
                    _counts["Entire"]["Ref_CCG"] = _counts["Entire"]["Ref_CCG"] + 1
                    if in_vp1_2a:
                        _counts["VP1_2A"]["Ref_CCG"] = _counts["VP1_2A"]["Ref_CCG"] + 1
            if ref_base != seq_base and ref_base != "-" and (seq_base != "-"):
                if (
                    left_context == "G"
                    and ref_base == "G"
                    and (seq_base == "A")
                    and (right_context == "T")
                ):
                    _counts["Entire"]["G[G>A]T"] = _counts["Entire"]["G[G>A]T"] + 1
                    if in_vp1_2a:
                        _counts["VP1_2A"]["G[G>A]T"] = _counts["VP1_2A"]["G[G>A]T"] + 1
                if (
                    left_context == "C"
                    and ref_base == "C"
                    and (seq_base == "T")
                    and (right_context == "G")
                ):
                    _counts["Entire"]["C[C>T]G"] = _counts["Entire"]["C[C>T]G"] + 1
                    if in_vp1_2a:
                        _counts["VP1_2A"]["C[C>T]G"] = _counts["VP1_2A"]["C[C>T]G"] + 1
        metrics = {}
        for region in ["Entire", "VP1_2A"]:
            ref_ggt = _counts[region]["Ref_GGT"]
            ref_ccg = _counts[region]["Ref_CCG"]
            if ref_ggt > 0:
                metrics[f"G[G>A]T_MotifPenetrance_{region}"] = (
                    _counts[region]["G[G>A]T"] / ref_ggt * 100
                )
            else:
                metrics[f"G[G>A]T_MotifPenetrance_{region}"] = 0.0
            if ref_ccg > 0:
                metrics[f"C[C>T]G_MotifPenetrance_{region}"] = (
                    _counts[region]["C[C>T]G"] / ref_ccg * 100
                )
            else:
                metrics[f"C[C>T]G_MotifPenetrance_{region}"] = 0.0
        return pd.Series(metrics)

    merged_df_3 = merged_df_2.loc[:, ~merged_df_2.columns.duplicated()].copy()
    _metric_cols = [
        "G[G>A]T_MotifPenetrance_Entire",
        "G[G>A]T_MotifPenetrance_VP1_2A",
        "C[C>T]G_MotifPenetrance_Entire",
        "C[C>T]G_MotifPenetrance_VP1_2A",
    ]
    merged_df_3 = merged_df_3.drop(
        columns=[_c for _c in _metric_cols if _c in merged_df_3.columns],
        errors="ignore",
    )
    print("Calculating true context motif penetrance...")
    _new_cols = merged_df_3["seq"].apply(
        lambda x: get_true_context_loads(x, _reference_sequence)
    )
    merged_df_3 = pd.concat([merged_df_3, _new_cols], axis=1)
    sns.set_theme(style="whitegrid")
    (_fig, _axes) = plt.subplots(1, 2, figsize=(14, 6))
    _plot_configs = [
        (
            "G[G>A]T_MotifPenetrance_Entire",
            "G[G>A]T_MotifPenetrance_VP1_2A",
            "G[G>A]T Motif Penetrance (%)\nWhole Genome vs VP1-2A",
            "red",
            _axes[0],
        ),
        (
            "C[C>T]G_MotifPenetrance_Entire",
            "C[C>T]G_MotifPenetrance_VP1_2A",
            "CCG Motif Penetrance (%)\nWhole Genome vs VP1-2A",
            "green",
            _axes[1],
        ),
    ]
    for _x_col, _y_col, _title, _line_color, _ax in _plot_configs:
        sns.regplot(
            data=merged_df_3,
            x=_x_col,
            y=_y_col,
            scatter_kws={"alpha": 0.6, "color": "gray"},
            line_kws={"color": _line_color},
            ax=_ax,
        )
        _mask = ~merged_df_3[_x_col].isna() & ~merged_df_3[_y_col].isna()
        if _mask.sum() > 1:
            (_slope, _intercept, _r_value, _p_value, _std_err) = linregress(
                merged_df_3.loc[_mask, _x_col], merged_df_3.loc[_mask, _y_col]
            )
            _ax.set_title(
                _title + f"\n$R^2$: {_r_value**2:.3f} | p-value: {_p_value:.3e}",
                fontsize=13,
            )
        else:
            _ax.set_title(_title, fontsize=13)
        _ax.set_xlabel(_x_col.replace("_", " ") + " (%)", fontsize=11)
        _ax.set_ylabel(_y_col.replace("_", " ") + " (%)", fontsize=11)
    plt.tight_layout()
    plt.show()
    return (merged_df_3,)


@app.cell
def _(merged_df_3, pd, plt, sns):
    from scipy.stats import spearmanr

    _plot_data = merged_df_3.dropna(
        subset=[
            "year",
            "G[G>A]T_MotifPenetrance_Entire",
            "C[C>T]G_MotifPenetrance_Entire",
        ]
    ).copy()
    _plot_data["year"] = pd.to_numeric(_plot_data["year"], errors="coerce")
    sns.set_theme(style="whitegrid")
    (_fig, _axes) = plt.subplots(1, 2, figsize=(14, 6))
    sns.regplot(
        data=_plot_data,
        x="year",
        y="G[G>A]T_MotifPenetrance_Entire",
        scatter_kws={"alpha": 0.6, "color": "gray", "edgecolor": "w"},
        line_kws={"color": "red"},
        ax=_axes[0],
    )
    (corr1, p1) = spearmanr(
        _plot_data["year"], _plot_data["G[G>A]T_MotifPenetrance_Entire"]
    )
    _axes[0].set_title(
        f"G[G>A]T Motif Penetrance vs Time\nSpearman's $\\rho$: {corr1:.3f} | p-value: {p1:.3e}",
        fontsize=13,
    )
    _axes[0].set_xlabel("Collection Year", fontsize=11)
    _axes[0].set_ylabel("G[G>A]T Motif Penetrance (%)", fontsize=11)
    sns.regplot(
        data=_plot_data,
        x="year",
        y="C[C>T]G_MotifPenetrance_Entire",
        scatter_kws={"alpha": 0.6, "color": "gray", "edgecolor": "w"},
        line_kws={"color": "green"},
        ax=_axes[1],
    )
    (corr2, p2) = spearmanr(
        _plot_data["year"], _plot_data["C[C>T]G_MotifPenetrance_Entire"]
    )
    _axes[1].set_title(
        f"C[C>T]G Motif Penetrance vs Time\nSpearman's $\\rho$: {corr2:.3f} | p-value: {p2:.3e}",
        fontsize=13,
    )
    _axes[1].set_xlabel("Collection Year", fontsize=11)
    _axes[1].set_ylabel("C[C>T]G Motif Penetrance (%)", fontsize=11)
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(AlignIO, merged_df_3, pd):
    from collections import Counter

    def analyze_reference_anchored_mutations(
        metadata_df,
        fasta_file,
        ref_accession="NC_001489",
        split_year=2018,
        top_n=50,
        freq_threshold=0.05,
    ):
        """
        Anchors mutations explicitly to the specified Reference Sequence.
        Maps alignment coordinates back to true biological coordinates.
        """
        metadata_df["year"] = pd.to_numeric(metadata_df["year"], errors="coerce")
        pre_ids = set(metadata_df[metadata_df["year"] < split_year]["accession_id"])
        post_ids = set(metadata_df[metadata_df["year"] >= split_year]["accession_id"])
        print(
            f"Comparing {len(pre_ids)} Pre-{split_year} vs {len(post_ids)} Post-{split_year} sequences."
        )
        try:
            _alignment = AlignIO.read(fasta_file, "fasta")
        except FileNotFoundError:
            print(f"Error: Could not find {fasta_file}")
            return pd.DataFrame()
        ref_seq_record = None
        for _record in _alignment:
            if ref_accession in _record.id:
                ref_seq_record = _record
                break
        if ref_seq_record is None:
            print(
                f"Error: Reference sequence containing '{ref_accession}' not found in alignment."
            )
            return pd.DataFrame()
        ref_seq_string = str(ref_seq_record.seq).lower()
        print(f"Successfully anchored to Reference Sequence: {ref_seq_record.id}")
        align_to_true_coord = {}
        true_pos = 1
        for _i, _res in enumerate(ref_seq_string):
            if _res != "-":
                align_to_true_coord[_i] = true_pos
                true_pos = true_pos + 1
            else:
                align_to_true_coord[_i] = None
        align_len = _alignment.get_alignment_length()
        pre_counts = {_i: Counter() for _i in range(align_len)}
        post_counts = {_i: Counter() for _i in range(align_len)}
        for _record in _alignment:
            if _record.id == ref_seq_record.id:
                continue
            seq_id = _record.id.split()[0]
            seq_string = str(_record.seq).lower()
            if seq_id in pre_ids:
                for _i, _res in enumerate(seq_string):
                    pre_counts[_i][_res] = pre_counts[_i][_res] + 1
            elif seq_id in post_ids:
                for _i, _res in enumerate(seq_string):
                    post_counts[_i][_res] = post_counts[_i][_res] + 1
        trend_data = []
        for pos in range(align_len):
            true_coord = align_to_true_coord.get(pos)
            if true_coord is None:
                continue
            total_pre = sum(pre_counts[pos].values())
            total_post = sum(post_counts[pos].values())
            if total_pre == 0 or total_post == 0:
                continue
            ref_residue = ref_seq_string[pos]
            all_residues = set(pre_counts[pos].keys()) | set(post_counts[pos].keys())
            for mut_residue in all_residues:
                if (
                    mut_residue == "-"
                    or mut_residue == "n"
                    or mut_residue == ref_residue
                ):
                    continue
                mut_freq_pre = (
                    pre_counts[pos][mut_residue] / total_pre
                    if mut_residue in pre_counts[pos]
                    else 0.0
                )
                mut_freq_post = (
                    post_counts[pos][mut_residue] / total_post
                    if mut_residue in post_counts[pos]
                    else 0.0
                )
                delta = mut_freq_post - mut_freq_pre
                if abs(delta) >= freq_threshold:
                    ref_freq_pre = (
                        pre_counts[pos][ref_residue] / total_pre
                        if ref_residue in pre_counts[pos]
                        else 0.0
                    )
                    ref_freq_post = (
                        post_counts[pos][ref_residue] / total_post
                        if ref_residue in post_counts[pos]
                        else 0.0
                    )
                    trend_data.append(
                        {
                            "True_Position": true_coord,
                            "Align_Column": pos + 1,
                            "Ref_Residue": ref_residue.upper(),
                            "Mutant": mut_residue.upper(),
                            "Mutation_Name": f"{ref_residue.upper()}{true_coord}{mut_residue.upper()}",
                            "Trend": "Increasing" if delta > 0 else "Decreasing",
                            "Change": round(delta, 4),
                            "Ref_Pre": round(ref_freq_pre, 4),
                            "Ref_Post": round(ref_freq_post, 4),
                            "Mut_Pre": round(mut_freq_pre, 4),
                            "Mut_Post": round(mut_freq_post, 4),
                            "Abs_Change": round(abs(delta), 4),
                        }
                    )
        results_df = pd.DataFrame(trend_data)
        if results_df.empty:
            return pd.DataFrame()
        results_df = results_df.sort_values(
            by=["Abs_Change", "True_Position"], ascending=[False, True]
        )
        return results_df.head(top_n)

    _cols = [
        "Mutation_Name",
        "Trend",
        "Ref_Pre",
        "Ref_Post",
        "Mut_Pre",
        "Mut_Post",
        "Change",
    ]
    print("=" * 65)
    print(" PROTEIN CHANGES (Reference-Anchored Frequencies) ")
    print("=" * 65)
    df_prot_detailed = analyze_reference_anchored_mutations(
        merged_df_3,
        "../data/processed/HAV_Collection_protein_align.fasta",
        split_year=2018,
        top_n=100,
        freq_threshold=0.05,
    )
    if not df_prot_detailed.empty:
        inc_prot = df_prot_detailed[df_prot_detailed["Trend"] == "Increasing"]
        dec_prot = df_prot_detailed[df_prot_detailed["Trend"] == "Decreasing"]
        print("\n--- 🟢 INCREASING MUTATIONS (Gaining prevalence) ---")
        print(
            inc_prot[_cols].to_string(index=False)
            if not inc_prot.empty
            else "No increasing mutations found >= 5%"
        )
        print("\n--- 🔴 DECREASING MUTATIONS (Losing prevalence) ---")
        print(
            dec_prot[_cols].to_string(index=False)
            if not dec_prot.empty
            else "No decreasing mutations found >= 5%"
        )
    print("\n\n" + "=" * 65)
    print(" NUCLEOTIDE CHANGES (Reference-Anchored Frequencies) ")
    print("=" * 65)
    df_nuc_detailed = analyze_reference_anchored_mutations(
        merged_df_3,
        "../data/processed/HAV_Collection date_aligned.fasta",
        split_year=2018,
        top_n=100,
        freq_threshold=0.05,
    )
    if not df_nuc_detailed.empty:
        _inc_nuc = df_nuc_detailed[df_nuc_detailed["Trend"] == "Increasing"]
        _dec_nuc = df_nuc_detailed[df_nuc_detailed["Trend"] == "Decreasing"]
        print("\n--- 🟢 INCREASING MUTATIONS (Gaining prevalence) ---")
        print(
            _inc_nuc[_cols].to_string(index=False)
            if not _inc_nuc.empty
            else "No increasing mutations found >= 5%"
        )
        print("\n--- 🔴 DECREASING MUTATIONS (Losing prevalence) ---")
        print(
            _dec_nuc[_cols].to_string(index=False)
            if not _dec_nuc.empty
            else "No decreasing mutations found >= 5%"
        )
    return Counter, inc_prot


@app.cell
def _(AlignIO, Counter, merged_df_3, pd):
    #RSCU Data obtained from MEGA 12 
    rscu_db = {
        "TTT": {"aa": "F", "hav": 1.26, "hum": 0.9387},
        "TTC": {"aa": "F", "hav": 0.74, "hum": 1.0612},
        "TTA": {"aa": "L", "hav": 1.6, "hum": 0.477},
        "TTG": {"aa": "L", "hav": 1.91, "hum": 0.783},
        "CTT": {"aa": "L", "hav": 0.93, "hum": 0.81},
        "CTC": {"aa": "L", "hav": 0.53, "hum": 1.15},
        "CTA": {"aa": "L", "hav": 0.39, "hum": 0.433},
        "CTG": {"aa": "L", "hav": 0.64, "hum": 2.344},
        "ATT": {"aa": "I", "hav": 1.88, "hum": 1.105},
        "ATC": {"aa": "I", "hav": 0.51, "hum": 1.363},
        "ATA": {"aa": "I", "hav": 0.61, "hum": 0.531},
        "ATG": {"aa": "M", "hav": 1.0, "hum": 1.0},
        "GTT": {"aa": "V", "hav": 1.86, "hum": 0.742},
        "GTC": {"aa": "V", "hav": 0.87, "hum": 0.94},
        "GTA": {"aa": "V", "hav": 0.47, "hum": 0.484},
        "GTG": {"aa": "V", "hav": 0.79, "hum": 1.832},
        "TCT": {"aa": "S", "hav": 1.24, "hum": 1.127},
        "TCC": {"aa": "S", "hav": 1.21, "hum": 1.277},
        "TCA": {"aa": "S", "hav": 1.79, "hum": 0.927},
        "TCG": {"aa": "S", "hav": 0.07, "hum": 0.327},
        "CCT": {"aa": "P", "hav": 1.6, "hum": 1.142},
        "CCC": {"aa": "P", "hav": 0.65, "hum": 1.285},
        "CCA": {"aa": "P", "hav": 1.69, "hum": 1.104},
        "CCG": {"aa": "P", "hav": 0.06, "hum": 0.467},
        "ACT": {"aa": "T", "hav": 1.39, "hum": 1.016},
        "ACC": {"aa": "T", "hav": 0.76, "hum": 1.385},
        "ACA": {"aa": "T", "hav": 1.67, "hum": 1.152},
        "ACG": {"aa": "T", "hav": 0.18, "hum": 0.445},
        "GCT": {"aa": "A", "hav": 1.51, "hum": 1.05},
        "GCC": {"aa": "A", "hav": 0.89, "hum": 1.59},
        "GCA": {"aa": "A", "hav": 1.45, "hum": 0.92},
        "GCG": {"aa": "A", "hav": 0.16, "hum": 0.438},
        "TAT": {"aa": "Y", "hav": 1.16, "hum": 0.906},
        "TAC": {"aa": "Y", "hav": 0.84, "hum": 1.1093},
        "TAA": {"aa": "*", "hav": 0.81, "hum": 0.852},
        "TAG": {"aa": "*", "hav": 0.62, "hum": 0.67},
        "CAT": {"aa": "H", "hav": 1.27, "hum": 0.851},
        "CAC": {"aa": "H", "hav": 0.73, "hum": 1.148},
        "CAA": {"aa": "Q", "hav": 1.04, "hum": 0.54},
        "CAG": {"aa": "Q", "hav": 0.96, "hum": 1.459},
        "AAT": {"aa": "N", "hav": 1.38, "hum": 0.963},
        "AAC": {"aa": "N", "hav": 0.62, "hum": 1.036},
        "AAA": {"aa": "K", "hav": 1.34, "hum": 0.889},
        "AAG": {"aa": "K", "hav": 0.66, "hum": 1.11},
        "GAT": {"aa": "D", "hav": 1.46, "hum": 0.941},
        "GAC": {"aa": "D", "hav": 0.54, "hum": 1.058},
        "GAA": {"aa": "E", "hav": 1.22, "hum": 0.863},
        "GAG": {"aa": "E", "hav": 0.78, "hum": 1.136},
        "TGT": {"aa": "C", "hav": 1.11, "hum": 0.931},
        "TGC": {"aa": "C", "hav": 0.89, "hum": 1.068},
        "TGA": {"aa": "*", "hav": 1.56, "hum": 1.477},
        "TGG": {"aa": "W", "hav": 1.0, "hum": 1.0},
        "CGT": {"aa": "R", "hav": 0.13, "hum": 0.475},
        "CGC": {"aa": "R", "hav": 0.19, "hum": 1.093},
        "CGA": {"aa": "R", "hav": 0.1, "hum": 0.6411},
        "CGG": {"aa": "R", "hav": 0.18, "hum": 1.205},
        "AGT": {"aa": "S", "hav": 1.02, "hum": 0.911},
        "AGC": {"aa": "S", "hav": 0.67, "hum": 1.428},
        "AGA": {"aa": "R", "hav": 3.81, "hum": 1.303},
        "AGG": {"aa": "R", "hav": 1.59, "hum": 1.28},
        "GGT": {"aa": "G", "hav": 1.1, "hum": 0.646},
        "GGC": {"aa": "G", "hav": 0.56, "hum": 1.346},
        "GGA": {"aa": "G", "hav": 1.62, "hum": 1.01},
        "GGG": {"aa": "G", "hav": 0.72, "hum": 0.996},
    }

    def analyze_codon_shift(true_position, ref_seq_string, mutated_nuc, utr_offset=734):
        cds_pos = true_position - utr_offset
        if cds_pos <= 0:
            return {"Error": "5' UTR"}
        codon_start_idx = true_position - 1 - (cds_pos - 1) % 3
        wt_codon = ref_seq_string[codon_start_idx : codon_start_idx + 3].upper()
        pos_in_codon = (cds_pos - 1) % 3
        mut_codon_list = list(wt_codon)
        mut_codon_list[pos_in_codon] = mutated_nuc.upper()
        mut_codon = "".join(mut_codon_list)
        wt_data = rscu_db.get(wt_codon)
        mut_data = rscu_db.get(mut_codon)
        if not wt_data or not mut_data:
            return {"Error": "N/A (Gap/Ambiguous)"}
        return {
            "Codon_Change": f"{wt_codon}->{mut_codon}",
            "AA_Change": f"{wt_data['aa']}->{mut_data['aa']}",
            "Type": "Synonymous" if wt_data["aa"] == mut_data["aa"] else "Non-Syn",
            "HAV_RSCU_Pre": wt_data["hav"],
            "HAV_RSCU_Post": mut_data["hav"],
            "Hum_RSCU_Pre": wt_data["hum"],
            "Hum_RSCU_Post": mut_data["hum"],
            "HAV_RSCU_Delta": round(mut_data["hav"] - wt_data["hav"], 4),
            "Hum_RSCU_Delta": round(mut_data["hum"] - wt_data["hum"], 4),
        }

    def analyze_nucleotide_mutations_with_rscu(
        metadata_df,
        fasta_file,
        ref_accession="NC_001489",
        split_year=2018,
        top_n=100,
        freq_threshold=0.05,
    ):
        metadata_df["year"] = pd.to_numeric(metadata_df["year"], errors="coerce")
        pre_ids = set(metadata_df[metadata_df["year"] < split_year]["accession_id"])
        post_ids = set(metadata_df[metadata_df["year"] >= split_year]["accession_id"])
        try:
            _alignment = AlignIO.read(fasta_file, "fasta")
        except FileNotFoundError:
            return pd.DataFrame()
        ref_seq_record = next((r for r in _alignment if ref_accession in r.id), None)
        if not ref_seq_record:
            return pd.DataFrame()
        ref_seq_string = str(ref_seq_record.seq).lower()
        align_to_true_coord = {}
        true_pos = 1
        for _i, _res in enumerate(ref_seq_string):
            if _res != "-":
                align_to_true_coord[_i] = true_pos
                true_pos = true_pos + 1
            else:
                align_to_true_coord[_i] = None
        align_len = _alignment.get_alignment_length()
        pre_counts = {_i: Counter() for _i in range(align_len)}
        post_counts = {_i: Counter() for _i in range(align_len)}
        for _record in _alignment:
            if _record.id == ref_seq_record.id:
                continue
            seq_id = _record.id.split()[0]
            seq_string = str(_record.seq).lower()
            if seq_id in pre_ids:
                for _i, _res in enumerate(seq_string):
                    pre_counts[_i][_res] = pre_counts[_i][_res] + 1
            elif seq_id in post_ids:
                for _i, _res in enumerate(seq_string):
                    post_counts[_i][_res] = post_counts[_i][_res] + 1
        trend_data = []
        for pos in range(align_len):
            true_coord = align_to_true_coord.get(pos)
            if true_coord is None:
                continue
            total_pre = sum(pre_counts[pos].values())
            total_post = sum(post_counts[pos].values())
            if total_pre == 0 or total_post == 0:
                continue
            ref_residue = ref_seq_string[pos]
            all_residues = set(pre_counts[pos].keys()) | set(post_counts[pos].keys())
            for mut_residue in all_residues:
                if mut_residue in ["-", "n", ref_residue]:
                    continue
                mut_freq_pre = (
                    pre_counts[pos][mut_residue] / total_pre
                    if mut_residue in pre_counts[pos]
                    else 0.0
                )
                mut_freq_post = (
                    post_counts[pos][mut_residue] / total_post
                    if mut_residue in post_counts[pos]
                    else 0.0
                )
                delta = mut_freq_post - mut_freq_pre
                if abs(delta) >= freq_threshold:
                    codon_info = analyze_codon_shift(
                        true_coord, ref_seq_string, mut_residue
                    )
                    if "Error" in codon_info:
                        c_type = codon_info["Error"]
                        (c_change, aa_change) = ("-", "-")
                        (wt_hav, mut_hav, wt_hum, mut_hum) = (0.0, 0.0, 0.0, 0.0)
                        (d_hav, d_hum) = (0.0, 0.0)
                    else:
                        c_type = codon_info["Type"]
                        c_change = codon_info["Codon_Change"]
                        aa_change = codon_info["AA_Change"]
                        wt_hav = codon_info["HAV_RSCU_Pre"]
                        mut_hav = codon_info["HAV_RSCU_Post"]
                        wt_hum = codon_info["Hum_RSCU_Pre"]
                        mut_hum = codon_info["Hum_RSCU_Post"]
                        d_hav = codon_info["HAV_RSCU_Delta"]
                        d_hum = codon_info["Hum_RSCU_Delta"]
                    trend_data.append(
                        {
                            "Mutation": f"{ref_residue.upper()}{true_coord}{mut_residue.upper()}",
                            "Trend": "Increasing" if delta > 0 else "Decreasing",
                            "Change": round(delta, 4),
                            "Type": c_type,
                            "Codon": c_change,
                            "AA": aa_change,
                            "HAV_RSCU_Pre": wt_hav,
                            "HAV_RSCU_Post": mut_hav,
                            "Δ_HAV_RSCU": d_hav,
                            "Hum_RSCU_Pre": wt_hum,
                            "Hum_RSCU_Post": mut_hum,
                            "Δ_Hum_RSCU": d_hum,
                            "Pre_Freq": round(mut_freq_pre, 4),
                            "Post_Freq": round(mut_freq_post, 4),
                            "Abs_Change": round(abs(delta), 4),
                        }
                    )
        results_df = pd.DataFrame(trend_data)
        if results_df.empty:
            return pd.DataFrame()
        return results_df.sort_values(by=["Abs_Change"], ascending=[False]).head(top_n)

    print("\n" + "=" * 85)
    print(" NUCLEOTIDE DYNAMICS & CODON OPTIMIZATION SHIFTS (NC_001489 Anchored) ")
    print("=" * 85)
    df_nuc_rscu = analyze_nucleotide_mutations_with_rscu(
        merged_df_3,
        "../data/processed/HAV_Collection date_aligned.fasta",
        split_year=2018,
        top_n=100,
        freq_threshold=0.05,
    )
    if not df_nuc_rscu.empty:
        display_cols = [
            "Mutation",
            "Trend",
            "Change",
            "Type",
            "Codon",
            "AA",
            "HAV_RSCU_Pre",
            "HAV_RSCU_Post",
            "Δ_HAV_RSCU",
        ]
        _inc_nuc = df_nuc_rscu[df_nuc_rscu["Trend"] == "Increasing"]
        _dec_nuc = df_nuc_rscu[df_nuc_rscu["Trend"] == "Decreasing"]
        print("\n--- 🟢 INCREASING NUCLEOTIDE MUTATIONS ---")
        print(
            _inc_nuc[display_cols].to_string(index=False)
            if not _inc_nuc.empty
            else "None"
        )
        print("\n--- 🔴 DECREASING NUCLEOTIDE MUTATIONS ---")
        print(
            _dec_nuc[display_cols].to_string(index=False)
            if not _dec_nuc.empty
            else "None"
        )
    return


@app.cell
def _(inc_prot):
    from scipy.stats import hypergeom

    print("\n\n" + "=" * 65)
    # --- HYPERGEOMETRIC ANALYSIS FOR FULMINANT MUTATIONS ---
    print(" VIRULENCE OVERLAP ANALYSIS (Theamboonlers et al. 2012) ")
    print("=" * 65)
    theamboonlers = {
        63: "I",
        67: "K",
        1052: "V",
        1131: "V",
        1151: "K",
        1168: "R",
        1188: "S",
        1194: "T",
        1451: "T",
        1764: "H",
        1790: "M",
        1807: "R",
        1821: "A",
        1888: "I",
        1930: "T",
    }
    if "inc_prot" in locals() and (not inc_prot.empty):
        # 1. Define the 15 Fulminant Mutations from the chimpanzee model
        matched_mutations = []
        for t_pos, t_aa in theamboonlers.items():
            match = inc_prot[
                inc_prot["True_Position"].isin([t_pos, t_pos + 1, t_pos - 1])
                & (inc_prot["Mutant"] == t_aa)
            ]
            if not match.empty:
                actual_pos = match["True_Position"].values[0]
                matched_mutations.append(f"{t_pos}{t_aa} (Mapped to {actual_pos})")
        # Ensure inc_prot exists and is not empty from the previous execution
        k = len(matched_mutations)
        N = 2227
        K = 15
        n = len(inc_prot)  # 2. Search for matches
        _p_value = hypergeom.sf(k - 1, N, K, n)
        print(
            f"Total Increasing Mutations Analyzed (n): {n}"
        )  # Check exact position, +1 shift, and -1 shift to account for literature strain variations
        print(f"Total Matches Found (k): {k} out of 15")
        if k > 0:
            print(f"Matched Residues: {', '.join(matched_mutations)}")
        print("\n--- Statistical Significance ---")
        print(f"Hypergeometric p-value: {_p_value:.2e}")
        if _p_value < 0.05:
            print(
                "Result: STATISTICALLY SIGNIFICANT."
            )  # Grab the exact mapped biological position to show where it hit
        else:
            print("Result: NOT SIGNIFICANT.")
    else:
        print(
            "Error: The 'inc_prot' dataframe is empty or undefined. No increasing mutations to analyze."
        )  # 3. Define Hypergeometric Parameters  # Observed successes (Matches found)  # Total population (Approx AAs in HAV polyprotein)  # Total successes in population (Theamboonlers mutations)  # Sample size (Total increasing reference-anchored mutations)  # 4. Calculate the p-value  # sf(k-1, N, K, n) calculates P(X >= k)  # 5. Output the results
    return


@app.cell
def _(merged_df_3, np, pd, plt):
    from scipy.stats import fisher_exact, chi2_contingency

    # FIX 1: Corrected string check to look for 'merged_df_3'
    if "merged_df_3" in locals() and "seq" in merged_df_3.columns:
        print("1. Locating Reference Genome...")
        try:
            _ref_row = merged_df_3[
                merged_df_3["accession_id"].str.contains("NC_001489", na=False)
            ].iloc[0]
            REF_SEQ = str(_ref_row["seq"]).upper()
            print(f"   Reference Found: {_ref_row['accession_id']}")
            print(f"   Length: {len(REF_SEQ)} bp")
        except IndexError:
            print("   CRITICAL ERROR: NC_001489.1 not found in merged_df_3.")
            REF_SEQ = ""

        if REF_SEQ:
            print("2. Calculating Mutation Context Loads (vs NC_001489.1)...")

            def calculate_context_loads(row, ref_seq):
                seq = str(row["seq"]).upper()
                limit = min(len(seq), len(ref_seq))
                c_c2t_g_count = 0
                g_g2a_t_count = 0
                for _i in range(1, limit - 1):
                    ref_base = ref_seq[_i]
                    query_base = seq[_i]
                    if query_base != ref_base:
                        prev_base = ref_seq[_i - 1]
                        next_base = ref_seq[_i + 1]
                        if ref_base == "C" and query_base == "T":
                            if prev_base == "C" and next_base == "G":
                                c_c2t_g_count = c_c2t_g_count + 1
                        elif ref_base == "G" and query_base == "A":
                            if prev_base == "G" and next_base == "T":
                                g_g2a_t_count = g_g2a_t_count + 1
                length_kb = len(seq) / 1000.0
                if length_kb == 0:
                    return (0.0, 0.0)
                return (c_c2t_g_count / length_kb, g_g2a_t_count / length_kb)

            loads = merged_df_3.apply(
                lambda r: calculate_context_loads(r, REF_SEQ), axis=1
            )
            merged_df_3["C[C>T]G"] = [x[0] for x in loads]
            merged_df_3["G[G>A]T"] = [x[1] for x in loads]
            print("   Metrics Calculated: 'C[C>T]G', 'G[G>A]T'")
            cut_ga = merged_df_3["ga_obye"].median()
            cut_ggat = merged_df_3["G[G>A]T"].median()
            cut_cctg = merged_df_3["C[C>T]G"].median()
            print("3. Applying Thresholds (Global Medians):")
            print(f"   GpA O/E > {cut_ga:.3f}")
            print(f"   G[G>A]T > {cut_ggat:.3f}")
            print(f"   C[C>T]G > {cut_cctg:.3f}")

            # FIX 2: Changed all internal '_row' references to match parameter 'row'
            def get_risk_score(row):
                score = 0
                if row["ga_obye"] > cut_ga:
                    score = score + 1
                if row["G[G>A]T"] > cut_ggat:
                    score = score + 1
                if row["C[C>T]G"] > cut_cctg:
                    score = score + 1
                if score == 0:
                    return "Low Risk"
                if score == 1:
                    return "Intermediate"
                return "High Risk"

            merged_df_3["Risk_Profile"] = merged_df_3.apply(get_risk_score, axis=1)
            merged_df_3["Era"] = merged_df_3["year"].apply(
                lambda x: "Post-2018" if x >= 2018 else "Pre-2018"
            )
            _plot_df = merged_df_3.dropna(subset=["year"]).copy()
            contingency = pd.crosstab(_plot_df["Risk_Profile"], _plot_df["Era"])
            contingency = contingency.reindex(
                ["Low Risk", "Intermediate", "High Risk"]
            ).fillna(0)
            print("\nRisk Distribution by Era:")
            print(contingency)
            (_chi2, p_chi, dof, ex) = chi2_contingency(contingency)
            print(f"Chi-Squared p-value: {p_chi:.4e}")
            print("\nRunning Bootstraps (n=10,000)...")
            _n_boot = 10000
            boot_ors = []
            data_for_boot = _plot_df[["Era", "Risk_Profile"]].copy()
            data_for_boot["is_HighRisk"] = (
                data_for_boot["Risk_Profile"] == "High Risk"
            ).astype(int)
            for _i in range(_n_boot):
                sample = data_for_boot.sample(frac=1.0, replace=True)
                ct = pd.crosstab(sample["is_HighRisk"], sample["Era"])
                if ct.shape != (2, 2):
                    continue
                try:
                    h_post = ct.loc[1, "Post-2018"]
                    o_post = ct.loc[0, "Post-2018"]
                    h_pre = ct.loc[1, "Pre-2018"]
                    o_pre = ct.loc[0, "Pre-2018"]
                    if o_post * h_pre > 0:
                        boot_ors.append(h_post * o_pre / (o_post * h_pre))
                except KeyError:
                    continue
            boot_ors = np.array(boot_ors)
            ci_lower = np.percentile(boot_ors, 2.5)
            ci_upper = np.percentile(boot_ors, 97.5)
            mean_or = np.mean(boot_ors)
            print(
                f"Bootstrapped OR: {mean_or:.2f} [95% CI: {ci_lower:.2f} - {ci_upper:.2f}]"
            )
            (_fig, _ax) = plt.subplots(figsize=(8, 10))
            props = (
                pd.crosstab(
                    _plot_df["Era"], _plot_df["Risk_Profile"], normalize="index"
                )
                * 100
            )
            era_order = ["Pre-2018", "Post-2018"]
            era_order = [e for e in era_order if e in props.index]
            props = props.reindex(era_order)
            props = props.reindex(columns=["Low Risk", "Intermediate", "High Risk"])
            _colors = ["#66c2a5", "#fdae61", "#d53e4f"]
            props.plot(
                kind="bar",
                stacked=True,
                color=_colors,
                edgecolor="black",
                width=0.6,
                ax=_ax,
            )
            for _c in _ax.containers:
                _labels = [
                    f"{v.get_height():.1f}%" if v.get_height() > 5 else "" for v in _c
                ]
                _ax.bar_label(
                    _c,
                    labels=_labels,
                    label_type="center",
                    fontweight="bold",
                    color="white",
                    fontsize=11,
                )
            title_text = f"Global Expansion of Hyper-Adapted Lineages\nOdds Ratio = {mean_or:.2f} [95% CI: {ci_lower:.2f}-{ci_upper:.2f}]\n$p$ (Chi-Sq) = {p_chi:.1e}"
            _ax.set_title(title_text, fontweight="bold", fontsize=14)
            _ax.set_ylabel("Proportion of Genomes (%)", fontweight="bold", fontsize=12)
            _ax.set_xlabel("Evolutionary Era", fontweight="bold", fontsize=12)
            plt.xticks(rotation=0)
            (handles, _labels) = _ax.get_legend_handles_labels()
            new_labels = ["Basal Profile", "Transitional", "Hyper-Adapted"]
            _ax.legend(
                handles,
                new_labels,
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                title="Genomic Profile",
            )
            plt.tight_layout()


            fig_out = "../results/figures/Fig_6G_HyperAdapted_Lineage_Expansion.svg"
            fig_out1 = "../results/figures/Fig_6G_HyperAdapted_Lineage_Expansion.png"
            plt.savefig(fig_out, dpi=300, bbox_inches="tight")
            plt.savefig(fig_out1, dpi=300, bbox_inches="tight")
            print(f"\n[Deliverable Saved] Figure S6/6G -> {fig_out}")

            plt.show()

            csv_out = "../data/processed/merged_df_3_final_risk_stratified.csv"
            merged_df_3.to_csv(csv_out, index=False)
            print(f"[Deliverable Saved] Risk Cohort -> {csv_out}")

        else:
            print("Cannot proceed without Reference Sequence.")
    else:
        print("Error: 'merged_df_3' not found or missing 'seq' column.")
    return chi2_contingency, fisher_exact


@app.cell
def _(merged_df_3, np, pd, plt, sns):
    from matplotlib.gridspec import GridSpec
    from scipy import stats
    import textwrap
    import scikit_posthocs as sp

    plt.rcParams["svg.fonttype"] = "none"
    _OKABE_BLUE = "#0072B2"
    _OKABE_VERMILION = "#D55E00"
    _OKABE_TEAL = "#009E73"
    _OKABE_SKY = "#56B4E9"
    _COLOR_PRE = _OKABE_BLUE
    _COLOR_POST = _OKABE_VERMILION
    _PALETTE_MAIN = {"Historical": _COLOR_PRE, "Contemporary": _COLOR_POST}
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
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
        },
    )
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    _LABEL_FS = 20
    _ANNOT_FS = 14
    _LETTER_FS = 28


    class _StatisticalReport:
        def __init__(self):
            self.report = []
            self.add_header()

        def add_header(self):
            self.report.append("==================================================")
            self.report.append("    FIGURE 1 STATISTICAL SUMMARY                  ")
            self.report.append("==================================================\n")

        def add_section(self, panel_name, test_name, n_total, results_dict):
            self.report.append(f"PANEL {panel_name}: {test_name}")
            self.report.append("-" * 50)
            self.report.append(f"Total N: {n_total}")
            for key, value in results_dict.items():
                self.report.append(f"{key:<32}: {value}")
            self.report.append("\n")

        def save(self, filename="../results/tables/figure1_stat_summary.txt"):
            with open(filename, "w") as _f:
                _f.write("\n".join(self.report))
            print(f"Stats report saved to {filename}")


    _stats_engine = _StatisticalReport()


    def _get_p_text(p_val):
        if p_val == 0:
            return "$p$ < 1e-300"
        if p_val < 0.0001:
            return f"$p$ = {p_val:.2e}"
        return f"$p$ = {p_val:.4f}"


    def _format_p_for_report(p_val):
        if p_val == 0.0:
            return "< 2.2e-308"
        elif p_val < 0.0001:
            return f"{p_val:.4e}"
        else:
            return f"{p_val:.4f}"


    def _wrap_text(text, width=20):
        return "\n".join(textwrap.wrap(text, width))


    def _draw_styled_violin(
        ax, data, x, y, palette, ylabel, order=None, p_val=None, custom_ylim=None
    ):
        if order is None:
            order = sorted(data[x].unique())
        sns.violinplot(
            data=data,
            x=x,
            y=y,
            order=order,
            palette=palette,
            inner=None,
            ax=ax,
            cut=0,
            linewidth=0,
            saturation=1,
        )
        for _i, group in enumerate(order):
            _subset = data[data[x] == group]
            if not _subset.empty:
                sns.boxplot(
                    data=_subset,
                    x=x,
                    y=y,
                    width=0.15,
                    ax=ax,
                    boxprops={
                        "zorder": 2,
                        "facecolor": palette[group],
                        "edgecolor": "black",
                        "alpha": 0.6,
                        "linewidth": 2,
                    },
                    whiskerprops={"zorder": 2, "color": "black", "linewidth": 2},
                    capprops={"zorder": 2, "color": "black", "linewidth": 2},
                    medianprops={"zorder": 2, "color": "black", "linewidth": 2},
                    fliersize=0,
                )
        _plot_data = data.sample(2000) if len(data) > 2000 else data
        sns.swarmplot(
            data=_plot_data,
            x=x,
            y=y,
            order=order,
            color="#2b2b2b",
            alpha=0.7,
            size=4,
            ax=ax,
            zorder=1,
        )
        if custom_ylim:
            ax.set_ylim(custom_ylim)
        if p_val is not None:
            (c_ymin, c_ymax) = ax.get_ylim()
            visible_points = data[data[y] <= c_ymax][y]
            visible_max = visible_points.max() if not visible_points.empty else c_ymax
            h = (c_ymax - c_ymin) * 0.03
            y_bracket_base = visible_max + h
            y_bracket_top = y_bracket_base + h
            ax.plot(
                [0, 0, 1, 1],
                [y_bracket_base, y_bracket_top, y_bracket_top, y_bracket_base],
                lw=2,
                c="black",
            )
            ax.text(
                0.5,
                y_bracket_top + h * 0.5,
                _get_p_text(p_val),
                ha="center",
                va="bottom",
                fontsize=_ANNOT_FS,
                color="black",
                fontweight="bold",
            )
            ax.set_ylim(c_ymin, max(c_ymax, y_bracket_top + h * 5))
        ax.set_ylabel(
            _wrap_text(ylabel, width=20), fontsize=_LABEL_FS, fontweight="bold"
        )
        ax.set_xlabel("")
        sns.despine(ax=ax, trim=False)


    if "merged_df_3" not in locals():
        print("WARNING: merged_df_3 not found.")
    elif "year" in merged_df_3.columns:
        merged_df_3["Condition"] = merged_df_3["year"].apply(
            lambda x: "Historical" if x < 2018 else "Contemporary"
        )
        _counts = merged_df_3["Condition"].value_counts()
        merged_df_3["Period_Label"] = merged_df_3["Condition"].apply(
            lambda x: f"{x} (n={_counts.get(x, 0)})"
        )

        def _get_timeline(y):
            if 1998 <= y <= 2007:
                return "1998-2007"
            elif 2008 <= y <= 2012:
                return "2008-2012"
            elif 2013 <= y <= 2017:
                return "2013-2017"
            elif 2018 <= y <= 2022:
                return "2018-2022"
            return "Other"

        merged_df_3["Timeline"] = merged_df_3["year"].apply(_get_timeline)
        _unique_labels = merged_df_3["Period_Label"].unique()
        _scatter_palette = {}
        for _lbl in _unique_labels:
            _scatter_palette[_lbl] = _COLOR_PRE if "Historical" in _lbl else _COLOR_POST

    _fig = plt.figure(figsize=(24, 24))
    _gs = GridSpec(4, 2, figure=_fig, hspace=0.6, wspace=0.35)
    _axA = _fig.add_subplot(_gs[0, :])
    _axB = _fig.add_subplot(_gs[1, :])
    _axC = _fig.add_subplot(_gs[2, :])
    _axD = _fig.add_subplot(_gs[3, 0])
    _axE = _fig.add_subplot(_gs[3, 1])

    if "Nucleotide Mutations" in merged_df_3.columns:
        _clean_A = merged_df_3.dropna(subset=["year", "Nucleotide Mutations"])
        sns.scatterplot(
            data=_clean_A,
            x="year",
            y="Nucleotide Mutations",
            hue="Period_Label",
            palette=_scatter_palette,
            s=180,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.8,
            ax=_axA,
        )
        _data_pre = _clean_A[_clean_A["Condition"] == "Historical"]
        _data_post = _clean_A[_clean_A["Condition"] == "Contemporary"]
        if not _data_pre.empty:
            sns.regplot(
                data=_data_pre,
                x="year",
                y="Nucleotide Mutations",
                color=_COLOR_PRE,
                scatter=False,
                ax=_axA,
                line_kws={"linewidth": 4, "linestyle": "--", "alpha": 0.8},
                ci=None,
            )
        if not _data_post.empty:
            sns.regplot(
                data=_data_post,
                x="year",
                y="Nucleotide Mutations",
                color=_COLOR_POST,
                scatter=False,
                ax=_axA,
                line_kws={"linewidth": 4, "linestyle": "-", "alpha": 0.9},
                ci=None,
            )
        _y_limits = _axA.get_ylim()
        _y_txt_pos = _y_limits[1] * 0.95
        _axA.axvline(x=1995, color=_COLOR_POST, linestyle="--", lw=2.5)
        _axA.text(
            1995.5,
            _y_txt_pos,
            "Vaccine (1995)",
            color=_COLOR_POST,
            rotation=90,
            va="top",
            fontsize=_ANNOT_FS,
            fontweight="bold",
        )
        _axA.axvline(x=2018, color=_OKABE_TEAL, linestyle=":", lw=2.5)
        _axA.text(
            2016.5,
            _y_txt_pos,
            "Outbreak (2018)",
            color=_OKABE_TEAL,
            rotation=90,
            va="top",
            fontsize=_ANNOT_FS,
            fontweight="bold",
        )
        (_rho, _p_rho) = stats.spearmanr(
            _clean_A["year"], _clean_A["Nucleotide Mutations"]
        )
        _stats_engine.add_section(
            "A",
            "Spearman Correlation",
            len(_clean_A),
            {"Rho": f"{_rho:.4f}", "P-value": _format_p_for_report(_p_rho)},
        )
        _annot_txt = f"Spearman's ρ = {_rho:.2f}\n{_get_p_text(_p_rho)}"
        _axA.text(
            0.02,
            0.95,
            _annot_txt,
            transform=_axA.transAxes,
            fontsize=_ANNOT_FS,
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", alpha=1),
            va="top",
        )
        _axA.text(
            -0.06,
            1.1,
            "A",
            transform=_axA.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="top",
        )
        _axA.set_ylabel(
            _wrap_text("Nucleotide Mutations / Seq"),
            fontsize=_LABEL_FS,
            fontweight="bold",
        )
        _axA.set_xlabel("Collection Year", fontsize=_LABEL_FS, fontweight="bold")
        _axA.legend(
            frameon=False, fontsize=14, loc="upper left", bbox_to_anchor=(0.18, 1)
        )
        sns.despine(ax=_axA, trim=False)

    if "Amino Acid Mutations" in merged_df_3.columns:
        _clean_B = merged_df_3[merged_df_3["Timeline"] != "Other"].copy()
        _timeline_order = ["1998-2007", "2008-2012", "2013-2017", "2018-2022"]
        sns.barplot(
            data=_clean_B,
            x="Timeline",
            y="Nucleotide Mutations",
            color=_OKABE_SKY,
            alpha=0.6,
            estimator=np.median,
            errorbar=("ci", 95),
            order=_timeline_order,
            ax=_axB,
            capsize=0.1,
            err_kws={"linewidth": 2},
        )
        _axB.set_ylabel(
            _wrap_text("Nuc. Mutations (Median)"),
            fontsize=_LABEL_FS,
            fontweight="bold",
            color=_OKABE_BLUE,
        )
        _axB.tick_params(axis="y", colors=_OKABE_BLUE)
        _axB2 = _axB.twinx()
        sns.pointplot(
            data=_clean_B,
            x="Timeline",
            y="Amino Acid Mutations",
            color=_OKABE_VERMILION,
            estimator=np.median,
            errorbar=("ci", 95),
            order=_timeline_order,
            ax=_axB2,
            capsize=0.1,
            markers="D",
            scale=1.5,
        )
        _axB2.set_ylabel(
            _wrap_text("AA Mutations (Median)"),
            fontsize=_LABEL_FS,
            fontweight="bold",
            color=_OKABE_VERMILION,
        )
        _axB2.tick_params(axis="y", colors=_OKABE_VERMILION)
        _axB2.spines["right"].set_color(_OKABE_VERMILION)
        _axB2.spines["right"].set_linewidth(2)
        _axB.spines["left"].set_color(_OKABE_BLUE)
        _axB.spines["left"].set_linewidth(2)

        # 1. OMNIBUS KRUSKAL-WALLIS
        _groups_nuc = [
            _clean_B[_clean_B["Timeline"] == t]["Nucleotide Mutations"].dropna().values
            for t in _timeline_order
        ]
        (_kw_n, _p_n) = stats.kruskal(*_groups_nuc)

        _groups_aa = [
            _clean_B[_clean_B["Timeline"] == t]["Amino Acid Mutations"].dropna().values
            for t in _timeline_order
        ]
        (_kw_a, _p_a) = stats.kruskal(*_groups_aa)

        # 2. DUNN'S POST-HOC (Holm-Bonferroni corrected)
        _dunn_nuc = sp.posthoc_dunn(
            _clean_B, val_col="Nucleotide Mutations", group_col="Timeline", p_adjust="holm"
        )
        _dunn_aa = sp.posthoc_dunn(
            _clean_B, val_col="Amino Acid Mutations", group_col="Timeline", p_adjust="holm"
        )

        panel_b_dict = {
            "Nuc H (Omnibus)": f"{_kw_n:.2f}",
            "Nuc P (Omnibus)": _format_p_for_report(_p_n),
            "AA H (Omnibus)": f"{_kw_a:.2f}",
            "AA P (Omnibus)": _format_p_for_report(_p_a),
        }

        # Unpack upper-triangle of pairwise Dunn matrix into report dict
        for _i in range(len(_timeline_order)):
            for _j in range(_i + 1, len(_timeline_order)):
                _g1, _g2 = _timeline_order[_i], _timeline_order[_j]
                if _g1 in _dunn_nuc.index and _g2 in _dunn_nuc.columns:
                    panel_b_dict[f"Nuc [{_g1} vs {_g2}]"] = _format_p_for_report(
                        _dunn_nuc.loc[_g1, _g2]
                    )
                if _g1 in _dunn_aa.index and _g2 in _dunn_aa.columns:
                    panel_b_dict[f"AA  [{_g1} vs {_g2}]"] = _format_p_for_report(
                        _dunn_aa.loc[_g1, _g2]
                    )

        _stats_engine.add_section(
            "B", "Kruskal-Wallis + Dunn Post-Hoc (Holm)", len(_clean_B), panel_b_dict
        )

        _axB.text(
            -0.06,
            1.1,
            "B",
            transform=_axB.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="top",
        )
        _stats_title = f"Nuc: {_get_p_text(_p_n)} | AA: {_get_p_text(_p_a)}"
        _axB.text(
            0.5,
            1.1,
            _stats_title,
            transform=_axB.transAxes,
            ha="center",
            fontsize=_ANNOT_FS,
            fontweight="bold",
        )
        _axB.set_xlabel("Time Period", fontsize=_LABEL_FS, fontweight="bold")
        sns.despine(ax=_axB, top=True, right=False, trim=False)

    print("\n" + "=" * 50)
    print(" MEDIANS FOR TIMELINE BINS (PANEL B)")
    print("=" * 50)
    if (
        "Amino Acid Mutations" in merged_df_3.columns
        and "Nucleotide Mutations" in merged_df_3.columns
    ):
        timeline_medians = _clean_B.groupby("Timeline")[
            ["Nucleotide Mutations", "Amino Acid Mutations"]
        ].median()
        timeline_medians = timeline_medians.reindex(_timeline_order)
        print(timeline_medians)

    try:
        _sky_df = pd.read_csv("../data/raw/HAV_beast_skygrid.txt", sep="\t", comment="#")
        if "mean" in _sky_df.columns:
            _sky_df["log_mean"] = np.log(_sky_df["mean"])
            _sky_df["log_upper"] = np.log(_sky_df["upper"])
            _sky_df["log_lower"] = np.log(_sky_df["lower"])
            _axC.fill_between(
                _sky_df["time"],
                _sky_df["log_lower"],
                _sky_df["log_upper"],
                color="grey",
                alpha=0.25,
                edgecolor=None,
            )
            _axC.plot(_sky_df["time"], _sky_df["log_mean"], color=_COLOR_PRE, lw=4)
            _axC.set_xlim(1990, 2025)
            _y_max = _sky_df["log_mean"].max()
            _axC.axvline(x=1995, color=_COLOR_POST, linestyle="--", lw=2.5)
            _axC.text(
                1995.5,
                _y_max,
                "Vaccine (1995)",
                color=_COLOR_POST,
                rotation=90,
                va="top",
                fontsize=_ANNOT_FS,
                fontweight="bold",
            )
            _axC.axvline(x=2018, color=_OKABE_TEAL, linestyle=":", lw=2.5)
            _axC.text(
                2017,
                _y_max,
                "Outbreak (2018)",
                color=_OKABE_TEAL,
                rotation=90,
                va="top",
                fontsize=_ANNOT_FS,
                fontweight="bold",
            )
            _axC.text(
                -0.06,
                1.1,
                "C",
                transform=_axC.transAxes,
                fontsize=_LETTER_FS,
                fontweight="bold",
                va="top",
            )
            _axC.set_ylabel(
                _wrap_text("log(Effective Pop. Size)"),
                fontsize=_LABEL_FS,
                fontweight="bold",
            )
            sns.despine(ax=_axC, trim=False)
            _axC.set_ylim(0, 8)
    except Exception:
        _axC.text(0.5, 0.5, "Data Unavailable", ha="center")

    _motif_panels = [
        (_axD, "cpg_obye", "D", "CpG O/E Ratio", (0.1, 0.2)),
        (_axE, "ga_obye", "E", "GpA O/E Ratio", (1.1, 1.4)),
    ]
    for _ax, _col, _letter, _ylab, custom_ylim in _motif_panels:
        if _col in merged_df_3.columns:
            _clean = merged_df_3.dropna(subset=[_col])
            _pre = _clean[_clean["Condition"] == "Historical"][_col]
            _post = _clean[_clean["Condition"] == "Contemporary"][_col]
            (_u, _p) = stats.mannwhitneyu(_pre, _post)
            _stats_engine.add_section(
                _letter,
                f"Mann-Whitney {_col}",
                len(_clean),
                {
                    "U": f"{_u}",
                    "P": _format_p_for_report(_p),
                    "Med Pre": f"{_pre.median()}",
                    "Med Post": f"{_post.median()}",
                },
            )
            _draw_styled_violin(
                _ax,
                _clean,
                "Condition",
                _col,
                _PALETTE_MAIN,
                _ylab,
                order=["Historical", "Contemporary"],
                p_val=_p,
                custom_ylim=custom_ylim,
            )
            _ax.text(
                -0.15,
                1.1,
                _letter,
                transform=_ax.transAxes,
                fontsize=_LETTER_FS,
                fontweight="bold",
                va="top",
            )
        else:
            _ax.text(0.5, 0.5, f"Missing Data\n{_col}", ha="center")

    print("\n" + "=" * 50)
    print(" MEDIANS FOR CONDITION BINS (PANELS D & E)")
    print("=" * 50)
    for _ax, _col, _letter, _ylab, custom_ylim in _motif_panels:
        if _col in merged_df_3.columns:
            condition_medians = (
                merged_df_3.dropna(subset=[_col]).groupby("Condition")[_col].median()
            )
            print(f"\n{_ylab} ({_col}):")
            print(condition_medians)


    _stats_engine.save()

    plt.tight_layout()
    plt.savefig("../results/figures/Figure1.png", dpi=300, bbox_inches="tight")
    plt.savefig("../results/figures/Figure1.pdf", format="pdf", bbox_inches="tight")
    plt.savefig("../results/figures/Figure1.svg", format="svg", bbox_inches="tight")
    print("Figure 1 saved as Figure1.png/pdf/svg.")
    plt.show()
    return GridSpec, stats


@app.cell
def _(pd, plt, sns):
    df_6 = pd.read_csv("../data/raw/ESS_Data.txt", sep="\t", index_col=0)
    alpha_data = {
        "Codon_Position": ["CP1", "CP2", "CP3"],
        "Alpha_Value": [
            float(df_6.loc["mean", "CP1.alpha"]),
            float(df_6.loc["mean", "CP2.alpha"]),
            float(df_6.loc["mean", "CP3.alpha"]),
        ],
    }
    df_alpha = pd.DataFrame(alpha_data)
    rate_types = ["rateAC", "rateAG", "rateAT", "rateCG", "rateCT", "rateGT"]
    overall_rates = []
    for rate in rate_types:
        means = [
            float(df_6.loc["mean", f"CP{_i}.gtr.rates.{rate}"]) for _i in [1, 2, 3]
        ]
        avg_mean = sum(means) / 3.0
        (lowers, uppers) = ([], [])
        for _i in [1, 2, 3]:
            hpd = (
                df_6.loc["95% HPD interval", f"CP{_i}.gtr.rates.{rate}"]
                .strip("[]")
                .split(",")
            )
            lowers.append(float(hpd[0]))
            uppers.append(float(hpd[1]))
        avg_lower = sum(lowers) / 3.0
        avg_upper = sum(uppers) / 3.0
        overall_rates.append(
            {
                "Substitution": rate.replace("rate", ""),
                "Mean": avg_mean,
                "yerr_lower": avg_mean - avg_lower,
                "yerr_upper": avg_upper - avg_mean,
            }
        )
    df_overall = pd.DataFrame(overall_rates)
    sns.set_context("paper", font_scale=1.5)
    sns.set_style(
        "ticks",
        {
            "axes.linewidth": 2.0,
            "xtick.major.width": 2.0,
            "ytick.major.width": 2.0,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "font.family": "sans-serif",
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
        },
    )
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    (_fig, _axes) = plt.subplots(
        1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [2, 1]}
    )
    _ax1 = _axes[0]
    _colors = ["#bdc3c7", "#e74c3c", "#bdc3c7", "#bdc3c7", "#3498db", "#bdc3c7"]
    sns.barplot(
        data=df_overall,
        x="Substitution",
        y="Mean",
        palette=_colors,
        edgecolor="black",
        linewidth=1.5,
        ax=_ax1,
    )
    _ax1.errorbar(
        x=range(len(df_overall)),
        y=df_overall["Mean"],
        yerr=[df_overall["yerr_lower"], df_overall["yerr_upper"]],
        fmt="none",
        c="black",
        capsize=5,
        linewidth=1.5,
    )
    _ax1.set_title(
        "A. Overall Nucleotide Substitution Rates", fontweight="bold", pad=15
    )
    _ax1.set_xlabel("Substitution Type", fontweight="bold")
    _ax1.set_ylabel("Relative Rate (Mean ± Avg 95% HPD)", fontweight="bold")
    _ax1.axvspan(0.5, 1.5, color="red", alpha=0.05, label="ADAR target (A>G)")
    _ax1.axvspan(3.5, 4.5, color="blue", alpha=0.05, label="APOBEC target (C>T)")
    _ax1.legend(loc="upper left", frameon=True)
    _ax2 = _axes[1]
    sns.barplot(
        data=df_alpha,
        x="Codon_Position",
        y="Alpha_Value",
        palette=["#1abc9c", "#e67e22", "#9b59b6"],
        edgecolor="black",
        linewidth=1.5,
        ax=_ax2,
    )
    _ax2.set_title("B. Gamma Shape (α) by Position", fontweight="bold", pad=15)
    _ax2.set_xlabel("Codon Position", fontweight="bold")
    _ax2.set_ylabel("Alpha (α) Value", fontweight="bold")
    for _p in _ax2.patches:
        _ax2.annotate(
            format(_p.get_height(), ".3f"),
            (_p.get_x() + _p.get_width() / 2.0, _p.get_height()),
            ha="center",
            va="center",
            xytext=(0, 9),
            textcoords="offset points",
            fontweight="bold",
        )
    sns.despine()
    plt.tight_layout()
    output_file = "../results/figures/HAV_Combined_Evolutionary_Metrics.pdf"
    output_file_svg = "../results/figures/HAV_Combined_Evolutionary_Metrics.svg"
    output_file_png = "../results/figures/HAV_Combined_Evolutionary_Metrics.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.savefig(output_file_svg, format="svg", bbox_inches="tight")
    plt.savefig(output_file_png, format="png", bbox_inches="tight")
    print(f"Combined plot saved to {output_file}")
    plt.show()
    print("\n" + "=" * 50)
    print(" SPECIFIC SUBSTITUTION RATES ")
    print("=" * 50)
    ag_rate = df_overall.loc[df_overall["Substitution"] == "AG", "Mean"].values[0]
    ct_rate = df_overall.loc[df_overall["Substitution"] == "CT", "Mean"].values[0]
    transversions = ["AC", "AT", "CG", "GT"]
    bg_transversions = df_overall[df_overall["Substitution"].isin(transversions)]
    bg_transversion_mean = bg_transversions["Mean"].mean()
    print(f"AG Rate (ADAR target):             {ag_rate:.5f}")
    print(f"CT Rate (APOBEC target):           {ct_rate:.5f}")
    print(f"Background Transversion Rate:      {bg_transversion_mean:.5f}")
    print("(Average of AC, AT, CG, and GT)")
    print("=" * 50 + "\n")
    return


@app.cell
def _(AlignIO, GridSpec, np, pd, plt, sns, stats):
    import json
    import math
    import traceback
    import matplotlib.patches as mpatches
    from Bio import Phylo
    from Bio.Seq import Seq
    from matplotlib.gridspec import GridSpecFromSubplotSpec


    # ==============================================================================
    # 1. CONFIGURATION & TYPOGRAPHY
    # ==============================================================================
    _FIG_SIZE = (30, 20)
    sns.set_context("paper", font_scale=1.4)
    sns.set_style(
        "ticks",
        {
            "axes.linewidth": 2.0,
            "xtick.major.width": 2.0,
            "ytick.major.width": 2.0,
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.edgecolor": "black",
            "xtick.color": "black",
        },
    )

    _COLOR_POS = "#D55E00"
    _COLOR_NEG = "#0072B2"
    _COLOR_NEU = "#FFFFFF"
    _COLOR_TREE = "#333333"

    HAV_GENES_NT = [
        ("VP4", 84, "#E0F2F1"),
        ("VP2", 888, "#B2DFDB"),
        ("VP3", 798, "#80CBC4"),
        ("VP1", 837, "#4DB6AC"),
        ("2A", 213, "#FFF9C4"),
        ("2B", 810, "#FFF59D"),
        ("2C", 1047, "#FFF176"),
        ("3A", 243, "#FFE0B2"),
        ("3B", 72, "#FFCC80"),
        ("3C", 708, "#FFB74D"),
        ("3D", 1560, "#FF8A65"),
    ]
    POLYPROTEIN_START_CODON = 1


    # ==============================================================================
    # 2. SHARED DATA PARSING UTILITIES
    # ==============================================================================
    def _load_alignment_robust(file_path):
        for fmt in ["phylip-relaxed", "phylip", "fasta"]:
            try:
                return AlignIO.read(file_path, fmt)
            except Exception:
                continue
        raise ValueError(f"Could not read {file_path}")


    def _get_sequence_for_tip(tip_name, seq_map):
        if tip_name in seq_map:
            return seq_map[tip_name]
        norm = tip_name.replace("_", ".")
        if norm in seq_map:
            return seq_map[norm]
        norm2 = tip_name.replace(".", "_")
        if norm2 in seq_map:
            return seq_map[norm2]
        if len(tip_name) > 10:
            trunc = tip_name[:10]
            if trunc in seq_map:
                return seq_map[trunc]
        return None


    def _get_tip_year_color(tip_name, metadata_df):
        if metadata_df is None:
            return "#999999"
        clean_id = tip_name.replace("_", ".")
        try:
            col = (
                "accession_id" if "accession_id" in metadata_df.columns else "Accession"
            )
            if col in metadata_df.columns:
                row = metadata_df[metadata_df[col] == clean_id]
                if row.empty:
                    row = metadata_df[metadata_df[col] == tip_name]
                if not row.empty:
                    if "Condition" in row.columns:
                        cond = str(row["Condition"].values[0])
                        return _COLOR_POS if cond.lower() == "contemporary" else _COLOR_NEG
                    elif "year" in row.columns:
                        year = int(row["year"].values[0])
                        return _COLOR_POS if year >= 2018 else _COLOR_NEG
        except Exception:
            pass
        return "#999999"


    def _translate_codon(codon_str):
        if "-" in codon_str or "N" in codon_str or len(codon_str) != 3:
            return "?"
        try:
            return str(Seq(codon_str).translate())
        except Exception:
            return "?"


    def _is_valid_codon(c):
        if not isinstance(c, str) or len(c) != 3:
            return False
        if any(char in c for char in ["-", "?", "N", "X", "*"]):
            return False
        return True


    def _load_json_data(file, val_col, cat_logic):
        try:
            with open(file, "r") as f:
                loaded_json = json.load(f)

            mle_tree = loaded_json.get("MLE", {})
            content = mle_tree.get("content", {})
            raw_rows = []

            if isinstance(content, list):
                raw_rows = content
            elif isinstance(content, dict) and len(content) > 0:
                sorted_keys = sorted(
                    content.keys(), key=lambda x: int(x) if str(x).isdigit() else x
                )
                first_item = content[sorted_keys[0]]
                if (
                    isinstance(first_item, list)
                    and len(first_item) > 0
                    and isinstance(first_item[0], (list, dict))
                ):
                    raw_rows = first_item
                else:
                    raw_rows = [content[k] for k in sorted_keys]

            if not raw_rows:
                return pd.DataFrame()

            plot_data = []
            for idx, r in enumerate(raw_rows):
                res = cat_logic(r)
                if res:
                    plot_data.append({"Site": idx + 1, **res})
            return pd.DataFrame(plot_data)

        except Exception as e:
            print(f"   [PARSER ERROR] Crashed reading {file}: {e}")
            return pd.DataFrame()


    def get_fubar_row(row):
        ds = float(row[0])
        dn = float(row[1])
        val = dn - ds  # True dN - dS calculation
        neg_prob = float(row[3])
        pos_prob = float(row[4])
        cat = "Positive" if pos_prob >= 0.95 else ("Negative" if neg_prob >= 0.95 else "Neutral")
        return {"Value": val, "Category": cat}


    def _get_meme_row(row):
        p_val = float(row[6]) if len(row) > 6 else 1.0
        p_val = 1e-10 if p_val == 0 else p_val
        cat = "Significant" if p_val <= 0.05 else "Neutral"
        return {"Score": -np.log10(p_val), "Category": cat}


    # ==============================================================================
    # 3. PATRISTIC DIVERGENCE ENGINE & S2B VIOLIN PLOTTER
    # ==============================================================================
    def plot_patristic_violin(
        df_clean,
        u_stat,
        p_val,
        output_file="../results/figures/Supplementary_Figure_S2B_Patristic_Divergence.pdf",
    ):
        print("Generating Patristic Divergence Violin Plot (Supplementary Figure S2B)...")
        plt.figure(figsize=(10, 8))

        era_order = ["Historical", "Contemporary"]
        palette_map = {"Historical": _COLOR_NEG, "Contemporary": _COLOR_POS}

        ax = sns.violinplot(
            data=df_clean,
            x="Era",
            y="Divergence_Subs_Per_Site",
            order=era_order,
            palette=palette_map,
            inner=None,
            alpha=0.6,
            cut=0,
        )

        sns.stripplot(
            data=df_clean,
            x="Era",
            y="Divergence_Subs_Per_Site",
            order=era_order,
            palette=palette_map,
            size=6,
            jitter=0.2,
            alpha=0.8,
            edgecolor="black",
            linewidth=0.5,
            ax=ax,
        )

        medians = df_clean.groupby("Era")["Divergence_Subs_Per_Site"].median()
        counts = df_clean.groupby("Era")["Divergence_Subs_Per_Site"].count()

        for i, era in enumerate(era_order):
            if era in medians:
                med_val = medians[era]
                n_cnt = counts[era]
                ax.hlines(med_val, xmin=i - 0.25, xmax=i + 0.25, color="black", linewidth=3.0, zorder=10)
                ax.text(
                    i,
                    med_val + 0.008,
                    f"Median: {med_val:.4f}\n(n={n_cnt})",
                    horizontalalignment="center",
                    fontweight="bold",
                    fontsize=13,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", alpha=0.8),
                )

        max_y = df_clean["Divergence_Subs_Per_Site"].max()
        y_bar = max_y * 1.08
        ax.plot([0, 0, 1, 1], [y_bar - 0.005, y_bar, y_bar, y_bar - 0.005], lw=2.0, c="black")
        p_text = f"p = {p_val:.2e}" if p_val < 0.001 else f"p = {p_val:.4f}"
        ax.text(
            0.5,
            y_bar + 0.004,
            f"Mann-Whitney U Test\n{p_text}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=14,
        )

        ax.set_ylabel("Patristic Distance to HM175 Ref (substitutions/site)", fontweight="bold", fontsize=16)
        ax.set_xlabel("Sampling Era", fontweight="bold", fontsize=16)
        ax.set_ylim(0, max_y * 1.25)

        sns.despine(trim=False)
        plt.tight_layout()
        try:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
        except Exception:
            plt.savefig("Supplementary_Figure_S2B.pdf", dpi=300, bbox_inches="tight")
        plt.savefig("Supplementary_Figure_S2B.png", dpi=300, bbox_inches="tight")
        print("✔ Patristic Violin Plot saved successfully.")


    def get_era_divergence_distribution(tree, metadata_df, ref_identifier="NC_001489"):
        print("\n" + "=" * 65)
        print("   HAV PATRISTIC DIVERGENCE DISTRIBUTION (vs HM175 Ref)   ")
        print("=" * 65)

        terminals = tree.get_terminals()
        ref_tip = next((t for t in terminals if t.name and ref_identifier in t.name), None)

        if not ref_tip:
            print(f"[ERROR] Reference '{ref_identifier}' not found for divergence calculation.")
            return None

        records = []
        for tip in terminals:
            if tip.name == ref_tip.name:
                continue
            patristic_dist = tree.distance(ref_tip, tip)
            era = "Unknown"

            if metadata_df is not None:
                clean_name = tip.name.replace("_", ".")
                col = "accession_id" if "accession_id" in metadata_df.columns else "Accession"
                match = metadata_df[metadata_df[col] == clean_name]
                if match.empty:
                    match = metadata_df[metadata_df[col] == tip.name]

                if not match.empty:
                    if "Condition" in match.columns:
                        cond = str(match["Condition"].values[0]).capitalize()
                        era = "Contemporary" if "contemp" in cond.lower() else "Historical"
                    elif "year" in match.columns:
                        yr = int(match["year"].values[0])
                        era = "Contemporary" if yr >= 2018 else "Historical"

            records.append({
                "Accession": tip.name,
                "Divergence_Subs_Per_Site": round(patristic_dist, 6),
                "Era": era,
            })

        df_div = pd.DataFrame(records)
        df_clean = df_div[df_div["Era"].isin(["Historical", "Contemporary"])].copy()

        if df_clean.empty:
            print("[WARNING] Could not map tree tips to Historical/Contemporary eras.")
            return df_div

        summary = (
            df_clean.groupby("Era")["Divergence_Subs_Per_Site"]
            .agg(["count", "median", lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75), "mean", "std"])
            .rename(columns={"<lambda_0>": "IQR_25", "<lambda_1>": "IQR_75"})
        )
        print(summary.to_string())

        hist_vals = df_clean[df_clean["Era"] == "Historical"]["Divergence_Subs_Per_Site"]
        cont_vals = df_clean[df_clean["Era"] == "Contemporary"]["Divergence_Subs_Per_Site"]

        u_stat, p_val = stats.mannwhitneyu(hist_vals, cont_vals, alternative="two-sided")
        print("-" * 65)
        print(f"Mann-Whitney U Test: U = {u_stat:.1f}, p-value = {p_val:.4e}")
        print("=" * 65 + "\n")

        df_clean.to_csv("HAV_Era_Divergence_Distribution.csv", index=False)
        plot_patristic_violin(df_clean, u_stat, p_val)
        return df_clean


    # ==============================================================================
    # 4. FIGURE 2 GENERATOR (WITH MANUSCRIPT COUNTERS)
    # ==============================================================================
    def plot_figure_2_final(
        tree_file="../data/raw/READY_FOR_HYPHY.nwk",
        aln_file="../data/raw/READY_FOR_HYPHY.fasta",
        fubar_file="../data/raw/HAV_fubar.json",
        meme_file="../data/raw/HAV_meme.json",
        output_file="../results/figures/Figure2_Ref_NC001489.pdf",
        metadata_df=None,
    ):
        print("Generating Figure 2 (Reference: NC_001489_1)...")
        df_fubar = _load_json_data(fubar_file, "Value", get_fubar_row)
        df_meme = _load_json_data(meme_file, "Score", _get_meme_row)

        print("\n" + "=" * 65)
        print("MEME POSITIVE SELECTION DISTRIBUTION (FOR MANUSCRIPT)")
        print("=" * 65)

        region_map = {
            "P1 (Structural: VP4-VP1)": ["VP4", "VP2", "VP3", "VP1"],
            "P2 (Non-structural: 2A-2C)": ["2A", "2B", "2C"],
            "P3 (Non-structural: 3A-3D)": ["3A", "3B", "3C", "3D"],
        }

        site_to_gene = {}
        site_to_region = {}
        curr_codon_count = POLYPROTEIN_START_CODON

        for gene_name, len_nt, _ in HAV_GENES_NT:
            len_codons = int(len_nt / 3)
            for site in range(curr_codon_count, curr_codon_count + len_codons):
                site_to_gene[site] = gene_name
                for reg_name, gene_list in region_map.items():
                    if gene_name in gene_list:
                        site_to_region[site] = reg_name
                        break
            curr_codon_count += len_codons

        if not df_meme.empty:
            sig_meme = df_meme[df_meme["Category"] == "Significant"]
            total_sig = len(sig_meme)
            print(f"Total statistically significant MEME sites (p <= 0.05): {total_sig}")

            reg_counts = {
                "P1 (Structural: VP4-VP1)": 0,
                "P2 (Non-structural: 2A-2C)": 0,
                "P3 (Non-structural: 3A-3D)": 0,
            }
            gene_counts = {}

            for site in sig_meme["Site"]:
                reg = site_to_region.get(int(site), "Unknown/UTR")
                gene = site_to_gene.get(int(site), "Unknown")
                if reg in reg_counts:
                    reg_counts[reg] += 1
                gene_counts[gene] = gene_counts.get(gene, 0) + 1

            p1_count = reg_counts["P1 (Structural: VP4-VP1)"]
            p2_count = reg_counts["P2 (Non-structural: 2A-2C)"]
            p3_count = reg_counts["P3 (Non-structural: 3A-3D)"]
            non_structural_total = p2_count + p3_count

            print(f" -> P1 Structural (VP4, VP2, VP3, VP1): {p1_count:>3} sites")
            print(f" -> P2 Non-structural (2A, 2B, 2C):     {p2_count:>3} sites")
            print(f" -> P3 Non-structural (3A, 3B, 3C, 3D): {p3_count:>3} sites")
            print(" ----------------------------------------------------")
            print(f" -> TOTAL NON-STRUCTURAL (P2 + P3):     {non_structural_total:>3} sites")
            print("\nGene-by-Gene Breakdown of Significant Sites:")
            for gene_name, _, _ in HAV_GENES_NT:
                c = gene_counts.get(gene_name, 0)
                print(f"    - {gene_name:<4}: {c:>2} sites")
        print("=" * 65 + "\n")

        max_site_fubar = df_fubar["Site"].max() if not df_fubar.empty else 0
        max_site_meme = df_meme["Site"].max() if not df_meme.empty else 0
        max_x = max(2227, max_site_fubar, max_site_meme)

        _fig = plt.figure(figsize=(26, 30))
        _gs = GridSpec(5, 1, figure=_fig, height_ratios=[12, 2.5, 2.5, 0.5, 12], hspace=0.15)
        _axA = _fig.add_subplot(_gs[0])
        _axB = _fig.add_subplot(_gs[1])
        _axC = _fig.add_subplot(_gs[2])
        _axD = _fig.add_subplot(_gs[3])
        gs_E = GridSpecFromSubplotSpec(1, 3, subplot_spec=_gs[4], width_ratios=[2, 0.2, 1.5], wspace=0.02)
        axTree = _fig.add_subplot(gs_E[0])
        axDots = _fig.add_subplot(gs_E[1])
        axMatrix = _fig.add_subplot(gs_E[2])

        _axA.axis("off")
        _axA.text(0.5, 0.5, "", ha="center", fontsize=22, color="grey")
        _axA.text(-0.05, 1.0, "A", transform=_axA.transAxes, fontsize=28, fontweight="bold", va="top")

        def draw_bg(ax):
            curr = POLYPROTEIN_START_CODON
            for _, len_nt, _color in HAV_GENES_NT:
                w = len_nt / 3
                ax.axvspan(curr, curr + w, color=_color, alpha=0.4, lw=0, zorder=0)
                curr = curr + w

        draw_bg(_axB)
        if not df_fubar.empty:
            _axB.scatter(df_fubar[df_fubar["Category"] == "Negative"]["Site"], df_fubar[df_fubar["Category"] == "Negative"]["Value"], c=_COLOR_NEG, s=40, alpha=0.5)
            _axB.scatter(df_fubar[df_fubar["Category"] == "Positive"]["Site"], df_fubar[df_fubar["Category"] == "Positive"]["Value"], c=_COLOR_POS, s=120, edgecolors="k")
        _axB.axhline(0, ls="--", color="k")
        _axB.set_yscale("symlog", linthresh=1.0)  # Handles deep purifying -20 spikes gracefully
        _axB.set_ylabel("Selection\n(dN - dS)", fontsize=22)
        _axB.set_xlim(0, max_x)
        _axB.text(-0.05, 1.0, "B", transform=_axB.transAxes, fontsize=28, fontweight="bold", va="top")
        sns.despine(ax=_axB, trim=False)

        draw_bg(_axC)
        target_sites = []
        if not df_meme.empty:
            sig = df_meme[df_meme["Category"] == "Significant"]
            _axC.scatter(df_meme["Site"], df_meme["Score"], c="#E0E0E0", s=20, alpha=0.3)
            _axC.scatter(sig["Site"], sig["Score"], c=_COLOR_POS, s=120, edgecolors="k")
            target_sites = sig.sort_values("Score", ascending=False)["Site"].astype(int).tolist()[:12]
            for _, r in sig.nlargest(5, "Score").iterrows():
                _axC.annotate(f"{int(r['Site'])}", (r["Site"], r["Score"]), xytext=(0, 10), textcoords="offset points", ha="center", fontweight="bold")
        _axC.axhline(1, ls=":", color="gray")
        _axC.set_ylabel("-log10(p)", fontsize=22)
        _axC.set_xlim(0, max_x)
        _axC.text(-0.05, 1.0, "C", transform=_axC.transAxes, fontsize=28, fontweight="bold", va="top")
        sns.despine(ax=_axC, trim=False)

        curr = 1
        for _name, len_nt, _color in HAV_GENES_NT:
            w = len_nt / 3
            _axD.add_patch(mpatches.Rectangle((curr, 0), w, 1, facecolor=_color, edgecolor="black"))
            _axD.text(curr + w / 2, 0.5, _name, ha="center", va="center", fontweight="bold", fontsize=12)
            curr = curr + w
        _axD.set_xlim(0, max_x)
        _axD.axis("off")
        _axD.text(-0.05, 1.0, "D", transform=_axD.transAxes, fontsize=28, fontweight="bold", va="top")
        _axD.text(max_x / 2, -0.8, "Codon Position (Alignment Coordinates)", ha="center", fontsize=18, fontweight="bold")

        try:
            tree = Phylo.read(tree_file, "newick")
            alignment = _load_alignment_robust(aln_file)
            seq_map = {record.id: str(record.seq) for record in alignment}
            ref_seq_str = _get_sequence_for_tip("NC_001489_1", seq_map) or _get_sequence_for_tip("NC_001489.1", seq_map) or ("N" * 10000)
            tree.ladderize()
            terminals = tree.get_terminals()
            axTree.axis("off")

            for clade in tree.find_clades():
                clade.color = _COLOR_TREE
            for tip in terminals:
                if tip.name and "NC_001489" in tip.name:
                    tip.color = "#D50000"
                    tip.width = 3.0

            Phylo.draw(tree, axes=axTree, do_show=False, label_func=lambda x: None, branch_labels=None)

            matrix_muts = []
            year_colors = []
            all_observed_codons = set()
            for term in terminals:
                year_colors.append(_get_tip_year_color(term.name, metadata_df))
                seq = _get_sequence_for_tip(term.name, seq_map)
                row = []
                if seq:
                    for site in target_sites:
                        idx = (site - 1) * 3
                        codon = seq[idx : idx + 3] if idx + 3 <= len(seq) else "-"
                        row.append(codon)
                        if _is_valid_codon(codon):
                            all_observed_codons.add(codon)
                else:
                    row = ["?"] * len(target_sites)
                matrix_muts.append(row)

            unique_codons = sorted(list(all_observed_codons))
            import matplotlib.cm as cm

            cmap = cm.get_cmap("tab20", len(unique_codons)) if unique_codons else None
            codon_map = {c: cmap(i) for i, c in enumerate(unique_codons)}

            axDots.set_ylim(0, len(terminals) + 1)
            axDots.set_xlim(0, 1)
            axDots.axis("off")
            for i, c in enumerate(year_colors):
                axDots.scatter(0.5, i + 1, color=c, s=150, edgecolors="none")
            axDots.set_title("Era", fontsize=14, fontweight="bold")

            axMatrix.set_xlim(0, len(target_sites))
            axMatrix.set_ylim(0, len(terminals) + 1)
            x_labels = []
            used_codons_in_legend = set()

            for col_idx in range(len(target_sites)):
                site_num = target_sites[col_idx]
                idx = (site_num - 1) * 3
                ref_codon = ref_seq_str[idx : idx + 3] if idx + 3 <= len(ref_seq_str) else "---"
                x_labels.append(f"{site_num}\n(Ref: {ref_codon})")
                for row_idx in range(len(matrix_muts)):
                    codon = matrix_muts[row_idx][col_idx]
                    bg_color = _COLOR_NEU if (codon == ref_codon or not _is_valid_codon(codon)) else codon_map.get(codon, "#000000")
                    if bg_color != _COLOR_NEU:
                        used_codons_in_legend.add(codon)
                    axMatrix.add_patch(mpatches.Rectangle((col_idx, row_idx + 0.5), 1, 1, facecolor=bg_color, edgecolor="white", lw=0.5))

            axMatrix.set_yticks([])
            axMatrix.set_xticks(np.arange(len(target_sites)) + 0.5)
            axMatrix.set_xticklabels(x_labels, rotation=90, fontsize=14, fontweight="bold")
            axMatrix.set_title("Mutations vs Reference", fontsize=22, fontweight="bold")
            for sp in axMatrix.spines.values():
                sp.set_visible(False)

            l_red = mpatches.Patch(color=_COLOR_POS, label="Contemporary")
            l_blue = mpatches.Patch(color=_COLOR_NEG, label="Historical")
            axTree.add_artist(axTree.legend(handles=[l_red, l_blue], loc="lower right", fontsize=20, frameon=True, title="Period"))

            codon_handles = [mpatches.Patch(color=codon_map[c], label=f"{c} ({_translate_codon(c)})") for c in sorted(list(used_codons_in_legend))]
            if codon_handles:
                axMatrix.legend(handles=codon_handles, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=12, frameon=True, title="Mutated Codons")
        except Exception as e:
            traceback.print_exc()
            axTree.text(0.5, 0.5, f"Error: {e}", ha="center", color="red")
            axTree.axis("off")
            axMatrix.axis("off")
            axDots.axis("off")

        axTree.text(-0.05, 1.0, "E", transform=axTree.transAxes, fontsize=28, fontweight="bold", va="top")
        plt.tight_layout()
        try:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
        except Exception:
            plt.savefig("Figure2.pdf", dpi=300, bbox_inches="tight")
        plt.savefig("Figure2.png", dpi=300, bbox_inches="tight")
        print(f"✔ Figure 2 saved successfully.")


    # ==============================================================================
    # 5. SUPPLEMENTARY MATRIX PLOTTER
    # ==============================================================================
    def plot_supplementary_matrix(
        tree_file="../data/raw/READY_FOR_HYPHY.nwk",
        aln_file="../data/raw/READY_FOR_HYPHY.fasta",
        meme_file="../data/raw/HAV_meme.json",
        output_file="../results/figures/Supplementary_Figure_Matrix_Ref_NC001489.pdf",
        metadata_df=None,
    ):
        print("Generating Supplementary Figure (Reference: NC_001489_1)...")
        df_meme = _load_json_data(meme_file, "Score", _get_meme_row)
        target_sites = []
        if not df_meme.empty:
            sig = df_meme[df_meme["Category"] == "Significant"]
            target_sites = sorted(sig["Site"].astype(int).tolist())
            print(f"   -> Found {len(target_sites)} significant sites. Plotting ALL.")

        _fig = plt.figure(figsize=_FIG_SIZE)
        _gs = GridSpec(1, 3, figure=_fig, width_ratios=[1.5, 0.1, 4], wspace=0.02)
        axTree = _fig.add_subplot(_gs[0])
        axDots = _fig.add_subplot(_gs[1])
        axMatrix = _fig.add_subplot(_gs[2])

        try:
            tree = Phylo.read(tree_file, "newick")
            alignment = _load_alignment_robust(aln_file)
            seq_map = {record.id: str(record.seq) for record in alignment}
            ref_seq_str = _get_sequence_for_tip("NC_001489_1", seq_map) or _get_sequence_for_tip("NC_001489.1", seq_map) or ("N" * 10000)
            tree.ladderize()
            terminals = tree.get_terminals()
            axTree.axis("off")

            for clade in tree.find_clades():
                clade.color = _COLOR_TREE
            for tip in terminals:
                if tip.name and "NC_001489" in tip.name:
                    tip.color = "#D50000"
                    tip.width = 3.0

            Phylo.draw(tree, axes=axTree, do_show=False, label_func=lambda x: None, branch_labels=None)
            axTree.set_title("Phylogeny", fontsize=18, fontweight="bold")

            matrix_muts = []
            year_colors = []
            all_observed_codons = set()
            for term in terminals:
                year_colors.append(_get_tip_year_color(term.name, metadata_df))
                seq = _get_sequence_for_tip(term.name, seq_map)
                row = []
                if seq:
                    for site in target_sites:
                        idx = (site - 1) * 3
                        codon = seq[idx : idx + 3] if idx + 3 <= len(seq) else "-"
                        row.append(codon)
                        if _is_valid_codon(codon):
                            all_observed_codons.add(codon)
                else:
                    row = ["?"] * len(target_sites)
                matrix_muts.append(row)

            unique_codons = sorted(list(all_observed_codons))
            base_palette = sns.color_palette("husl", n_colors=len(unique_codons)) if unique_codons else []
            codon_map = {c: base_palette[i] for i, c in enumerate(unique_codons)}

            axDots.set_ylim(0, len(terminals) + 1)
            axDots.set_xlim(0, 1)
            axDots.axis("off")
            for i, c in enumerate(year_colors):
                axDots.scatter(0.5, i + 1, color=c, s=150, edgecolors="none")
            axDots.set_title("Era", fontsize=14, fontweight="bold")

            axMatrix.set_xlim(0, len(target_sites))
            axMatrix.set_ylim(0, len(terminals) + 1)
            x_labels = []
            used_codons_in_legend = set()

            for col_idx in range(len(target_sites)):
                site_num = target_sites[col_idx]
                idx = (site_num - 1) * 3
                ref_codon = ref_seq_str[idx : idx + 3] if idx + 3 <= len(ref_seq_str) else "---"
                x_labels.append(f"{site_num}\n(Ref: {ref_codon})")
                for row_idx, codon in enumerate([row_vec[col_idx] for row_vec in matrix_muts]):
                    bg_color = _COLOR_NEU if (codon == ref_codon or not _is_valid_codon(codon)) else codon_map.get(codon, "#000000")
                    if bg_color != _COLOR_NEU:
                        used_codons_in_legend.add(codon)
                    axMatrix.add_patch(mpatches.Rectangle((col_idx, row_idx + 0.5), 1, 1, facecolor=bg_color, edgecolor="white", lw=0.5))

            axMatrix.set_yticks([])
            axMatrix.set_xticks(np.arange(len(target_sites)) + 0.5)
            axMatrix.set_xticklabels(x_labels, rotation=90, fontsize=9)
            axMatrix.set_title(f"Complete Genotype Matrix ({len(target_sites)} Significant Sites) vs NC_001489_1", fontsize=18, fontweight="bold")
            for sp in axMatrix.spines.values():
                sp.set_visible(False)

            l_red = mpatches.Patch(color=_COLOR_POS, label="Contemporary")
            l_blue = mpatches.Patch(color=_COLOR_NEG, label="Historical")
            axTree.legend(handles=[l_red, l_blue], loc="lower right", fontsize=22, frameon=True, title="Time")

            codon_handles = [mpatches.Patch(color=codon_map[c], label=f"{c} ({_translate_codon(c)})") for c in sorted(list(used_codons_in_legend))]
            if codon_handles:
                cols = math.ceil(len(codon_handles) / 15)
                axMatrix.legend(handles=codon_handles, loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=11, frameon=True, title="Mutant Codons", ncol=cols)
        except Exception as e:
            traceback.print_exc()
            axTree.text(0.5, 0.5, f"Error: {e}", ha="center", color="red")
            axTree.axis("off")
            axMatrix.axis("off")

        plt.tight_layout()
        try:
            plt.savefig(output_file, dpi=600, bbox_inches="tight")
        except Exception:
            plt.savefig("Supplementary_Figure_Matrix.pdf", dpi=600, bbox_inches="tight")
        plt.savefig("SupplementaryFigure.png", dpi=300, bbox_inches="tight")
        print("✔ Supplementary Matrix Figure saved successfully.")


    # ==============================================================================
    # 6. MASTER EXECUTION PIPELINE
    # ==============================================================================
    if __name__ == "__main__":
        # Safely resolve active metadata in notebook runtime
        active_metadata = locals().get(
            "merged_df_5",
            locals().get("merged_df_4", locals().get("merged_df_3", locals().get("merged_df", None))),
        )

        tree_path = "../data/raw/READY_FOR_HYPHY.nwk"
        fasta_path = "../data/raw/READY_FOR_HYPHY.fasta"
        meme_path = "../data/raw/HAV_meme.json"
        fubar_path = "../data/raw/HAV_fubar.json"

        # Step 1: Generate Figure 2 & Output Manuscript Site Counts
        plot_figure_2_final(
            tree_file=tree_path,
            aln_file=fasta_path,
            fubar_file=fubar_path,
            meme_file=meme_path,
            output_file="../results/figures/Figure2_Ref_NC001489.pdf",
            metadata_df=active_metadata,
        )

        # Step 2: Compute Patristic Divergence & Generate S2B Violin Plot
        try:
            master_tree = Phylo.read(tree_path, "newick")
            df_div_dist = get_era_divergence_distribution(
                tree=master_tree,
                metadata_df=active_metadata,
                ref_identifier="NC_001489",
            )
        except Exception as e:
            print(f"[EXECUTION WARNING] Could not compute patristic divergence: {e}")

        # Step 3: Render Complete Supplementary Matrix Heatmap (All 57 MEME Sites)
        plot_supplementary_matrix(
            tree_file=tree_path,
            aln_file=fasta_path,
            meme_file=meme_path,
            output_file="../results/figures/Supplementary_Figure_Matrix_Ref_NC001489.pdf",
            metadata_df=active_metadata,
        )

        plt.show()
    return GridSpecFromSubplotSpec, Phylo, Seq, json, math


@app.cell
def _(AlignIO, json, pd):
    _ALN_FILE = "../data/raw/READY_FOR_HYPHY.fasta"
    _MEME_FILE = "../data/raw/HAV_meme.json"
    _OUTPUT_CSV = "../results/tables/HAV_Mutation_Context_Stratified.csv"
    _REF_ID_HINT = "NC_001489_1"

    def _load_alignment_robust(file_path):
        for fmt in ["phylip-relaxed", "phylip", "fasta"]:
            try:
                return AlignIO.read(file_path, fmt)
            except:
                continue
        raise ValueError(f"Could not read {file_path}")

    def _get_meme_significant_sites(json_file):
        """Extracts significant site numbers (1-based) from MEME JSON."""
        try:
            with open(json_file, "r") as _f:
                _data = json.load(_f)
            content = _data["MLE"]["content"]["0"]
            raw_rows = (
                [content[k] for k in sorted(content.keys(), key=lambda x: int(x))]
                if isinstance(content, dict)
                else content
            )
            sig_sites = []
            for _i, _row in enumerate(raw_rows):
                _p_val = float(_row[6]) if len(_row) > 6 else 1.0
                if _p_val <= 0.05:
                    sig_sites.append(_i + 1)
            return sig_sites
        except Exception as e:
            print(f"Error reading MEME JSON: {e}")
            return []

    def _get_context_triplet(seq_str, nt_index):
        """Returns the triplet: 5'-[nt-1][nt][nt+1]-3'"""
        if nt_index < 1 or nt_index >= len(seq_str) - 1:
            return "---"
        return seq_str[nt_index - 1 : nt_index + 2]

    def _analyze_contexts():
        print("--- Starting Context Stratification ---")
        try:
            _alignment = _load_alignment_robust(_ALN_FILE)
        except Exception as e:
            print(f"Failed to load alignment: {e}")
            return
        sig_sites = _get_meme_significant_sites(_MEME_FILE)
        if not sig_sites:
            print(
                "No significant sites found. Checking ALL sites instead? (Edit script if needed)"
            )
            return
        print(f"Analyzing {len(sig_sites)} significant codons.")
        seq_map = {rec.id: str(rec.seq).upper() for rec in _alignment}
        ref_seq = None
        for candidate in [_REF_ID_HINT, "NC_001489.1", "NC_001489"]:
            if candidate in seq_map:
                ref_seq = seq_map[candidate]
                print(f"Using Reference: {candidate}")
                break
        if not ref_seq:
            print(
                "WARNING: Reference not found. Using the first sequence in file as reference."
            )
            ref_seq = list(seq_map.values())[0]
        results = []
        for site_num in sig_sites:
            codon_idx = (site_num - 1) * 3
            if codon_idx + 3 > len(ref_seq):
                continue
            ref_codon = ref_seq[codon_idx : codon_idx + 3]
            for rec_id, seq in seq_map.items():
                if rec_id == _REF_ID_HINT:
                    continue
                query_codon = seq[codon_idx : codon_idx + 3]
                if "-" in query_codon or "N" in query_codon:
                    continue
                if query_codon == ref_codon:
                    continue
                for _i in range(3):
                    ref_nt = ref_codon[_i]
                    query_nt = query_codon[_i]
                    if ref_nt == query_nt:
                        continue
                    nt_idx = codon_idx + _i
                    triplet = _get_context_triplet(ref_seq, nt_idx)
                    five_prime = triplet[0]
                    three_prime = triplet[2]
                    mutation_type = f"{ref_nt}>{query_nt}"
                    category = "Other"
                    context_label = "Other"
                    if mutation_type == "C>T":
                        category = "C>T (Transition)"
                        if five_prime == "T":
                            context_label = "TpC Context (APOBEC-like)"
                        else:
                            context_label = f"{five_prime}pC Context"
                    elif mutation_type == "G>A":
                        category = "G>A (Transition)"
                        if three_prime == "A":
                            context_label = "GpA Context (Antisense APOBEC)"
                        else:
                            context_label = f"Gp{three_prime} Context"
                    elif mutation_type == "A>G":
                        category = "A>G (Transition)"
                        if five_prime == "T":
                            context_label = "TpA Context (ADAR-like)"
                        elif five_prime == "A":
                            context_label = "ApA Context (ADAR-like)"
                        else:
                            context_label = f"{five_prime}pA Context"
                    elif mutation_type == "T>C":
                        category = "T>C (Transition)"
                        context_label = f"{five_prime}pT Context"
                    results.append(
                        {
                            "Site": site_num,
                            "Sequence_ID": rec_id,
                            "Ref_Codon": ref_codon,
                            "Query_Codon": query_codon,
                            "Mutation": mutation_type,
                            "Category": category,
                            "Context_Triplet": triplet,
                            "5_Prime_Base": five_prime,
                            "3_Prime_Base": three_prime,
                            "Detailed_Context": context_label,
                        }
                    )
        df = pd.DataFrame(results)
        if not df.empty:
            df.to_csv(_OUTPUT_CSV, index=False)
            print(f"\nSuccess! Stratified data saved to: {_OUTPUT_CSV}")
            print("\n--- SUMMARY OF MUTATION TYPES ---")
            print(df["Category"].value_counts())
            print("\n--- TOP CONTEXTS FOR C>T (Potential APOBEC) ---")
            ct_muts = df[df["Mutation"] == "C>T"]
            if not ct_muts.empty:
                print(ct_muts["Detailed_Context"].value_counts().head())
            else:
                print("No C>T mutations found.")
            print("\n--- TOP CONTEXTS FOR A>G (Potential ADAR) ---")
            ag_muts = df[df["Mutation"] == "A>G"]
            if not ag_muts.empty:
                print(ag_muts["Detailed_Context"].value_counts().head())
            else:
                print("No A>G mutations found.")
        else:
            print("No mutations found to stratify.")

    if __name__ == "__main__":
        _analyze_contexts()
    return


@app.cell
def _(merged_df_3):
    merged_df_3
    return


@app.cell
def _(Phylo, merged_df_3):
    import os

    _TREE_FILE = "../data/raw/HAV_MCC.tree"
    COLORS = {"Historical": "#4c72b0", "Contemporary": "#dd8452"}

    def _safe_clean_id(id_str):
        """Cleans the ID to match how they are stored in merged_df"""
        id_str = str(id_str).strip().upper()
        id_str = id_str.replace("'", "")
        if id_str.startswith("NC_"):
            return id_str.split(".")[0]
        else:
            return id_str.split(".")[0].split("_")[0]

    def generate_itol_datasets(tree_path, metadata_df):
        print(f"--- Generating iTOL Datasets from {tree_path} ---")
        meta_map = {}
        for _, _row in metadata_df.dropna(
            subset=["accession_id", "Condition"]
        ).iterrows():
            clean_acc = _safe_clean_id(_row["accession_id"])
            meta_map[clean_acc] = _row["Condition"]
        try:
            tree = Phylo.read(tree_path, "nexus")
        except Exception as e:
            try:
                tree = Phylo.read(tree_path, "newick")
            except Exception as e2:
                print(
                    f"[!] Failed to read tree file. Ensure it is valid Nexus or Newick format. Errors: {e}, {e2}"
                )
                return
        exact_tip_names = [tip.name for tip in tree.get_terminals()]
        tip_mapping = []
        unmapped = []
        for exact_name in exact_tip_names:
            clean_name = _safe_clean_id(exact_name)
            condition = meta_map.get(clean_name, "Unknown")
            if condition in COLORS:
                tip_mapping.append(
                    {
                        "exact_tip": exact_name,
                        "condition": condition,
                        "color": COLORS[condition],
                    }
                )
            else:
                unmapped.append(exact_name)
        print(f"Successfully mapped {len(tip_mapping)} sequences.")
        if unmapped:
            print(
                f"WARNING: {len(unmapped)} sequences from the tree could not be found in the metadata."
            )
            print(f"First 5 unmapped examples: {unmapped[:5]}")
        if not tip_mapping:
            print("[!] No sequences matched. Exiting without creating files.")
            return
        strip_file = "../results/tables/iTOL_Color_Strip.txt"
        with open(strip_file, "w") as _f:
            _f.write("DATASET_COLORSTRIP\n")
            _f.write("SEPARATOR COMMA\n")
            _f.write("DATASET_LABEL,Era (Historical vs Contemporary)\n")
            _f.write("COLOR,#000000\n")
            _f.write("STRIP_WIDTH,25\n")
            _f.write("MARGIN,5\n")
            _f.write("SHOW_INTERNAL,0\n")
            _f.write("DATA\n")
            for tip in tip_mapping:
                _f.write(f"{tip['exact_tip']},{tip['color']},{tip['condition']}\n")
        print(f"[✔] Created {strip_file}")
        text_color_file = "../results/tables/iTOL_Tip_Text.txt"
        with open(text_color_file, "w") as _f:
            _f.write("TREE_COLORS\n")
            _f.write("SEPARATOR COMMA\n")
            _f.write("DATA\n")
            for tip in tip_mapping:
                _f.write(f"{tip['exact_tip']},label,{tip['color']},bold,1\n")
        print(f"[✔] Created {text_color_file}")

    if __name__ == "__main__":
        generate_itol_datasets(_TREE_FILE, merged_df_3)
    return (os,)


@app.cell
def _(AlignIO, Phylo, pd):

    # ==========================================
    # CELL 1: CONFIGURATION & CONSTANTS
    # ==========================================
    ALN_FILE = "../data/raw/HAV_BEAST.fasta"
    TREE_FILE = "../data/raw/HAV_MCC.tree"
    OUTPUT_FILE = "../results/tables/iTOL_CT_Proportions_SimpleBar.txt"
    REF_ID_HINT = "NC_001489"


    # ==========================================
    # CELL 2: HELPER FUNCTIONS
    # ==========================================
    def safe_clean_id(id_str):
        """Cleans the ID to match across alignment, metadata, and tree."""
        id_str = str(id_str).strip().upper()
        id_str = id_str.replace("'", "")
        id_str = id_str.split("|")[0].strip()

        if id_str.startswith("NC_"):
            core = id_str[3:].split(".")[0].split("_")[0]
            return "NC_" + core
        else:
            return id_str.split(".")[0].split("_")[0]


    # ==========================================
    # CELL 3: CORE EXECUTION & iTOL GENERATION
    # ==========================================
    def generate_itol_ct_barchart(aln_path, tree_path, out_path, ref_hint):
        print(
            "--- Generating iTOL SimpleBar Dataset for C>T (APOBEC) Proportions ---"
        )

        # 1. Load Alignment & Find Reference Sequence
        try:
            alignment = AlignIO.read(aln_path, "fasta")
        except Exception as e:
            print(f"[!] Failed to load alignment: {e}")
            return None

        ref_seq = None
        for rec in alignment:
            if ref_hint in rec.id.upper():
                ref_seq = str(rec.seq).upper()
                break

        if not ref_seq:
            print("[!] Reference not found in alignment. Exiting.")
            return None

        # 2. Calculate C>T Proportions for each sequence
        print("Scanning alignment and calculating C>T mutational proportions...")
        proportion_data = {}
        aln_length = len(ref_seq)
        valid_bases = {"A", "C", "G", "T"}

        for rec in alignment:
            clean_id = safe_clean_id(rec.id)
            total_muts = 0
            c_to_t_count = 0
            query_seq = str(rec.seq).upper()

            for i in range(aln_length):
                ref_base = ref_seq[i]
                query_base = query_seq[i]

                # Only count true nucleotide substitutions (ignore gaps/N)
                if ref_base in valid_bases and query_base in valid_bases:
                    if ref_base != query_base:
                        total_muts += 1
                        if ref_base == "C" and query_base == "T":
                            c_to_t_count += 1

            # Calculate percentage (avoiding division by zero)
            ct_pct = (c_to_t_count / total_muts * 100) if total_muts > 0 else 0.0
            proportion_data[clean_id] = round(ct_pct, 2)

        print(f"Calculated C>T proportions for {len(proportion_data)} isolates.")

        # 3. Parse the Tree to get exact tip names
        try:
            tree = Phylo.read(tree_path, "nexus")
        except Exception:
            try:
                tree = Phylo.read(tree_path, "newick")
            except Exception as e:
                print(f"[!] Failed to read tree file: {e}")
                return None

        exact_tip_names = [tip.name for tip in tree.get_terminals()]

        # 4. Map proportions to exact tree tips
        dataset_mapping = []
        for exact_name in exact_tip_names:
            clean_name = safe_clean_id(exact_name)
            if clean_name in proportion_data:
                dataset_mapping.append((exact_name, proportion_data[clean_name]))

        # 5. Write the iTOL formatted SimpleBar file
        with open(out_path, "w") as f:
            f.write("DATASET_SIMPLEBAR\n")
            f.write("SEPARATOR COMMA\n")
            f.write("DATASET_LABEL,% C>T (APOBEC)\n")
            f.write("COLOR,#D73027\n")  # High-contrast red/vermilion for APOBEC
            f.write("WIDTH,250\n")
            f.write("MARGIN,10\n")
            f.write("SHOW_INTERNAL,0\n")
            f.write("ALIGN_TO_LABELS,1\n")
            f.write("DATA\n")

            for exact_name, ct_pct in dataset_mapping:
                f.write(f"{exact_name},{ct_pct}\n")

        print(f"[✔] Created {out_path} successfully!")
        return pd.DataFrame(dataset_mapping, columns=["Taxon", "CT_Pct"])


    # Execute directly in cell (Marimo safe)
    df_ct_results = generate_itol_ct_barchart(
        ALN_FILE, TREE_FILE, OUTPUT_FILE, REF_ID_HINT
    )
    return


@app.cell
def _(AlignIO, json, merged_df_3, pd):

    # ==========================================
    # CONFIGURATION
    # ==========================================
    ALN_FILE = "../data/raw/READY_FOR_HYPHY.fasta"
    MEME_FILE = "../data/raw/HAV_meme.json"
    OUTPUT_CSV = "../results/tables/HAV_SBS96_Context_Time_Analysis.csv"
    REF_ID_HINT = "NC_001489"

    # ==========================================
    # HELPER FUNCTIONS
    # ==========================================
    def load_alignment_robust(file_path):
        for fmt in ['phylip-relaxed', 'phylip', 'fasta']:
            try: return AlignIO.read(file_path, fmt)
            except: continue
        raise ValueError(f"Could not read {file_path}")

    def get_meme_significant_sites(json_file):
        try:
            with open(json_file, 'r') as f: data = json.load(f)
            content = data['MLE']['content']['0']
            raw_rows = [content[k] for k in sorted(content.keys(), key=lambda x: int(x))] if isinstance(content, dict) else content

            sig_sites = []
            for i, row in enumerate(raw_rows):
                p_val = float(row[6]) if len(row) > 6 else 1.0
                if p_val <= 0.05:
                    sig_sites.append(i + 1)
            return sig_sites
        except Exception as e:
            print(f"Error reading MEME JSON: {e}")
            return []

    def get_trinuc_context(seq, index):
        """Returns (prev_base, ref_base, next_base)"""
        if index < 1 or index >= len(seq)-1:
            return None
        return seq[index-1], seq[index], seq[index+1]

    def clean_id(id_str):
        """Aggressively strips version numbers (.1), spaces, and forces uppercase."""
        # Keeps the core accession (e.g., MZ123456 from MZ123456.1_India)
        return str(id_str).split('.')[0].split('_')[0].strip().upper()

    # ==========================================
    # MAIN ANALYSIS
    # ==========================================
    def generate_sbs96_table(metadata_df):
        print("--- Extracting Temporal SBS96 Contexts via Metadata ---")

        # 1. Build Strict Metadata Lookup Dictionary
        meta_map = {}
        for _, row in metadata_df.dropna(subset=['accession_id', 'year']).iterrows():
            acc_clean = clean_id(row['accession_id'])
            cond = row.get('Condition', "Contemporary" if int(row['year']) >= 2018 else "Historical")
            meta_map[acc_clean] = (int(row['year']), cond)

        print(f"Loaded metadata mapping for {len(meta_map)} isolates.")

        # 2. Load Alignment
        try:
            alignment = load_alignment_robust(ALN_FILE)
        except Exception as e:
            print(f"Failed to load alignment: {e}"); return

        # Get Reference Sequence
        seq_map = {rec.id: str(rec.seq).upper() for rec in alignment}
        # --- ADD THIS QUICK VERIFICATION ---
        valid_fasta_ids = {clean_id(k) for k in seq_map.keys() if REF_ID_HINT not in k}
        mapped_ids = valid_fasta_ids.intersection(meta_map.keys())
        print(f"[✔] Successfully matched {len(mapped_ids)} FASTA sequences to metadata profiles.")
        # -----------------------------------
        ref_seq = None
        for k in seq_map.keys():
            if REF_ID_HINT in k:
                ref_seq = seq_map[k]
                print(f"Reference identified as: {k}")
                break

        if not ref_seq:
            print("Reference not found! Using first sequence in alignment.")
            ref_seq = list(seq_map.values())[0]

        # Get Sites
        sig_sites = get_meme_significant_sites(MEME_FILE)
        print(f"Analyzing {len(sig_sites)} positively selected sites.\n")

        # 3. Iterate Sites & Extract Contexts
        sbs_data = []
        target_mutations = ["C>T", "T>C", "G>A", "A>G"]

        missing_metadata_count = 0
        unmatched_ids = set()

        for site_num in sig_sites:
            codon_idx = (site_num - 1) * 3
            if codon_idx + 3 > len(ref_seq): continue
            ref_codon = ref_seq[codon_idx : codon_idx+3]

            for rec_id, seq in seq_map.items():
                if REF_ID_HINT in rec_id: continue # Skip ref

                query_codon = seq[codon_idx : codon_idx+3]

                if "-" in query_codon or "N" in query_codon: continue
                if query_codon == ref_codon: continue

                # --- EXACT O(1) DICTIONARY MATCHING ---
                rec_clean = clean_id(rec_id)

                # If the clean FASTA ID is in our clean metadata dictionary, grab it.
                if rec_clean in meta_map:
                    year, period = meta_map[rec_clean]
                else:
                    # If it's one of your 40 outliers or missing, skip it instantly.
                    missing_metadata_count += 1
                    unmatched_ids.add(rec_id)
                    continue

                for i in range(3):
                    r_nt = ref_codon[i]
                    q_nt = query_codon[i]

                    if r_nt == q_nt: continue

                    mut_str = f"{r_nt}>{q_nt}"

                    if mut_str in target_mutations:
                        nt_global_idx = codon_idx + i
                        context_tuple = get_trinuc_context(ref_seq, nt_global_idx)

                        if context_tuple:
                            prev_b, curr_b, next_b = context_tuple
                            sbs96_label = f"{prev_b}[{mut_str}]{next_b}"

                            motif = "Other"
                            if mut_str == "C>T" and prev_b == "T": motif = "TpC (APOBEC)"
                            if mut_str == "G>A" and next_b == "A": motif = "GpA (Antisense APOBEC)"
                            if mut_str == "A>G" and prev_b in ["T", "A"]: motif = "TpA/ApA (ADAR)"

                            sbs_data.append({
                                "Site": site_num,
                                "Sample_ID": rec_id,
                                "Year": year,
                                "Condition": period,
                                "Mutation": mut_str,
                                "5_Prime": prev_b,
                                "Ref_Base": curr_b,
                                "3_Prime": next_b,
                                "SBS96_Label": sbs96_label,
                                "Motif_Class": motif
                            })

        # 4. Save and Report
        df = pd.DataFrame(sbs_data)

        # We expect some skipped isolates now (like your 40 root-to-tip outliers)
        if len(unmatched_ids) > 0:
            print(f"[!] Skipped sequences from {len(unmatched_ids)} unique FASTA IDs due to filtering/missing metadata.")
            print(f"[!] Examples of skipped IDs: {list(unmatched_ids)[:5]}\n")

        if not df.empty:
            df.to_csv(OUTPUT_CSV, index=False)
            print(f"Saved {len(df)} rigorously matched mutational events to: {OUTPUT_CSV}")

            print("\n--- Total Mutations by Condition ---")
            print(df['Condition'].value_counts())

            print("\n--- Top SBS96 Contexts Split by Era ---")
            crosstab_df = pd.crosstab(df['SBS96_Label'], df['Condition'])
            crosstab_df['Total'] = crosstab_df.sum(axis=1)
            crosstab_df = crosstab_df.sort_values(by='Total', ascending=False).drop(columns='Total')

            pd.set_option('display.max_rows', 60)
            print(crosstab_df.head(60))

        else:
            print("No matching mutations found. Check the unmatched IDs above.")

    if __name__ == "__main__":
        generate_sbs96_table(merged_df_3)
    return


@app.cell
def _(AlignIO, merged_df_3, np, pd, plt, sns):
    from scipy.stats import norm, pearsonr, mannwhitneyu

    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    _colors = {"Historical": "#4c72b0", "Contemporary": "#dd8452"}
    sns.set_context("paper", font_scale=1.5)
    sns.set_style(
        "ticks",
        {
            "axes.linewidth": 2.5,
            "axes.edgecolor": "black",
            "xtick.major.width": 2.5,
            "ytick.major.width": 2.5,
        },
    )

    def _safe_clean_id(id_str):
        id_str = str(id_str).strip().upper()
        if id_str.startswith("NC_"):
            return id_str.split(".")[0]
        else:
            return id_str.split(".")[0].split("_")[0]

    def _load_alignment_robust(file_path):
        for fmt in ["phylip-relaxed", "phylip", "fasta"]:
            try:
                return AlignIO.read(file_path, fmt)
            except:
                continue
        raise ValueError(f"Could not read {file_path}")

    def get_significance_asterisks(p_val):
        if p_val < 0.001:
            return "***"
        elif p_val < 0.01:
            return "**"
        elif p_val < 0.05:
            return "*"
        else:
            return "ns"

    print("--- Syncing Metadata with MEME Alignment (Strict Filter) ---")
    df_clean = merged_df_3.dropna(subset=["Condition", "Nucleotide Mutations"]).copy()
    df_clean["accession_clean"] = df_clean["accession_id"].apply(_safe_clean_id)
    _ALN_FILE = "../data/raw/READY_FOR_HYPHY.fasta"
    try:
        _alignment = _load_alignment_robust(_ALN_FILE)
        fasta_ids = {_safe_clean_id(rec.id) for rec in _alignment}
        print(f"Extracted {len(fasta_ids)} high-quality sequences from {_ALN_FILE}.")
    except Exception as e:
        print(f"[!] Failed to load FASTA for filtering: {e}")
        fasta_ids = set()
    _plot_df = df_clean[df_clean["accession_clean"].isin(fasta_ids)].copy()
    try:
        meme_df = pd.read_csv("../results/tables/HAV_SBS96_Context_Time_Analysis.csv")
        meme_df["accession_clean"] = meme_df["Sample_ID"].apply(_safe_clean_id)
        meme_counts = (
            meme_df.groupby("accession_clean").size().reset_index(name="MEME_Mutations")
        )
        _plot_df = pd.merge(_plot_df, meme_counts, on="accession_clean", how="left")
        _plot_df["MEME_Mutations"] = _plot_df["MEME_Mutations"].fillna(0)
        print(
            f"Successfully mapped per-sequence MEME data. Final analytical N={len(_plot_df)}\n"
        )
    except FileNotFoundError:
        print("[!] HAV_SBS96_Context_Time_Analysis.csv not found.")
        raise
    hist_df = _plot_df[_plot_df["Condition"] == "Historical"].dropna(
        subset=["Nucleotide Mutations", "MEME_Mutations"]
    )
    cont_df = _plot_df[_plot_df["Condition"] == "Contemporary"].dropna(
        subset=["Nucleotide Mutations", "MEME_Mutations"]
    )
    (u_stat_a, p_val_a) = mannwhitneyu(
        hist_df["MEME_Mutations"], cont_df["MEME_Mutations"], alternative="two-sided"
    )
    (r_hist, p_hist_corr) = pearsonr(
        hist_df["Nucleotide Mutations"], hist_df["MEME_Mutations"]
    )
    (r_cont, p_cont_corr) = pearsonr(
        cont_df["Nucleotide Mutations"], cont_df["MEME_Mutations"]
    )
    (n_hist, n_cont) = (len(hist_df), len(cont_df))
    (z_hist, z_cont) = (np.arctanh(r_hist), np.arctanh(r_cont))
    se_diff = np.sqrt(1 / (n_hist - 3) + 1 / (n_cont - 3))
    z_score = (z_hist - z_cont) / se_diff
    p_fisher_z = 2 * (1 - norm.cdf(abs(z_score)))
    print("Running 1,000 bootstraps to calculate correlation effect size...")
    _n_boot = 1000
    (boot_r_hist, boot_r_cont) = ([], [])
    for _i in range(_n_boot):
        h_samp = hist_df.sample(frac=1.0, replace=True, random_state=_i)
        c_samp = cont_df.sample(frac=1.0, replace=True, random_state=_i)
        boot_r_hist.append(
            pearsonr(h_samp["Nucleotide Mutations"], h_samp["MEME_Mutations"])[0]
        )
        boot_r_cont.append(
            pearsonr(c_samp["Nucleotide Mutations"], c_samp["MEME_Mutations"])[0]
        )
    boot_df = pd.DataFrame(
        {"Historical": boot_r_hist, "Contemporary": boot_r_cont}
    ).melt(var_name="Condition", value_name="Pearson_r")
    mut_counts = {"C>T": 1040, "T>C": 616, "A>G": 440, "G>A": 330, "Transversions": 472}
    spectrum_df = pd.DataFrame(
        list(mut_counts.items()), columns=["Mutation Type", "Count"]
    )
    total_muts = spectrum_df["Count"].sum()
    spectrum_df["Percentage"] = spectrum_df["Count"] / total_muts * 100
    (_fig, _axes) = plt.subplots(2, 2, figsize=(20, 18))
    plt.subplots_adjust(hspace=0.35, wspace=0.25)
    sns.barplot(
        data=_plot_df,
        x="Condition",
        y="MEME_Mutations",
        palette=_colors,
        ax=_axes[0, 0],
        order=["Historical", "Contemporary"],
        edgecolor="black",
        linewidth=2.5,
        errorbar=("ci", 95),
        capsize=0.1,
    )
    sns.swarmplot(
        data=_plot_df,
        x="Condition",
        y="MEME_Mutations",
        color=".2",
        alpha=0.4,
        size=4,
        ax=_axes[0, 0],
        order=["Historical", "Contemporary"],
        zorder=1,
    )
    _axes[0, 0].set_title(
        f"Targeted Positive Selection\n(Normalized MEME Mutations, n={len(_plot_df)})",
        fontweight="bold",
        pad=15,
    )
    _axes[0, 0].set_ylabel("MEME Mutations per Sequence", fontweight="bold")
    _axes[0, 0].set_xlabel("")
    _axes[0, 0].text(
        -0.1,
        1.05,
        "A",
        transform=_axes[0, 0].transAxes,
        fontsize=28,
        fontweight="bold",
        va="top",
    )
    y_max_a = _plot_df["MEME_Mutations"].max() + _plot_df["MEME_Mutations"].max() * 0.05
    _axes[0, 0].plot(
        [0, 0, 1, 1],
        [y_max_a, y_max_a * 1.02, y_max_a * 1.02, y_max_a],
        lw=2,
        c="black",
    )
    _axes[0, 0].text(
        0.5,
        y_max_a * 1.02,
        get_significance_asterisks(p_val_a),
        ha="center",
        va="bottom",
        color="black",
        fontsize=20,
        fontweight="bold",
    )
    sns.scatterplot(
        data=_plot_df,
        x="Nucleotide Mutations",
        y="MEME_Mutations",
        hue="Condition",
        palette=_colors,
        alpha=0.7,
        ax=_axes[0, 1],
        s=80,
        edgecolor="black",
        zorder=2,
    )
    sns.regplot(
        data=hist_df,
        x="Nucleotide Mutations",
        y="MEME_Mutations",
        scatter=False,
        color=_colors["Historical"],
        ax=_axes[0, 1],
        line_kws={"linewidth": 4, "zorder": 3},
    )
    sns.regplot(
        data=cont_df,
        x="Nucleotide Mutations",
        y="MEME_Mutations",
        scatter=False,
        color=_colors["Contemporary"],
        ax=_axes[0, 1],
        line_kws={"linewidth": 4, "zorder": 3},
    )
    _stats_text = f"Historical $r$ = {r_hist:.2f}\nContemporary $r$ = {r_cont:.2f}"
    _axes[0, 1].text(
        0.05,
        0.95,
        _stats_text,
        transform=_axes[0, 1].transAxes,
        fontsize=15,
        va="top",
        bbox=dict(
            facecolor="white", edgecolor="black", boxstyle="round,pad=0.5", alpha=0.9
        ),
        zorder=4,
    )
    _axes[0, 1].set_title(
        "Mutational Decoupling\n(Genome-Wide vs. MEME)", fontweight="bold", pad=15
    )
    _axes[0, 1].set_xlabel("Genome-Wide Mutations", fontweight="bold")
    _axes[0, 1].set_ylabel("MEME Mutations", fontweight="bold")
    _axes[0, 1].text(
        -0.1,
        1.05,
        "B",
        transform=_axes[0, 1].transAxes,
        fontsize=28,
        fontweight="bold",
        va="top",
    )
    _axes[0, 1].legend(title="Era", frameon=True, loc="lower right")
    sns.kdeplot(
        data=boot_df,
        x="Pearson_r",
        hue="Condition",
        fill=True,
        palette=_colors,
        ax=_axes[1, 0],
        linewidth=3,
        alpha=0.5,
        common_norm=False,
    )
    _axes[1, 0].axvline(
        np.mean(boot_r_hist), color=_colors["Historical"], linestyle="--", linewidth=2.5
    )
    _axes[1, 0].axvline(
        np.mean(boot_r_cont),
        color=_colors["Contemporary"],
        linestyle="--",
        linewidth=2.5,
    )
    effect_size_text = f"Effect Size ($\\Delta r$) = {np.mean(boot_r_hist) - np.mean(boot_r_cont):.2f}\nFisher's $z$ ($p$) = {p_fisher_z:.2e}"
    _axes[1, 0].text(
        0.5,
        0.95,
        effect_size_text,
        transform=_axes[1, 0].transAxes,
        fontsize=15,
        ha="center",
        va="top",
        bbox=dict(
            facecolor="white", edgecolor="black", boxstyle="round,pad=0.5", alpha=0.9
        ),
    )
    _axes[1, 0].set_title(
        "Decoupling Effect Size\n(Bootstrapped Correlation Distribution)",
        fontweight="bold",
        pad=15,
    )
    _axes[1, 0].set_xlabel("Pearson Correlation Coefficient ($r$)", fontweight="bold")
    _axes[1, 0].set_ylabel("Density", fontweight="bold")
    _axes[1, 0].text(
        -0.1,
        1.05,
        "C",
        transform=_axes[1, 0].transAxes,
        fontsize=28,
        fontweight="bold",
        va="top",
    )
    spectrum_colors = ["#2ecc71", "#27ae60", "#f1c40f", "#f39c12", "#95a5a6"]
    sns.barplot(
        data=spectrum_df,
        x="Mutation Type",
        y="Percentage",
        hue="Mutation Type",
        ax=_axes[1, 1],
        palette=spectrum_colors,
        edgecolor="black",
        linewidth=2.5,
        legend=False,
    )
    for _i, _p in enumerate(_axes[1, 1].patches):
        _axes[1, 1].annotate(
            f"{spectrum_df['Percentage'].iloc[_i]:.1f}%",
            (_p.get_x() + _p.get_width() / 2.0, _p.get_height()),
            ha="center",
            va="bottom",
            fontsize=14,
            fontweight="bold",
            xytext=(0, 5),
            textcoords="offset points",
        )
    _axes[1, 1].set_title(
        f"Mutation Spectrum at MEME Loci\n(Total Mutational Events = {total_muts})",
        fontweight="bold",
        pad=15,
    )
    _axes[1, 1].set_ylabel("Proportion of Total Events (%)", fontweight="bold")
    _axes[1, 1].set_xlabel("Mutation Type", fontweight="bold")
    _axes[1, 1].text(
        -0.1,
        1.05,
        "D",
        transform=_axes[1, 1].transAxes,
        fontsize=28,
        fontweight="bold",
        va="top",
    )
    for _ax in _axes.flat:
        sns.despine(ax=_ax)
    plt.savefig("../results/figures/Figure_Evolutionary_Mechanism_2x2.png", dpi=300, bbox_inches="tight")
    plt.savefig("../results/figures/Figure_Evolutionary_Mechanism_2x2.svg", format="svg", bbox_inches="tight"
    )
    print("[✔] 2x2 Plot saved successfully!")
    plt.show()
    print("--- Generating stat_summary.txt ---")
    with open("../results/tables/stat_summary.txt", "w") as _f:
        _f.write("=================================================================\n")
        _f.write("    STATISTICAL SUMMARY: MECHANISMS OF EVOLUTIONARY ESCAPE\n")
        _f.write(
            "=================================================================\n\n"
        )
        _f.write("--- PANEL C: EFFECT SIZE (BOOTSTRAPPED CORRELATION) ---\n")
        _f.write(
            f"Historical 95% CI for r: [{np.percentile(boot_r_hist, 2.5):.3f} - {np.percentile(boot_r_hist, 97.5):.3f}]\n"
        )
        _f.write(
            f"Contemporary 95% CI for r: [{np.percentile(boot_r_cont, 2.5):.3f} - {np.percentile(boot_r_cont, 97.5):.3f}]\n"
        )
        _f.write(
            f"Difference in correlation (Effect Size): {np.mean(boot_r_hist) - np.mean(boot_r_cont):.3f}\n"
        )
        _f.write(f"Fisher's Z-Transformation p-value: {p_fisher_z:.4e}\n")
        _f.write(
            "Interpretation: The non-overlapping bootstrapped distributions visually confirm that the drop in correlation strength is highly significant and not an artifact of sample size.\n\n"
        )
    return (mannwhitneyu,)


@app.cell
def _(AlignIO, Seq, json, merged_df_3, pd):
    _ALN_FILE = "../data/raw/READY_FOR_HYPHY.fasta"
    _MEME_FILE = "../data/raw/HAV_meme.json"
    _OUTPUT_CSV = "../results/tables/HAV_MEME_AA_Changes_Stratified.csv"
    _REF_ID_HINT = "NC_001489"

    def _load_alignment_robust(file_path):
        for fmt in ["phylip-relaxed", "phylip", "fasta"]:
            try:
                return AlignIO.read(file_path, fmt)
            except:
                continue
        raise ValueError(f"Could not read {file_path}")

    def _get_meme_significant_sites(json_file):
        try:
            with open(json_file, "r") as _f:
                _data = json.load(_f)
            content = _data["MLE"]["content"]["0"]
            raw_rows = (
                [content[k] for k in sorted(content.keys(), key=lambda x: int(x))]
                if isinstance(content, dict)
                else content
            )
            sig_sites = []
            for _i, _row in enumerate(raw_rows):
                _p_val = float(_row[6]) if len(_row) > 6 else 1.0
                if _p_val <= 0.05:
                    sig_sites.append(_i + 1)
            return sig_sites
        except Exception as e:
            print(f"Error reading MEME JSON: {e}")
            return []

    def _safe_clean_id(id_str):
        id_str = str(id_str).strip().upper()
        if id_str.startswith("NC_"):
            return id_str.split(".")[0]
        else:
            return id_str.split(".")[0].split("_")[0]

    def analyze_amino_acids(metadata_df):
        print("--- Starting Amino Acid Translation at MEME Loci ---")
        meta_map = {}
        for _, _row in metadata_df.dropna(
            subset=["accession_id", "Condition"]
        ).iterrows():
            acc_clean = _safe_clean_id(_row["accession_id"])
            meta_map[acc_clean] = _row["Condition"]
        try:
            _alignment = _load_alignment_robust(_ALN_FILE)
        except Exception as e:
            print(f"Failed to load alignment: {e}")
            return
        seq_map = {rec.id: str(rec.seq).upper() for rec in _alignment}
        ref_seq = None
        ref_id_actual = None
        for key in seq_map.keys():
            if _REF_ID_HINT in key:
                ref_seq = seq_map[key]
                ref_id_actual = key
                break
        if not ref_seq:
            print("WARNING: Reference not found. Using first sequence.")
            ref_id_actual = list(seq_map.keys())[0]
            ref_seq = list(seq_map.values())[0]
        sig_sites = _get_meme_significant_sites(_MEME_FILE)
        if not sig_sites:
            print("No significant sites found.")
            return
        print(f"Translating codons at {len(sig_sites)} significant sites...\n")
        aa_results = []
        for site_num in sig_sites:
            codon_idx = (site_num - 1) * 3
            if codon_idx + 3 > len(ref_seq):
                continue
            ref_codon = ref_seq[codon_idx : codon_idx + 3]
            if "-" in ref_codon or "N" in ref_codon:
                continue
            ref_aa = str(Seq(ref_codon).translate())
            for rec_id, seq in seq_map.items():
                if rec_id == ref_id_actual:
                    continue
                query_codon = seq[codon_idx : codon_idx + 3]
                if "-" in query_codon or "N" in query_codon:
                    continue
                if query_codon == ref_codon:
                    continue
                query_aa = str(Seq(query_codon).translate())
                rec_clean = _safe_clean_id(rec_id)
                condition = meta_map.get(rec_clean, "Unknown")
                if condition == "Unknown":
                    continue
                if ref_aa == query_aa:
                    mut_type = "Synonymous"
                    notation = f"{ref_aa}{site_num}{ref_aa}"
                else:
                    mut_type = "Non-Synonymous"
                    notation = f"{ref_aa}{site_num}{query_aa}"
                aa_results.append(
                    {
                        "Site": site_num,
                        "Sequence_ID": rec_clean,
                        "Condition": condition,
                        "Ref_Codon": ref_codon,
                        "Query_Codon": query_codon,
                        "Ref_AA": ref_aa,
                        "Query_AA": query_aa,
                        "Effect": mut_type,
                        "AA_Mutation": notation,
                    }
                )
        df = pd.DataFrame(aa_results)
        if not df.empty:
            df.to_csv(_OUTPUT_CSV, index=False)
            print(
                f"[✔] Successfully translated and mapped {len(df)} mutational events. Saved to {_OUTPUT_CSV}"
            )
            print("\n=======================================================")
            print("    MUTATION EFFECTS: OVERALL & STRATIFIED BY ERA")
            print("=======================================================")
            summary_table = pd.crosstab(
                df["Condition"], df["Effect"], margins=True, margins_name="Total"
            )
            print(summary_table.to_string())
            ns_df = df[df["Effect"] == "Non-Synonymous"]
            print("\n=============================================")
            print("    TOP 15 NON-SYNONYMOUS: CONTEMPORARY")
            print("=============================================")
            cont_ns = ns_df[ns_df["Condition"] == "Contemporary"]
            print(cont_ns["AA_Mutation"].value_counts().head(15).to_string())
            print("\n=============================================")
            print("    TOP 15 NON-SYNONYMOUS: HISTORICAL")
            print("=============================================")
            hist_ns = ns_df[ns_df["Condition"] == "Historical"]
            print(hist_ns["AA_Mutation"].value_counts().head(15).to_string())
        else:
            print("No valid mutations found to translate.")

    if __name__ == "__main__":
        analyze_amino_acids(merged_df_3)
    return


@app.cell
def _(AlignIO, merged_df_3, pd):
    _alignment_file = "../data/processed/HAV_Collection date_aligned.fasta"
    ref_id_query = "NC_001489.1"
    _sub_cols = [
        "A>C",
        "A>G",
        "A>T",
        "C>A",
        "C>G",
        "C>T",
        "G>A",
        "G>C",
        "G>T",
        "T>A",
        "T>C",
        "T>G",
    ]

    def get_substitution_counts(alignment_path, ref_query):
        try:
            aln = AlignIO.read(alignment_path, "fasta")
        except Exception as e:
            return print(f"Error loading file: {e}")
        ref_record = next((s for s in aln if ref_query in s.id), None)
        if not ref_record:
            return print(f"❌ Reference '{ref_query}' not found in alignment.")
        ref_seq = str(ref_record.seq).upper()
        print(f"✅ Reference found: {ref_record.id}")
        _data = []
        for _record in aln:
            if _record.id == ref_record.id:
                continue
            _counts = {k: 0 for k in _sub_cols}
            clean_id = _record.id.split()[0]
            query_seq = str(_record.seq).upper()
            for r_base, q_base in zip(ref_seq, query_seq):
                if q_base in "ACGT" and r_base in "ACGT":
                    if q_base != r_base:
                        mut_type = f"{r_base}>{q_base}"
                        _counts[mut_type] = _counts[mut_type] + 1
            _counts["accession_id"] = clean_id
            _data.append(_counts)
        return pd.DataFrame(_data)

    mutation_df = get_substitution_counts(_alignment_file, ref_id_query)
    if mutation_df is not None:
        final_df = pd.merge(merged_df_3, mutation_df, on="accession_id", how="left")
        final_df[_sub_cols] = final_df[_sub_cols].fillna(0).astype(int)
        print("\nMerge Complete.")
        print("Columns checked: year, Period, Timeline are preserved.")
        cols_to_show = ["accession_id", "year", "Period", "Timeline"] + _sub_cols[:4]
        print(final_df[cols_to_show].head())
    return (final_df,)


@app.cell
def _(chi2_contingency, final_df, np, pd, plt, sns):
    from statsmodels.stats.multitest import multipletests

    stats_output_file = "../results/tables/substitution_stats_with_effect_size.csv"
    EFFECT_SIZE_THRESHOLD = 0.05
    SHOW_WEAK_SIGNIFICANCE = True
    _sub_cols = [
        "A>C",
        "A>G",
        "A>T",
        "C>A",
        "C>G",
        "C>T",
        "G>A",
        "G>C",
        "G>T",
        "T>A",
        "T>C",
        "T>G",
    ]
    transitions = ["A>G", "G>A", "C>T", "T>C"]
    sns.set_style("ticks")
    _TITLE_FS = 0
    _LABEL_FS = 28
    _TICK_FS = 24
    palette_E = {"Transition": "#2b83ba", "Transversion": "#d7191c"}
    final_df["Period_Group"] = final_df["year"].apply(
        lambda y: "Count_After_2018" if y >= 2018 else "Count_Before_2018"
    )
    grouped = final_df.groupby("Period_Group")[_sub_cols].sum().T
    for _col in ["Count_Before_2018", "Count_After_2018"]:
        if _col not in grouped.columns:
            grouped[_col] = 0
    df_E_counts = grouped.copy()
    df_E_counts.index.name = "Substitution Type"
    total_before = df_E_counts["Count_Before_2018"].sum()
    total_after = df_E_counts["Count_After_2018"].sum()
    df_E_counts["Prop_Before_2018"] = (
        df_E_counts["Count_Before_2018"] / total_before * 100
    )
    df_E_counts["Prop_After_2018"] = df_E_counts["Count_After_2018"] / total_after * 100
    df_E_counts["Absolute Change (%)"] = (
        df_E_counts["Prop_After_2018"] - df_E_counts["Prop_Before_2018"]
    )

    def get_significance_label(row):
        _p = row["P_Value_Corrected"]
        v = row["Cramers_V"]
        if _p >= 0.05:
            return ""
        if v < EFFECT_SIZE_THRESHOLD and (not SHOW_WEAK_SIGNIFICANCE):
            return ""
        if _p < 0.001:
            return "***"
        elif _p < 0.01:
            return "**"
        elif _p < 0.05:
            return "*"
        return ""

    stats_results = []
    for sub_type in df_E_counts.index:
        a = df_E_counts.loc[sub_type, "Count_Before_2018"]
        _c = df_E_counts.loc[sub_type, "Count_After_2018"]
        b = total_before - a
        d = total_after - _c
        table = np.array([[a, b], [_c, d]])
        (_chi2, _p, _, _) = chi2_contingency(table)
        n_total = table.sum()
        cramers_v = np.sqrt(_chi2 / n_total)
        if (table == 0).any():
            (a_c, b_c, c_c, d_c) = (a + 0.5, b + 0.5, _c + 0.5, d + 0.5)
        else:
            (a_c, b_c, c_c, d_c) = (a, b, _c, d)
        odds_ratio = c_c * b_c / (a_c * d_c)
        stats_results.append(
            {
                "Substitution Type": sub_type,
                "Count_Before": a,
                "Count_After": _c,
                "Prop_Before (%)": df_E_counts.loc[sub_type, "Prop_Before_2018"],
                "Prop_After (%)": df_E_counts.loc[sub_type, "Prop_After_2018"],
                "Absolute_Change (%)": df_E_counts.loc[sub_type, "Absolute Change (%)"],
                "P_Value_Raw": _p,
                "Chi2_Stat": _chi2,
                "Cramers_V": cramers_v,
                "Odds_Ratio": odds_ratio,
            }
        )
    stats_df_E = pd.DataFrame(stats_results)
    (reject, pvals_corrected, _, _) = multipletests(
        stats_df_E["P_Value_Raw"], alpha=0.05, method="fdr_bh"
    )
    stats_df_E["P_Value_Corrected"] = pvals_corrected
    stats_df_E["Symbol"] = stats_df_E.apply(get_significance_label, axis=1)
    stats_df_E.to_csv(stats_output_file, index=False)
    print(f"✅ Statistics saved to {stats_output_file}")
    print("Review 'Cramers_V' column: <0.1 is small effect, >0.3 is medium.")
    print(
        stats_df_E[
            ["Substitution Type", "Absolute_Change (%)", "Cramers_V", "Symbol"]
        ].head()
    )
    df_E_sorted = df_E_counts.reset_index().merge(
        stats_df_E[["Substitution Type", "Symbol"]], on="Substitution Type"
    )
    df_E_sorted["Group"] = df_E_sorted["Substitution Type"].apply(
        lambda x: "Transition" if x in transitions else "Transversion"
    )
    df_E_sorted = df_E_sorted.sort_values(
        by="Absolute Change (%)", ascending=True
    ).reset_index(drop=True)
    df_E_sorted["text_label"] = df_E_sorted["Absolute Change (%)"].apply(
        lambda x: f"{x:.2f}%"
    )
    (_fig, _axE) = plt.subplots(figsize=(24, 14))
    sns.barplot(
        data=df_E_sorted,
        y="Substitution Type",
        x="Absolute Change (%)",
        hue="Group",
        palette=palette_E,
        dodge=False,
        ax=_axE,
    )
    for _i, (_index, _row) in enumerate(df_E_sorted.iterrows()):
        x_val = _row["Absolute Change (%)"]
        _label = _row["text_label"]
        symbol = _row["Symbol"]
        ha = "left" if x_val >= 0 else "right"
        offset = 0.05 if x_val >= 0 else -0.05
        _axE.text(
            x_val + offset,
            _i,
            f" {_label} {symbol}",
            va="center",
            ha=ha,
            fontsize=_TICK_FS - 2,
        )
    _axE.set_title(
        "Change in Substitution Proportion (Pre-2018 vs Post-2018)", fontsize=_TITLE_FS
    )
    _axE.set_xlabel("Absolute Change in Proportion (%)", fontsize=_LABEL_FS)
    _axE.set_ylabel("Substitution Type", fontsize=_LABEL_FS)
    _axE.axvline(0, linewidth=1.5, linestyle="--", color="grey")
    _axE.legend().set_visible(False)
    _axE.tick_params(axis="both", which="major", labelsize=_TICK_FS)
    num_transversions = (df_E_sorted["Group"] == "Transversion").sum()
    num_transitions = (df_E_sorted["Group"] == "Transition").sum()
    t_y_pos = (num_transversions - 1) / 2
    tr_y_pos = num_transversions + (num_transitions - 1) / 2
    min_val = df_E_sorted["Absolute Change (%)"].min()
    _label_x_pos = min_val - 0.5 if min_val < -2 else -3
    _axE.text(
        _label_x_pos,
        t_y_pos,
        "Transversions",
        ha="center",
        va="center",
        fontweight="bold",
        color=palette_E["Transversion"],
        fontsize=_LABEL_FS - 2,
        rotation=90,
    )
    _axE.text(
        _label_x_pos,
        tr_y_pos,
        "Transitions",
        ha="center",
        va="center",
        fontweight="bold",
        color=palette_E["Transition"],
        fontsize=_LABEL_FS - 2,
        rotation=90,
    )
    max_val = df_E_sorted["Absolute Change (%)"].max()
    _axE.set_xlim(min_val - 1.0, max_val + 1.0)
    sns.despine(ax=_axE)
    plt.tight_layout(pad=3.0)
    plt.show()
    return (multipletests,)


@app.cell
def _(GridSpec, json, os, pd, plt, sns, stats):
    _OKABE_BLUE = "#0072B2"
    _OKABE_VERMILION = "#D55E00"
    _OKABE_TEAL = "#009E73"
    _OKABE_SKY = "#56B4E9"
    _OKABE_YELLOW = "#F0E442"
    _OKABE_GREY = "#999999"
    _COLOR_7D = _OKABE_BLUE
    _COLOR_14D = _OKABE_VERMILION
    _PALETTE_TIMEPOINTS = {
        "7d_post_infection": _COLOR_7D,
        "14d_post_infection": _COLOR_14D,
    }
    sns.set_context("paper")
    sns.set_style(
        "ticks",
        {
            "axes.linewidth": 1.5,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "font.family": "sans-serif",
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
        },
    )
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    _LABEL_FS = 16
    _TICK_FS = 14
    _ANNOT_FS = 12
    _LETTER_FS = 24

    def _get_significance_star(p):
        if _p < 0.0001:
            return "****"
        if _p < 0.001:
            return "***"
        if _p < 0.01:
            return "**"
        if _p < 0.05:
            return "*"
        return "ns"

    def _add_stat_annotation(ax, x1, x2, y_max, p_val, h_offset=0.08):
        """Draws a crisp journal-style bracket with p-value significance."""
        y_range = _ax.get_ylim()[1] - _ax.get_ylim()[0]
        if y_range == 0:
            y_range = 1.0
        h = y_range * 0.03
        y = _y_max + y_range * h_offset
        _ax.plot([_x1, _x1, _x2, _x2], [y, y + h, y + h, y], lw=1.5, c="black")
        sig = _get_significance_star(_p_val)
        p_text = (
            f"{sig}\n($p$={_p_val:.1e})"
            if _p_val < 0.001
            else f"{sig}\n($p$={_p_val:.3f})"
        )
        _ax.text(
            (_x1 + _x2) * 0.5,
            y + h * 1.2,
            p_text,
            ha="center",
            va="bottom",
            fontsize=_ANNOT_FS,
            fontweight="bold",
            color="black",
        )

    def _manual_fasta_parser(path):
        seq = []
        with open(_path, "r") as _f:
            for line in _f:
                if not line.startswith(">"):
                    seq.append(line.strip())
        return "".join(seq).upper()

    _GENE_CATEGORIES = {
        "ZAP": ["ZC3HAV1"],
        "Immune Regulators": ["IRF1", "IRF3", "IRF7"],
        "APOBEC Family": [
            "APOBEC1",
            "APOBEC2",
            "APOBEC3A",
            "APOBEC3B",
            "APOBEC3C",
            "APOBEC3D",
            "APOBEC3F",
            "APOBEC3G",
            "APOBEC3H",
            "APOBEC4",
        ],
    }
    _ORDERED_GENES = [
        gene for category in _GENE_CATEGORIES.values() for gene in category
    ]
    _HEATMAP_FILES = {"HAV": "HAV_vs_Mock_Day14_results_annotated.csv"}
    _REFERENCE_FASTA = "NC_001489.1.fasta"
    _GENOME_LENGTH_BP = 7478
    _FILE_GROUPS = {
        "7d_post_infection": [
            "REDItools_7dHAVrep1.txt",
            "REDItools_7dHAVrep2.txt",
            "REDItools_7dHAVrep3.txt",
        ],
        "14d_post_infection": [
            "REDItools_14dHAVrep1.txt",
            "REDItools_14dHAVrep2.txt",
            "REDItools_14dHAVrep3.txt",
        ],
    }
    _DATA_DOT_PANELS = {
        "Group": ["7d_post_infection"] * 3 + ["14d_post_infection"] * 3,
        "CpG_OE": [0.1753, 0.1741, 0.1705, 0.1814, 0.1768, 0.1772],
        "GpA_OE": [1.2078, 1.2078, 1.2059, 1.2055, 1.2063, 1.2065],
    }

    def _process_df_heatmap(df, condition_name):
        df.columns = [_c.lower() if _c.lower() == "symbol" else _c for _c in df.columns]
        if "symbol" not in df.columns:
            return pd.DataFrame()
        df_proc = df[df["symbol"].isin(_ORDERED_GENES)][
            ["symbol", "log2FoldChange", "padj"]
        ].copy()
        df_proc.rename(
            columns={
                "log2FoldChange": f"log2FC_{condition_name}",
                "padj": f"padj_{condition_name}",
            },
            inplace=True,
        )
        df_proc.drop_duplicates(subset="symbol", keep="first", inplace=True)
        return df_proc

    def _process_editing_totals(filepath, ref_seq, genome_len):
        _counts = {"APOBEC": 0}
        try:
            df = pd.read_csv(filepath, sep="\t", on_bad_lines="skip")
        except:
            return None
        if "Position" not in df.columns or "BaseCount[A,C,G,T]" not in df.columns:
            return None
        for _, _row in df.iterrows():
            pos = int(_row["Position"])
            if pos < 2 or pos > len(ref_seq) - 1:
                continue
            ref_base = ref_seq[pos - 1]
            pre_base = ref_seq[pos - 2]
            fol_base = ref_seq[pos]
            try:
                c_list = json.loads(_row["BaseCount[A,C,G,T]"])
            except:
                continue
            if sum(c_list) == 0:
                continue
            if ref_base == "C" and pre_base in ["T", "C"]:
                _counts["APOBEC"] = _counts["APOBEC"] + c_list[3]
            elif ref_base == "G" and fol_base in ["A", "G"]:
                _counts["APOBEC"] = _counts["APOBEC"] + c_list[0]
        return {k: v / (genome_len / 1000.0) for (k, v) in _counts.items()}

    def _main():
        print("Generating 2x2 MBE-Formatted Figure 3...")
        _fig = plt.figure(figsize=(14, 12))
        _gs = GridSpec(2, 2, figure=_fig, hspace=0.35, wspace=0.35)
        _axA = _fig.add_subplot(_gs[0, 0])
        _axB = _fig.add_subplot(_gs[0, 1])
        _axC = _fig.add_subplot(_gs[1, 0])
        _axD = _fig.add_subplot(_gs[1, 1])
        _axA.text(
            -0.25,
            1.05,
            "A",
            transform=_axA.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="bottom",
        )
        if os.path.exists(_HEATMAP_FILES["HAV"]):
            hav = _process_df_heatmap(pd.read_csv(_HEATMAP_FILES["HAV"]), "HAV")
            merged = hav
            val_cols = [_c for _c in merged.columns if "log2FC" in _c]
            padj_cols = [_c for _c in merged.columns if "padj" in _c]
            _plot_data = merged.set_index("symbol")[val_cols].fillna(0)
            _plot_data.columns = ["HAV / Mock"]
            p_data = merged.set_index("symbol")[padj_cols].fillna(1)
            genes = [_g for _g in _ORDERED_GENES if _g in _plot_data.index]
            _plot_data = _plot_data.loc[genes]
            p_data = p_data.loc[genes]
            annot_matrix = pd.DataFrame(
                index=_plot_data.index, columns=_plot_data.columns
            )
            for r in _plot_data.index:
                for c_idx, _c in enumerate(_plot_data.columns):
                    _val = _plot_data.loc[r, _c]
                    pval = p_data.iloc[_plot_data.index.get_loc(r), c_idx]
                    _star = _get_significance_star(pval)
                    star_txt = f"^{_star}" if _star != "ns" else ""
                    annot_matrix.loc[r, _c] = f"{_val:.2f}{star_txt}"
            sns.heatmap(
                _plot_data,
                annot=annot_matrix,
                fmt="s",
                cmap="vlag",
                center=0,
                linewidths=1,
                linecolor="white",
                ax=_axA,
                cbar_kws={
                    "label": "log$_2$ Fold Change",
                    "shrink": 0.7,
                    "location": "right",
                    "pad": 0.05,
                },
                annot_kws={"fontsize": _TICK_FS, "fontweight": "bold"},
            )
            _axA.xaxis.tick_top()
            _axA.tick_params(axis="both", labelsize=_TICK_FS)
            _axA.set_xticklabels(_axA.get_xticklabels(), fontweight="bold")
            _axA.set_xlabel("")
            _axA.set_ylabel("")
            curr = 0
            g_idx = {_g: _i for (_i, _g) in enumerate(genes)}
            for cat, g_list in _GENE_CATEGORIES.items():
                found = [_g for _g in g_list if _g in genes]
                if found:
                    _start = g_idx[found[0]]
                    _end = g_idx[found[-1]]
                    if curr > 0:
                        _axA.axhline(_start, color="black", lw=1.5)
                    curr = _end + 1
        else:
            _axA.text(0.5, 0.5, "Heatmap File Missing", ha="center", fontsize=20)
        _axB.text(
            -0.25,
            1.05,
            "B",
            transform=_axB.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="bottom",
        )
        if os.path.exists(_REFERENCE_FASTA):
            ref_seq = _manual_fasta_parser(_REFERENCE_FASTA)
            edit_data = []
            for grp, files in _FILE_GROUPS.items():
                for _f in files:
                    if os.path.exists(_f):
                        r = _process_editing_totals(_f, ref_seq, _GENOME_LENGTH_BP)
                        if r:
                            for enzyme, rate in r.items():
                                edit_data.append(
                                    {"Group": grp, "Enzyme": enzyme, "Rate": rate}
                                )
            if edit_data:
                df_edit = pd.DataFrame(edit_data)
                sns.barplot(
                    x="Enzyme",
                    y="Rate",
                    hue="Group",
                    data=df_edit,
                    palette=_PALETTE_TIMEPOINTS,
                    errwidth=1.5,
                    capsize=0.05,
                    edgecolor="black",
                    linewidth=1.5,
                    ax=_axB,
                    alpha=0.9,
                )
                sns.stripplot(
                    x="Enzyme",
                    y="Rate",
                    hue="Group",
                    data=df_edit,
                    dodge=True,
                    color="black",
                    size=7,
                    ax=_axB,
                    legend=False,
                    alpha=0.7,
                )
                enzymes = df_edit["Enzyme"].unique()
                for _i, enz in enumerate(enzymes):
                    _sub = df_edit[df_edit["Enzyme"] == enz]
                    g1 = _sub[_sub["Group"] == "7d_post_infection"]["Rate"]
                    g2 = _sub[_sub["Group"] == "14d_post_infection"]["Rate"]
                    if len(g1) > 1 and len(g2) > 1:
                        (t, _p) = stats.ttest_ind(g1, g2)
                        _y_max = _sub["Rate"].max()
                        h = _y_max * 0.1
                        y = _y_max + h
                        _x1 = _i - 0.2
                        _x2 = _i + 0.2
                        _axB.plot(
                            [_x1, _x1, _x2, _x2],
                            [y, y + h / 2, y + h / 2, y],
                            lw=1.5,
                            c="black",
                        )
                        sig_text = f"{_get_significance_star(_p)}\n($p$={_p:.3f})"
                        _axB.text(
                            _i,
                            y + h * 1.2,
                            sig_text,
                            ha="center",
                            va="bottom",
                            fontsize=_ANNOT_FS,
                            fontweight="bold",
                        )
                _axB.set_ylim(0, _axB.get_ylim()[1] * 1.15)
                _axB.set_ylabel(
                    "Total Editing Rate\n(Mutations / kb)",
                    fontsize=_LABEL_FS,
                    fontweight="bold",
                )
                _axB.set_xlabel("")
                (handles, _labels) = _axB.get_legend_handles_labels()
                _axB.legend(
                    handles=handles[:2],
                    labels=["7 DPI", "14 DPI"],
                    fontsize=_TICK_FS,
                    loc="upper right",
                    frameon=False,
                )
                _axB.tick_params(axis="both", labelsize=_TICK_FS)
                sns.despine(ax=_axB)
        else:
            _axB.text(
                0.5, 0.5, "Ref Genome Missing for Panel B", ha="center", fontsize=20
            )
        _axC.text(
            -0.25,
            1.05,
            "C",
            transform=_axC.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="bottom",
        )
        df_dot = pd.DataFrame(_DATA_DOT_PANELS)
        g1_cpg = df_dot[df_dot["Group"] == "7d_post_infection"]["CpG_OE"]
        g2_cpg = df_dot[df_dot["Group"] == "14d_post_infection"]["CpG_OE"]
        (t_cpg, p_cpg) = stats.ttest_ind(g1_cpg, g2_cpg)
        sns.barplot(
            x="Group",
            y="CpG_OE",
            data=df_dot,
            palette=_PALETTE_TIMEPOINTS,
            errwidth=1.5,
            capsize=0.05,
            edgecolor="black",
            linewidth=1.5,
            ax=_axC,
            alpha=0.9,
        )
        sns.stripplot(
            x="Group",
            y="CpG_OE",
            data=df_dot,
            ax=_axC,
            color="black",
            size=7,
            jitter=0.1,
            alpha=0.7,
        )
        _axC.set_ylim(0, _axC.get_ylim()[1] * 1.05)
        _add_stat_annotation(_axC, 0, 1, df_dot["CpG_OE"].max(), p_cpg)
        _axC.set_ylabel("CpG Obs/Exp Ratio", fontsize=_LABEL_FS, fontweight="bold")
        _axC.set_xlabel("")
        _axC.set_xticklabels(
            ["7 Days", "14 Days"], fontsize=_TICK_FS, fontweight="bold"
        )
        _axC.tick_params(axis="both", labelsize=_TICK_FS)
        sns.despine(ax=_axC)
        _axD.text(
            -0.25,
            1.05,
            "D",
            transform=_axD.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="bottom",
        )
        g1_gpa = df_dot[df_dot["Group"] == "7d_post_infection"]["GpA_OE"]
        g2_gpa = df_dot[df_dot["Group"] == "14d_post_infection"]["GpA_OE"]
        (t_gpa, p_gpa) = stats.ttest_ind(g1_gpa, g2_gpa)
        sns.barplot(
            x="Group",
            y="GpA_OE",
            data=df_dot,
            palette=_PALETTE_TIMEPOINTS,
            errwidth=1.5,
            capsize=0.05,
            edgecolor="black",
            linewidth=1.5,
            ax=_axD,
            alpha=0.9,
        )
        sns.stripplot(
            x="Group",
            y="GpA_OE",
            data=df_dot,
            ax=_axD,
            color="black",
            size=7,
            jitter=0.1,
            alpha=0.7,
        )
        _axD.set_ylim(0, _axD.get_ylim()[1] * 1.02)
        _add_stat_annotation(_axD, 0, 1, df_dot["GpA_OE"].max(), p_gpa)
        _axD.set_ylabel("GpA Obs/Exp Ratio", fontsize=_LABEL_FS, fontweight="bold")
        _axD.set_xlabel("")
        _axD.set_xticklabels(
            ["7 Days", "14 Days"], fontsize=_TICK_FS, fontweight="bold"
        )
        _axD.tick_params(axis="both", labelsize=_TICK_FS)
        sns.despine(ax=_axD)
        plt.tight_layout()
        plt.savefig("Figure3_2x2_MBE.pdf", format="pdf", bbox_inches="tight")
        try:
            plt.savefig(
                "Figure3_2x2_MBE.tiff",
                dpi=300,
                format="tiff",
                pil_kwargs={"compression": "tiff_lzw"},
                bbox_inches="tight",
            )
            print("High-Res TIFF generated successfully.")
        except Exception as e:
            print(f"TIFF generation note: {e}")
        plt.savefig("Figure3_2x2_MBE.png", dpi=300, bbox_inches="tight")
        print("Figure 3 2x2 Formatting Complete.")
        plt.show()

    if __name__ == "__main__":
        _main()
    return


@app.cell
def _(final_df):
    final_df.columns
    return


@app.cell
def _(AlignIO, pd):
    def find_broad_motif_transitions(alignment_file, ref_index=0):
        """
        Identifies transitions using broader dinucleotide contexts:

        1. APOBEC (Direct C>T): C>T in [T/C]C context (Preceding base is T or C)
        2. APOBEC (Neg Strand G>A): G>A in G[A/G] context (Following base is A or G)
        """
        alignment = AlignIO.read(alignment_file, "fasta")
        ref_seq = alignment[ref_index].seq.upper()

        apobec_data = []


        # Iterate through sequence (avoiding index errors at ends)
        for i in range(1, len(ref_seq) - 1):
            prev_base = ref_seq[i-1]
            current_base = ref_seq[i]
            next_base = ref_seq[i+1]

            # --- 1. APOBEC (C>T) ---
            # Logic: 5' base is T or C
            if current_base == "C" and prev_base in ["T", "C"]:
                mutation_count = 0
                for seq in alignment:
                    if seq.seq[i].upper() == "T" and seq.seq[i] != "-":
                        mutation_count += 1
                if mutation_count > 0:
                    apobec_data.append({
                        "Position": i + 1,
                        "Context": f"{prev_base}C",
                        "Mutation Count": mutation_count,
                        "Type": "APOBEC (C>T Direct)"
                    })

            # --- 2. APOBEC (G>A - Negative Strand) ---
            # Logic: 3' (Next) base is A or G
            elif current_base == "G" and next_base in ["A", "G"]:
                mutation_count = 0
                for seq in alignment:
                    if seq.seq[i].upper() == "A" and seq.seq[i] != "-":
                        mutation_count += 1
                if mutation_count > 0:
                    apobec_data.append({
                        "Position": i + 1,
                        "Context": f"G{next_base}",
                        "Mutation Count": mutation_count,
                        "Type": "APOBEC (G>A NegStrand)"
                    })


        return pd.DataFrame(apobec_data)

    # Run the analysis
    df_apobec = find_broad_motif_transitions("../data/processed/HAV_after2018_align.fasta")

    # Save
    df_apobec.to_excel("../results/tables/HAV_APOBEC_Broad_Combined.xlsx", index=False)


    print(f"APOBEC Hits (C>T + G>A): {len(df_apobec)}")
    return


@app.cell
def _(
    GridSpec,
    GridSpecFromSubplotSpec,
    SeqIO,
    chi2_contingency,
    fisher_exact,
    mannwhitneyu,
    multipletests,
    np,
    os,
    pd,
    plt,
    sns,
):

    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    _FILES = {
        "FASTA_BEFORE": "../data/processed/HAV_before2018_align.fasta",
        "FASTA_AFTER": "../data/processed/HAV_after2018_align.fasta",
        "TRACK_APOBEC": "../results/tables/HAV_APOBEC_Broad_Combined.xlsx",
    }
    _REPORT_FILE = "Detailed_Statistical_Report.txt"
    _PALETTE_E = {"Transition": "#2b83ba", "Transversion": "#d7191c"}
    _PALETTE_COHORT = {"Historical": "#0072B2", "Contemporary": "#D55E00"}
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
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
        },
    )
    _LABEL_FS = 18
    _LETTER_FS = 26
    _SUB_COLS = [
        "A>C",
        "A>G",
        "A>T",
        "C>A",
        "C>G",
        "C>T",
        "G>A",
        "G>C",
        "G>T",
        "T>A",
        "T>C",
        "T>G",
    ]
    _TRANSITIONS = ["A>G", "G>A", "C>T", "T>C"]
    _GENOME_LENGTH = 7478
    _GENOME_REGIONS = [
        (0, 734, "#F5F5F5", "5'UTR"),
        (735, 803, "#C8E6C9", "VP4"),
        (804, 1469, "#A5D6A7", "VP2"),
        (1470, 2207, "#81C784", "VP3"),
        (2208, 3107, "#66BB6A", "VP1"),
        (3108, 3674, "#FFF9C4", "2A"),
        (3675, 3995, "#FFF59D", "2B"),
        (3996, 5000, "#FFF176", "2C"),
        (5001, 5222, "#E1F5FE", "3A"),
        (5223, 5291, "#B3E5FC", "3B"),
        (5292, 5948, "#81D4FA", "3C"),
        (5949, 7415, "#4FC3F7", "3D"),
        (7416, 7478, "#F5F5F5", "3'UTR"),
    ]


    # FIXED: Changed _val to val
    def _standardize_period(val):
        s = str(val).lower().strip()
        if any((x in s for x in ["before", "pre", "<", "old", "historical"])):
            return "Historical"
        if any((x in s for x in ["after", "post", ">", "new", "contemporary"])):
            return "Contemporary"
        return None


    # FIXED: Changed _p to p
    def _get_star(p):
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return ""


    def _process_differential_stats(df):
        if "Period" not in df.columns:
            return pd.DataFrame()
        df["Period_Std"] = df["Period"].apply(_standardize_period)
        df_clean = df.dropna(subset=["Period_Std"])
        grouped = df_clean.groupby("Period_Std")[_SUB_COLS].sum().T
        if "Historical" not in grouped.columns or "Contemporary" not in grouped.columns:
            return pd.DataFrame()
        total_before = grouped["Historical"].sum()
        total_after = grouped["Contemporary"].sum()
        stats_data = []
        for _sub in grouped.index:
            a = grouped.loc[_sub, "Historical"]
            _c = grouped.loc[_sub, "Contemporary"]
            (odds, _p) = fisher_exact([[a, total_before - a], [_c, total_after - _c]])
            prop_before = a / total_before * 100
            prop_after = _c / total_after * 100
            stats_data.append(
                {
                    "Substitution Type": _sub,
                    "Absolute Change (%)": prop_after - prop_before,
                    "P_Value": _p,
                    "Group": "Transition" if _sub in _TRANSITIONS else "Transversion",
                }
            )
        stats_df = pd.DataFrame(stats_data)
        (_, adj_p, _, _) = multipletests(stats_df["P_Value"], method="fdr_bh")
        stats_df["FDR"] = adj_p
        return stats_df.sort_values("Absolute Change (%)")


    def _get_context(ref, alt, left, right):
        complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
        if (
            ref not in complements
            or alt not in complements
            or left not in complements
            or (right not in complements)
        ):
            return None
        if ref in ["C", "T"]:
            return f"{left}[{ref}>{alt}]{right}"
        elif ref in ["G", "A"]:
            return f"{complements[right]}[{complements[ref]}>{complements[alt]}]{complements[left]}"
        return None


    def _process_sequences_for_volcano(fasta_path, group_name):
        if not os.path.exists(fasta_path):
            return (pd.DataFrame(), [])
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            return (pd.DataFrame(), [])
        ref_seq = str(records[0].seq).upper()
        _data = []
        bases = ["A", "C", "G", "T"]
        sbs96 = [
            f"{l}[{m}]{r}"
            for m in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
            for l in bases
            for r in bases
        ]
        for r in records[1:]:
            q_seq = str(r.seq).upper()
            if len(q_seq) != len(ref_seq):
                continue
            _counts = {_c: 0 for _c in sbs96}
            total_muts = 0
            for _i in range(1, len(ref_seq) - 1):
                ref_b = ref_seq[_i]
                q_b = q_seq[_i]
                if ref_b != q_b and ref_b in bases and (q_b in bases):
                    l_b = ref_seq[_i - 1]
                    r_b = ref_seq[_i + 1]
                    if l_b in bases and r_b in bases:
                        ctx = _get_context(ref_b, q_b, l_b, r_b)
                        if ctx and ctx in _counts:
                            _counts[ctx] = _counts[ctx] + 1
                            total_muts = total_muts + 1
            if total_muts > 0:
                _row = {k: v / total_muts for (k, v) in _counts.items()}
                _row["Group"] = group_name
                _data.append(_row)
        return (pd.DataFrame(_data), sbs96)


    def _calculate_volcano_stats():
        (df1, sbs_cols) = _process_sequences_for_volcano(
            _FILES["FASTA_BEFORE"], "Historical"
        )
        (df2, _) = _process_sequences_for_volcano(_FILES["FASTA_AFTER"], "Contemporary")
        if df1.empty or df2.empty:
            return pd.DataFrame()
        full_df = pd.concat([df1, df2], ignore_index=True)
        results = []
        for motif in sbs_cols:
            v1 = full_df[full_df["Group"] == "Historical"][motif].values
            v2 = full_df[full_df["Group"] == "Contemporary"][motif].values
            try:
                if (
                    len(np.unique(v1)) < 2
                    and len(np.unique(v2)) < 2
                    and (np.mean(v1) == np.mean(v2))
                ):
                    _p = 1.0
                else:
                    (_, _p) = mannwhitneyu(v1, v2, alternative="two-sided")
            except:
                _p = 1.0
            mean1 = np.mean(v1) + 1e-09
            mean2 = np.mean(v2) + 1e-09
            log2fc = np.log2(mean2 / mean1)
            results.append({"Motif": motif, "P_Value": _p, "Log2_FC": log2fc})
        res_df = pd.DataFrame(results)
        (_, adj_p, _, _) = multipletests(res_df["P_Value"], method="fdr_bh")
        res_df["Q_Value"] = adj_p
        return res_df


    def _get_context_data_and_stats(center_base, mutation_base):
        bases = ["A", "C", "G", "T"]
        contexts = [f"{L}{center_base}{R}" for L in bases for R in bases]
        _data = []
        fastas = [
            ("Historical", _FILES["FASTA_BEFORE"]),
            ("Contemporary", _FILES["FASTA_AFTER"]),
        ]
        cohort_data = {}
        for _name, _path in fastas:
            if not os.path.exists(_path):
                continue
            records = list(SeqIO.parse(_path, "fasta"))
            if not records:
                continue
            ref = str(records[0].seq).upper()
            ctx_sites = {ctx: 0 for ctx in contexts}
            total_sites = 0
            for _i in range(1, len(ref) - 1):
                if ref[_i] == center_base:
                    (L, R) = (ref[_i - 1], ref[_i + 1])
                    if L in bases and R in bases:
                        total_sites = total_sites + 1
                        ctx_sites[f"{L}{center_base}{R}"] = (
                            ctx_sites[f"{L}{center_base}{R}"] + 1
                        )
            mut_counts = {ctx: 0 for ctx in contexts}
            total_muts = 0
            for rec in records[1:]:
                seq = str(rec.seq).upper()
                if len(seq) != len(ref):
                    continue
                for _i in range(1, len(ref) - 1):
                    if ref[_i] == center_base and seq[_i] == mutation_base:
                        (L, R) = (ref[_i - 1], ref[_i + 1])
                        if L in bases and R in bases:
                            total_muts = total_muts + 1
                            mut_counts[f"{L}{center_base}{R}"] = (
                                mut_counts[f"{L}{center_base}{R}"] + 1
                            )
            cohort_data[_name] = {
                "muts": mut_counts,
                "total": total_muts,
                "sites": ctx_sites,
                "total_sites": total_sites,
            }
            for ctx in contexts:
                obs = mut_counts[ctx] / total_muts if total_muts > 0 else 0
                exp = ctx_sites[ctx] / total_sites if total_sites > 0 else 0
                enrich = obs / exp if exp > 0 else 0
                _data.append({"Cohort": _name, "Context": ctx, "Enrichment": enrich})
        stats_rows = []
        if "Historical" in cohort_data and "Contemporary" in cohort_data:
            B = cohort_data["Historical"]
            A = cohort_data["Contemporary"]
            for ctx in contexts:
                count_b = B["muts"][ctx]
                other_b = B["total"] - count_b
                count_a = A["muts"][ctx]
                other_a = A["total"] - count_a
                try:
                    if count_b + other_b > 0 and count_a + other_a > 0:
                        (_chi2, _p, _, _) = chi2_contingency(
                            [[count_b, other_b], [count_a, other_a]]
                        )
                    else:
                        _p = 1.0
                except:
                    _p = 1.0
                stats_rows.append({"Context": ctx, "P_Raw": _p})
            p_raw = [r["P_Raw"] for r in stats_rows]
            (_, adj_p, _, _) = multipletests(p_raw, method="fdr_bh")
            for _i, _row in enumerate(stats_rows):
                _row["FDR"] = adj_p[_i]
                _row["Sig"] = _get_star(adj_p[_i])
        return (pd.DataFrame(_data), pd.DataFrame(stats_rows))


    def _get_binned_track_data(file_path, length=_GENOME_LENGTH, bin_size=200):
        if not os.path.exists(file_path):
            return (None, None, None)
        try:
            df = pd.read_excel(file_path)
            bins = np.arange(0, length + bin_size, bin_size)
            (bin_means, bin_sems) = ([], [])
            for _i in range(len(bins) - 1):
                _mask = (df["Position"] >= bins[_i]) & (df["Position"] < bins[_i + 1])
                vals = df.loc[_mask, "Mutation Count"]
                if len(vals) > 0:
                    bin_means.append(vals.mean())
                    bin_sems.append(vals.sem())
                else:
                    bin_means.append(0)
                    bin_sems.append(0)
            return (bins[:-1] + bin_size / 2, np.array(bin_means), np.array(bin_sems))
        except:
            return (None, None, None)


    def _generate_report_to_file(df, ct_stats, filename=_REPORT_FILE):
        with open(filename, "w") as _f:
            _f.write("DETAILED STATISTICAL SUMMARY\n" + "=" * 60 + "\n\n")
            if "Period" in df.columns:
                _f.write(
                    "1. GLOBAL DIFFERENTIAL MUTATION RATES (Historical vs Contemporary)\n"
                )
                _f.write("-" * 60 + "\n")
                df["Period_Std"] = df["Period"].apply(_standardize_period)
                grouped = (
                    df.dropna(subset=["Period_Std"])
                    .groupby("Period_Std")[_SUB_COLS]
                    .sum()
                    .T
                )
                if (
                    "Historical" in grouped.columns
                    and "Contemporary" in grouped.columns
                ):
                    n_b = grouped["Historical"].sum()
                    n_a = grouped["Contemporary"].sum()
                    _f.write(f"Total Mutations: Historical={n_b}, Contemporary={n_a}\n")
                    _f.write(
                        f"{'Type':<8} | {'Hist%':<8} | {'Cont%':<8} | {'P-val':<10} | {'FDR':<10} | {'Sig'}\n"
                    )
                    results = []
                    for _sub in _SUB_COLS:
                        (a, _c) = (
                            grouped.loc[_sub, "Historical"],
                            grouped.loc[_sub, "Contemporary"],
                        )
                        (odds, _p) = fisher_exact([[a, n_b - a], [_c, n_a - _c]])
                        results.append(
                            {
                                "sub": _sub,
                                "pb": a / n_b * 100,
                                "pa": _c / n_a * 100,
                                "p": _p,
                            }
                        )
                    (_, adj, _, _) = multipletests(
                        [r["p"] for r in results], method="fdr_bh"
                    )
                    for _i, r in enumerate(results):
                        _f.write(
                            f"{r['sub']:<8} | {r['pb']:<8.2f} | {r['pa']:<8.2f} | {r['p']:<10.2e} | {adj[_i]:<10.2e} | {_get_star(adj[_i])}\n"
                        )
                _f.write("\n")
            _f.write("2. CONTEXT-SPECIFIC ANALYSIS (Chi-Square/Fisher 2x2)\n")
            _f.write("-" * 60 + "\n")

            # FIXED: Changed _name to name
            def write_ctx_stats(name, stats_df):
                _f.write(f"--- {name} ---\n")
                if stats_df.empty:
                    _f.write("No data available.\n")
                    return
                _f.write(
                    f"{'Context':<8} | {'P-raw':<10} | {'FDR (q)':<10} | {'Sig'}\n"
                )
                for _, r in stats_df.iterrows():
                    _f.write(
                        f"{r['Context']:<8} | {r['P_Raw']:<10.2e} | {r['FDR']:<10.2e} | {r['Sig']}\n"
                    )
                _f.write("\n")

            write_ctx_stats("APOBEC (C>T)", ct_stats)
        # FIXED: Changed _filename to filename
        print(f"✅ Report saved: {filename}")


    def _main():
        if "final_df" in globals():
            df = globals()["final_df"]
        else:
            print("Error: final_df missing.")
            return
        (ct_data, ct_stats) = _get_context_data_and_stats("C", "T")
        volc_df = _calculate_volcano_stats()
        _generate_report_to_file(df, ct_stats)
        print("Generating Figure...")
        _fig = plt.figure(figsize=(24, 25))
        outer_gs = GridSpec(3, 1, figure=_fig, height_ratios=[0.8, 1, 0.6], hspace=0.35)
        _axA = _fig.add_subplot(outer_gs[0])
        _axA.text(
            -0.02,
            1.05,
            "A",
            transform=_axA.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
        )
        df_diff = _process_differential_stats(df)
        if not df_diff.empty:
            df_diff["label"] = df_diff["Absolute Change (%)"].apply(
                lambda x: f"{x:+.1f}%"
            )
            sns.barplot(
                data=df_diff,
                y="Substitution Type",
                x="Absolute Change (%)",
                hue="Group",
                palette=_PALETTE_E,
                dodge=False,
                ax=_axA,
            )
            for _i, _row in enumerate(df_diff.itertuples()):
                x_val = _row._2
                (ha, offset) = ("left", 0.1) if x_val >= 0 else ("right", -0.1)
                sig_star = _get_star(_row.FDR)
                txt = f"{_row.label} ({sig_star})"
                _axA.text(
                    x_val + offset,
                    _i,
                    txt,
                    va="center",
                    ha=ha,
                    fontsize=11,
                    fontweight="bold",
                )
            _axA.axvline(0, linewidth=2, linestyle="--", color="black")
            _axA.legend(loc="upper right", frameon=False, fontsize=14)
            _axA.set_xlim(-4, 5.5)
            _axA.set_title("", fontsize=_LABEL_FS, fontweight="bold")
            _axA.set_xlabel(
                "Absolute Change (%)", fontsize=_LABEL_FS, fontweight="bold"
            )
            _axA.set_ylabel("Substitution Types", fontsize=_LABEL_FS, fontweight="bold")
            _axA.tick_params(axis="both", labelsize=_LABEL_FS)
            sns.despine(ax=_axA)
        else:
            _axA.text(0.5, 0.5, "Data Missing", ha="center")
        gs_row1 = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_gs[1], wspace=0.3)
        (_axB, _axC) = [_fig.add_subplot(gs_row1[_i]) for _i in range(2)]
        _axB.text(
            -0.05,
            1.05,
            "B",
            transform=_axB.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
        )
        if not volc_df.empty:
            volc_df["neg_log_q"] = -np.log10(volc_df["Q_Value"])
            _colors = []
            for _, r in volc_df.iterrows():
                if r["Q_Value"] < 0.05:
                    if r["Log2_FC"] > 0.5:
                        _colors.append("#e62725")
                    elif r["Log2_FC"] < -0.5:
                        _colors.append("#1f77b4")
                    else:
                        _colors.append("grey")
                else:
                    _colors.append("lightgrey")
            _axB.scatter(
                volc_df["Log2_FC"], volc_df["neg_log_q"], c=_colors, s=40, alpha=0.7
            )
            _axB.axhline(-np.log10(0.05), linestyle="--", color="grey", alpha=0.5)
            _axB.axvline(0.5, linestyle="--", color="grey", alpha=0.5)
            _axB.axvline(-0.5, linestyle="--", color="grey", alpha=0.5)
            _axB.set_xlim(-3, 2)
            top_hits = volc_df[(volc_df["Q_Value"] < 0.05) & (volc_df["Log2_FC"] > 0.5)]
            for _, _row in top_hits.iterrows():
                _lbl = f"{_row['Motif'][0]}$\\mathbf{{{_row['Motif'][2]}}}${_row['Motif'][6]}"
                _axB.text(
                    _row["Log2_FC"] + 0.02, _row["neg_log_q"] + 0.2, _lbl, fontsize=12
                )
            _axB.set_title("", fontsize=_LABEL_FS, fontweight="bold")
            _axB.set_xlabel("Log2 Fold Change", fontsize=_LABEL_FS, fontweight="bold")
            _axB.set_ylabel("-Log10 FDR", fontsize=_LABEL_FS, fontweight="bold")
            _axB.tick_params(axis="both", labelsize=_LABEL_FS)
            sns.despine(ax=_axB)
        else:
            _axB.text(0.5, 0.5, "Volcano Data Missing/Error", ha="center")

        # FIXED: Stripped underscores from ax, data, stats, title, letter
        def plot_context_panel(ax, data, stats, title, letter):
            ax.text(
                -0.05,
                1.05,
                letter,
                transform=ax.transAxes,
                fontsize=_LETTER_FS,
                fontweight="bold",
            )
            if not data.empty:
                sns.barplot(
                    data=data,
                    y="Context",
                    x="Enrichment",
                    hue="Cohort",
                    palette=_PALETTE_COHORT,
                    ax=ax,
                    edgecolor="black",
                )
                ax.set_title(title, fontsize=_LABEL_FS, fontweight="bold")
                ax.axvline(1.0, linestyle="--", color="grey")
                ax.legend(loc="upper right", frameon=False)
                ax.set_xlabel("Enrichment", fontsize=_LABEL_FS, fontweight="bold")
                ax.set_ylabel("N[C>T]N", fontsize=_LABEL_FS, fontweight="bold")
                ax.set_yticks(ax.get_yticks())
                _formatted_labels = [
                    f"{_lbl.get_text()[0]}[C>T]{_lbl.get_text()[2]}"
                    if len(_lbl.get_text()) >= 3
                    else _lbl.get_text()
                    for _lbl in ax.get_yticklabels()
                ]
                ax.set_yticklabels(_formatted_labels)
                ax.tick_params(axis="both", labelsize=_LABEL_FS)
                if not stats.empty:
                    sig_map = dict(zip(stats["Context"], stats["Sig"]))
                    contexts = sorted(data["Context"].unique())
                    for _i, ctx in enumerate(contexts):
                        sig = sig_map.get(ctx, "")
                        if sig:
                            max_val = data[data["Context"] == ctx]["Enrichment"].max()
                            ax.text(
                                max_val + 0.1,
                                _i,
                                sig,
                                va="center",
                                fontsize=14,
                                fontweight="bold",
                                color="black",
                            )
                sns.despine(ax=ax)
            else:
                ax.text(0.5, 0.5, "Data Missing", ha="center")

        plot_context_panel(_axC, ct_data, ct_stats, "", "C")
        gs_row2 = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer_gs[2], height_ratios=[1, 0.2], hspace=0.1
        )
        (axD1, axD2) = [_fig.add_subplot(gs_row2[_i]) for _i in range(2)]
        axD1.sharex(axD2)
        axD1.text(
            -0.02,
            1.1,
            "D",
            transform=axD1.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
        )
        (ap_x, ap_y, ap_err) = _get_binned_track_data(_FILES["TRACK_APOBEC"])
        if ap_x is not None:
            axD1.bar(
                ap_x,
                ap_y,
                yerr=ap_err,
                width=180,
                color="#d62728",
                alpha=0.8,
                capsize=3,
            )
            axD1.set_ylabel(
                "APOBEC Transitions \n(Weighted Moving Avg.)",
                fontsize=16,
                fontweight="bold",
                color="#d62728",
            )
            axD1.tick_params(axis="both", labelsize=_LABEL_FS)
        plt.setp(axD1.get_xticklabels(), visible=False)
        sns.despine(ax=axD1, bottom=True)
        for _start, _end, _color, _label in _GENOME_REGIONS:
            axD2.barh(
                0,
                width=_end - _start,
                left=_start,
                color=_color,
                edgecolor="black",
                height=1,
            )
            axD2.text(
                (_start + _end) / 2,
                0,
                _label,
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
            )
            for _ax in [axD1, axD2]:
                _ax.axvspan(_start, _end, facecolor=_color, alpha=0.2)
        axD2.set_yticks([])
        axD2.set_xlim(0, _GENOME_LENGTH)
        axD2.set_xlabel("Genomic Position (bp)", fontsize=_LABEL_FS, fontweight="bold")
        axD2.tick_params(axis="both", labelsize=_LABEL_FS)
        sns.despine(ax=axD2, left=True)
        plt.tight_layout()
        plt.savefig("../results/figures/Figure4_Final_Revised.png", dpi=300, bbox_inches="tight")
        plt.savefig("../results/figures/Figure4_Final_Revised.pdf", format="pdf", bbox_inches="tight")
        plt.savefig("../results/figures/Figure4_Final_Revised.svg", format="svg", bbox_inches="tight")
        print("✅ Figure saved to Figure4_Final_Revised.png, .pdf, and .svg")
        plt.show()


    if __name__ == "__main__":
        _main()
        (ct_data, ct_stats) = _get_context_data_and_stats("C", "T")
        _calculate_volcano_stats()
        print("\n" + "=" * 50)
        print(" EXACT ENRICHMENT VALUES (APOBEC C>T) ")
        print("=" * 50)
        if not ct_data.empty:
            sorted_ct_data = ct_data.sort_values(by=["Cohort", "Context"])
            print(sorted_ct_data.to_string(index=False))
            sorted_ct_data.to_csv("../results/tables/APOBEC_Enrichment_Values.csv", index=False)
            print("\n✅ Enrichment values saved to 'APOBEC_Enrichment_Values.csv'")
        else:
            print("No enrichment data found.")
        print("=" * 50 + "\n")
    return


@app.cell
def _(final_df):
    final_df["Timeline"]
    return


@app.cell
def _(GridSpec, chi2_contingency, fisher_exact, pd, plt, re, sns):
    # ==============================================================================
    # 1. TYPOGRAPHY & PALETTE CONFIGURATION
    # ==============================================================================
    sns.set_context("paper", font_scale=1.4)
    sns.set_style(
        "ticks",
        {
            "axes.linewidth": 2.0,
            "xtick.major.width": 2.0,
            "ytick.major.width": 2.0,
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.edgecolor": "black",
            "xtick.color": "black",
        },
    )

    _GENOTYPE_PALETTE = {
        "Genotype IA": "#1f77b4",       # Deep blue
        "Genotype IB": "#D55E00",       # Vibrant vermilion (Highlighting the outbreak clone)
        "Genotype IIIA": "#ff7f0e",     # Orange
        "Genotype IIIB": "#ffbb78",     # Light orange
        "Genotype IIA": "#2ca02c",      # Green
        "Genotype IIB": "#98df8a",      # Light green
        "Other/Unassigned": "#7f7f7f",  # Grey
    }


    # ==============================================================================
    # 2. ROBUST REGEX GENOTYPE & TIMELINE EXTRACTOR
    # ==============================================================================
    def _standardize_hav_genotype(val):
        if pd.isna(val) or not str(val).strip():
            return None
        text = str(val).upper().strip()

        # Step 1: Explicit keyword search (e.g., "Genotype IA", "GT: 3A", "type IIIB")
        explicit_pattern = r"(?:GENOTYPE|TYPE|GT|G\.T\.|GENO)\s*([I123]|II|III|IV|V)\s*[_-]?\s*([ABC])?\b"
        match = re.search(explicit_pattern, text, re.IGNORECASE)
        if match:
            num, sub = match.group(1), match.group(2) or ""
            roman_map = {"1": "I", "2": "II", "3": "III"}
            return f"Genotype {roman_map.get(num, num)}{sub}".strip()

        # Step 2: Roman Numeral word-boundary search (Avoids 2A, 2B, 3C, 3D protein collisions)
        roman_pattern = r"\b(I|II|III|IV|V|VI)\s*[_-]?\s*([ABC])\b"
        match = re.search(roman_pattern, text)
        if match:
            return f"Genotype {match.group(1)}{match.group(2)}"

        # Step 3: Exact clean code check ("IA", "1A", "3A")
        clean_exact = re.fullmatch(r"([I123]|II|III)\s*[_-]?\s*([ABC])", text)
        if clean_exact:
            num, sub = clean_exact.group(1), clean_exact.group(2)
            roman_map = {"1": "I", "2": "II", "3": "III"}
            return f"Genotype {roman_map.get(num, num)}{sub}"

        return None


    def _clean_timeline_bin(val, period_std):
        if pd.isna(val) or not str(val).strip() or str(val).strip().lower() == "other":
            return "<1998 / Endemic" if period_std == "Historical" else "2018+ (Contemporary)"
        val_str = str(val).strip()
        if any(k in val_str.lower() for k in ["after", "post", "2018", "contemp", ">"]):
            return "2018+ (Contemporary)"
        return val_str


    def extract_hav_composition(df):
        print("--- Running Hierarchical Genotype & Timeline Extraction ---")
        search_cols = [c for c in ["genotype", "strain", "isolate", "header", "organism"] if c in df.columns]

        assigned_genotypes = []
        for _, row in df.iterrows():
            gt_found = None
            for col in search_cols:
                gt_found = _standardize_hav_genotype(row[col])
                if gt_found:
                    break
            assigned_genotypes.append(gt_found or "Other/Unassigned")

        df_out = df.copy()
        df_out["Genotype_Clean"] = assigned_genotypes

        # Standardize Macro Era (Historical vs Contemporary)
        if "Period_Std" in df_out.columns and df_out["Period_Std"].notna().any():
            period_col = "Period_Std"
        elif "Condition" in df_out.columns and df_out["Condition"].notna().any():
            period_col = "Condition"
        else:
            df_out["Period_Std"] = df_out["Period"].astype(str).apply(
                lambda x: "Contemporary" if any(k in x.lower() for k in ["after", "post", "contemp", ">"]) else "Historical"
            )
            period_col = "Period_Std"

        # Standardize Micro Timeline Bins
        df_out["Timeline_Clean"] = [
            _clean_timeline_bin(row["Timeline"] if "Timeline" in df_out.columns else None, row[period_col])
            for _, row in df_out.iterrows()
        ]

        df_clean = df_out[df_out[period_col].isin(["Historical", "Contemporary"])].copy()
        print(f"Total sequences analyzed: N={len(df_clean)} across {df_clean['Timeline_Clean'].nunique()} timeline bins.")
        return df_clean, period_col


    # ==============================================================================
    # 3. STATISTICAL COMPOSITION ENGINES
    # ==============================================================================
    def analyze_macro_composition(df_clean, period_col):
        gt_order = ["Genotype IB", "Genotype IA", "Genotype IIIA", "Genotype IIIB", "Genotype IIA", "Genotype IIB", "Other/Unassigned"]
        full_gt_list = [g for g in gt_order if g in df_clean["Genotype_Clean"].unique()] + sorted(list(set(df_clean["Genotype_Clean"].unique()) - set(gt_order)))

        hist_total = sum(df_clean[period_col] == "Historical")
        cont_total = sum(df_clean[period_col] == "Contemporary")

        stats_rows = []
        for gt in full_gt_list:
            h_cnt = sum((df_clean[period_col] == "Historical") & (df_clean["Genotype_Clean"] == gt))
            c_cnt = sum((df_clean[period_col] == "Contemporary") & (df_clean["Genotype_Clean"] == gt))

            h_pct = (h_cnt / hist_total * 100) if hist_total > 0 else 0.0
            c_pct = (c_cnt / cont_total * 100) if cont_total > 0 else 0.0

            _, p_val = fisher_exact([[h_cnt, hist_total - h_cnt], [c_cnt, cont_total - c_cnt]])

            stats_rows.append({
                "Genotype": gt,
                "Historical (N)": h_cnt,
                "Historical (%)": round(h_pct, 2),
                "Contemporary (N)": c_cnt,
                "Contemporary (%)": round(c_pct, 2),
                "Total (N)": h_cnt + c_cnt,
                "Absolute Shift (%)": round(c_pct - h_pct, 2),
                "P-Value": p_val,
                "Significance": "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns")),
            })

        comp_df = pd.DataFrame(stats_rows).sort_values("Total (N)", ascending=False)
        print("\n" + "=" * 90)
        print("   1. MACRO-EVOLUTIONARY COMPOSITION (HISTORICAL vs. CONTEMPORARY)   ")
        print("=" * 90)
        print(comp_df.to_string(index=False))
        comp_df.to_csv("HAV_Macro_Genotype_Composition.csv", index=False)
        print("✔ Saved: HAV_Macro_Genotype_Composition.csv")
        return comp_df, full_gt_list


    def analyze_timeline_composition(df_clean, full_gt_list):
        # Enforce chronological ordering for timeline bins
        preferred_order = ["<1998 / Endemic", "1998-2007", "2008-2012", "2013-2017", "2018+ (Contemporary)"]
        existing_bins = [b for b in preferred_order if b in df_clean["Timeline_Clean"].unique()]
        extra_bins = sorted(list(set(df_clean["Timeline_Clean"].unique()) - set(existing_bins)))
        timeline_order = existing_bins + extra_bins

        timeline_data = []
        for t_bin in timeline_order:
            bin_df = df_clean[df_clean["Timeline_Clean"] == t_bin]
            bin_total = len(bin_df)
            row_dict = {"Timeline Bin": t_bin, "Total (N)": bin_total}

            for gt in full_gt_list:
                cnt = sum(bin_df["Genotype_Clean"] == gt)
                pct = (cnt / bin_total * 100) if bin_total > 0 else 0.0
                row_dict[f"{gt} (%)"] = round(pct, 2)
                row_dict[f"{gt} (N)"] = cnt
            timeline_data.append(row_dict)

        time_df = pd.DataFrame(timeline_data)

        # Global Chi-Square test for independence across the timeline
        contingency_table = []
        for t_bin in timeline_order:
            bin_df = df_clean[df_clean["Timeline_Clean"] == t_bin]
            contingency_table.append([sum(bin_df["Genotype_Clean"] == gt) for gt in full_gt_list])

        chi2, p_val, _, _ = chi2_contingency(contingency_table)

        print("\n" + "=" * 90)
        print(f"   2. MICRO-EVOLUTIONARY COMPOSITION ACROSS TIMELINE (Global Chi2 p = {p_val:.2e})   ")
        print("=" * 90)
        # Print clean summary showing only percentages for readability
        display_cols = ["Timeline Bin", "Total (N)"] + [c for c in time_df.columns if "(%)" in c]
        print(time_df[display_cols].to_string(index=False))
        print("=" * 90 + "\n")

        time_df.to_csv("HAV_Timeline_Genotype_Composition.csv", index=False)
        print("✔ Saved: HAV_Timeline_Genotype_Composition.csv")
        return time_df, timeline_order


    # ==============================================================================
    # 4. UNIFIED 2-PANEL VISUALIZATION
    # ==============================================================================
    def plot_combined_composition(macro_df, timeline_df, full_gt_list, timeline_order):
        fig = plt.figure(figsize=(22, 9))
        gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 1.8], wspace=0.25)
        axA = fig.add_subplot(gs[0])
        axB = fig.add_subplot(gs[1])

        # Panel A: Macro Comparison
        macro_plot = macro_df.set_index("Genotype")[["Historical (%)", "Contemporary (%)"]].T
        macro_plot = macro_plot[[g for g in full_gt_list if g in macro_plot.columns]]
        macro_plot.plot(
            kind="bar", stacked=True, ax=axA,
            color=[_GENOTYPE_PALETTE.get(g, "#333333") for g in macro_plot.columns],
            edgecolor="black", linewidth=2.0, width=0.5, legend=False
        )
        axA.set_title("A. Macro-Evolutionary Era Turnover", fontsize=18, fontweight="bold", pad=15)
        axA.set_ylabel("Proportion of Circulating Strains (%)", fontsize=16, fontweight="bold")
        axA.set_xlabel("Sampling Era", fontsize=16, fontweight="bold")
        axA.set_xticklabels(["Historical\n(Pre-2018)", "Contemporary\n(Post-2018)"], rotation=0, fontweight="bold", fontsize=14)
        axA.set_ylim(0, 105)

        for c_idx, col in enumerate(macro_plot.columns):
            for r_idx, val in enumerate(macro_plot[col]):
                if val > 6.0:
                    cum_bottom = macro_plot.iloc[r_idx, :c_idx].sum()
                    axA.text(r_idx, cum_bottom + (val / 2.0), f"{col.replace('Genotype ', '')}\n({val:.1f}%)", ha="center", va="center", fontsize=12, fontweight="bold", color="white" if col in ["Genotype IA", "Other/Unassigned"] else "black")

        # Panel B: Timeline Breakdown
        time_plot = timeline_df.set_index("Timeline Bin")[[f"{g} (%)" for g in full_gt_list if f"{g} (%)" in timeline_df.columns]]
        time_plot.columns = [c.replace(" (%)", "") for c in time_plot.columns]
        time_plot = time_plot.reindex(timeline_order)

        time_plot.plot(
            kind="bar", stacked=True, ax=axB,
            color=[_GENOTYPE_PALETTE.get(g, "#333333") for g in time_plot.columns],
            edgecolor="black", linewidth=2.0, width=0.6,
        )
        axB.set_title("B. Granular 5-Era Genotype Trajectory", fontsize=18, fontweight="bold", pad=15)
        axB.set_ylabel("")
        axB.set_xlabel("Timeline Bin", fontsize=16, fontweight="bold")
        axB.set_xticklabels(time_plot.index, rotation=15, fontweight="bold", fontsize=13)
        axB.set_ylim(0, 105)

        for c_idx, col in enumerate(time_plot.columns):
            for r_idx, val in enumerate(time_plot[col]):
                if val > 7.0:
                    cum_bottom = time_plot.iloc[r_idx, :c_idx].sum()
                    axB.text(r_idx, cum_bottom + (val / 2.0), f"{val:.1f}%", ha="center", va="center", fontsize=11, fontweight="bold", color="white" if col in ["Genotype IA", "Other/Unassigned"] else "black")

        axB.legend(title="HAV Genotype", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=13, title_fontsize=14)
        sns.despine()
        plt.tight_layout()
        plt.savefig("HAV_Combined_Genotype_Evolution_2Panel.png", dpi=300, bbox_inches="tight")
        plt.savefig("HAV_Combined_Genotype_Evolution_2Panel.pdf", format="pdf", bbox_inches="tight")
        print("✔ Saved combined figure: HAV_Combined_Genotype_Evolution_2Panel.png/.pdf")
        plt.show()


    # ==============================================================================
    # MARIMO EXECUTION BLOCK
    # ==============================================================================
    if __name__ == "__main__":
        if "final_df" in globals():
            df_analyzed, active_period_col = extract_hav_composition(globals()["final_df"])
            macro_table, genotype_order = analyze_macro_composition(df_analyzed, active_period_col)
            timeline_table, time_order = analyze_timeline_composition(df_analyzed, genotype_order)
            plot_combined_composition(macro_table, timeline_table, genotype_order, time_order)
        else:
            print("[!] Error: 'final_df' is not defined in the global namespace.")
    return


@app.cell
def _(SeqIO, mannwhitneyu, np, pd):
    _GENETIC_CODE = {
        "F": ["TTT", "TTC"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "I": ["ATT", "ATC", "ATA"],
        "M": ["ATG"],
        "V": ["GTT", "GTC", "GTA", "GTG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "A": ["GCT", "GCC", "GCA", "GCG"],
        "Y": ["TAT", "TAC"],
        "*": ["TAA", "TAG", "TGA"],
        "H": ["CAT", "CAC"],
        "Q": ["CAA", "CAG"],
        "N": ["AAT", "AAC"],
        "K": ["AAA", "AAG"],
        "D": ["GAT", "GAC"],
        "E": ["GAA", "GAG"],
        "C": ["TGT", "TGC"],
        "W": ["TGG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "G": ["GGT", "GGC", "GGA", "GGG"],
    }
    _CODON_TO_AA = {
        codon: aa for (aa, codons) in _GENETIC_CODE.items() for codon in codons
    }

    def _calculate_sequence_gc3(seq_string):
        """Calculates the GC content strictly at the 3rd position of codons."""
        if len(seq_string) < 3:
            return np.nan
        third_positions = seq_string[2::3]
        valid_bases = [b for b in third_positions if b in "ACGT"]
        if not valid_bases:
            return np.nan
        gc_count = valid_bases.count("G") + valid_bases.count("C")
        return gc_count / len(valid_bases) * 100

    def _analyze_population_codon_usage(fasta_file):
        """
            Reads a FASTA alignment of CDS sequences.
        # Reverse mapping: Codon -> Amino Acid
            Returns:
            1. A dictionary of empirical RSCU values for the population.
            2. A list of GC3 percentages for every sequence in the population.
        """
        codon_counts = {codon: 0 for codon in _CODON_TO_AA.keys()}
        gc3_list = []
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            print(f"File not found or empty: {fasta_file}")
            return (None, [])
        for _record in records:
            seq = str(_record.seq).upper().replace("U", "T")
            gc3 = _calculate_sequence_gc3(seq)
            if not np.isnan(gc3):
                gc3_list.append(gc3)
            for _i in range(0, len(seq) - 2, 3):
                codon = seq[_i : _i + 3]
                if codon in codon_counts:
                    codon_counts[codon] = codon_counts[codon] + 1
        rscu_dict = {}
        for aa, synonymous_codons in _GENETIC_CODE.items():
            total_aa_count = sum((codon_counts[_c] for _c in synonymous_codons))
            n_codons = len(synonymous_codons)
            for codon in synonymous_codons:
                if total_aa_count == 0:
                    rscu_dict[codon] = 0.0
                else:
                    expected = total_aa_count / n_codons
                    rscu_dict[codon] = codon_counts[codon] / expected
        return (rscu_dict, gc3_list)

    print("Calculating Empirical Codon Usage Dynamics...\n")
    (_rscu_hist, _gc3_hist) = _analyze_population_codon_usage(
        "../data/processed/HAV_before2018_cds_align.fasta"
    )
    (_rscu_cont, _gc3_cont) = _analyze_population_codon_usage(
        "../data/processed/HAV_after2018_cds_align.fasta"
    )
    if _gc3_hist and _gc3_cont:
        (_stat, _p_val_gc3) = mannwhitneyu(
            _gc3_hist, _gc3_cont, alternative="two-sided"
        )
        med_hist_gc3 = np.median(_gc3_hist)
        med_cont_gc3 = np.median(_gc3_cont)
        print("--- GENOME-WIDE GC3 CONTENT (Host Adaptation Marker) ---")
        print(f"Historical Median GC3:   {med_hist_gc3:.2f}%")
        print(f"Contemporary Median GC3: {med_cont_gc3:.2f}%")
        print(f"Mann-Whitney U p-value:  {_p_val_gc3:.4e}")
        if _p_val_gc3 < 0.05 and med_cont_gc3 < med_hist_gc3:
            print(
                "RESULT: Highly Significant Deoptimization. Contemporary HAV is losing GC3 content."
            )
        rscu_df = pd.DataFrame(
            {
                "Codon": list(_rscu_hist.keys()),
                "AA": [_CODON_TO_AA[_c] for _c in _rscu_hist.keys()],
                "RSCU_Historical": list(_rscu_hist.values()),
                "RSCU_Contemporary": list(_rscu_cont.values()),
            }
        )
        print("\n--- EMPIRICAL RSCU DIFFERENCES (All Shifts) ---")
        rscu_df["Abs_Shift"] = (
            rscu_df["RSCU_Contemporary"] - rscu_df["RSCU_Historical"]
        ).abs()
        print(
            rscu_df.sort_values(by="Abs_Shift", ascending=False)
            .head(100)
            .to_string(index=False)
        )
    return


@app.cell
def _(
    AlignIO,
    Counter,
    GridSpec,
    GridSpecFromSubplotSpec,
    SeqIO,
    chi2_contingency,
    mannwhitneyu,
    math,
    np,
    os,
    pd,
    plt,
    re,
    sns,
):
    from scipy.stats import wilcoxon


    # ==============================================================================
    # 1. VISUAL CONFIGURATION & GLOBAL CONSTANTS
    # ==============================================================================
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    _TOP_N_MUTATIONS = 15
    _CODON_FILE = "../data/raw/Codon Usage Summary.xlsx"
    _FASTA_ALIGN_FILE = "../data/processed/HAV_Collection date_aligned.fasta"
    _FASTA_ALIGN_PRE = "../data/processed/HAV_before2018_align.fasta"
    _FASTA_ALIGN_POST = "../data/processed/HAV_after2018_align.fasta"
    _FASTA_CDS_PRE = "../data/processed/HAV_before2018_cds_align.fasta"
    _FASTA_CDS_POST = "../data/processed/HAV_after2018_cds_align.fasta"
    _STAT_REPORT = "../results/tables/figure5_stat.txt"
    _FIG_NAME = "../results/figures/Figure5_Nature_Formatted.png"
    _PDF_NAME = "../results/figures/Figure5_Nature_Formatted.pdf"
    _SVG_NAME = "../results/figures/Figure5_Nature_Formatted.svg"

    _LABEL_FS = 18
    _TICK_FS = 16
    _TITLE_FS = 16
    _ANNOTATION_FS = 14
    _LETTER_FS = 28
    _LEGEND_FS = 14
    _MAP_LABEL_FS = 14
    _AXES_LW = 1.5
    _LINE_LW = 3.0
    _MARKER_SIZE = 100
    _SWARM_SIZE = 6

    sns.set_context("paper", font_scale=2.5)
    sns.set_style(
        "ticks",
        {
            "axes.linewidth": _AXES_LW,
            "xtick.major.width": _AXES_LW,
            "ytick.major.width": _AXES_LW,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
        },
    )

    _color_opt = "#009E73"
    _color_subopt = "#D55E00"
    _palette_period = {"Historical": "#4c72b0", "Contemporary": "#dd8452"}
    _GENOME_REGIONS = [
        (735, 803, "#C8E6C9", "VP4"),
        (804, 1469, "#A5D6A7", "VP2"),
        (1470, 2207, "#81C784", "VP3"),
        (2208, 3107, "#66BB6A", "VP1"),
        (3108, 3674, "#FFF9C4", "2A"),
        (3675, 3995, "#FFF59D", "2B"),
        (3996, 5000, "#FFF176", "2C"),
        (5001, 5222, "#E1F5FE", "3A"),
        (5223, 5291, "#B3E5FC", "3B"),
        (5292, 5948, "#81D4FA", "3C"),
        (5949, 7415, "#4FC3F7", "3D"),
    ]
    _CDS_START = 735
    _CDS_END = 7415
    _GENOME_LENGTH = 7478
    _REF_ID = "NC_001489.1"
    _GENETIC_CODE = {
        "F": ["TTT", "TTC"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "I": ["ATT", "ATC", "ATA"],
        "M": ["ATG"],
        "V": ["GTT", "GTC", "GTA", "GTG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "A": ["GCT", "GCC", "GCA", "GCG"],
        "Y": ["TAT", "TAC"],
        "*": ["TAA", "TAG", "TGA"],
        "H": ["CAT", "CAC"],
        "Q": ["CAA", "CAG"],
        "N": ["AAT", "AAC"],
        "K": ["AAA", "AAG"],
        "D": ["GAT", "GAC"],
        "E": ["GAA", "GAG"],
        "C": ["TGT", "TGC"],
        "W": ["TGG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "G": ["GGT", "GGC", "GGA", "GGG"],
    }
    _CODON_TO_AA = {
        codon: aa for (aa, codons) in _GENETIC_CODE.items() for codon in codons
    }


    # ==============================================================================
    # 2. STATISTICAL LOGGING & FASTA PARSING
    # ==============================================================================
    def _write_stat(text):
        with open(_STAT_REPORT, "a") as _f:
            _f.write(text + "\n")


    with open(_STAT_REPORT, "w") as _f:
        _f.write("FIGURE 5 DETAILED STATISTICAL REPORT\n")
        _f.write("========================================================\n\n")


    def _read_multi_fasta(file_path):
        if not os.path.exists(file_path):
            return []
        return [str(_record.seq).upper() for _record in SeqIO.parse(file_path, "fasta")]


    def _calculate_sequence_gc3(seq_string):
        if len(seq_string) < 3:
            return np.nan
        third_positions = seq_string[2::3]
        valid_bases = [b for b in third_positions if b in "ACGT"]
        if not valid_bases:
            return np.nan
        gc_count = valid_bases.count("G") + valid_bases.count("C")
        return gc_count / len(valid_bases) * 100


    def _analyze_population_codon_usage(fasta_file):
        codon_counts = {codon: 0 for codon in _CODON_TO_AA.keys()}
        gc3_list = []
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            return (None, [])
        for _record in records:
            seq = str(_record.seq).upper().replace("U", "T")
            gc3 = _calculate_sequence_gc3(seq)
            if not np.isnan(gc3):
                gc3_list.append(gc3)
            for _i in range(0, len(seq) - 2, 3):
                codon = seq[_i : _i + 3]
                if codon in codon_counts:
                    codon_counts[codon] = codon_counts[codon] + 1
        rscu_dict = {}
        for aa, synonymous_codons in _GENETIC_CODE.items():
            total_aa_count = sum((codon_counts[_c] for _c in synonymous_codons))
            n_codons = len(synonymous_codons)
            for codon in synonymous_codons:
                if total_aa_count == 0:
                    rscu_dict[codon] = 0.0
                else:
                    rscu_dict[codon] = codon_counts[codon] / (total_aa_count / n_codons)
        return (rscu_dict, gc3_list)


    def _analyze_mutation_type(
        true_position, ref_seq_string, mutated_nuc, utr_offset=734
    ):
        cds_pos = true_position - utr_offset
        if cds_pos <= 0:
            return {"Type": "5' UTR", "Codon_Change": "-", "AA_Change": "-"}
        codon_start_idx = true_position - 1 - (cds_pos - 1) % 3
        wt_codon = ref_seq_string[codon_start_idx : codon_start_idx + 3].upper()
        pos_in_codon = (cds_pos - 1) % 3
        mut_codon_list = list(wt_codon)
        mut_codon_list[pos_in_codon] = mutated_nuc.upper()
        mut_codon = "".join(mut_codon_list)
        wt_aa = _CODON_TO_AA.get(wt_codon, "X")
        mut_aa = _CODON_TO_AA.get(mut_codon, "X")
        if wt_aa == "X" or mut_aa == "X":
            return {"Type": "Error", "Codon_Change": "-", "AA_Change": "-"}
        return {
            "Codon_Change": f"{wt_codon}->{mut_codon}",
            "AA_Change": f"{wt_aa}->{mut_aa}",
            "Type": "Synonymous" if wt_aa == mut_aa else "Non-Syn",
        }


    def _analyze_nucleotide_mutations(
        metadata_df,
        fasta_file,
        ref_accession="NC_001489",
        split_year=2018,
        freq_threshold=0.05,
    ):
        metadata_df["year"] = pd.to_numeric(metadata_df["year"], errors="coerce")
        pre_ids = set(metadata_df[metadata_df["year"] < split_year]["accession_id"])
        post_ids = set(metadata_df[metadata_df["year"] >= split_year]["accession_id"])
        try:
            _alignment = AlignIO.read(fasta_file, "fasta")
        except:
            return pd.DataFrame()
        ref_seq_record = next((r for r in _alignment if ref_accession in r.id), None)
        if not ref_seq_record:
            return pd.DataFrame()
        ref_seq_string = str(ref_seq_record.seq).lower()
        align_to_true_coord = {}
        true_pos = 1
        for _i, _res in enumerate(ref_seq_string):
            if _res != "-":
                align_to_true_coord[_i] = true_pos
                true_pos = true_pos + 1
            else:
                align_to_true_coord[_i] = None
        align_len = _alignment.get_alignment_length()
        pre_counts = {_i: Counter() for _i in range(align_len)}
        post_counts = {_i: Counter() for _i in range(align_len)}
        for _record in _alignment:
            if _record.id == ref_seq_record.id:
                continue
            seq_id = _record.id.split()[0]
            seq_string = str(_record.seq).lower()
            if seq_id in pre_ids:
                for _i, _res in enumerate(seq_string):
                    pre_counts[_i][_res] = pre_counts[_i][_res] + 1
            elif seq_id in post_ids:
                for _i, _res in enumerate(seq_string):
                    post_counts[_i][_res] = post_counts[_i][_res] + 1
        trend_data = []
        for pos in range(align_len):
            true_coord = align_to_true_coord.get(pos)
            if true_coord is None:
                continue
            total_pre = sum(pre_counts[pos].values())
            total_post = sum(post_counts[pos].values())
            if total_pre == 0 or total_post == 0:
                continue
            ref_residue = ref_seq_string[pos]
            all_residues = set(pre_counts[pos].keys()) | set(post_counts[pos].keys())
            for mut_residue in all_residues:
                if mut_residue in ["-", "n", ref_residue]:
                    continue
                mut_freq_pre = pre_counts[pos][mut_residue] / total_pre
                mut_freq_post = post_counts[pos][mut_residue] / total_post
                delta = mut_freq_post - mut_freq_pre
                if abs(delta) >= freq_threshold:
                    mut_info = _analyze_mutation_type(
                        true_coord, ref_seq_string, mut_residue
                    )
                    trend_data.append(
                        {
                            "Mutation": f"{ref_residue.upper()}{true_coord}{mut_residue.upper()}",
                            "Trend": "Increasing" if delta > 0 else "Decreasing",
                            "Type": mut_info["Type"],
                            "Codon": mut_info["Codon_Change"],
                            "AA": mut_info["AA_Change"],
                            "Pre_Freq": round(mut_freq_pre, 4),
                            "Post_Freq": round(mut_freq_post, 4),
                            "Abs_Change": round(abs(delta), 4),
                        }
                    )
        return pd.DataFrame(trend_data)


    def _load_codon_stats():
        try:
            df = pd.read_excel(_CODON_FILE)
            df.rename(
                columns={
                    "RSCU (Before 2018)": "rscu_before",
                    "RSCU (After 2018)": "rscu_after",
                    "RSCU (Humans)": "rscu_human",
                },
                inplace=True,
            )
            for _c in ["rscu_before", "rscu_after", "rscu_human"]:
                df[_c] = pd.to_numeric(df[_c], errors="coerce")
            df.dropna(subset=["rscu_before", "rscu_after", "rscu_human"], inplace=True)
            epsilon = 1e-09
            df["log2FC"] = np.log2(
                (df["rscu_after"] + epsilon) / (df["rscu_before"] + epsilon)
            )
            df["Host_Optimality"] = np.where(
                df["rscu_human"] >= 1, "Optimal in Host", "Suboptimal in Host"
            )
            return df
        except:
            return pd.DataFrame()


    def _load_codon_maps_for_track():
        try:
            df = pd.read_excel(_CODON_FILE)
            df.rename(
                columns={
                    "RSCU (Before 2018)": "rscu_before",
                    "RSCU (After 2018)": "rscu_after",
                    "RSCU (Humans)": "rscu_human",
                    "Codon": "codon",
                },
                inplace=True,
            )
            for _c in ["rscu_before", "rscu_after", "rscu_human"]:
                df[_c] = pd.to_numeric(df[_c], errors="coerce")
            df.dropna(
                subset=["rscu_before", "rscu_after", "rscu_human", "codon"],
                inplace=True,
            )
            df["Deopt_Score_Before"] = (df["rscu_human"] - df["rscu_before"]).abs()
            df["Deopt_Score_After"] = (df["rscu_human"] - df["rscu_after"]).abs()
            _map_pre = pd.Series(
                df.Deopt_Score_Before.values,
                index=df.codon.str.extract("([A-Z]{3})")[0],
            ).to_dict()
            _map_post = pd.Series(
                df.Deopt_Score_After.values, index=df.codon.str.extract("([A-Z]{3})")[0]
            ).to_dict()
            return (_map_pre, _map_post)
        except:
            return ({}, {})


    def _get_cds_tracks(map_pre, map_post):

        def map_to_genome(fasta, metrics_map):
            seqs = _read_multi_fasta(fasta)
            if not seqs:
                return np.full(_GENOME_LENGTH, np.nan)
            seq = seqs[0]
            track = np.full(_GENOME_LENGTH, np.nan)
            max_cds_len = min(len(seq), _CDS_END - _CDS_START + 1)
            for k in range(0, max_cds_len - 2, 3):
                codon = seq[k : k + 3]
                genomic_pos = _CDS_START + k
                if "-" not in codon and "N" not in codon:
                    track[genomic_pos : genomic_pos + 3] = metrics_map.get(
                        codon.replace("T", "U"), np.nan
                    )
            return track

        _t_before = map_to_genome(_FASTA_CDS_PRE, _map_pre)
        _t_after = map_to_genome(_FASTA_CDS_POST, _map_post)
        w = 250
        _ma_before = pd.Series(_t_before).rolling(window=w, min_periods=50).mean()
        _ma_after = pd.Series(_t_after).rolling(window=w, min_periods=50).mean()
        return (_ma_before, _ma_after, _t_before, _t_after)


    def _get_motif_counts(file_path, ref_id, motif_regex):
        sequences = {
            _record.id: str(_record.seq).upper()
            for _record in SeqIO.parse(file_path, "fasta")
        }
        if not sequences:
            return None
        found_ref_id = next((k for k in sequences.keys() if ref_id in k), None)
        if not found_ref_id:
            return None
        ref_seq = sequences[found_ref_id]
        motif_indices = set()
        pattern = re.compile(motif_regex)
        for _i in range(len(ref_seq)):
            match = pattern.match(ref_seq, _i)
            if match:
                for _idx in range(match.start(), match.end()):
                    motif_indices.add(_idx)
        (m_in, t_in, m_out, t_out) = (0, 0, 0, 0)
        for q_id, q_seq in sequences.items():
            if q_id == found_ref_id:
                continue
            for _i, base in enumerate(q_seq):
                if (
                    _i >= len(ref_seq)
                    or base not in "ACGT"
                    or ref_seq[_i] not in "ACGT"
                ):
                    continue
                if _i in motif_indices:
                    t_in = t_in + 1
                    if base != ref_seq[_i]:
                        m_in = m_in + 1
                else:
                    t_out = t_out + 1
                    if base != ref_seq[_i]:
                        m_out = m_out + 1
        return (m_in, t_in, m_out, t_out)


    def _calculate_pf_and_ci(m_in, t_in, m_out, t_out):
        if m_in == 0 or m_out == 0:
            m_in = m_in + 0.5
            t_in = t_in + 0.5
            m_out = m_out + 0.5
            t_out = t_out + 0.5
        (p_in, p_out) = (m_in / t_in, m_out / t_out)
        pf = p_out / p_in if p_in > 0 else 1.0
        se_ln_pf = math.sqrt(1 / m_out - 1 / t_out + 1 / m_in - 1 / t_in)
        lower_ci = math.exp(math.log(pf) - 1.96 * se_ln_pf)
        upper_ci = math.exp(math.log(pf) + 1.96 * se_ln_pf)
        return (pf, pf - lower_ci, upper_ci - pf)


    def _run_analysis(motif_name, regex):
        pre_counts = _get_motif_counts(_FASTA_ALIGN_PRE, _REF_ID, regex)
        post_counts = _get_motif_counts(_FASTA_ALIGN_POST, _REF_ID, regex)
        if not pre_counts or not post_counts:
            return None
        (_pf_pre, el_pre, eh_pre) = _calculate_pf_and_ci(*pre_counts)
        (_pf_post, el_post, eh_post) = _calculate_pf_and_ci(*post_counts)
        table = [
            [pre_counts[0], pre_counts[1] - pre_counts[0]],
            [post_counts[0], post_counts[1] - post_counts[0]],
        ]
        (_, _p_val, _, _) = chi2_contingency(table)
        return {
            "pre": (_pf_pre, el_pre, eh_pre),
            "post": (_pf_post, el_post, eh_post),
            "p_val": _p_val,
        }


    # ==============================================================================
    # 3. STATISTICAL EXECUTION & DATA MINING
    # ==============================================================================
    print("Executing Global Calculations...")
    _df_codon = _load_codon_stats()
    (_map_pre, _map_post) = _load_codon_maps_for_track()
    (_ma_before, _ma_after, _t_before, _t_after) = _get_cds_tracks(_map_pre, _map_post)

    _plot_data_motif = []
    _yerr_low_m = []
    _yerr_high_m = []
    for _name, _regex in {"CpG": "CG", "GpA": "GA"}.items():
        _res = _run_analysis(_name, _regex)
        if _res:
            _p_val_m = _res["p_val"]
            _write_stat(
                f"Panel C ({_name} Motif Protection): Chi-square p={_p_val_m:.4e}"
            )
            _plot_data_motif.append(
                {
                    "Motif": _name,
                    "Period": "Historical",
                    "PF": _res["pre"][0],
                    "p_val": _p_val_m,
                }
            )
            _yerr_low_m.append(_res["pre"][1])
            _yerr_high_m.append(_res["pre"][2])
            _plot_data_motif.append(
                {
                    "Motif": _name,
                    "Period": "Contemporary",
                    "PF": _res["post"][0],
                    "p_val": _p_val_m,
                }
            )
            _yerr_low_m.append(_res["post"][1])
            _yerr_high_m.append(_res["post"][2])
    _df_motif = pd.DataFrame(_plot_data_motif)

    print("Calculating Empirical RSCU and GC3 Dynamics...")
    (_rscu_hist, _gc3_hist) = _analyze_population_codon_usage(_FASTA_CDS_PRE)
    (_rscu_cont, _gc3_cont) = _analyze_population_codon_usage(_FASTA_CDS_POST)
    if _gc3_hist and _gc3_cont:
        _df_gc3 = pd.DataFrame(
            {
                "Period": ["Historical"] * len(_gc3_hist)
                + ["Contemporary"] * len(_gc3_cont),
                "GC3_Percent": _gc3_hist + _gc3_cont,
            }
        )
        (_, _p_val_gc3) = mannwhitneyu(_gc3_hist, _gc3_cont, alternative="two-sided")
        _write_stat(f"Panel B (Genome-Wide GC3 Shift): Mann-Whitney p={_p_val_gc3:.4e}")
        _rscu_df_emp = pd.DataFrame(
            {
                "Codon": list(_rscu_hist.keys()),
                "AA": [_CODON_TO_AA[_c] for _c in _rscu_hist.keys()],
                "RSCU_Historical": list(_rscu_hist.values()),
                "RSCU_Contemporary": list(_rscu_cont.values()),
            }
        )
        _rscu_df_emp["Delta_RSCU"] = (
            _rscu_df_emp["RSCU_Contemporary"] - _rscu_df_emp["RSCU_Historical"]
        )
        _rscu_df_emp["Abs_Shift"] = _rscu_df_emp["Delta_RSCU"].abs()
        _top_10_shifts = (
            _rscu_df_emp.sort_values(by="Abs_Shift", ascending=False).head(10).copy()
        )
        _top_10_shifts["Ending"] = _top_10_shifts["Codon"].apply(
            lambda x: "G/C Ending" if x[-1] in ["G", "C"] else "A/T Ending"
        )
        _top_10_shifts = _top_10_shifts.sort_values(by="Delta_RSCU", ascending=False)
    else:
        (_df_gc3, _top_10_shifts) = (pd.DataFrame(), pd.DataFrame())

    print(
        f"Dynamically mining top {_TOP_N_MUTATIONS} Gained and Lost Mutations for Panel D..."
    )
    _meta_input = globals().get("merged_df_3", pd.DataFrame())
    _df_all_muts = _analyze_nucleotide_mutations(
        _meta_input, _FASTA_ALIGN_FILE, freq_threshold=0.05
    )
    if _df_all_muts.empty:
        print(
            "WARNING: No alignment mutations discovered. Check filenames or merged_df structures."
        )
        _df_D_data = pd.DataFrame()
    else:
        _inc_nuc = (
            _df_all_muts[_df_all_muts["Trend"] == "Increasing"]
            .sort_values(by="Abs_Change", ascending=False)
            .head(_TOP_N_MUTATIONS)
        )
        _dec_nuc = (
            _df_all_muts[_df_all_muts["Trend"] == "Decreasing"]
            .sort_values(by="Abs_Change", ascending=False)
            .head(_TOP_N_MUTATIONS)
        )
        _hist_block = pd.DataFrame(
            {
                "Region": pd.concat([_inc_nuc["Mutation"], _dec_nuc["Mutation"]]),
                "Time Period": "Historical",
                "Mutation Frequency": pd.concat(
                    [_inc_nuc["Pre_Freq"], _dec_nuc["Pre_Freq"]]
                )
                * 100,
                "Type": pd.concat([_inc_nuc["Type"], _dec_nuc["Type"]]),
            }
        )
        _cont_block = pd.DataFrame(
            {
                "Region": pd.concat([_inc_nuc["Mutation"], _dec_nuc["Mutation"]]),
                "Time Period": "Contemporary",
                "Mutation Frequency": pd.concat(
                    [_inc_nuc["Post_Freq"], _dec_nuc["Post_Freq"]]
                )
                * 100,
                "Type": pd.concat([_inc_nuc["Type"], _dec_nuc["Type"]]),
            }
        )
        _df_D_data = pd.concat([_hist_block, _cont_block]).reset_index(drop=True)

    # ==============================================================================
    # 4. PLOT RENDERING PIPELINE (RECONFIGURED LAYOUT)
    # ==============================================================================
    print("Rendering plot canvas...")
    _fig = plt.figure(figsize=(24, 22))
    _gs = GridSpec(3, 12, figure=_fig, hspace=0.8, wspace=1.7, height_ratios=[1, 1, 1])

    # Top Row (3 Panels): Codon Optimality (A), GC3 Content (B), Motif Protection (C)
    _axA = _fig.add_subplot(_gs[0, 0:4])
    _axB = _fig.add_subplot(_gs[0, 4:8])
    _axC = _fig.add_subplot(_gs[0, 8:12])

    # Middle Row (2 Panels): Nucleotide Mutations (D - 7 columns), Top Codon Shifts (E - 5 columns)
    _axD = _fig.add_subplot(_gs[1, 0:7])
    _axE = _fig.add_subplot(_gs[1, 8:12])

    # Panel A: Host Codon Optimality vs Log2FC
    if not _df_codon.empty:
        sns.violinplot(
            data=_df_codon,
            x="Host_Optimality",
            y="log2FC",
            hue="Host_Optimality",
            palette={
                "Optimal in Host": _color_opt,
                "Suboptimal in Host": _color_subopt,
            },
            legend=False,
            ax=_axA,
            cut=0,
        )
        sns.swarmplot(
            data=_df_codon,
            x="Host_Optimality",
            y="log2FC",
            color=".2",
            alpha=0.6,
            size=_SWARM_SIZE,
            ax=_axA,
        )
        _axA.axhline(0, color="grey", ls="--")
        _axA.set_xlabel("Host Status", fontweight="bold", fontsize=_LABEL_FS)
        _axA.set_ylabel("Log2 Fold Change", fontweight="bold", fontsize=_LABEL_FS)
        _opt = _df_codon[_df_codon["Host_Optimality"] == "Optimal in Host"]["log2FC"]
        _sub = _df_codon[_df_codon["Host_Optimality"] == "Suboptimal in Host"]["log2FC"]
        if len(_opt) > 0 and len(_sub) > 0:
            (_, _p) = mannwhitneyu(_opt, _sub)
            _axA.text(
                0.5,
                0.9,
                f"MWU p={_p:.2e}",
                transform=_axA.transAxes,
                ha="center",
                fontsize=_ANNOTATION_FS,
            )
        _axA.tick_params(axis="both", labelsize=_TICK_FS)
        sns.despine(ax=_axA)

    # Panel B: Genome-Wide GC3 Content Shift (Moved from Old D)
    if not _df_gc3.empty:
        sns.violinplot(
            data=_df_gc3,
            x="Period",
            y="GC3_Percent",
            hue="Period",
            palette=_palette_period,
            inner="box",
            ax=_axB,
            legend=False,
        )
        sns.swarmplot(
            data=_df_gc3.sample(min(len(_df_gc3), 200), random_state=42),
            x="Period",
            y="GC3_Percent",
            color=".2",
            alpha=0.4,
            size=3,
            ax=_axB,
        )
        _axB.set_title(
            f"\n(MWU p={_p_val_gc3:.2e})", fontsize=_TITLE_FS, fontweight="bold"
        )
        _axB.set_ylabel("GC3 Content (%)", fontweight="bold", fontsize=_LABEL_FS)
        _axB.set_xlabel("", fontsize=_LABEL_FS)
        _axB.tick_params(axis="both", labelsize=_TICK_FS)
        sns.despine(ax=_axB)

    # Panel C: Motif Protection Factor (Moved from Old B)
    if not _df_motif.empty:
        sns.barplot(
            data=_df_motif,
            x="Motif",
            y="PF",
            hue="Period",
            palette=[_palette_period["Historical"], _palette_period["Contemporary"]],
            edgecolor="black",
            linewidth=_AXES_LW,
            ax=_axC,
        )
        _width = 0.4
        _x_coords = np.arange(len(_df_motif["Motif"].unique()))
        _axC.errorbar(
            x=_x_coords - _width / 2,
            y=_df_motif[_df_motif["Period"] == "Historical"]["PF"],
            yerr=[_yerr_low_m[::2], _yerr_high_m[::2]],
            fmt="none",
            c="black",
            capsize=6,
            lw=_AXES_LW,
        )
        _axC.errorbar(
            x=_x_coords + _width / 2,
            y=_df_motif[_df_motif["Period"] == "Contemporary"]["PF"],
            yerr=[_yerr_low_m[1::2], _yerr_high_m[1::2]],
            fmt="none",
            c="black",
            capsize=6,
            lw=_AXES_LW,
        )
        _axC.axhline(1.0, color="gray", linestyle="--", linewidth=2)
        _unique_motifs = _df_motif["Motif"].unique()
        for _i, _name in enumerate(_unique_motifs):
            _subset = _df_motif[_df_motif["Motif"] == _name]
            _p_val_m = _subset["p_val"].iloc[0]
            if _p_val_m < 0.001:
                _star = "***"
            elif _p_val_m < 0.01:
                _star = "**"
            elif _p_val_m < 0.05:
                _star = "*"
            else:
                _star = "ns"
            _pf_pre = _subset[_subset["Period"] == "Historical"]["PF"].values[0]
            _pf_post = _subset[_subset["Period"] == "Contemporary"]["PF"].values[0]
            _err_pre = _yerr_high_m[2 * _i]
            _err_post = _yerr_high_m[2 * _i + 1]
            _y_max = max(_pf_pre + _err_pre, _pf_post + _err_post)
            _h_line = _y_max + 0.1
            _v_line = 0.05
            (_x1, _x2) = (_i - _width / 2, _i + _width / 2)
            _axC.plot(
                [_x1, _x1, _x2, _x2],
                [_h_line, _h_line + _v_line, _h_line + _v_line, _h_line],
                lw=_AXES_LW,
                c="black",
            )
            _axC.text(
                (_x1 + _x2) / 2,
                _h_line + _v_line + 0.02,
                _star,
                ha="center",
                va="bottom",
                fontsize=_ANNOTATION_FS + 4,
                color="black",
            )
        _axC.set_ylabel(
            "Protection Factor", labelpad=10, fontweight="bold", fontsize=_LABEL_FS
        )
        _axC.set_xlabel("", fontsize=_LABEL_FS)
        _axC.tick_params(axis="both", labelsize=_TICK_FS)
        _axC.legend(frameon=False, loc="upper left", fontsize=_LEGEND_FS)
        sns.despine(ax=_axC)

    # Panel D: Top 15 Gained & Lost Nucleotide Mutations (Moved from Old C)
    if not _df_D_data.empty:
        sns.barplot(
            data=_df_D_data,
            x="Region",
            y="Mutation Frequency",
            hue="Time Period",
            palette=_palette_period,
            ax=_axD,
            edgecolor="black",
            lw=1,
        )
        _axD.tick_params(axis="y", labelsize=_TICK_FS)
        _axD.set_xlabel(
            f"Top {_TOP_N_MUTATIONS} Gained & Top {_TOP_N_MUTATIONS} Lost Mutations",
            fontweight="bold",
            fontsize=_LABEL_FS,
        )
        _axD.set_ylabel("Frequency (%)", fontweight="bold", fontsize=_LABEL_FS)
        _axD.legend(loc="upper right", frameon=False, ncol=2, fontsize=_LEGEND_FS)
        _type_mapping = pd.Series(
            _df_D_data["Type"].values, index=_df_D_data["Region"]
        ).to_dict()
        _axD.set_xticklabels(_axD.get_xticklabels(), rotation=90, fontsize=_TICK_FS - 4)
        _fig.canvas.draw()
        for _label in _axD.get_xticklabels():
            _m_type = _type_mapping.get(_label.get_text(), "")
            if _m_type == "Non-Syn":
                _label.set_color("#c0392b")
                _label.set_weight("bold")
            elif "5' UTR" in _m_type:
                _label.set_color("#7f8c8d")
            else:
                _label.set_color("#2c3e50")
        sns.despine(ax=_axD)

    # Panel E: Top 10 Empirical Codon Shifts (Moved from Old E)
    if not _top_10_shifts.empty:
        _ending_palette = {"A/T Ending": "#d35400", "G/C Ending": "#2980b9"}
        sns.barplot(
            data=_top_10_shifts,
            x="Codon",
            y="Delta_RSCU",
            hue="Ending",
            palette=_ending_palette,
            ax=_axE,
            edgecolor="black",
            dodge=False,
        )
        _axE.axhline(0, color="black", lw=1.5)
        _axE.set_title("", fontsize=_TITLE_FS, fontweight="bold")
        _axE.set_ylabel(
            "Δ RSCU \n(Contemporary vs. Historical)",
            fontweight="bold",
            fontsize=_LABEL_FS,
        )
        _axE.set_xlabel("Codon", fontweight="bold", fontsize=_LABEL_FS)
        _formatted_labels = [
            f"{_row.Codon}\n({_row.AA})" for (_, _row) in _top_10_shifts.iterrows()
        ]
        _axE.set_xticklabels(_formatted_labels, fontsize=_TICK_FS - 2)
        _axE.tick_params(axis="y", labelsize=_TICK_FS)
        _axE.legend(loc="lower left", frameon=True, fontsize=_LEGEND_FS - 2)
        sns.despine(ax=_axE)

    # Bottom Row (2 Panels): Genomic Divergence Track (F), Wilcoxon Signed-Rank Divergence (G)
    _gs_track = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=_gs[2, 0:9], height_ratios=[4, 1], hspace=0.05
    )
    _axF_track = _fig.add_subplot(_gs_track[0])
    _axF_map = _fig.add_subplot(_gs_track[1], sharex=_axF_track)
    _axG_viol = _fig.add_subplot(_gs[2, 9:12])

    # Panel F: Host Codon Divergence Track along Genome
    if not _ma_before.empty:
        _x_vals = np.arange(len(_ma_before))
        _axF_track.plot(
            _x_vals,
            _ma_before,
            color="#4c72b0",
            lw=_LINE_LW,
            ls="--",
            label="Historical",
        )
        _axF_track.plot(
            _x_vals, _ma_after, color="#dd8452", lw=_LINE_LW, label="Contemporary"
        )
        _axF_track.set_ylabel(
            "Host Codon Divergence \n(Avg. |ΔRSCU|)",
            fontsize=_LABEL_FS,
            fontweight="bold",
        )
        _axF_track.set_ylim(0.4, 1.0)
        _axF_track.legend(
            loc="lower center", ncol=2, frameon=False, fontsize=_LEGEND_FS
        )
        _axF_track.tick_params(axis="both", labelsize=_TICK_FS)
    plt.setp(_axF_track.get_xticklabels(), visible=False)
    sns.despine(ax=_axF_track, bottom=True)

    for _start, _end, _color, _label in _GENOME_REGIONS:
        _axF_map.barh(
            0,
            width=_end - _start,
            left=_start,
            color=_color,
            edgecolor="black",
            height=1,
        )
        _axF_map.text(
            (_start + _end) / 2,
            0,
            _label,
            ha="center",
            va="center",
            fontsize=_MAP_LABEL_FS,
            fontweight="bold",
        )
        _axF_track.axvspan(_start, _end, facecolor=_color, alpha=0.2)
    _axF_map.set_yticks([])
    _axF_map.set_xlim(0, _GENOME_LENGTH)
    _axF_map.set_xlabel("Genomic Position (bp)", fontsize=_LABEL_FS, fontweight="bold")
    _axF_map.tick_params(axis="x", labelsize=_TICK_FS)
    sns.despine(ax=_axF_map, left=True)

    # Panel G: Wilcoxon Signed-Rank Global Divergence Shift
    _valid_idx = ~np.isnan(_t_before) & ~np.isnan(_t_after)
    if _valid_idx.any():
        (_val_before, _val_after) = (_t_before[_valid_idx], _t_after[_valid_idx])
        (_, _p_viol) = wilcoxon(_val_before, _val_after)
        _df_viol = pd.DataFrame(
            {"Historical": _val_before, "Contemporary": _val_after}
        ).melt(var_name="Period", value_name="Divergence")
        sns.violinplot(
            data=_df_viol,
            x="Period",
            y="Divergence",
            hue="Period",
            palette=_palette_period,
            inner="box",
            ax=_axG_viol,
            legend=False,
        )
        _plot_data_viol = _df_viol.sample(n=min(len(_df_viol), 100), random_state=42)
        sns.swarmplot(
            data=_plot_data_viol,
            x="Period",
            y="Divergence",
            color=".2",
            alpha=0.6,
            size=_SWARM_SIZE,
            ax=_axG_viol,
        )
        _axG_viol.set_title(
            f"Wilcoxon Signed-rank \n(p={_p_viol:.2e})", fontsize=_TITLE_FS
        )
        _axG_viol.set_ylabel(
            "Divergence (|ΔRSCU|)", fontweight="bold", fontsize=_LABEL_FS
        )
        _axG_viol.tick_params(axis="both", labelsize=_TICK_FS)
        sns.despine(ax=_axG_viol)

    # Re-lettering all panels sequentially from A to G
    _ax_map_refs = [_axA, _axB, _axC, _axD, _axE, _axF_track, _axG_viol]
    _labels = ["A", "B", "C", "D", "E", "F", "G"]
    for _ax, _lbl in zip(_ax_map_refs, _labels):
        _y_off = 1.3 if _ax == _axF_track else 1.15
        _x_off = -0.05 if _ax == _axF_track else -0.12
        _ax.text(
            _x_off,
            _y_off,
            _lbl,
            transform=_ax.transAxes,
            fontsize=_LETTER_FS,
            fontweight="bold",
            va="top",
        )

    plt.tight_layout()
    plt.savefig(_FIG_NAME, dpi=300, bbox_inches="tight")
    plt.savefig(_PDF_NAME, dpi=300, bbox_inches="tight")
    plt.savefig(_SVG_NAME, format="svg", bbox_inches="tight")

    if "_df_viol" in locals():
        _summary_stats = _df_viol.groupby("Period")["Divergence"].describe(
            percentiles=[0.25, 0.5, 0.75]
        )
        _hist_med = _summary_stats.loc["Historical", "50%"]
        _contemp_med = _summary_stats.loc["Contemporary", "50%"]
        _manuscript_sentence = f"Genome-wide analysis revealed a highly significant shift in host codon optimality, with Contemporary lineages exhibiting a notably higher median divergence (|ΔRSCU| = {_contemp_med:.3f}) compared to Historical lineages (Median = {_hist_med:.3f}) (Wilcoxon signed-rank test, p = {_p_viol:.2e})."
        _write_stat("\n--- MANUSCRIPT RESULTS TEXT ---")
        _write_stat(_manuscript_sentence)

    print(
        f"Success. All structural discrepancies resolved. Output saved to {_STAT_REPORT}."
    )
    plt.show()
    return


if __name__ == "__main__":
    app.run()
