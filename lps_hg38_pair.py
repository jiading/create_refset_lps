import os
import sys
import subprocess
import pandas as pd


def load_exec_map(map_file_path):
    """
    Reads the CSV/TSV mapping file containing 'sample_key' and 'exec_id'.
    Returns a dictionary mapping sample_key -> exec_id (and base gccode -> exec_id as fallback).
    """
    if not os.path.exists(map_file_path):
        raise FileNotFoundError(f"Mapping file not found: {map_file_path}")

    # Auto-detect tab/space/comma separators
    df = pd.read_csv(map_file_path, sep=r'\s+|,', engine='python')

    # Ensure column names are stripped of whitespace
    df.columns = df.columns.str.strip()

    exec_map = {}
    for _, row in df.iterrows():
        sample_key = str(row['sample_key']).strip()
        exec_id = str(row['exec_id']).strip()

        # Key by exact sample_key (e.g. GC184872-I4A02)
        exec_map[sample_key] = exec_id

        # Also key by base gccode prefix (e.g. GC184872) as fallback
        base_gc = sample_key.split('-')[0]
        if base_gc not in exec_map:
            exec_map[base_gc] = exec_id

    return exec_map


def download_missing_file(sample_key, local_dir, exec_map):
    """
    Downloads a single sample file from GCS if it does not exist locally.
    Uses exec_map to determine the dynamic exec_id in GCS path.
    """
    gz_file = os.path.join(local_dir, f"{sample_key}.partition.bin.counts.gcc.mean.tsv.gz")
    tsv_file = os.path.join(local_dir, f"{sample_key}.partition.bin.counts.gcc.mean.tsv")

    # Return early if file already exists locally
    if os.path.exists(gz_file):
        return gz_file
    if os.path.exists(tsv_file):
        return tsv_file

    # Lookup exec_id from mapping
    exec_id = exec_map.get(sample_key)
    if not exec_id:
        print(f"Error: Could not find exec_id for sample '{sample_key}' in mapping file.")
        return None

    # Parse base gccode for GCS path directory structure
    base_gccode = sample_key.split('-')[0]  # e.g. GC184872-I4A02 -> GC184872
    gc_prefix_5 = base_gccode[:5]           # e.g. GC184
    gc_prefix_8 = base_gccode[:8]           # e.g. GC184872

    # Construct GCS path with target exec_id from previous script
    gcs_path = (
        f"gs://dev-gcpi-ducn/{gc_prefix_5}/{gc_prefix_8}/{exec_id}/"
        f"partition.bins/bin.counts/gcc.mean/CM_N/{sample_key}.partition.bin.counts.gcc.mean.tsv.gz"
    )

    print(f"Downloading missing file from GCS: {gcs_path}")
    try:
        subprocess.run(["gsutil", "cp", gcs_path, gz_file], check=True)
        if os.path.exists(gz_file):
            return gz_file
    except Exception as e:
        print(f"Error downloading {sample_key} from GCS: {e}")

    return None


def combine_gcc_per_sex(csv_file, input_arg):
    """
    1. Loads exec_map from CSV mapping file.
    2. Reads sample ID list files for XX and XY.
    3. Downloads files using mapping-based exec_ids.
    4. Merges TSVs on BINDEX and exports AVG_COUNT files.
    """
    base_dir = os.getcwd()
    local_dir = os.path.join(base_dir, "inputs/bin_counts/gcc_mean/CM_N/gccode_genome")
    output_dir = os.path.join(base_dir, "inputs/bin_counts/gcc_mean/CM_N/sex_genome")
    os.makedirs(local_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # Load exec_map mapping
    exec_map = load_exec_map(csv_file)

    # Determine filename prefix and suffix based on input_arg
    if "XX" in input_arg:
        prefix_str, suffix_str = input_arg.split("XX", 1)
    elif "XY" in input_arg:
        prefix_str, suffix_str = input_arg.split("XY", 1)
    else:
        prefix_str = ""
        suffix_str = f".{input_arg}.hg38"

    for sex in ['XX', 'XY']:
        print(f'\nProcessing {sex}')
        list_filename = f"{prefix_str}{sex}{suffix_str}"
        file_gccodes = os.path.join(base_dir, list_filename)

        if not os.path.exists(file_gccodes):
            print(f"Warning: List file '{file_gccodes}' not found. Skipping {sex}.")
            continue

        dfs = []
        with open(file_gccodes) as file:
            for line in file:
                sample_key = line.strip()
                if not sample_key:
                    continue

                target_file = download_missing_file(sample_key, local_dir, exec_map)
                if not target_file:
                    print(f"Skipping {sample_key}: File not found locally or on GCS.")
                    continue

                print(f"Processing {sex} input: {target_file}")
                try:
                    df = pd.read_csv(target_file, sep="\t", compression='infer')
                    dfs.append(df)
                except Exception as e:
                    print(f"Error reading {target_file}: {e}")

        if not dfs:
            print(f"No valid data loaded for {sex}, skipping.")
            continue

        print(f"Loaded {len(dfs)} files for {sex}. Merging...")

        # Merge DataFrames on 'BINDEX'
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = merged_df.merge(df, on="BINDEX", how="outer", suffixes=('', '_drop'))

        # Calculate average count across all sample columns (excluding BINDEX)
        count_columns = [col for col in merged_df.columns if col != "BINDEX" and not col.endswith('_drop')]
        merged_df["AVG_COUNT"] = merged_df[count_columns].mean(axis=1)

        # Export final BINDEX and AVG_COUNT
        final_df = merged_df[["BINDEX", "AVG_COUNT"]]
        out_base_name = list_filename.rstrip(".hg38") if list_filename.endswith(".hg38") else list_filename
        output_file = os.path.join(output_dir, f"{out_base_name}.hg38.genome")
        final_df.to_csv(output_file, sep="\t", index=False)
        print(f"Averaged file saved as '{output_file}'")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python3 lps_hg38_pair.py <mapping_csv> [sex_list_arg]")
        print("\nExample:")
        print("  python3 lps_hg38_pair.py gc_exec.csv lpspost.XX.novaseqxplus.hg38")
        sys.exit(1)

    csv_file = sys.argv[1]
    input_arg = sys.argv[2] if len(sys.argv) > 2 else "lpspost.XX.novaseqxplus.hg38"

    combine_gcc_per_sex(csv_file, input_arg)
