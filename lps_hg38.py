import os
import sys
import subprocess
import pandas as pd
def download_missing_file(gccode, local_dir):
    """Downloads a single sample file from GCS if it does not exist locally."""
    gz_file = os.path.join(local_dir, f"{gccode}.partition.bin.counts.gcc.mean.tsv.gz")
    tsv_file = os.path.join(local_dir, f"{gccode}.partition.bin.counts.gcc.mean.tsv")
    if os.path.exists(gz_file):
        return gz_file
    if os.path.exists(tsv_file):
        return tsv_file
    # GCS source path (without extra gccode_genome folder)
    gcs_path = f"gs://gipseq-clusters/pools/36bp-hg38/{gccode[:5]}/{gccode[:8]}/partition.bins/bin.counts/gcc.mean/CM_N/{gccode}.partition.bin.counts.gcc.mean.tsv.gz"
    print(f"Downloading missing file from GCS: {gcs_path}")
    try:
        # Pass full file destination path (gz_file) to prevent gsutil directory errors
        subprocess.run(["gsutil", "cp", gcs_path, gz_file], check=True)
        if os.path.exists(gz_file):
            return gz_file
    except Exception as e:
        print(f"Error downloading {gccode} from GCS: {e}")
    return None
def combine_gcc_per_sex(input_arg):
    base_dir = os.getcwd()
    local_dir = os.path.join(base_dir, "inputs/bin_counts/gcc_mean/CM_N/gccode_genome")
    prefix = os.path.join(base_dir, "inputs/bin_counts/gcc_mean/CM_N/sex_genome")
    os.makedirs(local_dir, exist_ok=True)
    os.makedirs(prefix, exist_ok=True)
    # Determine filename prefix and suffix based on input argument
    if "XX" in input_arg:
        prefix_str, suffix_str = input_arg.split("XX", 1)
    elif "XY" in input_arg:
        prefix_str, suffix_str = input_arg.split("XY", 1)
    else:
        prefix_str = ""
        suffix_str = f".{input_arg}.hg38"

    for sex in ['XX', 'XY']:
        print(f'Processing {sex}')
        list_filename = f"{prefix_str}{sex}{suffix_str}"
        file_gccodes = os.path.join(base_dir, list_filename)
        if not os.path.exists(file_gccodes):
            print(f"Warning: List file '{file_gccodes}' not found. Skipping {sex}.")
            continue
        dfs = []
        with open(file_gccodes) as file:
            for line in file:
                gccode = line.strip()
                if not gccode:
                    continue
                target_file = download_missing_file(gccode, local_dir)
                if not target_file:
                    print(f"Skipping {gccode}: File not found locally or on GCS.")
                    continue
                print(f'Processing {sex} input {target_file}')
                try:
                    df = pd.read_csv(target_file, sep="\t", compression='infer')
                    dfs.append(df)
                except Exception as e:
                    print(f"Error reading {target_file}: {e}")
        if not dfs:
            print(f"No valid data loaded for {sex}, skipping.")
            continue
        print(f"Loaded {len(dfs)} files for {sex}")
        # Merge DataFrames on 'BINDEX'
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = merged_df.merge(df, on="BINDEX", how="outer", suffixes=('', '_drop'))
        # Calculate average count across all samples
        count_columns = [col for col in merged_df.columns if col != "BINDEX"]
        merged_df["AVG_COUNT"] = merged_df[count_columns].mean(axis=1)
        final_df = merged_df[["BINDEX", "AVG_COUNT"]]
        out_base_name = list_filename.rstrip(".hg38") if list_filename.endswith(".hg38") else list_filename
        output_file = os.path.join(prefix, f"{out_base_name}.hg38.genome")
        final_df.to_csv(output_file, sep="\t", index=False)
        print(f"Averaged file saved as '{output_file}'")
if __name__ == "__main__":
    input_arg = sys.argv[1] if len(sys.argv) > 1 else "lpspost.XX.novaseqxplus.hg38"
    combine_gcc_per_sex(input_arg)

