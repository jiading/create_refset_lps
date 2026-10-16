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
    
    Parameters:
    - sample_key: Sample identifier, e.g. 'GC184872-I4A02' or 'GC184872'
    - local_dir: Destination directory path
    - exec_map: Dictionary returned by load_exec_map()
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

    # Construct GCS path with target exec_id
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


if __name__ == "__main__":
    import sys

    # Default settings
    OUTPUT_DIR = "./inputs"
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Check command-line arguments
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Process all samples in CSV:   python3 lps_hg38_pair.py gc_exec.csv")
        print("  Process a single sample:      python3 lps_hg38_pair.py gc_exec.csv GC184872-I4A02")
        sys.exit(1)

    csv_file = sys.argv[1]

    # Load mapping from the CSV passed as arg 1
    exec_map = load_exec_map(csv_file)

    # Case 1: Specific sample key provided as 2nd argument
    if len(sys.argv) >= 3:
        target_sample = sys.argv[2]
        print(f"Processing single sample: {target_sample}")
        download_missing_file(target_sample, local_dir=OUTPUT_DIR, exec_map=exec_map)

    # Case 2: Only CSV provided -> process all samples in the CSV
    else:
        print(f"Processing all samples from '{csv_file}'...")
        for sample_key in exec_map.keys():
            # Process actual sample keys containing hyphens (e.g. GC184872-I4A02)
            if '-' in sample_key:
                download_missing_file(sample_key, local_dir=OUTPUT_DIR, exec_map=exec_map)
