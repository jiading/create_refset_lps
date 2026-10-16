
import os
import sys
import subprocess
import pandas as pd


def sort_integer_file(file_to_sort):
  # Read the integers from the file
  with open(file_to_sort, 'r') as infile:
    numbers = [int(line.strip()) for line in infile if line.strip()]  # Ignore empty lines

  # Sort the numbers numerically
  numbers.sort()

  # Write the sorted numbers to a new file
  with open(file_to_sort, 'w') as outfile:
    for number in numbers:
        outfile.write(f"{number}\n")


# USED FOR BOTH REFSET BUILDING AND VALIDATION SAMPLE PREPARING
def localize_gccode_bin_counts_gcc_mean_cm_n(file_gccodes):
  # Open a file with one gccode per line (e.g., GC157745-I4E08)
  with open(file_gccodes) as file:
    for line in file:
      gccode = line.strip()
      # Construct the desired GCS path
      gcs_path = f"gs://gipseq-clusters/pools/36bp-hg38/{gccode[:5]}/{gccode[:8]}/partition.bins/bin.counts/gcc.mean/CM_N/gccode_genome/{gccode}.partition.bin.counts.gcc.mean.tsv.gz"
      local_dir = "/work/LPS/hg38/inputs/bin_counts/gcc_mean/CM_N/gccode_genome"
      local_file = f"{local_dir}/{gccode}.partition.bin.counts.gcc.mean.tsv.gz"
      unzipped_file = local_file.rstrip(".gz")  # Remove .gz extension

      # Ensure local directory exists
      os.makedirs(local_dir, exist_ok=True)

      # Copy file from GCS
      subprocess.run(["gsutil", "-m", "cp", gcs_path, local_dir], check=True)

      # Unzip the file
      subprocess.run(["gzip", "-d", local_file], check=True)  # or use "gunzip"

      print(f"Copied and extracted: {unzipped_file}")


def combine_gcc_per_sex(machine):
    local_dir = "/work/LPS/hg38/inputs/bin_counts/gcc_mean/CM_N/gccode_genome"
    prefix = "/work/LPS/hg38/inputs/bin_counts/gcc_mean/CM_N/sex_genome"
    
    for sex in ['XX', 'XY']:
        print(f'Processing {sex}')
        file_gccodes = f'/work/LPS/hg38/{sex}.{machine}.hg38'   # XX.novaseq6000=27; XY=33
        dfs = []

        with open(file_gccodes) as file:
            for line in file:
                gccode = line.strip()
                local_file = f"{local_dir}/{gccode}.partition.bin.counts.gcc.mean.tsv"  # 186 files
                print(f'Processing {sex} input {local_file}')

                try:
                    df = pd.read_csv(local_file, sep="\t")  
                    dfs.append(df)
                except Exception as e:
                    print(f"Error reading {local_file}: {e}")

        if not dfs:  # Skip if no data was loaded
            print(f"No data for {sex}, skipping.")
            continue

        print(f"Loaded {len(dfs)} files for {sex}")

        # Step 3: Merge all DataFrames on 'BINDEX' using outer join
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = merged_df.merge(df, on="BINDEX", how="outer", suffixes=('', '_drop'))

        print("Merged DataFrame before dropping '_drop' columns:")
        print(merged_df.head())

        # Step 4: Compute the average COUNT across all COUNT columns (excluding 'BINDEX')
        count_columns = [col for col in merged_df.columns if col != "BINDEX"]
        print(f"Columns used for averaging: {count_columns}")

        merged_df["AVG_COUNT"] = merged_df[count_columns].mean(axis=1)

        # Step 5: Drop the extra "_drop" columns after averaging
        merged_df = merged_df.loc[:, ~merged_df.columns.str.contains('_drop')]

        # Keep only 'BINDEX' and 'AVG_COUNT' columns
        final_df = merged_df[["BINDEX", "AVG_COUNT"]]
        print("Final DataFrame:")
        print(final_df.head())

        # Step 6: Save the result to a new TSV file
        output_file = f"{prefix}/{sex}.{machine}.hg38.genome"
        final_df.to_csv(output_file, sep="\t", index=False)

        print(f"Averaged file saved as '{output_file}'")


def localize_gccode_read_starts(file_gccodes):
  # Open a file with one gccode per line (eg, GC157745-I4E08)
  with open(file_gccodes) as file:
    for line in file:
      gccode = line.strip()
      # Construct the desired GCS path
      prefix = f"gs://gipseq-clusters/pools/36bp-hg38/{gccode[:5]}/{gccode[:8]}/read.starts/seqCBS.inputs/{gccode}/{gccode}.seqCBS.inputs.chr"

      # Localize all chr files
      for chr in range(1,25):
        chrStr = str(chr) if chr < 23 else 'X' if chr==23 else 'Y'
        path=f"{prefix}{chrStr}.txt.gz"
        command = ["gsutil", "-m", "cp",  path , f"/work/LPS/hg38/inputs/gccode_chr/."]
        print(command)
        print(path)

        try:
          result = subprocess.run(command, check=True, text=True, capture_output=True)
          print(result.stdout)
        except subprocess.CalledProcessError as e:
          print(f"An error occurred: {e.stderr}")


# USED FOR REFSET BUILDING  --> 
def combine_per_chr(machine):
  for sex in ['XX', 'XY']:
    file_gccodes = f'/work/LPS/hg38/{sex}.{machine}.hg38'

    for chr in range(1, 25):
        chrStr = str(chr) if chr < 23 else 'X' if chr == 23 else 'Y'
        print(f"Processing chromosome: {chrStr}")
        chr_file = f'/work/LPS/hg38/inputs/sex_chr/{machine}/{sex}.hg38.chr{chrStr}'

        # Open the output file in write mode
        with open(chr_file, 'w') as outfile:
            with open(file_gccodes) as file:
                for line in file:
                    gccode = line.strip()
                    gz_file = f"/work/LPS/hg38/inputs/gccode_chr/{gccode}.seqCBS.inputs.chr{chrStr}.txt.gz"    # 4464/24=186 gccode file
                    
                    # Use zcat to decompress the current gzipped file
                    command = ['zcat', gz_file]
                    
                    # Capture the output of zcat
                    result = subprocess.run(command, capture_output=True, text=True, check=True)
                    
                    # Process each line of the gzipped file
                    for gz_line in result.stdout.splitlines():
                        # Split the line into components and write each on a new line
                        for component in gz_line.split():
                            outfile.write(f"{component}\n")

        sort_integer_file(chr_file)

def combine_autosome_per_gccode(file_gccodes):
    with open(file_gccodes) as file:
      for line in file:
        gccode = line.strip()

        with open(f'/work/LPS/hg38/inputs/gccode_autosome/{gccode}.hg38.autosome', 'w') as outfile:
          for chr in range(1, 23):
            chrStr = str(chr) 
            print(f"Processing chromosome: {gccode} {chrStr}")
            gz_file = f"/work/LPS/hg38/inputs/gccode_chr/{gccode}.seqCBS.inputs.chr{chrStr}.txt.gz"
        
            # Open the output file in write mode
            # Use zcat to decompress the current gzipped file
            command = ['zcat', gz_file]
                 
            # Capture the output of zcat
            result = subprocess.run(command, capture_output=True, text=True, check=True)
                    
            # Process each line of the gzipped file
            for gz_line in result.stdout.splitlines():
                # Split the line into components and write each on a new line
                for component in gz_line.split():
                    outfile.write(f"{int(component)+chr*1000000000}\n")

def combine_genome_per_gccode(file_gccodes):
    with open(file_gccodes) as file:
      for line in file:
        gccode = line.strip()

        with open(f'/work/LPS/hg38/inputs/gccode_genome/{gccode}.hg38.genome', 'w') as outfile:
          for chr in range(1, 25):
            chrStr = str(chr) if chr < 23 else 'X' if chr == 23 else 'Y'
            print(f"Processing chromosome: {gccode} {chrStr}")
            gz_file = f"/work/LPS/hg38/inputs/gccode_chr/{gccode}.seqCBS.inputs.chr{chrStr}.txt.gz"
        
            # Open the output file in write mode
            # Use zcat to decompress the current gzipped file
            command = ['zcat', gz_file]
                 
            # Capture the output of zcat
            result = subprocess.run(command, capture_output=True, text=True, check=True)
                    
            # Process each line of the gzipped file
            for gz_line in result.stdout.splitlines():
                # Split the line into components and write each on a new line
                for component in gz_line.split():
                    outfile.write(f"{int(component)+chr*1000000000}\n")

# assumes combine_autosome_per_gccode done
def combine_autosome_per_type():
  for sample_type in ['amnios','cvs']:
    print(f'Processing {sample_type}')
    file_gccodes = f"/work/LPS/hg38/{sample_type}.hg38"
    type_file = f'/work/LPS/hg38/inputs/type_autosome/{sample_type}.hg38.autosome'
  
    # Open the output file in write mode
    with open(type_file, 'w') as outfile:
      with open(file_gccodes) as file:
        for line in file:
          gccode = line.strip()
          print(f'   Processing {sample_type} {gccode}')
          gccode_autosome_file = f'/work/LPS/hg38/inputs/gccode_autosome/{gccode}.hg38.autosome'

          with open(gccode_autosome_file, "r") as infile:
            outfile.write(infile.read())  # Copy content
            outfile.write("\n")  # Ensure separation between files

    #print(f'   Sorting...')
    #sort_integer_file(type_file)

def combine_genome_per_sex():
  for sex in ['female','male']:
    print(f'Processing {sex}')
    prefix = f'/work/LPS/hg38/inputs/sex_genome/{sex}.hg38'
    sex_file = f'{prefix}.genome'
  
    # Open the output file in write mode
    with open(sex_file, 'w') as outfile:
      for chr in range(1, 25):
        chrStr = str(chr) if chr < 23 else 'X' if chr == 23 else 'Y'
        print(f'    {chrStr}')
        chr_file = f"{prefix}.chr{chrStr}"

        command = ['cat', chr_file]
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        # Process each line 
        for line in result.stdout.splitlines():
          outfile.write(f"{int(line)+chr*1000000000}\n")

    sort_integer_file(sex_file)

def copy_to_cloud(file_gccodes):
  # Open a file with one gccode per line (e.g., GC157745-I4E08)
  with open(file_gccodes) as file:
    for line in file:
      gccode = line.strip()
      print(gccode)

      search_directories = "/work/LPS/hg38/cytosure_input"    # double check here
      command_string = f"find {search_directories} -name '{gccode}*' | grep tar"
      print(command_string)

      try:
        result = subprocess.run(command_string, shell=True, capture_output=True, text=True, check=True)
        found_files = result.stdout.strip().split('\n')

        print("Found files (tar):", [f for f in found_files if f]) # Filter out empty strings
        for file in found_files:
          command_string = f"gsutil cp {file} gs://gcap-resources/tmp/{file_gccodes}/{gccode}/"    #TODO: change the folder name
          subprocess.run(command_string, shell=True, capture_output=True, text=True, check=True)

      except Exception as e:
        print(f"Error reading {file}: {e}")

if __name__ == "__main__":
  file_gccodes = sys.argv[1] #
  #localize_gccode_read_starts(file_gccodes)
  #localize_gccode_bin_counts_gcc_mean_cm_n(file_gccodes)
  combine_gcc_per_sex('novaseqxplus')    # todo: this need to be a parameter in nf
  #combine_per_chr('novaseq6000')
  
  # copy_to_cloud(file_gccodes)  # todo: to fix


  #combine_autosome_per_gccode(file_gccodes)
  #combine_genome_per_gccode(file_gccodes)
  #combine_genome_per_sex()
  #combine_autosome_per_type()
   
  #GENERATE VALIDATION SET
  #combine_genome_per_gccode(file_gccodes)

#SELECTING 100M-50M-25M-10M random read starts
#shuf -n 100000000 /work/LPS/hg38/inputs/sex_genome/male.hg38.genome | sort -n > /work/LPS/hg38/inputs/sex_genome/male.hg38.genome.100M
#shuf -n 100000000 /work/LPS/hg38/inputs/sex_genome/female.hg38.genome | sort -n > /work/LPS/hg38/inputs/sex_genome/female.hg38.genome.100M
#shuf -n 50000000 /work/LPS/hg38/inputs/sex_genome/male.hg38.genome | sort -n > /work/LPS/hg38/inputs/sex_genome/male.hg38.genome.50M
#shuf -n 50000000 /work/LPS/hg38/inputs/sex_genome/female.hg38.genome | sort -n > /work/LPS/hg38/inputs/sex_genome/female.hg38.genome.50M
#shuf -n 10000000 /work/LPS/hg38/inputs/sex_genome/female.hg38.genome.50M | sort -n > /work/LPS/hg38/inputs/sex_genome/female.hg38.genome.10M
#shuf -n 50000000 /work/LPS/hg38/inputs/type_autosome/amnios.hg38.autosome | sort -n > /work/LPS/hg38/inputs/type_autosome/amnios.hg38.autosome.50M
#shuf -n 50000000 /work/LPS/hg38/inputs/type_autosome/cvs.hg38.autosome | sort -n > /work/LPS/hg38/inputs/type_autosome/cvs.hg38.autosome.50M
#cat /work/LPS/hg38/inputs/type_autosome/amnios.hg38.autosome  /work/LPS/hg38/inputs/type_autosome/cvs.hg38.autosome | shuf -n 50000000 | sort -n > /work/LPS/hg38/inputs/type_autosome/mixed.hg38.autosome.50M

