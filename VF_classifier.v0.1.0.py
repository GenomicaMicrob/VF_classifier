#!/usr/bin/env python3

import argparse
import subprocess
import re
import pandas as pd
from Bio import SeqIO
import sys
import os
import gzip
import shutil
import urllib.request
import urllib.error
import datetime
import tarfile
from pathlib import Path

# Try importing plotting libraries
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False

# --- Colors & Icons ---
class Msg:
    HEADER = '\033[95m\033[1m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    END = '\033[0m'

    @staticmethod
    def info(text):
        print(f"\n{Msg.BLUE}{text}{Msg.END}")
    
    @staticmethod
    def success(text):
        print(f"{Msg.GREEN}{text}{Msg.END}")
    
    @staticmethod
    def processing(text):
        print(f"{Msg.CYAN}{text}{Msg.END}")
    
    @staticmethod
    def file(text):
        print(f"  {Msg.GREEN}{text}{Msg.END}")
    
    @staticmethod
    def warn(text):
        print(f"{Msg.YELLOW}[WARNING] {text}{Msg.END}")
    
    @staticmethod
    def error(text):
        print(f"\n{Msg.RED}{Msg.BOLD}[ERROR] {text}{Msg.END}")

    @staticmethod
    def spinner_hook(count, block_size, total_size):
        """Callback for urllib.request.urlretrieve to show a rotating spinner."""
        if total_size <= 0: return
        
        # Determine progress spin
        chars = ['/', '-', '\\', '|']
        idx = count % len(chars)
        
        percent = min(100, int(count * block_size * 100 / total_size))
        # Clear line and print spinner
        sys.stdout.write(f"\r  {Msg.CYAN}{chars[idx]} Downloading... {percent}% complete {Msg.END}")
        sys.stdout.flush()
        
        if percent >= 100:
            sys.stdout.write("\n")
            sys.stdout.flush()

# --- 0. Constants ---
BASE_URL = "https://www.mgc.ac.cn/VFs"
PAGE_URL = f"{BASE_URL}/download.htm"
DNA_URL = f"{BASE_URL}/Down/VFDB_setB_nt.fas.gz"
DEFAULT_DB_DIR = "VF_database"
DEFAULT_DB_NAME = "VFDB_db"
METADATA_FILE = "VFDB_metadata.txt"

# --- 1. VFDB Downloader & Manager ---
class VFDBManager:
    def __init__(self, db_dir=DEFAULT_DB_DIR):
        self.db_dir = Path(db_dir)
        self.metadata_path = self.db_dir / METADATA_FILE
        self.fasta_gz = self.db_dir / "VFDB_setB_nt.fas.gz"
        self.fasta_unzipped = self.db_dir / "VFDB_setB_nt.fas"
        self.db_prefix = self.db_dir / DEFAULT_DB_NAME

    def get_online_version(self):
        Msg.info(f"Checking for updates at {PAGE_URL}...")
        try:
            with urllib.request.urlopen(PAGE_URL, timeout=10) as response:
                html = response.read().decode('utf-8', errors='ignore')
                # Look for "Last update" string
                match = re.search(r'Last update[: ]*([^<]+)', html, re.IGNORECASE)
                if match:
                    return match.group(1).strip()
        except Exception as e:
            print(f"Warning: Could not check for online updates: {e}")
        return None

    def get_local_version(self):
        if self.metadata_path.exists():
            with open(self.metadata_path, 'r') as f:
                content = f.read()
                match = re.search(r'VFDB Last Update:\s*(.+)', content)
                if match:
                    return match.group(1).strip()
        return None

    def setup_database(self, force=False, local_fasta=None, skip_curation=False):
        # 0. Safety Check: Prevent accidental overwrite
        if self.db_dir.exists() and not force:
            Msg.warn(f"The database directory '{self.db_dir}' already exists.")
            Msg.info("If you want to re-run the setup and overwrite it, please use the --force flag.")
            Msg.info(f"Currently installed version: {self.get_local_version() or 'Unknown'}")
            return False

        online_ver = None
        local_ver = self.get_local_version()

        # Handle Local File Mode
        if local_fasta:
            Msg.info(f"Using local source file for setup: {local_fasta}")
            self.db_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy/Extract based on extension
            if local_fasta.lower().endswith(('.tar.gz', '.tgz', '.tar')):
                Msg.processing(f"Extracting archive {os.path.basename(local_fasta)}...")
                try:
                    with tarfile.open(local_fasta, 'r:*') as tar:
                        # Look for FASTA files inside
                        members = [m for m in tar.getmembers() if m.isfile() and 
                                  (m.name.lower().endswith('.fas') or m.name.lower().endswith('.fasta') or m.name.lower().endswith('.fna'))]
                        
                        if not members:
                            Msg.error("No FASTA files (.fas, .fasta, .fna) found inside the archive.")
                            return False
                        
                        # Use the first one found
                        fasta_member = members[0]
                        Msg.info(f"Extracting: {fasta_member.name}")
                        with tar.extractfile(fasta_member) as f_in, open(self.fasta_unzipped, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                except Exception as e:
                    Msg.error(f"Error extracting archive: {e}")
                    return False
            elif local_fasta.lower().endswith('.gz'):
                Msg.processing(f"Extracting gzipped file {os.path.basename(local_fasta)}...")
                try:
                    with gzip.open(local_fasta, 'rb') as f_in:
                        with open(self.fasta_unzipped, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                except Exception as e:
                    Msg.error(f"Error extracting local file: {e}")
                    return False
            else:
                Msg.processing(f"Copying {os.path.basename(local_fasta)}...")
                try:
                    shutil.copy2(local_fasta, self.fasta_unzipped)
                except Exception as e:
                    Msg.error(f"Error copying local file: {e}")
                    return False
            
            # Simple validation: is it really a FASTA?
            try:
                with open(self.fasta_unzipped, 'rb') as f:
                    first_byte = f.read(1)
                    if first_byte != b'>':
                        Msg.error("The extracted file does not appear to be a valid FASTA (missing '>' at start).")
                        Msg.info("Please ensure you are providing the VFDB nucleotide sequences.")
                        return False
            except Exception:
                pass
            
            online_ver = "Local_Import" # Placeholder for metadata
            
        else:
            # Standard Download Mode
            online_ver = self.get_online_version()
            
            if online_ver:
                Msg.info(f"Downloading VFDB version: {online_ver}")
            else:
                Msg.warn("Proceeding with download (could not verify online version).")

            self.db_dir.mkdir(parents=True, exist_ok=True)

            # 1. Download
            Msg.info("Downloading VFDB source files...")
            Msg.info("Be patient, this might take some time...")
            try:
                urllib.request.urlretrieve(DNA_URL, self.fasta_gz, reporthook=Msg.spinner_hook)
            except Exception as e:
                Msg.error(f"Error downloading VFDB: {e}")
                return False

            # 2. Extract
            Msg.processing(f"Extracting {self.fasta_gz.name}...")
            try:
                with gzip.open(self.fasta_gz, 'rb') as f_in:
                    with open(self.fasta_unzipped, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            except Exception as e:
                Msg.error(f"Error extracting database: {e}")
                return False

        # --- Following steps are common to both Local and Download mods ---

        # 3. Source-Level Curation
        if not skip_curation:
            self.curate_fasta(self.fasta_unzipped)
        else:
            Msg.warn("Skipping automated curation as requested. Using original VFDB headers.")

        # 4. Create BLAST DB
        title_suffix = "(Curated)" if not skip_curation else "(Original Headers)"
        Msg.processing(f"Creating BLAST database from {'curated ' if not skip_curation else ''}source...")
        cmd = [
            "makeblastdb",
            "-in", str(self.fasta_unzipped),
            "-dbtype", "nucl",
            "-out", str(self.db_prefix),
            "-parse_seqids",
            "-title", f"VFDB_{online_ver or 'Unknown'} {title_suffix}"
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except FileNotFoundError:
            Msg.error("'makeblastdb' command not found. Please install NCBI BLAST+.")
            return False
        except subprocess.CalledProcessError as e:
            Msg.error(f"Error running makeblastdb: {e.stderr.decode()}")
            return False

        # 5. Write Metadata
        with open(self.metadata_path, 'w') as f:
            f.write("VFDB Database Metadata\n")
            f.write("====================\n")
            f.write(f"Created Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"VFDB Last Update: {online_ver or 'Unknown'}\n")
            if not local_fasta:
                f.write(f"Source URL: {DNA_URL}\n")
            else:
                f.write(f"Source: Local File ({local_fasta})\n")
            
            status_str = "Automated Curation Applied (Longest Name Strategy)" if not skip_curation else "Source Headers Uncurated"
            f.write(f"Status: {status_str}\n")

        Msg.success(f"Database ready in '{self.db_dir}'")
        return True

    def curate_fasta(self, fasta_path):
        """Rewrites the FASTA file headers using the Longest Name Strategy for Category/Specific Name only."""
        Msg.processing("Curating FASTA headers (Longest Name Strategy for Specific Names)...")
        
        raw_records = []
        curated_cats = {}
        
        # Pass 1: Identify longest category names per VF_ID
        for record in SeqIO.parse(fasta_path, "fasta"):
            blast_id = record.id
            full_desc = record.description
            
            # Extract basic info
            gene_match = re.search(r'\s\(([^)]+)\)\s', full_desc)
            gene = gene_match.group(1) if gene_match else "N/A"
            
            cat_match = re.search(r'\[(.*?VF.*? - .*?)\]', full_desc)
            full_cat = cat_match.group(1) if cat_match else "N/A"
            
            # Extract internal VF ID (the common key)
            vf_id_match = re.search(r'\(VF\d+\)', full_cat)
            vf_id = vf_id_match.group(0) if vf_id_match else blast_id
            
            # Track longest category strings ONLY
            if vf_id not in curated_cats or len(full_cat) > len(curated_cats[vf_id]):
                curated_cats[vf_id] = full_cat
            
            raw_records.append((record, vf_id, full_cat))
            
        # Pass 2: Rewrite FASTA and log changes
        curation_log = []
        curated_records = []
        
        for record, vf_id, orig_cat in raw_records:
            new_cat = curated_cats[vf_id]
            
            # Log and replace ONLY Category/Specific Name
            if orig_cat != new_cat and orig_cat != "N/A":
                curation_log.append((vf_id, "VF_Specific_Name/Category", orig_cat, new_cat))
                record.description = record.description.replace(orig_cat, new_cat)
            
            curated_records.append(record)
            
        # Write curated FASTA
        temp_curated = f"{fasta_path}.curated"
        SeqIO.write(curated_records, temp_curated, "fasta")
        os.replace(temp_curated, fasta_path)
        
        # Write Log
        if curation_log:
            unique_log = sorted(list(set(curation_log)))
            log_path = Path(self.db_dir) / "curation_log.tsv"
            try:
                with open(log_path, 'w') as f:
                    f.write("VF_ID\tField\tOriginal_Value\tCurated_Value\n")
                    for item in unique_log:
                        f.write(f"{item[0]}\t{item[1]}\t{clean_for_excel(item[2])}\t{clean_for_excel(item[3])}\n")
                Msg.success(f"Applied automated curation to {len(unique_log)} unique category entries.")
                Msg.file(f"Curation log saved to: {log_path}")
            except Exception as e:
                Msg.warn(f"Could not save curation log: {e}")

# --- 2. Database Check ---
def get_vfdb_fasta_path(db_path, user_provided_path=None):
    if not (os.path.exists(db_path + ".nin") or os.path.exists(db_path + ".nhr")):
        Msg.error(f"BLAST Database not found at: {db_path}")
        print("-" * 60)
        print("To create this database automatically, run:")
        print(f"  {sys.argv[0]} --setup")
        print("\nOr download VFDB_setB_nt.fas and run manually:")
        print(f"  makeblastdb -in VFDB_setB_nt.fas -dbtype nucl -out {db_path}")
        print("-" * 60)
        sys.exit(1)

    if user_provided_path:
        if os.path.exists(user_provided_path):
            return user_provided_path
        else:
            sys.exit(f"Error: The provided FASTA file was not found: {user_provided_path}")

    db_dir = os.path.dirname(db_path)
    if not db_dir: db_dir = "."
        
    default_filename = "VFDB_setB_nt.fas"
    auto_path = os.path.join(db_dir, default_filename)

    if os.path.exists(auto_path):
        Msg.info(f"Auto-detected VFDB FASTA: {auto_path}")
        return auto_path
    else:
        Msg.error(f"Could not find the original FASTA file at: {auto_path}")
        sys.exit(1)

# --- 2. Indexing VFDB Headers ---
def parse_vfdb_headers(vfdb_fasta):
    Msg.processing(f"Indexing curated headers from {vfdb_fasta}...")
    
    header_map = {}
    
    # Handle both compressed and decompressed files
    if vfdb_fasta.endswith('.gz'):
        handle = gzip.open(vfdb_fasta, "rt")
    else:
        handle = open(vfdb_fasta, "r")
    
    try:
        for record in SeqIO.parse(handle, "fasta"):
            blast_id = record.id 
            full_desc = record.description
            
            # Extract Gene
            gene_match = re.search(r'\s\(([^)]+)\)\s', full_desc)
            gene = gene_match.group(1) if gene_match else "N/A"
            
            # Extract Categories
            cat_match = re.search(r'\[(.*?VF.*? - .*?)\]', full_desc)
            full_category_string = cat_match.group(1) if cat_match else "N/A"
            
            # Extract Organism
            org_match = re.findall(r'\[([^\]]+)\]', full_desc)
            organism = org_match[-1] if org_match else "N/A"
            
            header_map[blast_id] = {
                'VF_Gene': gene,
                'VF_Full_Category_String': full_category_string,
                'VF_Organism': organism
            }
            
    finally:
        handle.close()

    return header_map

# --- 3. Running BLAST ---
def run_blast(genome_file, db_path, out_file="blast_results.tsv", threads=4):
    Msg.info(f"Running BLASTN (Genome: {genome_file} vs DB: {db_path})...")
    
    # Handle compressed input files
    query_file = genome_file
    temp_query = None
    
    if genome_file.endswith('.gz'):
        temp_query = genome_file.replace('.gz', '_temp.fna')
        Msg.processing(f"Decompressing {genome_file} to {temp_query}...")
        with gzip.open(genome_file, 'rt') as gz_file:
            with open(temp_query, 'w') as temp_f:
                temp_f.write(gz_file.read())
        query_file = temp_query
    
    cmd_blast = [
        "blastn", 
        "-query", query_file,
        "-db", db_path,
        "-out", out_file,
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sstrand",
        "-evalue", "1e-5",
        "-num_threads", str(threads)
    ]
    
    try:
        subprocess.run(cmd_blast, check=True)
        Msg.success("BLAST finished.")
    except subprocess.CalledProcessError as e:
        Msg.error(f"BLAST Error: {e}")
        sys.exit(1)
    finally:
        # Clean up temporary decompressed file
        if temp_query and os.path.exists(temp_query):
            os.remove(temp_query)

# --- 4. Plotting Function ---
def create_summary_plot(df, output_img_path):
    if not PLOTTING_AVAILABLE:
        Msg.warn("matplotlib/seaborn not installed. Skipping plot generation.")
        return

    Msg.processing(f"Generating summary plot: {output_img_path}...")
    
    category_counts = df['VF_Broad_Category'].value_counts().reset_index()
    category_counts.columns = ['Category', 'Count']
    
    # Smaller figure size and optimized settings
    plt.figure(figsize=(8, 5))
    sns.set_style("whitegrid")
    
    ax = sns.barplot(
        data=category_counts, 
        x='Count', 
        y='Category',
        hue='Category',
        palette='viridis', 
        edgecolor='black',
        legend=False
    )
    
    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=10, fontweight='bold')
    
    plt.title("Virulence Factors by Category (Best Hits)", fontsize=14, fontweight='bold')
    plt.xlabel("Number of Unique Genes Detected", fontsize=12)
    plt.ylabel("", fontsize=12)
    
    if not category_counts.empty:
        max_count = category_counts['Count'].max()
        plt.xlim(0, max_count * 1.15) 
    
    plt.tight_layout()
    # Optimized save with smaller DPI and compression
    plt.savefig(output_img_path, dpi=200, bbox_inches='tight')
    plt.close()

# --- 5. Krona Generation ---
def create_krona_chart(df, krona_txt, krona_html):
    Msg.processing(f"Generating Krona chart: {krona_html}...")
    
    # Matches: cut -f 5-6 (VF_Broad_Category, VF_Specific_Name)
    hierarchy_cols = ["VF_Broad_Category", "VF_Specific_Name"]
    
    krona_df = df.groupby(hierarchy_cols).size().reset_index(name='Count')
    krona_export = krona_df[['Count'] + hierarchy_cols]
    
    krona_export.to_csv(krona_txt, sep='\t', index=False, header=True)
    
    cmd_krona = [
        "ktImportText",
        krona_txt,
        "-o", krona_html,
        "-n", "Main virulence factor"
    ]
    
    try:
        subprocess.run(cmd_krona, check=True)
        Msg.file(f"Krona HTML saved: {krona_html}")
        # Note: data4plot.tsv file is preserved for user use
    except FileNotFoundError:
        Msg.warn("'ktImportText' not found. Is Krona installed? Skipping HTML generation.")
    except subprocess.CalledProcessError as e:
        Msg.error(f"Krona generation failed: {e}")

def get_sample_base_name(filename):
    """Strip common extensions to get a clean sample name."""
    base_name = os.path.basename(filename)
    for ext in ['.fna.gz', '.fasta.gz', '.fa.gz', '.fna', '.fasta', '.fa']:
        if base_name.lower().endswith(ext):
            return base_name[:-len(ext)]
    return os.path.splitext(base_name)[0]

def clean_for_excel(text):
    if not text or pd.isna(text): return text
    text = str(text).strip()
    # Excel triggers formulas on =, +, -, @
    if text.startswith(('=', '+', '-', '@')):
        return f"'{text}"
    return text

# --- 6. Processing & Saving ---
def process_results(blast_output, header_map, input_filename):
    Msg.processing(f"Processing results for: {os.path.basename(input_filename)}...")
    
    # 1. Get clean sample name
    base_name = get_sample_base_name(input_filename)

    # 2. Create Sample-Specific Output Directory
    master_results_dir = "VF_results"
    sample_output_dir = os.path.join(master_results_dir, base_name)
    os.makedirs(sample_output_dir, exist_ok=True)
    
    # 3. Define Output Paths inside the subfolder
    out_all = os.path.join(sample_output_dir, f"{base_name}_all_hits.tsv")
    out_best = os.path.join(sample_output_dir, f"{base_name}_best_hits.tsv")
    out_plot = os.path.join(sample_output_dir, f"{base_name}_summary.png")
    out_krona_txt = os.path.join(sample_output_dir, f"{base_name}_data4plot.tsv")
    out_krona_html = os.path.join(sample_output_dir, f"{base_name}_VF.html")

    cols = ["Query_Sequence", "VF_ID", "Identity_Pct", "Align_Len", "Query_Start", "Query_End", "Subj_Start", "Subj_End", "E_Value", "Bit_Score", "Strand"]
    
    try:
        df = pd.read_csv(blast_output, sep="\t", names=cols)
    except pd.errors.EmptyDataError:
        Msg.warn(f"No hits found for {base_name}! Skipping file generation.")
        return None
    except FileNotFoundError:
        Msg.error("BLAST output file was not created.")
        return None

    df['VF_Gene_Symbol'] = df['VF_ID'].map(lambda x: header_map.get(x, {}).get('VF_Gene', 'Unknown'))
    df['VF_Origin_Organism'] = df['VF_ID'].map(lambda x: header_map.get(x, {}).get('VF_Organism', 'Unknown'))
    df['Temp_Category'] = df['VF_ID'].map(lambda x: header_map.get(x, {}).get('VF_Full_Category_String', 'Unknown'))
    
    split_data = df['Temp_Category'].str.split(' - ', n=1, expand=True)
    df['VF_Specific_Name'] = split_data[0]
    df['VF_Broad_Category'] = split_data[1] if len(split_data.columns) > 1 else "Unspecified"

    df['VF_Specific_Name'] = df['VF_Specific_Name'].apply(clean_for_excel)
    df['VF_Broad_Category'] = df['VF_Broad_Category'].apply(clean_for_excel)
    df['VF_Gene_Symbol'] = df['VF_Gene_Symbol'].apply(clean_for_excel)

    final_cols = ["Query_Sequence", "Query_Start", "Query_End", "Strand", "VF_Broad_Category", "VF_Specific_Name", "VF_Gene_Symbol", "Identity_Pct", "E_Value", "VF_Origin_Organism", "VF_ID"]

    df_all = df.sort_values(by=["Query_Sequence", "Bit_Score"], ascending=[True, False])
    df_all[final_cols].to_csv(out_all, sep='\t', index=False, encoding='utf-8-sig')
    Msg.file(f"Saved All Hits:  {out_all}")

    # Filter for Best Hits
    df_best = df_all.drop_duplicates(subset=["Query_Sequence", "VF_Gene_Symbol"], keep="first").copy()
    df_best[final_cols].to_csv(out_best, sep='\t', index=False, encoding='utf-8-sig')
    Msg.file(f"Saved Best Hits: {out_best}")

    # Visualizations
    create_summary_plot(df_best, out_plot)
    
    # Generate Krona
    create_krona_chart(df_best, out_krona_txt, out_krona_html)
    
    return out_best

# --- 7. Merging Results ---
def merge_results(file_list=None, output_path=None):
    output_dir = "VF_results"
    
    if not file_list:
        # Auto-scan VF_results RECURSIVELY
        if not os.path.exists(output_dir):
            Msg.error(f"Directory '{output_dir}' not found. Run analysis first.")
            return False
        
        # Use pathlib for easy recursive globbing
        from pathlib import Path
        file_list = [str(p) for p in Path(output_dir).rglob("*_best_hits.tsv")]
        
    if not file_list:
        Msg.error("No '*_best_hits.tsv' files found to merge.")
        return False

    if not output_path:
        output_path = os.path.join(output_dir, "multiple_samples_best_hits.csv")

    Msg.processing(f"Merging {len(file_list)} sample(s) into matrix...")
    
    all_data = []
    for f in file_list:
        # Explicit check to ignore .csv files (safety)
        if f.lower().endswith('.csv'):
            continue
            
        try:
            # Get sample name from filename or parent folder
            base = os.path.basename(f)
            sample_name = base.replace("_best_hits.tsv", "")
            
            df = pd.read_csv(f, sep='\t')
            
            if df.empty:
                Msg.warn(f"Skipping empty file: {f}")
                continue
                
            # Keep Gene, Identity, and Categories
            cols_to_keep = ['VF_Broad_Category', 'VF_Specific_Name', 'VF_Gene_Symbol', 'Identity_Pct']
            # Safety check if columns exist
            existing_cols = [c for c in cols_to_keep if c in df.columns]
            temp_df = df[existing_cols].copy()
            temp_df['Sample'] = sample_name
            all_data.append(temp_df)
            Msg.file(f"Loaded {sample_name}")
        except Exception as e:
            Msg.error(f"Failed to read {f}: {e}")

    if not all_data:
        Msg.error("No valid data found to merge.")
        return False

    # Combine
    combined = pd.concat(all_data, ignore_index=True)
    
    # Pivot: Index with categories to facilitate clustered heatmaps
    index_cols = ['VF_Broad_Category', 'VF_Specific_Name', 'VF_Gene_Symbol']
    active_index = [c for c in index_cols if c in combined.columns]
    
    matrix = combined.pivot_table(
        index=active_index, 
        columns='Sample', 
        values='Identity_Pct', 
        aggfunc='max'
    )
    
    # Fill NAs with 0
    matrix = matrix.fillna(0)
    
    # Save
    matrix.to_csv(output_path)
    Msg.success(f"Successfully created clustered matrix: {output_path}")
    Msg.info(f"Matrix shape: {matrix.shape[0]} genes x {matrix.shape[1]} samples")
    
    return True

# --- 8. Main ---
def main():
    description_text = """
    Virulence Factor BLAST Pipeline (v19)
    -------------------------------------
    1. Runs BLASTN against VFDB.
    2. Outputs results into 'VF_results/'.
    
    To initialize or update the database:
      %(prog)s --setup

    To compile multiple results into a heatmap matrix:
      %(prog)s --merge
    """
    
    parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
    
    # Setup Group
    setup_group = parser.add_argument_group("Database Setup")
    setup_group.add_argument("--setup", action="store_true", help="Download VFDB, curate it, and create BLAST database.")
    setup_group.add_argument("--input-fasta", help="Optional: Use a locally downloaded VFDB FASTA file instead of downloading.")
    setup_group.add_argument("--no-curation", action="store_true", help="Skip automated curation and use original VFDB headers.")
    setup_group.add_argument("--db-dir", default=DEFAULT_DB_DIR, help=f"Directory for database (default: {DEFAULT_DB_DIR})")

    # Analysis Group
    input_group = parser.add_argument_group("Analysis")
    input_group.add_argument("-i", "--input", nargs='+', help="Path to input genome FASTA file(s). Supports multiple files.")
    input_group.add_argument("--force", action="store_true", help="Force database re-download OR overwrite existing sample results.")
    input_group.add_argument("-db", "--database", help=f"Path to BLAST DB prefix (default: {DEFAULT_DB_DIR}/{DEFAULT_DB_NAME})")
    input_group.add_argument("-v", "--vfdb_fasta", help="Optional: Path to VFDB .fas file.")
    input_group.add_argument("-t", "--threads", default=4, type=int, help="Number of CPU threads.")
    
    # Merge Group
    merge_group = parser.add_argument_group("Results Consolidation")
    merge_group.add_argument("--merge", nargs='*', help="Consolidate multiple '_best_hits.tsv' files into a single matrix. If no files are provided, it scans the 'VF_results/' directory recursively.")
    merge_group.add_argument("--out-matrix", help="Custom output path for the merged CSV matrix (default: VF_results/multiple_samples_best_hits.csv)")

    args = parser.parse_args()

    # Handle Merge (Manual Call)
    if args.merge is not None:
        merge_results(file_list=args.merge, output_path=args.out_matrix)
        sys.exit(0)

    # Handle Setup
    if args.setup:
        manager = VFDBManager(db_dir=args.db_dir)
        success = manager.setup_database(force=args.force, local_fasta=args.input_fasta, skip_curation=args.no_curation)
        sys.exit(0 if success else 1)

    # Validate Analysis Inputs
    if not args.input:
        parser.print_help()
        sys.exit("\nError: --input is required for analysis.")
    
    db_path = args.database if args.database else os.path.join(args.db_dir, DEFAULT_DB_NAME)
    vfdb_fasta_path = get_vfdb_fasta_path(db_path, args.vfdb_fasta)
    header_data = parse_vfdb_headers(vfdb_fasta_path)
    
    # Batch Processing
    processed_count = 0
    successful_best_hits = []

    for genome_file in args.input:
        try:
            sample_name = get_sample_base_name(genome_file)
            sample_out_dir = os.path.join("VF_results", sample_name)
            
            # Safety Check: Prevent accidental result overwrite
            if os.path.exists(sample_out_dir) and not args.force:
                print("\n" + "-"*60)
                Msg.warn(f"Skipping sample '{sample_name}': Output directory already exists.")
                Msg.info(f"Directory: {sample_out_dir}")
                Msg.info("Use the --force flag if you wish to overwrite these results.")
                print("-"*60)
                continue

            print("\n" + "="*60)
            Msg.info(f"ANALYZING SAMPLE: {os.path.basename(genome_file)}")
            print("="*60)
            
            # Save temp file in current dir
            blast_temp_file = f"temp_blast_{os.path.basename(genome_file)}.tsv"
            
            run_blast(genome_file, db_path, out_file=blast_temp_file, threads=args.threads)
            
            best_hits_path = process_results(blast_temp_file, header_data, genome_file)
            
            if best_hits_path:
                successful_best_hits.append(best_hits_path)
                processed_count += 1
            
            if os.path.exists(blast_temp_file):
                os.remove(blast_temp_file)
                
        except Exception as e:
            Msg.error(f"Failed to process {genome_file}: {e}")

    # Auto-Merge if more than one sample processed successfully
    if processed_count > 1:
        Msg.info(f"Batch complete. Auto-merging {processed_count} samples...")
        merge_results(file_list=successful_best_hits, output_path=args.out_matrix)
    elif processed_count == 1:
        Msg.success(f"Analysis complete for 1 sample.")
    else:
        Msg.error("No samples were successfully processed.")

if __name__ == "__main__":
    main()