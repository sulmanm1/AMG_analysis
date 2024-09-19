from utils import *
from mmseqs_utils import *
import glob, os, logging, time
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np

def concat_fasta(fastas: list, output_file):
    records = []
    for fasta in fastas:
        file_basename = os.path.basename(fasta).split('.')[0]
        og_file_fasta = read_fasta_no_progress(fasta)

        i = 1
        for record in og_file_fasta:
            record.id = f"{file_basename}||{i}"
            records.append(record)

    write_fasta(output_file, records)

def find_closest_distances(go_cyano_genome_df, vogdb_cyano_df):
    distances = []
    
    for _, go_hit in go_cyano_genome_df.iterrows():
        go_tend = go_hit['tend']
        
        # Calculate the absolute distances between the current go hit and all vogdb hits
        vog_distances = np.abs(vogdb_cyano_df['tstart'] - go_tend)
        
        # Find the minimum distance
        closest_distance = vog_distances.min()
        distances.append(closest_distance)

    return distances

# Function to find distances between consecutive hits in VOGdb
def find_vogdb_distances(vogdb_cyano_df):
    # Sort the hits by tstart to ensure they are in order
    vogdb_cyano_df = vogdb_cyano_df.sort_values(by='tstart')
    
    # Calculate the differences between consecutive tstart values
    distances = np.diff(vogdb_cyano_df['tstart'])
    
    return distances

# Function to plot and save histograms
def plot_histogram(distances, filename, title):
    # Filter the distances to only include values between 0 and 10,000
    filtered_distances = [d for d in distances if 0 <= d <= 10000]
    
    plt.figure(figsize=(10, 6))
    
    # Plot the histogram with density=True to normalize it
    counts, bins, _ = plt.hist(filtered_distances, bins=50, density=True, color='skyblue', edgecolor='black')
    
    # Convert counts to percentages
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y * 100:.0f}%'))
    
    plt.title(title)
    plt.xlabel('Distance')
    plt.ylabel('Percentage')
    plt.xlim(0, 10000)  # Set x-axis limits between 0 and 10,000
    
    # Save the plot as an image
    plt.savefig(filename)
    plt.close()

def prep_amg_data(data_dir, temp_dir):
    bacterial_amg = glob.glob(f"{data_dir}/amgs/bacteria/*")
    viral_amg = glob.glob(f"{data_dir}/amgs/virus/*")
    logging.info(f"Found {len(bacterial_amg)} bacterial AMGs and {len(viral_amg)} viral AMGs")
    bacterial_amg_fasta = f"{temp_dir}/bacterial_amg.fasta"
    viral_amg_fasta = f"{temp_dir}/viral_amg.fasta"
    concat_fasta(bacterial_amg, bacterial_amg_fasta)
    concat_fasta(viral_amg, viral_amg_fasta)

    logging.info("Creating mmseqs databases")
    bacterial_db_path = f"{temp_dir}/bacterial_db"
    viral_db_path = f"{temp_dir}/viral_db"
    mmseqs_makedb(bacterial_amg_fasta, bacterial_db_path)
    mmseqs_makedb(viral_amg_fasta, viral_db_path)
    return bacterial_db_path, viral_db_path

def prep_reference_data(data_dir, temp_dir):
    vogdb_path = f"{data_dir}/vogdb.fa"
    vogdb_db_path = f"{temp_dir}/vogdb_db"
    mmseqs_makedb(vogdb_path, vogdb_db_path)
    
    go_cyano_path = f"{data_dir}/go_cyano_rep.fasta"
    go_cyano_db_path = f"{temp_dir}/go_cyano_db"
    mmseqs_makedb(go_cyano_path, go_cyano_db_path)
    return go_cyano_db_path, vogdb_db_path

def prep_genome_data(data_dir, temp_dir):
    ref_cyanobacteria_genome = f"{data_dir}/ref_cyanobacteria.fna"
    ref_cyanobacteria_genome_db = f"{temp_dir}/ref_cyanobacteria_genome_db"
    mmseqs_makedb(ref_cyanobacteria_genome, ref_cyanobacteria_genome_db)

    ref_cyanophage_genome = f"{data_dir}/ref_cyanophage.fna"
    ref_cyanophage_genome_db = f"{temp_dir}/ref_cyanophage_genome_db"
    mmseqs_makedb(ref_cyanophage_genome, ref_cyanophage_genome_db)
    return ref_cyanobacteria_genome_db, ref_cyanophage_genome_db


def search_to_df(search_result):
    columns = ["query", "target", "pident", "aln_len", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore"]
    df = pd.read_csv(search_result, sep="\t", header=None, names=columns)
    return df

def run(temp_dir, args):
    # Prepare the AMG data
    bacterial_db_path, viral_db_path = prep_amg_data(args.data, temp_dir)

    #prep reference data
    go_cyano_db_path, vogdb_db_path = prep_reference_data(args.data, temp_dir)

    #prep genome data
    ref_cyanobacteria_genome_db, ref_cyanophage_genome_db = prep_genome_data(args.data, temp_dir)


    # Search bacterial and viral AMGs against each other
    search_results_dir = f"{args.output}/bacterial_viral_search"
    os.makedirs(search_results_dir, exist_ok=True)
    df = search_to_df(mmseqs_search(bacterial_db_path, viral_db_path, search_results_dir, f"{temp_dir}/tmp"))

    # Search viral protein db for metabolic genes, remove metabolic genes from vog
    go_vog_search_dir = f"{args.output}/go_vog_search"
    os.makedirs(go_vog_search_dir, exist_ok=True)
    go_vog_df = search_to_df(mmseqs_search(go_cyano_db_path, vogdb_db_path, go_vog_search_dir, f"{temp_dir}/tmp"))
    go_vog_df=go_vog_df[go_vog_df["bitscore"]>50]
    to_remove_vogdb = go_vog_df["target"].unique()


    #search cyanobacteria for metabolic genes
    go_cyanobacteria_search_dir = f"{args.output}/go_cyanobacteria_genome_search"
    os.makedirs(go_cyanobacteria_search_dir, exist_ok=True)
    go_cyano_genome_df = search_to_df(mmseqs_search(go_cyano_db_path, ref_cyanobacteria_genome_db, go_cyanobacteria_search_dir, f"{temp_dir}/tmp"))

    #search cyanophage for metabolic genes
    go_cyanophage_search_dir = f"{args.output}/go_cyanophage_genome_search"
    os.makedirs(go_cyanophage_search_dir, exist_ok=True)
    go_cyanophage_df = search_to_df(mmseqs_search(go_cyano_db_path, ref_cyanophage_genome_db, go_cyanophage_search_dir, f"{temp_dir}/tmp"))
    go_cyanophage_df=go_cyanophage_df[go_cyanophage_df["bitscore"]>50]
    potential_amg=go_cyanophage_df["query"].unique()

    vogdb_cyano_search_dir = f"{args.output}/vogdb_cyano_search"
    os.makedirs(vogdb_cyano_search_dir, exist_ok=True)
    vogdb_cyano_df = search_to_df(mmseqs_search(vogdb_db_path, ref_cyanobacteria_genome_db, vogdb_cyano_search_dir, f"{temp_dir}/tmp", extra_params="--search-type 3 "))

    vogdb_cyano_df=vogdb_cyano_df[vogdb_cyano_df["bitscore"]>50]
    vogdb_cyano_df=vogdb_cyano_df[~vogdb_cyano_df["target"].isin(to_remove_vogdb)]
    
    go_vog_distances = find_closest_distances(go_cyano_genome_df, vogdb_cyano_df)
    vogdb_distances = find_vogdb_distances(vogdb_cyano_df)

    # Plot and save both histograms
    plot_histogram(go_vog_distances, 'go_vog_distances_histogram.png', 'Closest Distances Between GO Cyano Hits and VOGdb Hits (0-10,000)')
    plot_histogram(vogdb_distances, 'vogdb_distances_histogram.png', 'Distances Between Consecutive VOGdb Hits (0-10,000)')





def main():
    temp_dir = "temp"
    args = arguments("data", "output", "rerunTruth", "keeptempTruth")

    if not args.rerun:
        logging.info("Doing from beginning, if Re-run is needed, please set the flag, but also use the --keeptemp flag")
        i=0
        if os.path.exists(temp_dir):
            while os.path.exists(temp_dir):
                i += 1
            temp_dir = f"{temp_dir}_{i}"

    try:
        #make directories for temporary files and output
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(args.output, exist_ok=True)

        # Initialize logging
        log_filename = f"{args.output}/amg_analysis.log"
        init_logging(log_filename)
        run(temp_dir=temp_dir, args = args)
    
    except Exception as e:
        logging.error(f"An error occurred: {e}", exc_info=True)
    
    finally:
        if not args.keeptemp:
            import shutil
            logging.info("Cleaning up temporary files")
            shutil.rmtree(temp_dir)
            logging.info("Done, cleaning up")

if __name__ == "__main__":
    main()