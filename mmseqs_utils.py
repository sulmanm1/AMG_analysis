import os, logging
from utils import running_message, run_command

@running_message
def mmseqs_makedb(input, outpath):
    if not os.path.exists(outpath):
        cmd = f"mmseqs createdb {input} {outpath}"
        run_command(cmd, shell=True)
    else:
        logging.info("Database already exists, use existing database")

@running_message
def mmseqs_cluster_cmd(db, outdir, threads, sensitivity, coverage=0.5, min_id = 0.3 ,e_value=1e5):
    tmp = f"{outdir}/tmp"
    os.makedirs(tmp, exist_ok=True)
    
    cluster_output = f"{outdir}/cluster_output"
    if not os.path.exists(f"{cluster_output}.index"):
        cmd = f"mmseqs cluster -s {sensitivity} -c {coverage} --min-seq-id {min_id} -e {e_value} --threads {threads} {db} {cluster_output} {tmp}"
        run_command(cmd, shell=True)
    else:
        print("Output file already exists, using existing file")
    return cluster_output

@running_message
def mmseqs_createtsv(db, cluster_output, outdir):
    tsv_path = f"{outdir}/cluster_output.tsv"
    if not os.path.exists(tsv_path):
        cmd = f"mmseqs createtsv {db} {db} {cluster_output} {tsv_path}"
        run_command(cmd, shell=True)
    else:
        print("Output file already exists, using existing file")
    return tsv_path

@running_message
def mmseqs_search(query_db, ref_db, outdir, tmp_path, extra_params=""):
    if not os.path.exists(f"{outdir}/network.m8"):
        cmd = f"mmseqs search {query_db} {ref_db} {outdir}/network_int {tmp_path} {extra_params}"
        run_command(cmd, shell=True, check=True)
        cmd2 = f"mmseqs convertalis {query_db} {ref_db} {outdir}/network_int {outdir}/network.m8"
        run_command(cmd2, shell=True, check=True)
    return f"{outdir}/network.m8"
