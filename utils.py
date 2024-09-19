from Bio import SeqIO
import os
from datetime import datetime, timedelta
import pandas as pd
from tqdm import tqdm
import io
import logging

def arguments(*args):
    import argparse
    parser = argparse.ArgumentParser()
    for arg in args:
        if not "Truth" in arg:
            parser.add_argument(f"--{arg}", required=True)
        else:
            print(arg.split("Truth")[0])
            parser.add_argument(f"--{arg.split("Truth")[0]}", action="store_true")
    args = parser.parse_args()
    return args

def get_unique_log_filename(base_log_filename):
    i = 1
    base_name, extension = os.path.splitext(base_log_filename)
    while os.path.exists(base_log_filename):
        base_log_filename = f"{base_name}_{i}{extension}"
        i += 1
    return base_log_filename

def init_logging(log_filename):
    from sys import stdout as sys_stdout
    class TeeHandler(logging.Handler):
        def __init__(self, filename, mode='a'):
            super().__init__()
            self.file = open(filename, mode)
            self.stream_handler = logging.StreamHandler(sys_stdout)

        def emit(self, record):
            log_entry = self.format(record)
            self.file.write(log_entry + '\n')
            self.file.flush()  # Ensure log entry is written to file immediately
            self.stream_handler.emit(record)

        def close(self):
            self.file.close()
            super().close()
    i=1
    while os.path.exists(log_filename):
        log_filename = get_unique_log_filename(log_filename)
    print(f"Logging to {log_filename}")
    # Configure logging to use the TeeHandler
    tee_handler = TeeHandler(log_filename)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    tee_handler.setFormatter(formatter)

    logging.basicConfig(level=logging.INFO, handlers=[tee_handler])

def main_logging(function):
    def wrapper(*args, **kwargs):
        try:
            result = function(*args, **kwargs)
            return result
        except Exception as e:
            logging.error(f"An error occurred: {e}", exc_info=True)
            print(f"An error occurred: {e}")
    return wrapper

def run_command(cmd, **kwargs):
    import shlex
    import subprocess
    if not 'stderr' in kwargs:
        kwargs['stderr'] = subprocess.PIPE
    if not 'stdout' in kwargs:
        kwargs['stdout'] = subprocess.PIPE
    if not 'text' in kwargs:
        kwargs['text'] = True
    if not 'bufsize' in kwargs:
        kwargs['bufsize'] = 1


    logging.info(f"Running command: {cmd}")
    
    cmd_list = shlex.split(cmd)

    with subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1) as process:
        for stdout_line in iter(process.stdout.readline, ""):
            logging.info(stdout_line.strip())

        for stderr_line in iter(process.stderr.readline, ""):
            logging.error(stderr_line.strip())

        process.stdout.close()
        process.stderr.close()
        return_code = process.wait()

        if return_code:
            logging.error(f"Command '{cmd}' failed with return code {return_code}")
            raise subprocess.CalledProcessError(return_code, cmd)

    logging.info(f"Command '{cmd}' completed successfully")
    return subprocess.CompletedProcess(cmd, return_code)

def running_message(function):
    def wrapper(*args, **kwargs):
        def format_argument(arg):
            if isinstance(arg, pd.DataFrame):
                return f"DataFrame({len(arg)} rows x {len(arg.columns)} columns)"
            elif isinstance(arg, (list, dict)) and len(arg) > 10:
                return f"{type(arg).__name__}({len(arg)} items)"
            return repr(arg)

        def format_timedelta(delta):
            seconds = delta.total_seconds()
            if seconds < 60:
                return f"{seconds:.2f} seconds"
            elif seconds < 3600:
                minutes = seconds / 60
                return f"{minutes:.2f} minutes"
            elif seconds < 86400:
                hours = seconds / 3600
                return f"{hours:.2f} hours"
            else:
                days = seconds / 86400
                return f"{days:.2f} days"

        T1 = datetime.now()
        arg_names = function.__code__.co_varnames[:function.__code__.co_argcount]
        args_repr = [f"{arg}={format_argument(a)}" for arg, a in zip(arg_names, args)]
        kwargs_repr = [f"{k}={format_argument(v)}" for k, v in kwargs.items()]
        signature = ", ".join(args_repr + kwargs_repr)
        logging.info('')
        logging.info("=======================================")
        logging.info(f"Running {function.__name__}")
        logging.info(f"Using inputs:{function.__name__}({signature})")

        try:
            result = function(*args, **kwargs)
        except Exception as e:
            logging.error(f"An error occurred in {function.__name__}: {e}", exc_info=True)
        else:
            T2 = datetime.now()
            total_time = format_timedelta(T2 - T1)
            logging.info(f"{function.__name__} Completed")
            logging.info(f"Total time taken: {total_time}")
            return result

    return wrapper

def read_fasta(fastafile):
    # Get the total size of the file
    total_size = os.path.getsize(fastafile)
    
    # Get the total number of records to parse with a progress bar for reading lines
    with open(fastafile) as f, tqdm(total=total_size, desc="Reading FASTA file", unit="B", unit_scale=True, unit_divisor=1024) as pbar:
        total_records = 0
        for line in f:
            pbar.update(len(line))
            if line.startswith(">"):
                total_records += 1
    
    # Parse the FASTA file with a progress bar
    with tqdm(total=total_records, desc="Parsing FASTA file", unit=" Records") as pbar:
        records = []
        for record in SeqIO.parse(fastafile, "fasta"):
            records.append(record)
            pbar.update(1)
    
    return records


def read_fasta_no_progress(fastafile):
    # Get the total number of records
    with open(fastafile) as f:
        total_records = 0
        for line in f:
            if line.startswith(">"):
                total_records += 1
    
    # Parse the FASTA file
    records = []
    for record in SeqIO.parse(fastafile, "fasta"):
        records.append(record)
    
    return records


def write_fasta(outpath: str, recordlist: list)->None:
    '''
    Writes a fasta file to a given location
    :param path: location to write fasta
    :param recordlist: list of fasta records
    :return: None
    '''
    with open(outpath, "w") as output_handle:
        SeqIO.write(recordlist, output_handle, "fasta")


def makedir(path: str, force_make: bool=False)->str:
    '''
    Makes a directory if the given direcotry doesn't exist. If force_make true, makes a directory with a number
    :param path: location of new directory
    :param force_make: To make a new directory with a number if a directory already exists at given path
    :return: path to new dir
    '''
    try:
        os.mkdir(path)
        return path
    except:
        if not force_make:
            return path
    i = 1
    if force_make:
        while True:
            new_name = f"{path}_{i}"
            try:
                os.mkdir(new_name)
                break
            except:
                i += 1
    
    
def read_lines(file_path):
    file_size = os.path.getsize(file_path)
    lines = []
    
    with open(file_path, 'r') as file:
        with tqdm(total=file_size, unit='B', unit_scale=True, desc=f"Reading {file_path}") as pbar:
            for line in file:
                lines.append(line)
                pbar.update(len(line))
    
    return lines

def pd_read_csv(file_path, chunksize=10000, **kwargs):
    file_size = os.path.getsize(file_path)
    
    with tqdm(total=file_size, unit='B', unit_scale=True, desc=f"Reading {file_path}") as pbar:
        chunks = []
        with open(file_path, 'r') as file:
            chunk_data = ""
            while True:
                chunk = file.read(chunksize)
                if not chunk:
                    break
                chunk_data += chunk
                pbar.update(len(chunk))
                
                if '\n' in chunk:
                    last_newline = chunk_data.rfind('\n')
                    to_process = chunk_data[:last_newline]
                    chunk_data = chunk_data[last_newline + 1:]
                    
                    chunks.append(pd.read_csv(io.StringIO(to_process), **kwargs))
        
        if chunk_data:
            chunks.append(pd.read_csv(io.StringIO(chunk_data), **kwargs))
    
    df = pd.concat(chunks, ignore_index=True)
    return df