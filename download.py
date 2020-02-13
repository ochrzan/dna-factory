
import bz2
import concurrent.futures
import getopt
import sys
from ftplib import FTP
import time
import re
import os
from multiprocessing import Queue

from common.snp import RefSNP, chromosome_from_filename
from common.db import db

DBSNP_FTP_DOMAIN = 'ftp.ncbi.nih.gov'
DBSNP_LATEST_DIR = '/snp/latest_release/JSON'
DBSNP_JSON_PATTERN = 'chr.*.json.bz2$'
#DBSNP_JSON_PATTERN = 'sample.*.json.bz2$'
DOWNLOAD_DIR = 'tmp_download'
OUTPUT_DIR = 'output'

# TODO support append mode and selecting the chromosomes to load (with cleanup of old data), min MAF filter

def fetch_snp_file(json_file, queue, min_maf=0):
    """
    Fetch a NIH refSNP file then open it and add RefSNP objects to the work queue.
    :param json_file: NIH file to download via FTP
    :param queue: work queue for RefSNP objects
    :return:
    """
    # Not sure if ftplib is threadsafe so use a ftp login per call
    ftp = ftp_login()
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(DOWNLOAD_DIR + '/' + json_file, 'wb') as f:
        ftp.retrbinary('RETR ' + json_file, f.write)
    with bz2.BZ2File(DOWNLOAD_DIR + '/' + json_file, 'rb') as f_in:
        chromosome = chromosome_from_filename(json_file)
        for line in f_in:
            snp = RefSNP.from_nih_json(line, chromosome)
            if snp.total_count > 0:
                if 0 < min_maf <= snp.get_maf():
                    queue.put(snp)


def write_snps_to_db(q):
    while not q.empty():
        line_num = 0
        ins_ref_snps = []
        ins_alleles = []
        for snp in q.get():
            ins_ref_snps.append(snp)
            ins_alleles.extend(snp.alleles)
            line_num += 1
            if line_num % 1000 == 0:
                db.bulk_insert(ins_ref_snps)
                db.bulk_insert(ins_alleles)
                ins_ref_snps = []
                ins_alleles = []
        db.bulk_insert(ins_ref_snps)
        db.bulk_insert(ins_alleles)


def ftp_login():
    ftp = FTP(DBSNP_FTP_DOMAIN)  # connect to host, default port
    ftp.login()  # user anonymous, passwd anonymous@
    ftp.cwd(DBSNP_LATEST_DIR)
    return ftp


def download_ref_snps(chromosomes, num_workers=2, append=False, min_maf=0):
    """ Downloads all RefSNP data from NIH's FTP site. Requires ~250 GB of disk space
    """
    ftp = ftp_login()
    file_list = []
    ftp.retrlines('NLST', file_list.append)
    json_for_dl = [f for f in file_list if re.search(DBSNP_JSON_PATTERN, f)]
    q = Queue(1000)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Start the load operations and mark each future with its URL
        future_to_file = {executor.submit(fetch_snp_file, json_file, q, min_maf): json_file for json_file in json_for_dl}
        time.sleep(10)
        while any(not f.done() for f in future_to_file.keys()):
            write_snps_to_db(q)
            time.sleep(1)  # Sleep since queue is empty but downloads may be happening
        for future in concurrent.futures.as_completed(future_to_file):
            filename = future_to_file[future]
            try:
                data = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (filename, exc))
            else:
                print('Successfully downloaded %s and loaded into db.' % filename)

def print_help():
    pass
    # TODO write help string

if __name__ == '__main__':
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "n:c:f:a", ["help"])
    except getopt.GetoptError as err:
        print(err.msg)
        print_help()
        sys.exit(2)
    min_freq = None
    num_processes = 2
    chromosomes = None
    append = False
    for opt, arg in opts:
        if opt in ('-h', "-?", "--help"):
            print_help()
            sys.exit()
        elif opt in "-c":
            chromosomes = arg.split(",")
        elif opt in "-f":
            min_freq = float(arg)
        elif opt in "-n":
            num_processes = int(arg)
        elif opt in "-a":
            append = True

    download_ref_snps(chromosomes, num_processes, append, min_freq)
