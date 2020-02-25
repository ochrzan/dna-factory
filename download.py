
import bz2
import concurrent.futures
import getopt
import multiprocessing
import sys
from ftplib import FTP
import time
import re
import os
import hashlib

from common.snp import RefSNP, chromosome_from_filename
from common.db import db
from definitions import ROOT_DIR

DBSNP_FTP_DOMAIN = 'ftp.ncbi.nih.gov'
DBSNP_LATEST_DIR = '/snp/latest_release/JSON'
DBSNP_JSON_PATTERN = 'chr%s\\.json.bz2$'
DOWNLOAD_DIR = os.path.join(ROOT_DIR, 'tmp_download')


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
    download_path = DOWNLOAD_DIR + '/' + json_file
    download_needed = True
    if os.path.exists(download_path):
        md5 = []
        ftp.retrlines('RETR ' + json_file + ".md5", md5.append)
        if md5:
            md5 = md5[0].split(" ")[0]
        print("FTP MD5 of %s: %s" % (json_file, md5))
        block_size = 65536
        hasher = hashlib.md5()
        with open(DOWNLOAD_DIR + '/' + json_file, 'rb') as afile:
            buf = afile.read(block_size)
            while len(buf) > 0:
                hasher.update(buf)
                buf = afile.read(block_size)
        local_md5 = hasher.hexdigest()
        print("Local File MD5: %s" % local_md5)
        if local_md5 == md5:
            print("MD5 matches local copy. Skipping download.")
            download_needed = False
    if download_needed:
        with open(DOWNLOAD_DIR + '/' + json_file, 'wb') as f:
            ftp.retrbinary('RETR ' + json_file, f.write)
    with bz2.BZ2File(DOWNLOAD_DIR + '/' + json_file, 'rb') as f_in:
        chromosome = chromosome_from_filename(json_file)
        for line in f_in:
            snp = RefSNP.from_nih_json(line, chromosome)
            if snp.total_count > 0 and snp.alleles:
                if 0 <= min_maf <= snp.maf:
                    queue.put(snp)
    return True


def write_snps_to_db(q):
    line_num = 0
    ins_ref_snps = []
    ins_alleles = []
    snps_inserted = 0
    while not q.empty():
        snp = q.get()
        ins_ref_snps.append(snp)
        ins_alleles.extend(snp.alleles)
        line_num += 1
        if line_num % 1000 == 0:
            db.bulk_insert(ins_ref_snps, db.ref_snps)
            db.bulk_insert(ins_alleles, db.alleles)
            ins_ref_snps = []
            ins_alleles = []
            snps_inserted += 1000
    if ins_ref_snps and ins_alleles:
        snps_inserted += len(ins_ref_snps)
        db.bulk_insert(ins_ref_snps, db.ref_snps)
        db.bulk_insert(ins_alleles, db.alleles)
    return snps_inserted


def ftp_login():
    ftp = FTP(DBSNP_FTP_DOMAIN)  # connect to host, default port
    ftp.login()  # user anonymous, passwd anonymous@
    ftp.cwd(DBSNP_LATEST_DIR)
    return ftp


def download_ref_snps(chromosome_list, num_workers=2, append=False, min_maf=0):
    """ Downloads all RefSNP data from NIH's FTP site. Requires ~250 GB of disk space
    """
    ftp = ftp_login()
    file_list = []
    ftp.retrlines('NLST', file_list.append)
    search_pattern = DBSNP_JSON_PATTERN % ".*"
    if chromosome_list:
        chromosome_match = "(" + "|".join(chromosome_list).replace(" ", "") + ")"
        search_pattern = DBSNP_JSON_PATTERN % chromosome_match
    json_for_dl = [f for f in file_list if re.search(search_pattern, f)]
    if not append:
        print("Removing old data from DB.")
        try:
            if not chromosome_list:
                print("No chromosome list specified. Clearing entire DB.")
                db.ref_snps.drop(db.engine)
                db.alleles.drop(db.engine)
            else:
                RefSNP.delete_chromosomes(chromosome_list, db.connection)
        except:
            print('Warning - Error droping tables before loading. Possible tables do not exist. Continuing')
    # Create schema if missing
    db.metadata.create_all(db.engine)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor, multiprocessing.Manager() as m:
        q = m.Queue(10000)
        # Start the load operations and mark each future with its URL
        future_to_file = {executor.submit(fetch_snp_file, json_file, q, min_maf): json_file for json_file in json_for_dl}
        time.sleep(10)
        count_inserted = 0
        while any(not f.done() for f in future_to_file.keys()):
            try:
                count_inserted += write_snps_to_db(q)
                print("Inserted %i refSNPs." % count_inserted)
                time.sleep(2)  # Wait for more items in the queue
            except Exception as e:
                print("Exception writing snps to DB.")
                print(e)
                for f in future_to_file.keys():
                    f.cancel()
                raise e
        for future in future_to_file.keys():
            filename = future_to_file[future]
            try:
                data = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (filename, exc))
            else:
                print('Successfully downloaded %s and loaded into db.' % filename)
    m.shutdown()


def print_help():
    print("""
    download is used to download data from the NIH refSNP ftp site and load it into a database
    Accepted Inputs are:
    -c             allows providing a comma delimited list of chromosomes to load (ie 3,5,6,X)
    --chromosomes  same as -c
    -f 0.n         min minor allele frequency for a SNP to be loaded, default is 0 (no filtering)
    -n n           number of worker processes to use
    -a             run in append mode. Will not delete any existing data. This can cause primary key insert errors if
                   a refSNP with the same ID is already in the database.
    This app uses a single writer process and multiple worker download processes that generate records for the writer. 
    If disk is slow the writer can bottleneck with a high worker process count (-n option).
    """)


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
        elif opt in ("-c", "--chromosomes"):
            chromosomes = arg.split(",")
            print("Chromosome list provided. Will only load/reload refSNPs on chromosomes %s" % arg)
        elif opt in "-f":
            min_freq = float(arg)
            print("Min minor allele frequency for loading set to %f" % min_freq)
        elif opt in "-n":
            num_processes = int(arg)
            print("Will use %i worker processes for file downloads." % num_processes)
        elif opt in "-a":
            append = True
            print("Running in append mode. Will not delete any existing data in database.")

    db.default_init()
    download_ref_snps(chromosomes, num_processes, append, min_freq)
