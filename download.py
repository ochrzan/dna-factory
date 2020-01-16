
import bz2
import concurrent.futures
from ftplib import FTP
import gzip
import re
import os
from common.snp import RefSNP

DBSNP_FTP_DOMAIN = 'ftp.ncbi.nih.gov'
DBSNP_LATEST_DIR = '/snp/latest_release/JSON'
DBSNP_JSON_PATTERN = 'chr2.*.json.bz2$'
#DBSNP_JSON_PATTERN = 'sample.*.json.bz2$'
DOWNLOAD_DIR = 'tmp_download'
OUTPUT_DIR = 'output'

def fetch_snp_file(json_file):
    # Not sure if ftplib is threadsafe so use a ftp login per call
    ftp = ftp_login()
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(DOWNLOAD_DIR + '/' + json_file, 'wb') as f:
        ftp.retrbinary('RETR ' + json_file, f.write)
    with bz2.BZ2File(DOWNLOAD_DIR + '/' + json_file, 'rb') as f_in, \
            gzip.open(OUTPUT_DIR + '/snp_' + json_file.replace(".json.bz2", ".yml.gz"), 'w') as snp_out:
        for line in f_in:
            snp = RefSNP.from_json(line)
            # Only write snps with frequency data
            if snp.total_allele_count() > 0:
                snp_out.write(str(snp) + "\n")


def ftp_login():
    ftp = FTP(DBSNP_FTP_DOMAIN)  # connect to host, default port
    ftp.login()  # user anonymous, passwd anonymous@
    ftp.cwd(DBSNP_LATEST_DIR)
    return ftp


def download_ref_snps():
    """ Downloads all RefSNP data from NIH's FTP site. Requires ~250 GB of disk space
    """
    ftp = ftp_login()
    file_list = []
    ftp.retrlines('NLST', file_list.append)
    json_for_dl = [f for f in file_list if re.search(DBSNP_JSON_PATTERN, f)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        # Start the load operations and mark each future with its URL
        future_to_file = {executor.submit(fetch_snp_file, json_file): json_file for json_file in json_for_dl}
        for future in concurrent.futures.as_completed(future_to_file):
            filename = future_to_file[future]
            try:
                data = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (filename, exc))
            else:
                print('Successfully downloaded %s' % filename)


if __name__ == '__main__':
    download_ref_snps()
