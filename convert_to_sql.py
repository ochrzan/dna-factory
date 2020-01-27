import concurrent.futures
import gzip
import glob
import re
from common.snp import RefSNP, Allele, chromosome_from_filename
from common import Session, db_engine
from sqlite3 import OperationalError
import getopt
import sys


def convert_json_to_db(append_mode):
    snp_file_list = glob.glob("output/*chr*.json*")
    skip_filter = "chr(1|10|11|12|15|Y|MT)\\."
    # Clean out tables to be loaded
    if not append_mode:
        try:
            Allele.__table__.drop(db_engine)
            RefSNP.__table__.drop(db_engine)
        except:
            print('Warning - Error droping tables before loading. Possible tables do not exist. Continuing')

    RefSNP.__base__.metadata.create_all(db_engine)
    filtered_list = list(filter(lambda x: not re.search(skip_filter, x), snp_file_list))
    for f in filtered_list:
        # Start the load operations and mark each future with its file
        load_file_into_db(f)
        print('Successfully loaded %s into DB' % f)
    print('DONE')


def load_file_into_db(snp_file):
    session = Session()
    chromosome = chromosome_from_filename(snp_file)
    with gzip.open(snp_file, 'rt') as f:
        line_num = 0
        for json_line in f:
            ref_snp = RefSNP.from_json(json_line, chromosome)
            session.add(ref_snp)
            line_num += 1
            if line_num % 10000 == 0:
                session.commit()
        session.commit()

        # perform queries to set total_count and MAF
    print('Updating total counts in ref_snps table')
    RefSNP.update_total_counts(chromosome, session)

    print('Updating maf in ref_snps table')
    RefSNP.update_maf(chromosome, session)
    session.close()


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a")
    except getopt.GetoptError as err:
        print(err.msg)
        sys.exit(2)
    append_data = False
    for opt, arg in opts:
        if opt == '-a':
            append_data = True
    convert_json_to_db(append_data)
