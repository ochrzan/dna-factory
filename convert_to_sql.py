import gzip
import glob
import re
import common.snp as common
from common.db import db
import getopt
import sys


def convert_json_to_db(append_mode):
    snp_file_list = glob.glob("output/*chr*.json*")
    skip_filter = "chr(1|10|11|12|Y|MT)\\."
    # Clean out tables to be loaded
    if not append_mode:
        try:
            db.ref_snps.drop(db.engine)
            db.alleles.drop(db.engine)
        except:
            print('Warning - Error droping tables before loading. Possible tables do not exist. Continuing')

    db.metadata.create_all(db.engine)
    filtered_list = list(filter(lambda x: not re.search(skip_filter, x), snp_file_list))
    for f in filtered_list:
        # Start the load operations and mark each future with its file
        print('Start loading %s into DB' % f)
        load_file_into_db(f)
        print('Successfully loaded %s into DB' % f)
    print('DONE')


def load_file_into_db(snp_file):
    chromosome = common.chromosome_from_filename(snp_file)
    with gzip.open(snp_file, 'rt') as f:
        line_num = 0
        ins_ref_snps = []
        ins_alleles = []
        for json_line in f:
            ref_snp = common.RefSNP.from_json(json_line, chromosome)
            ins_ref_snps.append(ref_snp)
            ins_alleles.extend(ref_snp.alleles)
            line_num += 1
            if line_num % 10000 == 0:
                db.bulk_insert(ins_ref_snps)
                db.bulk_insert(ins_alleles)
                ins_ref_snps = []
                ins_alleles = []
        db.bulk_insert(ins_ref_snps)
        db.bulk_insert(ins_alleles)

        # perform queries to set total_count and MAF
    print('Updating total counts in ref_snps table')
    common.RefSNP.update_total_counts(chromosome, db.connection)

    print('Updating maf in ref_snps table')
    common.RefSNP.update_maf(chromosome, db.connection)


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
    db.default_init()
    convert_json_to_db(append_data)
    db.connection.close()
