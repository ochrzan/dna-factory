import concurrent.futures
import gzip
import glob
import re
from common.snp import RefSNP, Allele, chromosome_from_filename
from common import Session, db_engine
from sqlite3 import OperationalError


def convert_json_to_db():
    snp_file_list = glob.glob("output/*chr*.json*")
    skip_filter = "chrAAA\\."
    # Clean out tables to be loaded
    try:
        Allele.__table__.drop(db_engine)
        RefSNP.__table__.drop(db_engine)
    except OperationalError as e:
        print('Warning - Error droping tables before loading. Possible tables do not exist. Continuing')
    RefSNP.__base__.metadata.create_all(db_engine)
    filtered_list = list(filter(lambda x: not re.search(skip_filter, x), snp_file_list))
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        # Start the load operations and mark each future with its file
        future_to_file = {executor.submit(load_file_into_db, json_file): json_file for json_file in filtered_list}
        for future in concurrent.futures.as_completed(future_to_file):
            filename = future_to_file[future]
            try:
                data = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (filename, exc))
            else:
                print('Successfully loaded %s into DB' % filename)


    # perform queries to set total_count and MAF
    update_total_sql = """
    update ref_snps set total_count = 
    (select total_count from 
        (select ref_snp_id, sum(allele_count) as total_count from alleles group by ref_snp_id)
     where
     id = ref_snp_id);
    """
    update_maf_sql = """
    --- MAF query
    update ref_snps as r set maf = 
    ( select sec_high * 1.0 / r.total_count  from 
        ( select a.ref_snp_id, max(a.allele_count) as sec_high from 
            (select ref_snp_id, max(allele_count) as highest from alleles group by ref_snp_id) as x
             join alleles a on x.ref_snp_id = a.ref_snp_id and a.allele_count != highest 
          group by a.ref_snp_id) as s
      where r.id = s.ref_snp_id
    ) 
    """
    session = Session()
    print('Updating total counts in ref_snps table')
    session.execute(update_total_sql)
    session.commit()
    print('Updating maf in ref_snps table')
    session.execute(update_maf_sql)
    session.commit()
    session.close()

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
    session.close()


if __name__ == '__main__':
    convert_json_to_db()
