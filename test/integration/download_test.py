import unittest
from multiprocessing import Queue

from common.db import db
import download

class DownloadTestCase(unittest.TestCase):

    def setUp(self):
        db.db_init('sqlite:///integration_test.db')
        db.metadata.drop_all(bind=db.connection)
        db.metadata.create_all(db.engine)

    def test_download_chromosome(self):
        download.download_ref_snps(chromosome_list=["Y"])
        select_query = db.ref_snps.select().where(db.ref_snps.c.chromosome == 'Y')
        one_row = db.connection.execute(select_query).fetchone()
        self.assertTrue(one_row is not None)
        self.assertEqual(one_row[db.ref_snps.c.chromosome], 'Y')

    def test_fetch_sample_file(self):
        q = Queue(1000)
        download.fetch_snp_file('refsnp-sample.json.bz2', q)

        snp = q.get(timeout=10)
        db.bulk_insert([snp], db.ref_snps)
        db.bulk_insert(snp.alleles, db.alleles)

        q.close()
        self.assertTrue(snp.id is not None)

    def test_append_and_maf(self):
        herman_cain = 999999999999
        insert_vals = {"id": herman_cain, "chromosome": "Y"}
        dummy_val = db.connection.execute(db.ref_snps.insert(), insert_vals)

        download.download_ref_snps(chromosome_list=["Y"], append=True, min_maf=0.2)
        select_query = db.ref_snps.select().where(db.ref_snps.c.id == herman_cain)
        one_row = db.connection.execute(select_query).fetchone()
        self.assertTrue(one_row is not None)
        self.assertEqual(one_row[db.ref_snps.c.chromosome], 'Y')
        select_query = db.ref_snps.select().where(db.ref_snps.c.maf < 0.2)
        one_row = db.connection.execute(select_query).fetchone()
        self.assertTrue(one_row is None)

if __name__ == '__main__':
    unittest.main()
