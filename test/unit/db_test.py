import unittest
from common.db import db
from common.snp import RefSNP, Allele


class DbTest(unittest.TestCase):
    def setUp(self):
        db.db_init('sqlite:///unit_test.db')
        db.metadata.drop_all(bind=db.connection)
        db.metadata.create_all(db.engine)
        self.r1 = RefSNP(1000, '1')
        self.r2 = RefSNP(1001, '2')
        a1 = Allele('A', 'A', 100190109)
        a1.allele_count = 5000
        a1.total_count = 11000
        self.r1.put_allele(a1)
        a2 = Allele('A', 'G', 100190109)
        a2.allele_count = 6000
        a2.total_count = 11000
        self.r1.put_allele(a2)

    def test_bulk_insert(self):
        result = db.bulk_insert([self.r1, self.r2], db.ref_snps)
        self.assertEqual(2, result.rowcount)
        result = db.bulk_insert(self.r1.alleles, db.alleles)
        self.assertEqual(2, result.rowcount)

    def test_delete_chromosome(self):
        result = db.bulk_insert([self.r1, self.r2], db.ref_snps)
        result = db.bulk_insert(self.r1.alleles, db.alleles)
        RefSNP.delete_chromosomes(["1", "2"], db.connection)
        select_query = db.ref_snps.select().where(db.ref_snps.c.chromosome == '1')
        one_row = db.connection.execute(select_query).fetchone()
        self.assertIsNone(one_row)


    def test_default_init(self):
        db.default_init()
        self.assertTrue(True)



if __name__ == '__main__':
    unittest.main()
