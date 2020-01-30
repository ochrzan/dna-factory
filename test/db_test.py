import unittest
from common.db import db
from common.snp import RefSNP, Allele


class DbTest(unittest.TestCase):
    def setUp(self):
        db.db_init('sqlite:///unit_test.db')
        db.metadata.drop_all(bind=db.connection)
        db.metadata.create_all(db.engine)

    def test_bulk_insert(self):
        r1 = RefSNP(1000, '1')
        r2 = RefSNP(1001, '2')
        result = db.bulk_insert([r1, r2], db.ref_snps)
        self.assertEqual(2, result.rowcount)
        a1 = Allele('A', 'A', 100190109)
        a1.allele_count = 5000
        a1.total_count = 11000
        r1.put_allele(a1)
        a2 = Allele('A', 'G', 100190109)
        a2.allele_count = 6000
        a2.total_count = 11000
        r1.put_allele(a2)
        result = db.bulk_insert([a1, a2], db.alleles)
        self.assertEqual(2, result.rowcount)


if __name__ == '__main__':
    unittest.main()
