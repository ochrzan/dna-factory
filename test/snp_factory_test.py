import unittest
from pop_factory import SnpFactory
from common.timer import Timer
from common.snp import CHROMOSOME_PROB
import numpy


class SnpFactoryTest(unittest.TestCase):

    def setUp(self) -> None:
        self.snp_factory = SnpFactory.init_from_cdf_file()

    @Timer(logger=print, text="Elapsed time test_pick_mafs: {:0.4f} sec")
    def test_gen_snps(self):
        snps_size = 100000
        min_maf = 0.16
        snps = self.snp_factory.random_snp_tuples(snps_size, min_maf=min_maf)
        self.assertEqual(len(snps), snps_size)
        self.assertEqual(len(snps[0].tuples), 2)
        count_chromosome_one = 0
        count_largest_maf = 0
        largest_maf = self.snp_factory.sorted_maf[len(self.snp_factory.sorted_maf) - 1]
        min_maf_works = True
        for s in snps:
            if s.chromosome == "1":
                count_chromosome_one += 1
            if s.tuples[0][1] == 1 - largest_maf:
                count_largest_maf += 1
            if (1 - s.tuples[0][1]) < min_maf:
                min_maf_works = False
        self.assertTrue(min_maf_works, "Found tuple with MAF under the min_maf filter")
        self.assertAlmostEqual(count_largest_maf * 1.0 / snps_size,
                               self.snp_factory.pdf[len(self.snp_factory.sorted_maf) - 1], delta=0.01)
        self.assertAlmostEqual(count_chromosome_one * 1.0 / snps_size, CHROMOSOME_PROB[0], delta=0.01)



if __name__ == '__main__':
    unittest.main()
