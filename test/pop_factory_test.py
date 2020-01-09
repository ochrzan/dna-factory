import unittest
import pop_factory

class PopulationFactoryTest(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)


class SNPTupleTest(unittest.TestCase):

    def setUp(self):
        self.snp_tuple = pop_factory.SNPTuples(100)
        self.snp_tuple.add_tuple("G", 0.20)
        self.snp_tuple.add_tuple("A", 0.70)
        self.snp_tuple.add_tuple("T", 1.0)

    def test_pick_nt(self):
        picked = self.snp_tuple.pick_snp_value(0.95)
        self.assertEqual("T", picked, "Should pick T")
        picked = self.snp_tuple.pick_snp_value(0.4)
        self.assertEqual("A", picked, "Should pick A")

    def test_pick_pathogen(self):
        picked = self.snp_tuple.pick_pathogen_value()
        self.assertEqual("T", picked, "Should pick T as pathogen (2nd most frequent)")


if __name__ == '__main__':
    unittest.main()
