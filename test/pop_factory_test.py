import unittest
import pop_factory

def test_snp_tuple(snp_id):
    snp_tuple = pop_factory.SNPTuples(snp_id)
    snp_tuple.add_tuple("G", 0.20)
    snp_tuple.add_tuple("A", 0.70)
    snp_tuple.add_tuple("T", 1.0)
    return snp_tuple

class PopulationFactoryTest(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)


class SNPTupleTest(unittest.TestCase):

    def setUp(self):
        self.snp_tuple = test_snp_tuple(100)

    def test_pick_nt(self):
        picked = self.snp_tuple.pick_snp_value(0.95)
        self.assertEqual("T", picked, "Should pick T")
        picked = self.snp_tuple.pick_snp_value(0.4)
        self.assertEqual("A", picked, "Should pick A")

    def test_pick_pathogen(self):
        picked = self.snp_tuple.pick_pathogen_value()
        self.assertEqual("T", picked, "Should pick T as pathogen (2nd most frequent)")

class PathogenGroupTest(unittest.TestCase):

    def setUp(self):
        self.snp_data = [(1, [test_snp_tuple(1), test_snp_tuple(2), test_snp_tuple(3), test_snp_tuple(4)])]

    def test_from_yml(self):
        pg = pop_factory.PathogenGroup.from_yml({"mutation_weights": [0.5, 0.5, 0.5],
                                                 "num_instances": 2,
                                                 "population_weight": 5,
                                                 "min_minor_allele_freq": 0.01}, self.snp_data)
        self.assertEqual(5, pg.population_weight, "Population Weight should be 5")
        self.assertEqual(3, len(pg.pathogens), "Should have picked 3 pathogens")
        self.assertEqual(2, len(pg.select_mutations()))

if __name__ == '__main__':
    unittest.main()
