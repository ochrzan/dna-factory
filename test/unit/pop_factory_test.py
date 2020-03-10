import unittest
import pop_factory


def test_snp_tuple(snp_id):
    snp_tuple = pop_factory.SNPTuples(snp_id, "1", 50000)
    snp_tuple.add_tuple("G", 0.70)
    snp_tuple.add_tuple("A", 0.90)
    snp_tuple.add_tuple("T", 1.0)
    return snp_tuple


class SNPTupleTest(unittest.TestCase):

    def setUp(self):
        self.snp_tuple = test_snp_tuple(100)

    def test_pick_nt(self):
        picked = self.snp_tuple.pick_snp_value(0.95)
        self.assertEqual("T", picked, "Should pick T")
        picked = self.snp_tuple.pick_snp_value(0.4)
        self.assertEqual("G", picked, "Should pick G")

    def test_pick_allele_index(self):
        picked = self.snp_tuple.pick_allele_index(0.95)
        self.assertEqual(2, picked, "Should pick index 2")
        picked = self.snp_tuple.pick_snp_value(0.4)
        self.assertEqual(0, picked, "Should pick 0")

    def test_pick_deleterious(self):
        picked = self.snp_tuple.pick_deleterious_value()
        self.assertEqual("A", picked, "Should pick A as deleterious (2nd most frequent)")


class DeleteriousGroupTest(unittest.TestCase):

    def setUp(self):
        self.snp_data = [(1, [test_snp_tuple(1), test_snp_tuple(2), test_snp_tuple(3), test_snp_tuple(4)])]

    def test_from_yml(self):
        pg = pop_factory.DeleteriousGroup.from_yml({"mutation_weights": [0.5, 0.5, 0.5],
                                                 "num_instances": 2,
                                                 "population_weight": 5,
                                                 "min_minor_allele_freq": 0.01}, self.snp_data, "groupA")
        self.assertEqual(5, pg.population_weight, "Population Weight should be 5")
        self.assertEqual(3, len(pg.deleterious), "Should have picked 3 deleterious")
        self.assertEqual(2, len(pg.select_mutations()))


class PopulationFactoryTest(unittest.TestCase):

    def setUp(self):
        self.snp_tuples = [test_snp_tuple(100), test_snp_tuple(101), test_snp_tuple(102), test_snp_tuple(103),
                           test_snp_tuple(104)]

    def test_deleterious_groups(self):
        groups = []
        for x in range(1, 3):
            groups.append(pop_factory.DeleteriousGroup.init_with_snps("group%i" % x, [0.5, 0.5], self.snp_tuples, 1))
        for x in range(3, 5):
            groups.append(pop_factory.DeleteriousGroup.init_with_snps("group%i" % x, [0.5, 0.5], self.snp_tuples, 10))

        selected_groups = pop_factory.PopulationFactory.pick_deleterious_groups(groups, 100)
        counts = {}
        for s in selected_groups:
            if not counts.get(s.name):
                counts[s.name] = 0
            counts[s.name] += 1
        self.assertGreater(counts["group3"], 2)
        self.assertLess(counts["group1"], 99)


class ArgParserTest(unittest.TestCase):

    def test_parse_args(self):
        cmd = "-s 10 -c 20 -n 5 -z 3 -p path_config.yml -f 0.1 " \
              + "-m 0.7 -x 2500 -l " \
              + "--deleterious_file /home/ochrzan/workspace/deleterious.json " \
              + "--offset 300 --snps_file my_snps.json --outdir myoutput/tuesday"
        cmds = cmd.split(" ")
        args = pop_factory.parse_cmd_args(cmds)
        self.assertEqual(10, args.size)
        self.assertEqual(20, args.control_size)
        self.assertEqual(5, args.num_processes)
        self.assertEqual(3, args.compression_level)
        self.assertEqual("path_config.yml", args.deleterious_config)
        self.assertEqual(0.1, args.min_freq)
        self.assertEqual(0.7, args.male_odds)
        self.assertEqual(2500, args.max_snps)
        self.assertEqual(True, args.generate_snps)
        self.assertEqual("/home/ochrzan/workspace/deleterious.json", args.deleterious_file)
        self.assertEqual(300, args.offset)
        self.assertEqual("my_snps.json", args.snps_file)
        self.assertEqual("myoutput/tuesday", args.outdir)

        cmd = "-s 10 -c 20 -x 2500"
        cmds = cmd.split(" ")
        args = pop_factory.parse_cmd_args(cmds)
        self.assertIsNone(args.deleterious_file)


if __name__ == '__main__':
    unittest.main()
