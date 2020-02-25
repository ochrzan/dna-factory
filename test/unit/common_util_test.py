import unittest

from common.snp import split_list


class SplitListTest(unittest.TestCase):
    def test_split_list(self):
        x = list(range(100))
        y = list(split_list(x, 3))
        self.assertEqual(len(y), 3)
        self.assertEqual(len(y[0]), 33)
        self.assertEqual(len(y[2]), 34)


if __name__ == '__main__':
    unittest.main()
