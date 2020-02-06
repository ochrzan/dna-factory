import time
from random import random

import numpy
import pop_factory
import gzip
from common.timer import Timer


r = numpy.random.rand(10000)
snp_tuple = pop_factory.SNPTuples(100, "1", 50000)
snp_tuple.add_tuple("G", 0.70)
snp_tuple.add_tuple("A", 0.90)
snp_tuple.add_tuple("T", 1.0)

x = [100, 200]
start = time.perf_counter()
for i in range(10000000):
    with Timer("in_list") as t:
        y = 100 not in x
end = time.perf_counter()
print("elapsed %s" % str(end - start))
print(str(Timer.report_all()))

x = {100: 1, 200: 2}
start = time.perf_counter()
for i in range(10000000):
    with Timer("in_dict") as t:
        y = 100 not in x
end = time.perf_counter()
print("elapsed %s" % str(end - start))
print(str(Timer.report_all()))
