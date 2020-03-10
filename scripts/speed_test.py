import subprocess

import numpy
import common.timer as timer
import gzip

from Bio import bgzf
from common.snp import RefSNP, Allele
from common.db import db
from datetime import datetime

MIN_FREQ = 0.005


def speed_test():
    load_via_sql()


def load_via_sql():
    start = datetime.now()
    db.default_init()
    print("%s SQL loading started" % start.strftime("%Y-%m-%d %H:%M:%S"))
    query = db.ref_snps.select().where(db.ref_snps.columns.maf >= MIN_FREQ)
    print(str(query.compile(db.engine, compile_kwargs={"literal_binds": True})))
    result = db.connection.execute(query)
    ref_snps = {}
    for row in result:
        snp = RefSNP.from_row_proxy(row)
        ref_snps[snp.id] = snp
    print("%s RefSNP Query complete" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    allele_query = "Select ref_snp_id, deleted, inserted, position, allele_count from ref_snps r  " \
                   "join alleles a on r.id = a.ref_snp_id and r.maf >= %f" \
                   " " % MIN_FREQ
    result = db.connection.execute(allele_query)
    for row in result:
        allele = Allele.from_row_proxy(row)
        ref_snps[row.ref_snp_id].put_allele(allele)
    end = datetime.now()
    print("%s DB loading finished. %s elapsed" % (end.strftime("%Y-%m-%d %H:%M:%S"), str(end - start)))
    print("%i snps loaded from the DB" % len(ref_snps))
    return ref_snps


def gzip_speed():

    gzip_pipe = subprocess.Popen(args="gzip -c > tmp_file.gz", shell=True, stdin=subprocess.PIPE)
    randos = []
    for i in range(20000):
        rands = numpy.random.rand(300)
        string = " ".join(map(lambda x: str(x), rands))
        randos.append(string)

    with timer.Timer(logger=print, name="OS GZip") as t:

        for r in randos:
            gzip_pipe.stdin.write(r.encode())
        gzip_pipe.stdin.close()
        gzip_pipe.wait()

    with timer.Timer(logger=print, name="GZipLib") as t, gzip.open("tmp2_file.gz", 'wt+', compresslevel=4) as f:
        for r in randos:
            f.write(r)

    with timer.Timer(logger=print, name="BGZipLib") as t, bgzf.BgzfWriter("tmp3_file.gz", 'wt+', compresslevel=4) as f:
        for r in randos:
            f.write(r)


if __name__ == '__main__':
    gzip_speed()
