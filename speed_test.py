import concurrent
import gzip
import glob
import re
import pop_factory
from common.snp import RefSNP, Allele, chromosome_from_filename
from common import Session, db_engine
from datetime import datetime

MIN_FREQ = 0.005

def speed_test():
    load_via_db()
    load_via_json()


def load_via_db():
    start = datetime.now()
    print("%s DB loading started" % start.strftime("%Y-%m-%d %H:%M:%S") )
    session = Session()
    snps = list(session.query(RefSNP).filter(RefSNP.maf >= MIN_FREQ))
    end = datetime.now()
    print("%s DB loading finished. %s elapsed" % (end.strftime("%Y-%m-%d %H:%M:%S"), str(end - start)))
    print("%i snps loaded from the DB" % len(snps))
    return snps

def load_via_json():
    start = datetime.now()
    print("%s JSON loading started" % start.strftime("%Y-%m-%d %H:%M:%S") )
    pf = pop_factory.PopulationFactory()
    pf.load_snps("test_snp", MIN_FREQ)
    end = datetime.now()
    print("%s JSON loading finished. %s elapsed" % (end.strftime("%Y-%m-%d %H:%M:%S"), str(end - start)))



if __name__ == '__main__':
    speed_test()