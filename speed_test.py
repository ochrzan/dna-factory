import pop_factory
from common.snp import RefSNP, Allele, obj_from_rowproxy
from common.db import db
from datetime import datetime

MIN_FREQ = 0.005


def speed_test():
    load_via_sql()
    load_via_json()


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


def load_via_json():
    start = datetime.now()
    print("%s JSON loading started" % start.strftime("%Y-%m-%d %H:%M:%S"))
    pf = pop_factory.PopulationFactory()
    pf.load_snps_json("test_snp", MIN_FREQ)
    end = datetime.now()
    print("%s JSON loading finished. %s elapsed" % (end.strftime("%Y-%m-%d %H:%M:%S"), str(end - start)))


if __name__ == '__main__':
    speed_test()
