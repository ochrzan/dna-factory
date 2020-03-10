import os

from sqlalchemy import Column, Integer, Float, String, ForeignKey, Table, MetaData, create_engine
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

from definitions import ROOT_DIR


class DbLayer:

    def __init__(self):
        self.connection = None
        self.engine = None
        self.conn_string = None
        self.metadata = MetaData()

        self.ref_snps = Table('ref_snps', self.metadata,
                              Column('id', Integer, primary_key=True),
                              Column('chromosome', String),
                              Column('maf', Float, index=True),
                              Column('total_count', Integer)
                              )

        self.alleles = Table('alleles', self.metadata,
                             Column('id', Integer, primary_key=True),
                             Column('deleted', String),
                             Column('inserted', String),
                             Column('position', Integer),
                             Column('allele_count', Integer, index=True),
                             Column('ref_snp_id', Integer, ForeignKey('ref_snps.id'), nullable=False)
                             )

    def db_init(self, conn_string):
        self.conn_string = conn_string
        self.engine = create_engine(conn_string or self.conn_string)
        self.metadata.create_all(self.engine)
        self.connection = self.engine.connect()

    def default_init(self):
        connect_string = 'sqlite:////' + os.path.join(ROOT_DIR, 'snps/refSNP.db')
        db_yaml_file = os.path.join(ROOT_DIR, "db.yml")
        if os.path.exists(db_yaml_file):
            with open(db_yaml_file, 'r') as p:
                db_yml = load(p, Loader=Loader)
                if db_yml.get("connection_string"):
                    connect_string = db_yml["connection_sring"]
        self.db_init(connect_string)

    def bulk_insert(self, objs, table):
        if not objs or table is None:
            return 0
        insert_vals = list(map(lambda x: vars(x), objs))
        return self.connection.execute(table.insert(), insert_vals)


db = DbLayer()
