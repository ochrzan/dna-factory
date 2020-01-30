from sqlalchemy import Column, Integer, Float, String, ForeignKey, Table, MetaData, create_engine


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
        self.db_init('sqlite:///snps/refSNP.db')

    def bulk_insert(self, objs, table):
        insert_vals = list(map(lambda x: vars(x), objs))
        return self.connection.execute(table.insert(), insert_vals)


db = DbLayer()
