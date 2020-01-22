from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

db_engine = create_engine('sqlite:///snps/refSNP.db')

# create a Session class
Session = sessionmaker(bind=db_engine)
