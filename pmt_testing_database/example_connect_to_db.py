import sqlalchemy
import sys

if __name__=='__main__':

    key = sys.argv[1]

    engine = sqlalchemy.create_engine('postgresql://%s:%s@%s:%i/%s' %
                                      ('postgres', 'b33feroni',
                                       'localhost', 5432, 'pmt_testing'),
                                       pool_recycle=3600)

    conn = engine.connect()

    command = ('SELECT * FROM pmt_information WHERE key=%d' % key)

    result = conn.execute(command)
    row = result.fetchone()
    keys = result.keys()

    print(row)

