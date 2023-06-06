import sqlalchemy
import sys

def connect():

    engine = sqlalchemy.create_engine('postgresql://%s:%s@%s:%i/%s' %
                                      ('postgres', 'b33feroni',
                                       'localhost', 5432, 'pmt_testing'),
                                       pool_recycle=3600)

    conn = engine.connect()
    return conn


def select_info(conn, key):

    command = ('SELECT tts_sigma, dark_rate FROM pmt_information WHERE key=%d' % (key,))

    result = conn.execute(command)
    row = result.fetchall()

    for tts_sigma, dark_rate in row:
        print "TTS (sigma):", tts_sigma
        print "TTS (FWHM):", tts_sigma*2.335
        print "Dark-rate (Hz):", dark_rate


if __name__=='__main__':

    key = int(sys.argv[1])

    conn = connect()

    select_info(conn, key)

