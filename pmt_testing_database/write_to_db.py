import psycopg2

if __name__=='__main__':

    DB_HOST = 'localhost'
    DB_NAME = 'pmt_testing'
    DB_USER = 'db_user'
    DB_PASS = 'pmt'
    DB_PORT = 5432

    conn = psycopg2.connect('host=%s dbname=%s user=%s password=%s port=%d' % \
                            (DB_HOST, DB_NAME, DB_USER, DB_PASS, DB_PORT))

    cursor = conn.cursor()

    command = 'SELECT * FROM pmt_information'

    cursor.execute(command)
    row = cursor.fetchall()

    print(row)
