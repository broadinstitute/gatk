import argparse
import pyodbc
import struct

# All taken from
# https://techcommunity.microsoft.com/t5/apps-on-azure-blog/how-to-connect-azure-sql-database-from-python-function-app-using/ba-p/3035595


def build_token_struct(database_access_token):
    with open(database_access_token) as database_token_file:
        database_token_string = database_token_file.read().rstrip()
        exp_token = b''
        for i in bytes(database_token_string, "UTF-8"):
            exp_token += bytes({i})
            exp_token += bytes(1)

        return struct.pack("=i", len(exp_token)) + exp_token


def build_connection_string(server, database):
    driver = "{ODBC Driver 18 for SQL Server}"
    return f'DRIVER={driver};SERVER={server}.database.windows.net;DATABASE={database}'


def connect(connection_string, token_struct):
    SQL_COPT_SS_ACCESS_TOKEN = 1256
    return pyodbc.connect(connection_string, attrs_before={SQL_COPT_SS_ACCESS_TOKEN: token_struct})


def query_and_print(connection):
    query = """

        select @@version as "Hello Azure SQL Database!"

    """

    cursor = connection.cursor()
    cursor.execute(query)
    row = cursor.fetchone()

    while row:
        print(row[0])
        row = cursor.fetchone()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Say Hello to Azure SQL Database from Python')
    parser.add_argument('--sql-server', type=str, help='Azure SQL Server name', required=True)
    parser.add_argument('--sql-database', type=str, help='Azure SQL Server database', required=True)
    parser.add_argument('--database-access-token', type=str, help='Azure SQL Database access token', required=True)
    args = parser.parse_args()

    connection_string = build_connection_string(args.sql_server, args.sql_database)
    token_struct = build_token_struct(args.database_access_token)

    connection = connect(connection_string, token_struct)
    query_and_print(connection)

