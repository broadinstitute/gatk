from azure.identity import DefaultAzureCredential

import argparse
import pyodbc
import struct


def read_token_from_file(token_file_name):
    # https://learn.microsoft.com/en-us/azure/app-service/tutorial-connect-msi-azure-database?tabs=sqldatabase%2Csystemassigned%2Cpython%2Cwindowsclient#3-modify-your-code
    with open(token_file_name) as token_file:
        token_str = token_file.read().rstrip().encode("UTF-16-LE")
        return token_str


def fetch_token():
    credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)
    token_str = credential.get_token("https://database.windows.net/.default").token.encode("UTF-16-LE")
    return token_str


def token_to_struct(token_str: bytes):
    token_struct = struct.pack(f'<I{len(token_str)}s', len(token_str), token_str)
    return token_struct


def build_connection_string(server, database):
    driver = "{ODBC Driver 18 for SQL Server}"
    return f'DRIVER={driver};SERVER={server}.database.windows.net;DATABASE={database}'


def connect_via_token(connection_string, token_struct):
    # PEP8 has no appreciation for the beauty of ALL_CAPS symbolic constants in functions.
    # noinspection PyPep8Naming
    SQL_COPT_SS_ACCESS_TOKEN = 1256
    return pyodbc.connect(connection_string, attrs_before={SQL_COPT_SS_ACCESS_TOKEN: token_struct})


def connect_via_msi(connection_string):
    return pyodbc.connect(connection_string + ";Authentication=ActiveDirectoryMsi")


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
    # All taken from
    # https://techcommunity.microsoft.com/t5/apps-on-azure-blog/how-to-connect-azure-sql-database-from-python-function-app-using/ba-p/3035595
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Say Hello to Azure SQL Database from Python')
    parser.add_argument('--server', type=str, help='Azure SQL Server name', required=True)
    parser.add_argument('--database', type=str, help='Azure SQL Server database', required=True)
    parser.add_argument('--token-file', type=str, help='Azure SQL Database access token', required=False)
    parser.add_argument('--msi-auth', type=bool,
                        help='Use MSI (Managed Service Identity) Authentication.',
                        required=False, default=False)
    args = parser.parse_args()

    connection_string = build_connection_string(args.server, args.database)

    if args.token_file:
        token = read_token_from_file(args.token_file)
        token_struct = token_to_struct(token)
        connection = connect_via_token(connection_string, token_struct)
    elif args.msi_auth:
        connection = connect_via_msi(connection_string)
    else:
        token = fetch_token()
        token_struct = token_to_struct(token)
        connection = connect_via_token(connection_string, token_struct)

    query_and_print(connection)
