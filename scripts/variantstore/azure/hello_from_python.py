import argparse
import pyodbc
import struct


def build_token_struct(access_token):
    with open(access_token) as database_token_file:
        database_token_string = database_token_file.read()
        exp_token = b''
        for i in bytes(database_token_string, "UTF-8"):
            exp_token += bytes({i})
            exp_token += bytes(1)

        return struct.pack("=i", len(exp_token)) + exp_token


def build_connection_string(server, database):
    driver = "{ODBC Driver 18 for SQL Server}"
    return f'DRIVER={driver};SERVER={server}.database.windows.net;DATABASE={database}'


def connect_access_token(connection_string, token_struct):
    # PEP8 has no appreciation for the beauty of ALL_CAPS symbolic constants in functions.
    # noinspection PyPep8Naming
    SQL_COPT_SS_ACCESS_TOKEN = 1256
    return pyodbc.connect(connection_string, attrs_before={SQL_COPT_SS_ACCESS_TOKEN: token_struct})


def connect_msi(connection_string):
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
    #
    # The default MSI (Managed System Identity) branch works fine when I manually allocate a VM and assign a UAMI (User
    # Assigned Managed Identity) to it. Unfortunately this doesn't work with the Azure Batch VMs spun up by CoA which
    # don't seem to be assigned managed identities.
    #
    # I have not had any luck with the access token branch either authed as myself or a managed identity and the error
    # messages are completely unhelpful. It's not a firewall issue as I can connect with username / password, and it's
    # not a token issue as I can connect with sqlcmd from the same machine using the same access token. I also don't
    # think this is an ODBC issue as both Python (via pyodbc via unixodbc) and sqlcmd seem to be using the same
    # underlying "{ODBC Driver 18 for SQL Server}" driver. So the problem would appear to be in the way the token is
    # being passed into or interpreted by the pyodbc / unixodbc layers.
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Say Hello to Azure SQL Database from Python')
    parser.add_argument('--sql-server', type=str, help='Azure SQL Server name', required=True)
    parser.add_argument('--sql-database', type=str, help='Azure SQL Server database', required=True)
    parser.add_argument('--msi-auth', type=bool,
                        help='Use MSI (Managed Service Identity) Authentication. If false, --access-token is required.',
                        required=False, default=True)
    parser.add_argument('--access-token', type=str, help='Azure SQL Database access token', required=False)
    args = parser.parse_args()

    connection_string = build_connection_string(args.sql_server, args.sql_database)
    print(connection_string)

    if args.access_token:
        token_struct = build_token_struct(args.access_token)
        connection = connect_access_token(connection_string, token_struct)
    elif args.msi_auth:
        connection = connect_msi(connection_string)
    else:
        raise ValueError("Must specify --access-token <token file> or --msi-auth True")

    query_and_print(connection)
