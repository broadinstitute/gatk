version 1.0

workflow HelloAzure {
    input {
        String sql_server
        String sql_database
        File database_access_token
        File python_script
    }
    parameter_meta {
        sql_server: {
            description: "Name of the Azure SQL Database Server without .database.windows.net suffix"
        }
        sql_database: {
            description: "Name of the Database within the Azure SQL Database Server"
        }
        database_access_token: {
            description: "An access token generated with auth that can access the Azure SQL Database Server. e.g. `az account get-access-token --resource https://database.windows.net --output tsv | cut -f 1 | tr -d '\n' | iconv -f ascii -t UTF-16LE > db_access_token.txt`"
        }
    }

    call HelloFromSqlcmd {
        input:
            sql_server = sql_server,
            sql_database = sql_database,
            database_access_token = database_access_token
    }

    call HelloFromPython {
        input:
            sql_server = sql_server,
            sql_database = sql_database,
            python_script = python_script
    }

    output {
        String sqlcmd_out = HelloFromSqlcmd.out
        String python_out = HelloFromPython.out
    }
}

task HelloFromSqlcmd {
    input {
        String sql_server
        String sql_database
        File database_access_token
    }
    meta {
        # Database token generated from a regular `az login` or `az login --identity` MSI plus:
        #
        # az account get-access-token --resource https://database.windows.net --output tsv | cut -f 1 | tr -d '\n' | iconv -f ascii -t UTF-16LE > db_access_token.txt
        #
        # These tokens have ~1 hour TTL and threading them through explicitly isn't great, but they are enough for
        # "Hello World" purposes.
        description: "Say hello to Azure SQL Database from sqlcmd using a database access token"
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # smh sqlcmd actually refuses to deal with the 141 character length of this input path... ln to the rescue!
        # Sqlcmd: '-P': Argument too long (maximum is 128 characters).
        ln -s ~{database_access_token} /tmp/db_access_token.txt

        # A hopefully temporary hack to allow Azure Batch VMs to connect to the Azure SQL Database Server. We don't know
        # the Batch VM IPs in advance and they don't currently seem to have the auth to allow-list only their own IPs.
        # Run this on your personal machine before attempting to run this workflow.
        # https://learn.microsoft.com/en-us/azure/azure-sql/database/firewall-configure?view=azuresql#connections-from-inside-azure
        #
        # az sql server firewall-rule create -n AllowAllWindowsAzureIps --start-ip-address 0.0.0.0 --end-ip-address 0.0.0.0

        sqlcmd -S tcp:~{sql_server}.database.windows.net,1433 -d ~{sql_database} -G -Q 'select @@version as "Hello Azure SQL Database!"' -P /tmp/db_access_token.txt
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-16"
    }
    output {
        String out = read_string(stdout())
    }
}

task HelloFromPython {
    input {
        String sql_server
        String sql_database
        File python_script
    }
    meta {
        # This task is not currently working under CoA as the Azure Batch VMs spun up by CoA do not appear to have any
        # assigned identity. There is an additional code path in the script to support the use of access tokens but that
        # isn't working either for unknown reasons.
        description: "Say hello to Azure SQL Database from Python -> pyodbc -> unixodbc -> MS ODBC driver"
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        python3 ~{python_script} --sql-server ~{sql_server} --sql-database ~{sql_database}
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-16"
    }
    output {
        String out = read_string(stdout())
    }
}
