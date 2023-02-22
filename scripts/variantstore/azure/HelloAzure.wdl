version 1.0

workflow HelloAzure {
    input {
        String sql_server
        String sql_database
        File utf8_token_file
        File python_script
        File ammonite_script
    }
    parameter_meta {
        sql_server: {
            description: "Name of the Azure SQL Database Server without .database.windows.net suffix"
        }
        sql_database: {
            description: "Name of the Database within the Azure SQL Database Server"
        }
        token_file: {
            description: "A file with a UTF-8 encoded access token generated with auth that can access the Azure SQL Database Server. e.g. `az account get-access-token --resource=https://database.windows.net/ --query accessToken --output tsv > db_access_token.txt"
        }
    }

    call HelloFromSqlcmd {
        input:
            sql_server = sql_server,
            sql_database = sql_database,
            token_file = utf8_token_file
    }

    call HelloFromPython {
        input:
            sql_server = sql_server,
            sql_database = sql_database,
            python_script = python_script,
            token_file = utf8_token_file
    }

    call HelloFromAmmonite {
        input:
            sql_server = sql_server,
            sql_database = sql_database,
            ammonite_script = ammonite_script,
            token_file = utf8_token_file
    }
}

task HelloFromSqlcmd {
    input {
        String sql_server
        String sql_database
        File token_file
    }
    meta {
        description: "Say hello to Azure SQL Database from sqlcmd using a database access token"
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # sqlcmd is particular about the formatting and encoding of its access token: no whitespace and UTF-16LE.
        # Python is particular too but these manipulations are sprinkled into the code. Java / Ammonite doesn't
        # seem to care about encoding or autodetects and adapts?
        cat ~{token_file} | cut -f 1 | tr -d '\n' | iconv -f ascii -t UTF-16LE > /tmp/db_access_token.txt

        # A hopefully temporary hack to allow Azure Batch VMs to connect to the Azure SQL Database Server. We don't know
        # the Batch VM IPs in advance and they don't currently seem to have the auth to allow-list only their own IPs.
        # Run this on your personal machine before attempting to run this workflow.
        # https://learn.microsoft.com/en-us/azure/azure-sql/database/firewall-configure?view=azuresql#connections-from-inside-azure
        #
        # az sql server firewall-rule create -n AllowAllWindowsAzureIps --start-ip-address 0.0.0.0 --end-ip-address 0.0.0.0

        sqlcmd -S tcp:~{sql_server}.database.windows.net,1433 -d ~{sql_database} -G -Q 'select @@version as "Hello Azure SQL Database!"' -P /tmp/db_access_token.txt
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-22"
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
        File token_file
    }
    meta {
        description: "Say hello to Azure SQL Database from Python -> pyodbc -> unixodbc -> MS ODBC driver"
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        python3 ~{python_script} --server ~{sql_server} --database ~{sql_database} --token-file ~{token_file}
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-22"
    }
    output {
        String out = read_string(stdout())
    }
}

task HelloFromAmmonite {
    input {
        String sql_server
        String sql_database
        File ammonite_script
        File token_file
    }
    meta {
        description: "Say hello to Azure SQL Database from Ammonite/Java -> JDBC -> MS JDBC driver"
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        amm ~{ammonite_script} --sql-server ~{sql_server} --sql-database ~{sql_database} --tokenFile ~{token_file}
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-22"
    }
    output {
        String out = read_string(stdout())
    }
}
