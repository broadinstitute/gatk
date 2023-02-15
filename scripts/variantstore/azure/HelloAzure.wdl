version 1.0

workflow HelloAzure {
    input {
        String sql_server
        String sql_database
        File access_token
    }
    parameter_meta {
        sql_server: {
            description: "Name of the Azure SQL Database Server"
        }
        sql_database: {
            description: "Name of the Database within the Azure SQL Database Server"
        }
        access_token: {
            description: "An access token generated with auth that can access the Azure SQL Database Server. e.g. `az account get-access-token --resource https://database.windows.net --output tsv | cut -f 1 | tr -d '\n' | iconv -f ascii -t UTF-16LE > db_access_token.txt`"
        }
    }

    call Hello {
        input:
            sql_server = sql_server,
            sql_database = sql_database,
            access_token = access_token
    }

    output {
        String out = Hello.out
    }
}

task Hello {
    input {
        String sql_server
        String sql_database
        File access_token
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # smh sqlcmd actually refuses to deal with the 141 character length of this input path... ln to the rescue!
        # Sqlcmd: '-P': Argument too long (maximum is 128 characters).
        ln -s ~{access_token} /tmp/db_access_token.txt

        # A hopefully temporary hack to allow Azure Batch VMs to connect to the Azure SQL Database Server. We don't know
        # the Batch VM IPs in advance and they don't currently seem to have the auth to allow-list only their own IPs.
        # Run this on your personal machine before attempting to run this workflow.
        # https://learn.microsoft.com/en-us/azure/azure-sql/database/firewall-configure?view=azuresql#connections-from-inside-azure
        #
        # az sql server firewall-rule create -n AllowAllWindowsAzureIps --start-ip-address 0.0.0.0 --end-ip-address 0.0.0.0

        sqlcmd -S tcp:~{sql_server}.database.windows.net,1433 -d ~{sql_database} -G -Q 'select @@version as "Hello Azure SQL Database!"' -P /tmp/db_access_token.txt
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-09"
    }
    output {
        String out = read_string(stdout())
    }
}
