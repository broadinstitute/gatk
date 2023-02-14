version 1.0

workflow HelloAzure {
    input {
        String sql_server
        String sql_database
        File access_token
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

        # sqlcmd actually cannot cope with the 141 character natural length of this input path... ln to the rescue
        # Sqlcmd: '-P': Argument too long (maximum is 128 characters).
        ln -s ~{access_token} /tmp/db_access_token.txt

        sqlcmd -S tcp:~{sql_server}.database.windows.net,1433 -d ~{sql_database} -G -Q 'select @@version as "Hello Azure SQL Database!"' -P /tmp/db_access_token.txt
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:coa-2023-02-09"
    }
    output {
        String out = read_string(stdout())
    }
}
