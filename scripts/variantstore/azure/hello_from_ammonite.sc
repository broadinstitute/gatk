// Ammonite script to say Hello to Azure SQL Database from the Java ecosystem.

// Package management
import $ivy.`com.microsoft.sqlserver:mssql-jdbc:12.2.0.jre11`
import $ivy.`com.azure:azure-identity:1.4.6`
// ANTLR 4 appears to be required by one of the above Microsoft packages but the dependency is not expressed explicitly
// so it is not imported automatically. Without adding this import ourselves Ammonite fails to compile this script.
import $ivy.`org.antlr:antlr4:4.12.0`

// Imports
import com.azure.core.credential.*
import com.azure.identity.*
import com.microsoft.sqlserver.jdbc.SQLServerDataSource
import java.nio.charset.StandardCharsets
import java.nio.file.*
import java.sql.*
import java.util.*


// Nearly everything taken from
// https://learn.microsoft.com/en-us/azure/app-service/tutorial-connect-msi-azure-database?tabs=sqldatabase%2Csystemassigned%2Cjava%2Cwindowsclient#3-modify-your-code
def getAccessTokenViaRequest(): String = {
    val creds = new DefaultAzureCredentialBuilder().build()
    val request = new TokenRequestContext()
    request.addScopes("https://database.windows.net//.default");
    val accessToken = creds.getToken(request).block()
    accessToken.getToken()
}

// Generate a token via
// az account get-access-token --resource=https://database.windows.net/ --query accessToken --output tsv > db_access_token.txt
// Note this produces a token file with a confounding trailing newline, the code below has to `trim()`.
// Also note that unlike the sqlcmd and Python contexts, there is nothing here about UTF-16LE encoding; the Java
// ecosystem seems to deal with the ASCII / UTF-8 access token just fine without the caller doing anything special.
def getAccessTokenViaFile(tokenFile: String): String = {
    // https://www.digitalocean.com/community/tutorials/java-read-file-to-string
    val token = new String(Files.readAllBytes(Paths.get(tokenFile)))
    token.trim()
}

@main
def main(server: String, database: String, tokenFile: Option[String] = None) = {
    val ds = new SQLServerDataSource()

    val token = tokenFile match {
      case Some(file) => getAccessTokenViaFile(file)
      case None => getAccessTokenViaRequest()
    }
    ds.setAccessToken(token)
    ds.setServerName(s"${server}.database.windows.net")
    ds.setDatabaseName(database)

    val connection = ds.getConnection()
    val statement = connection.createStatement()

    val resultSet = statement.executeQuery("""

    select @@version as "Hello Azure SQL Database!"

    """)

    resultSet.next()
    val result = resultSet.getString(1)
    print(result)
}
