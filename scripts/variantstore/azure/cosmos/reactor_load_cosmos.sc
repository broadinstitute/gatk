import $ivy.`com.azure:azure-cosmos:4.41.0`
import $ivy.`com.google.code.gson:gson:2.10.1`
import $ivy.`com.lihaoyi:ammonite-ops_2.13:2.4.1`
import $ivy.`org.apache.avro:avro:1.11.1`
import $ivy.`org.xerial.snappy:snappy-java:1.1.8.4`


import ammonite.ops._
import com.azure.cosmos._
import com.azure.cosmos.models._
import com.google.gson._
import java.io.File
import java.util.NoSuchElementException
import org.apache.avro.file._
import org.apache.avro.generic._
import reactor.core.publisher.Flux
import reactor.core.publisher.Mono
import scala.jdk.CollectionConverters._


// Bulk Cosmos writes from Java
// https://learn.microsoft.com/en-us/azure/cosmos-db/nosql/bulk-executor-java
@main
def main(database: String, container: String, avrodir: String): Unit = {
  val (endpoint, key) = extractCosmosEndpointAndKey()
  val client = buildClient(endpoint, key)
  val cosmosContainer = client.getDatabase(database).getContainer(container)

  val avro_paths = determineAvroPaths(avrodir)

  processAvros(cosmosContainer, avro_paths)
}


def extractCosmosEndpointAndKey(): (String, String) = {
  try {
    (sys.env("COSMOS_ENDPOINT"), sys.env("COSMOS_KEY"))
  } catch {
    case e: NoSuchElementException =>
      println(f"Error: ${e.getMessage()} environment variable not set.")
      sys.exit(1)
  }
}

def buildClient(endpoint: String, key: String): CosmosAsyncClient = {
  new CosmosClientBuilder().
    endpoint(endpoint).
    key(key).
    preferredRegions(List("East US").asJava).
    contentResponseOnWriteEnabled(true).
    consistencyLevel(ConsistencyLevel.SESSION).
    buildAsyncClient()
}


def determineAvroPaths(avrodir: String): Seq[Path] = {
  val path = pwd / RelPath(avrodir)
  val listing = ls ! path

  listing.filter(_.isFile)
}


// https://stackoverflow.com/a/45648136/21269164
def processAvros(container: CosmosAsyncContainer, avro_paths: Iterable[Path]): Unit = {
  val gson = new Gson()
  var id = 1

  val itemOperations = Flux.fromIterable(avro_paths.asJava).
    flatMap(path => {
      val file = new File(path.toString())
      val reader = new GenericDatumReader()
      val dataFileReader = new DataFileReader(file, reader)
      Flux.fromIterable(dataFileReader).
        // I had to put this `map` here inside the `flatMap` and not directly ahead of the `take` otherwise nothing
        // happened; I don't know why.
        map(r => gson.fromJson(r.toString(), classOf[Vet]))
    }).
    take(10000).
    map(r => { r.asInstanceOf[Vet].id = "" + id; id = id + 1; r }).
    map(r => CosmosBulkOperations.getCreateItemOperation(r, new PartitionKey(r.getSample_id())))

  executeItemOperationsWithErrorHandling(container, itemOperations)

  // flux.subscribe(e => System.out.println(e))
}


def executeItemOperationsWithErrorHandling(container: CosmosAsyncContainer, itemOperations: Flux[CosmosItemOperation]): Unit = {
  container.executeBulkOperations(itemOperations).flatMap(operationResponse => {
    val itemResponse = operationResponse.getResponse()
    val itemOperation = operationResponse.getOperation()

    if (operationResponse.getException() != null) {
      println(f"Bulk operation failed: ${operationResponse.getException()}")
    } else if (itemResponse == null || !itemResponse.isSuccessStatusCode()) {
      println(String.format(
        "The operation for Item ID: [%s]  Item PartitionKey Value: [%s] did not complete " +
          "successfully with " + "a" + " %s/%s response code, response: %s.",
        itemOperation.getItem().asInstanceOf[Vet].id,
        itemOperation.getItem().asInstanceOf[Vet].sample_id,
        if (itemResponse != null) itemResponse.getStatusCode() else "n/a",
        if (itemResponse != null) itemResponse.getSubStatusCode() else "n/a",
        itemResponse.getCosmosDiagnostics()));
    }

    if (itemResponse == null) {
      Mono.error(new IllegalStateException("No response retrieved."));
    } else {
      Mono.just(itemResponse);
    }

  }).blockLast()
}


case class Vet(
                val sample_id: Long                 = 0   ,
                val location: Long                  = 0   ,
                val ref: String                     = null,
                val alt: String                     = null,
                val AS_RAW_MQ: String               = null,
                val AS_RAW_MQRankSum: String        = null,
                val QUALapprox: String              = null,
                val AS_QUALapprox: String           = null,
                val AS_RAW_ReadPosRankSum: String   = null,
                val AS_SB_TABLE: String             = null,
                val AS_VarDP: String                = null,
                val call_GT: String                 = null,
                val call_AD: String                 = null,
                val call_GQ: String                 = null,
                val call_PGT: String                = null,
                val call_PID: String                = null,
                val call_PL: String                 = null
              ) {
  var id: String = null

  def getId(): String = id
  def getSample_id(): Long = this.sample_id
  def getLocation(): Long = this.location
  def getRef(): String = this.ref
  def getAlt(): String = this.alt
  def getAS_RAW_MQ(): String = this.AS_RAW_MQ
  def getAS_RAW_MQRankSum(): String = this.AS_RAW_MQRankSum
  def getQUALapprox(): String = this.QUALapprox
  def getAS_RAW_ReadPosRankSum(): String = this.AS_RAW_ReadPosRankSum
  def getAS_SB_TABLE(): String = this.AS_SB_TABLE
  def getAS_VarDP(): String = this.AS_VarDP
  def getCall_GT(): String = this.call_GT
  def getCall_AD(): String = this.call_AD
  def getCall_GQ(): String = this.call_GQ
  def getCall_PGT(): String = this.call_PGT
  def getCall_PID(): String = this.call_PID
  def getCall_PL(): String = this.call_PL
}