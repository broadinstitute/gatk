import $ivy.`ch.qos.logback:logback-core:1.4.5`
import $ivy.`ch.qos.logback:logback-classic:1.4.5`
import $ivy.`com.azure:azure-cosmos:4.41.0`
import $ivy.`com.lihaoyi:ammonite-ops_2.13:2.4.1`
import $ivy.`org.apache.avro:avro:1.11.1`
import $ivy.`org.slf4j:slf4j-api:2.0.6`
import $ivy.`org.xerial.snappy:snappy-java:1.1.8.4`


import ammonite.ops._
import com.azure.cosmos._
import com.azure.cosmos.models._
import com.fasterxml.jackson.databind.{JsonNode, ObjectMapper}
import com.fasterxml.jackson.databind.node.ObjectNode
import java.io.File
import java.util.NoSuchElementException
import java.util.concurrent.atomic.AtomicLong
import org.apache.avro.file._
import org.apache.avro.generic._
import org.slf4j.LoggerFactory
import reactor.core.publisher.Flux
import reactor.core.publisher.Mono
import scala.jdk.CollectionConverters._


val logger = LoggerFactory.getLogger("cosmos_load_data")


// Script to load vet and ref_ranges Avro files into Cosmos at the specified database / container coordinates.
// This script expects the environment variables `COSMOS_ENDPOINT` and `COSMOS_KEY` to be set per the instructions in
// `intro_to_azure_cosmos.md`.
//
// The Avro records are assumed to contain a String field `sample_id` which will be used as the Cosmos partition key; no
// other assumptions are made about the shape of the input data. The target container should be empty as this script
// automatically generates an `id` field which is a stringified Long starting from 1.
//
// Bulk Cosmos writes from Java
// https://learn.microsoft.com/en-us/azure/cosmos-db/nosql/bulk-executor-java
//
// Much of the code here is borrowed from / inspired by the sample application described in the link above.
@main
def main(database: String, container: String, avro_dir: String, num_records: Long=100000L, num_progress: Long=10000L): Unit = {
  configureLogging()
  val (endpoint, key) = extractCosmosEndpointAndKey()
  val client = buildClient(endpoint, key)
  val cosmosContainer = client.getDatabase(database).getContainer(container)

  val avroPaths = determineAvroPaths(avro_dir)

  loadAvros(cosmosContainer, avroPaths, num_records, num_progress)
}


def extractCosmosEndpointAndKey(): (String, String) = {
  try {
    (sys.env("COSMOS_ENDPOINT"), sys.env("COSMOS_KEY"))
  } catch {
    case e: NoSuchElementException =>
      logger.error(f"Error: ${e.getMessage()} environment variable not set.")
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


def determineAvroPaths(avroDirectory: String): Seq[Path] = {
  val path = pwd / RelPath(avroDirectory)
  val listing = ls ! path

  listing.filter(_.isFile)
}


// https://stackoverflow.com/a/45648136/21269164
def loadAvros(container: CosmosAsyncContainer, avroPaths: Iterable[Path], numRecordsToLoad: Long, numRecordsProgress: Long): Unit = {
  // There can occasionally be 409 race conditions on `id` if using a regular `Long` so go atomic!
  val id = AtomicLong(0L)

  // Where the Cosmos JSON serialization magic happens:
  // https://github.com/Azure/azure-sdk-for-java/blob/80b12e48aeb6ad2f49e86643dfd7223bde7a9a0c/sdk/cosmos/azure-cosmos/src/main/java/com/azure/cosmos/implementation/JsonSerializable.java#L255
  //
  // Making a Jackson `ObjectNode` (subtype of `JsonNode`) from a `GenericRecord#toString` JSON String seems like the
  // least bad option for Cosmos serialization (`JsonSerializable` looks like a Jackson type but is actually an
  // internal Cosmos type). This approach still necessitates an unsavory amount of bit shuffling from
  // Avro => Stringified JSON => Jackson object => Stringified JSON, but it certainly beats using brittle, clunky POJOs.
  // And despite all the shuffling this script generates Cosmos items far faster than we can actually push them to
  // Cosmos with 100K RU/s container bandwidth configuration.
  val objectMapper = new ObjectMapper()

  def itemFluxFromAvroPath(path: Path): Flux[CosmosItemOperation] = {
    val file = new File(path.toString())
    val reader = new GenericDatumReader()
    val dataFileReader = new DataFileReader(file, reader)
    Flux.fromIterable(dataFileReader).
      map(record => {
        val (objectNode, idLong) = objectNodeFromAvroRecord(objectMapper, record.toString(), id)
        val itemOperation =
          CosmosBulkOperations.getCreateItemOperation(objectNode, new PartitionKey(objectNode.get("sample_id").asLong))
        if (idLong % numRecordsProgress == 0L) logger.info(f"$idLong...");
        itemOperation
      })
  }

  for {
    avroPath <- avroPaths
    _ = logger.info(s"Processing Avro file '${avroPath}'...")
    itemFlux: Flux[CosmosItemOperation] = itemFluxFromAvroPath(avroPath)
    _ = executeItemOperationsWithErrorHandling(container, itemFlux)
    _ = logger.info(s"Avro file '${avroPath}' processing complete.")
  } yield ()
}


def executeItemOperationsWithErrorHandling(container: CosmosAsyncContainer, itemOperations: Flux[CosmosItemOperation]): Unit = {
  // Only the first and last few lines are the "execute" bits, all the rest is error handling iff something goes wrong.
  container.executeBulkOperations(itemOperations).flatMap(operationResponse => {
    val itemResponse = operationResponse.getResponse()
    val itemOperation = operationResponse.getOperation()

    if (operationResponse.getException() != null) {
      logger.error(f"Bulk operation failed: ${operationResponse.getException()}")
    } else if (itemResponse == null || !itemResponse.isSuccessStatusCode()) {
      val objectNode = itemOperation.getItem().asInstanceOf[ObjectNode]
      logger.error(String.format(
        "The operation for Item ID: [%s]  Item PartitionKey Value: [%s] did not complete " +
          "successfully with a %s/%s response code.",
        objectNode.get("id"),
        objectNode.get("sample_id"),
        if (itemResponse != null) itemResponse.getStatusCode() else "n/a",
        if (itemResponse != null) itemResponse.getSubStatusCode() else "n/a"))
    }

    if (itemResponse == null) {
      Mono.error(new IllegalStateException("No response retrieved."));
    } else {
      Mono.just(itemResponse);
    }

  }).blockLast()
}


def configureLogging(): Unit = {
  import ch.qos.logback.classic.Level
  val root = LoggerFactory.getLogger(org.slf4j.Logger.ROOT_LOGGER_NAME).asInstanceOf[ch.qos.logback.classic.Logger]
  root.setLevel(Level.INFO)
}


def objectNodeFromAvroRecord(objectMapper: ObjectMapper, jsonString: String, id: AtomicLong): (ObjectNode, Long) = {
  val jsonNode = objectMapper.readTree(jsonString)
  val objectNode = jsonNode.asInstanceOf[ObjectNode]
  val idLong = id.addAndGet(1L).longValue
  objectNode.put("id", objectMapper.convertValue("" + idLong, classOf[JsonNode]))
  (objectNode, idLong)
}
