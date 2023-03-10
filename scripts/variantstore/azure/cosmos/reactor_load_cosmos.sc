import $ivy.`com.azure:azure-cosmos:4.41.0`
import $ivy.`com.lihaoyi:ammonite-ops_2.13:2.4.1`
import $ivy.`org.apache.avro:avro:1.11.1`
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
  // There can occasionally be 409 race conditions on id if using a regular Long so go atomic.
  val id = AtomicLong(0L)

  // Where the Cosmos JSON serialization magic happens:
  // https://github.com/Azure/azure-sdk-for-java/blob/80b12e48aeb6ad2f49e86643dfd7223bde7a9a0c/sdk/cosmos/azure-cosmos/src/main/java/com/azure/cosmos/implementation/JsonSerializable.java#L255
  // Making a Jackson `ObjectNode` from `GenericRecord#toString` JSON seems like the least bad option
  // (`JsonSerializable` looks like a Jackson thing but is actually an internal Cosmos thing).
  val objectMapper = new ObjectMapper()

  val numDocsToLoad = 100000
  val numDocsProgressIncrement = 10000

  // On a Standard E4-2ads v5 Azure VM this Avro processing easily saturates the default 400 / 4000 RU/s maximum
  // throughput that Cosmos DB containers have when created in the Azure Portal. For 100K items, all the items print out
  // in the progress printlns, but then the script sits there waiting for the documents to be written to Cosmos.
  // This behavior was something of a surprise given the "reactor" structure; I would have expected backpressure to
  // limit the number of documents that were created if they couldn't actually be written to Cosmos. Is Cosmos silently
  // buffering documents up the RU/s limit?
  //
  // To see RU consumption:
  //
  // Azure Portal -> CosmosDB -> Metrics
  // Metric Namespace = 'Cosmos DB standard metrics'
  // Metric = 'Normalized RU Consumption'
  //
  // To see the number of documents written:
  //
  // CONTAINER_NAME=vets
  // az cosmosdb sql container show --database-name cosmos_gvs --name ${CONTAINER_NAME} --account-name ${COSMOS_DB_NAME}  | jq -r '..|.documentCount? //empty'
  val itemOperations = Flux.fromIterable(avro_paths.asJava).
    flatMap(path => {
      val file = new File(path.toString())
      val reader = new GenericDatumReader()
      val dataFileReader = new DataFileReader(file, reader)
      Flux.fromIterable(dataFileReader).
        // I had to put this `map` here inside the `flatMap` and not directly ahead of the `take` otherwise nothing
        // happened; I don't know why.
        map(record => {
          val jsonNode = objectMapper.readTree(record.toString())
          val objectNode = jsonNode.asInstanceOf[ObjectNode]
          val idLong = id.addAndGet(1L).longValue
          objectNode.put("id", objectMapper.convertValue("" + idLong, classOf[JsonNode]))
          if (idLong % numDocsProgressIncrement == 0L) println(f"${idLong}...");
          objectNode
        })
    }).
    take(numDocsToLoad).
    map(r => CosmosBulkOperations.getCreateItemOperation(r, new PartitionKey(r.get("sample_id").asLong)))

  executeItemOperationsWithErrorHandling(container, itemOperations)

  // itemOperations.subscribe(e => e)
}


def executeItemOperationsWithErrorHandling(container: CosmosAsyncContainer, itemOperations: Flux[CosmosItemOperation]): Unit = {
  // Only the first and last lines are the "execute" bits, the rest is error handling if something goes wrong.
  container.executeBulkOperations(itemOperations).flatMap(operationResponse => {
    val itemResponse = operationResponse.getResponse()
    val itemOperation = operationResponse.getOperation()

    if (operationResponse.getException() != null) {
      println(f"Bulk operation failed: ${operationResponse.getException()}")
    } else if (itemResponse == null || !itemResponse.isSuccessStatusCode()) {
      val objectNode = itemOperation.getItem().asInstanceOf[ObjectNode]
      println(String.format(
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
