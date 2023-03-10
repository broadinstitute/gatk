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
  var id: Long = 1L
  // Where the Cosmos JSON serialization magic happens:
  // https://github.com/Azure/azure-sdk-for-java/blob/80b12e48aeb6ad2f49e86643dfd7223bde7a9a0c/sdk/cosmos/azure-cosmos/src/main/java/com/azure/cosmos/implementation/JsonSerializable.java#L255
  // Making a Jackson `ObjectNode` from `GenericRecord#toString` JSON seems like the least bad option
  // (`JsonSerializable` looks like a Jackson thing but is actually an internal Cosmos thing).
  val objectMapper = new ObjectMapper()

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
          objectNode.put("id", objectMapper.convertValue("" + id, classOf[JsonNode]))
          id = id + 1
          objectNode
        })
    }).
    take(100000).
    // map(r => { if (id % 10000 == 0) println(f"${id}..."); r }).
    map(r => CosmosBulkOperations.getCreateItemOperation(r, new PartitionKey(r.get("sample_id").asLong)))

  executeItemOperationsWithErrorHandling(container, itemOperations)

  // itemOperations.subscribe(e => System.out.println(e))
}


def executeItemOperationsWithErrorHandling(container: CosmosAsyncContainer, itemOperations: Flux[CosmosItemOperation]): Unit = {
  container.executeBulkOperations(itemOperations).flatMap(operationResponse => {
    val itemResponse = operationResponse.getResponse()
    val itemOperation = operationResponse.getOperation()

    if (operationResponse.getException() != null) {
      println(f"Bulk operation failed: ${operationResponse.getException()}")
    } else if (itemResponse == null || !itemResponse.isSuccessStatusCode()) {
      val objectNode = itemOperation.getItem().asInstanceOf[ObjectNode]
      println(String.format(
        "The operation for Item ID: [%s]  Item PartitionKey Value: [%s] did not complete " +
          "successfully with " + "a" + " %s/%s response code, response: %s.",
        objectNode.get("id"),
        objectNode.get("sample_id"),
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
