import $ivy.`org.apache.avro:avro:1.11.1`
import $ivy.`com.azure:azure-cosmos:4.41.0`
import $ivy.`com.lihaoyi:ammonite-ops_2.13:2.4.1`
import $ivy.`org.xerial.snappy:snappy-java:1.1.8.4`


import ammonite.ops._
import com.azure.cosmos._
import java.io.File
import java.util.NoSuchElementException
import org.apache.avro.file._
import org.apache.avro.generic._
import scala.jdk.CollectionConverters._


@main
def main(database: String, container: String, vets_dir: String, refs_dir: String): Unit = {
  val (endpoint, key) = extractEndpointAndKey()
  // val client = buildClient(endpoint, key)
  // println(f"Client is $client")

  val (vet_paths, ref_paths) = determineVetAndRefPaths(vets_dir, refs_dir)
  print(f"vets: $vet_paths\n")
  print(f"refs: $ref_paths\n")

  processFile(vet_paths.head)
}


def extractEndpointAndKey(): (String, String) = {
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


def determineVetAndRefPaths(vets_dir: String, refs_dir: String): (Seq[Path], Seq[Path]) = {
  val vets_path = pwd / RelPath(vets_dir)
  val refs_path = pwd / RelPath(refs_dir)
  val vet_listing = ls ! vets_path
  val refs_listing = ls ! refs_path

  (vet_listing.filter(_.isFile), refs_listing.filter(_.isFile))
}


// https://stackoverflow.com/a/45648136/21269164
def processFile(path: Path): Unit = {
  val file = new File(path.toString())
  val reader = new GenericDatumReader()
  val dataFileReader = new DataFileReader(file, reader)

  while (dataFileReader.hasNext) {
    val record: AnyRef = dataFileReader.next()
    System.out.println(record)
  }
}
