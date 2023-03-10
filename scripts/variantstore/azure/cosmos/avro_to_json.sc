import $ivy.`org.apache.avro:avro:1.11.1`
import $ivy.`com.lihaoyi:ammonite-ops_2.13:2.4.1`
import $ivy.`org.xerial.snappy:snappy-java:1.1.8.4`


import ammonite.ops._
import java.io.File
import org.apache.avro.file._
import org.apache.avro.generic._
import com.fasterxml.jackson.databind.{JsonNode, ObjectMapper}
import com.fasterxml.jackson.databind.node.ObjectNode

@main
def main(avro: String): Unit = {
  val path = pwd / avro
  val file = new File(path.toString())

  val genericReader = new GenericDatumReader()
  val genericFileReader = new DataFileReader(file, genericReader)

  val genericRecord: GenericRecord = genericFileReader.next()
  // Where the Cosmos JSON serialization magic happens:
  // https://github.com/Azure/azure-sdk-for-java/blob/80b12e48aeb6ad2f49e86643dfd7223bde7a9a0c/sdk/cosmos/azure-cosmos/src/main/java/com/azure/cosmos/implementation/JsonSerializable.java#L255
  // Making a Jackson `JsonNode` seems like the least bad option (`JsonSerializable` is an internal Cosmos thing).
  val objectMapper = new ObjectMapper()
  val jsonNode = objectMapper.readTree(genericRecord.toString())
  val objectNode = jsonNode.asInstanceOf[ObjectNode]
  // println(objectNode.getClass)
  val id: Long = 1L
  objectNode.put("id", objectMapper.convertValue("" + id, classOf[JsonNode]))
  val sample_id: Long = objectNode.get("sample_id").asLong
  println(sample_id)
  println(objectNode)
}
