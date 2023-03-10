import $ivy.`org.apache.avro:avro:1.11.1`
import $ivy.`com.lihaoyi:ammonite-ops_2.13:2.4.1`
import $ivy.`org.xerial.snappy:snappy-java:1.1.8.4`
import $ivy.`com.google.code.gson:gson:2.10.1`


import ammonite.ops._
import com.google.gson._
import java.io.File
import org.apache.avro._
import org.apache.avro.file._
import org.apache.avro.generic._
import org.apache.avro.specific._


class Vet() {
  var id: String = null
  var sample_id: Long = 0
  var location: Long = 0
  var ref: String = null
  var alt: String = null
  var AS_RAW_MQ: String = null
  var AS_RAW_MQRankSum: String = null
  var QUALapprox: String = null
  var AS_QUALapprox: String = null
  var AS_RAW_ReadPosRankSum: String = null
  var AS_SB_TABLE: String = null
  var AS_VarDP: String = null
  var call_GT: String = null
  var call_AD: String = null
  var call_GQ: String = null
  var call_PGT: String = null
  var call_PID: String = null
  var call_PL: String = null

  def getId(): String = this.id

  def setId(id: String): Unit = {
    this.id = id
  }
}

@main
def main(avro: String): Unit = {
  val path = pwd / avro
  val file = new File(path.toString())

  val genericReader = new GenericDatumReader()
  val genericFileReader = new DataFileReader(file, genericReader)

  val genericRecord: GenericRecord = genericFileReader.next()
  val gson = new Gson()
  val vet = gson.fromJson(genericRecord.toString(), classOf[Vet])
  println(vet.call_AD)
}
