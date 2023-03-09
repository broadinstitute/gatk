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
   val id: String                      = null
   val sample_id: Long                 = 0
   val location: Long                  = 0
   val ref: String                     = null
   val alt: String                     = null
   val AS_RAW_MQ: String               = null
   val AS_RAW_MQRankSum: String        = null
   val QUALapprox: String              = null
   val AS_QUALapprox: String           = null
   val AS_RAW_ReadPosRankSum: String   = null
   val AS_SB_TABLE: String             = null
   val AS_VarDP: String                = null
   val call_GT: String                 = null
   val call_AD: String                 = null
   val call_GQ: String                 = null
   val call_PGT: String                = null
   val call_PID: String                = null
   val call_PL: String                 = null
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
