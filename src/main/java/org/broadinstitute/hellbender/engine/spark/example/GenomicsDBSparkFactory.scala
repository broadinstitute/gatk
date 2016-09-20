/**
  * The MIT License (MIT)
  * Copyright (c) 2016 Intel Corporation
  *
  * Permission is hereby granted, free of charge, to any person obtaining a copy of
  * this software and associated documentation files (the "Software"), to deal in
  * the Software without restriction, including without limitation the rights to
  * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
  * the Software, and to permit persons to whom the Software is furnished to do so,
  * subject to the following conditions:
  *
  * The above copyright notice and this permission notice shall be included in all
  * copies or substantial portions of the Software.
  *
  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
  * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
  * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
  * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
  */

package org.broadinstitute.hellbender.engine.spark

import htsjdk.tribble.readers.PositionalBufferedStream
import htsjdk.variant.variantcontext.VariantContext
import org.apache.hadoop.io.LongWritable
import org.apache.spark.rdd.RDD
import org.apache.spark.{SparkConf, SparkContext}

object GenomicsDBSparkFactory {

  def usingGenomicsRDD(args: Array[String]): Unit = {

    val master = args(0)
    val port = args(1)
    val loaderJsonFile = args(2)
    val queryJsonFile = args(3)
    val hostfile = args(4)

    val conf = new SparkConf()
      .setMaster("spark://" + master + ":" + port)
      .setAppName("GenomicsDBTest using GenomicsDBRDD")
    val sc = new SparkContext(conf)
    val hadoopConf = sc.hadoopConfiguration
    hadoopConf.set(GenomicsDBConf.LOADERJSON, loaderJsonFile)
    hadoopConf.set(GenomicsDBConf.QUERYJSON, queryJsonFile)
    hadoopConf.set(GenomicsDBConf.MPIHOSTFILE, hostfile)

    val gc = new GenomicsDBContext(hadoopConf, sc)
    val myrdd = gc.getVariantContexts
    System.out.println(myrdd.count())
    myrdd.collect().foreach(println)
  }

  def usingNewAPIHadoopRDD(args: Array[String]): Unit = {
    val master = args(0)
    val port = args(1)
    val loaderJsonFile = args(2)
    val queryJsonFile = args(3)
    val hostfile = args(4)
    val conf = new SparkConf()
      .setMaster("spark://" + master + ":" + port)
      .setAppName("GenomicsDBTest using newAPIHadoopRDD")

    val classObjects: Array[Class[_]] = Array(Class.forName("org.apache.hadoop.io.LongWritable"),
      Class.forName("org.apache.hadoop.io.Text"))
    conf.registerKryoClasses(classObjects)

    val sc = new SparkContext(conf)
    val hadoopConf = sc.hadoopConfiguration
    hadoopConf.set(GenomicsDBConf.LOADERJSON, loaderJsonFile)
    hadoopConf.set(GenomicsDBConf.QUERYJSON, queryJsonFile)
    hadoopConf.set(GenomicsDBConf.MPIHOSTFILE, hostfile)

    val myrdd: RDD[(LongWritable, VariantContext)] =
      sc.newAPIHadoopRDD[LongWritable, VariantContext,
      GenomicsDBInputFormat[VariantContext, PositionalBufferedStream]](sc.hadoopConfiguration,
      classOf[GenomicsDBInputFormat[VariantContext, PositionalBufferedStream]],
      classOf[LongWritable], classOf[VariantContext])

    System.out.println(myrdd.count())
  }

  def main(args: Array[String]): Unit = {
//    usingNewAPIHadoopRDD(args)
    usingGenomicsRDD(args)
  }
}
