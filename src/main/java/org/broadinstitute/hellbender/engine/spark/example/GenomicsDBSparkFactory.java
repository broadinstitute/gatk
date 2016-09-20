/*
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

package org.broadinstitute.hellbender.engine.spark.example;

import com.intel.genomicsdb.GenomicsDBConf;
import com.intel.genomicsdb.GenomicsDBContext;
import com.intel.genomicsdb.GenomicsDBInputFormat;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.util.ReflectionUtils;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import scala.Tuple2;

import java.util.List;

/**
 * This factory class exposes the two ways JavaRDD of variant contexts (htsjdk)
 * can be retrieved from GenomicsDB. In case of the newAPIHadoopRDD(), GenomicsDB
 * returns a JavaPairRDD where the genomics positions are the key. However, this
 * is seldom used in the variant contexts as downstream applications in HellBender
 * code uses only the values and ignores the key
 */
public final class GenomicsDBSparkFactory {

  private static void usingGenomicsRDD(String[] args) {

    String master = args[0];
    int port = Integer.parseInt(args[1]);
    String loaderJsonFile = args[2];
    String queryJsonFile = args[3];
    String hostfile = args[4];

    SparkConf conf = new SparkConf();
    conf.setMaster("spark://" + master + ":" + port);
    conf.setAppName("GenomicsDBTest using GenomicsDBRDD");
    JavaSparkContext sc = new JavaSparkContext(conf);

    Configuration hadoopConf = sc.hadoopConfiguration();
    hadoopConf.set(GenomicsDBConf.LOADERJSON, loaderJsonFile);
    hadoopConf.set(GenomicsDBConf.QUERYJSON, queryJsonFile);
    hadoopConf.set(GenomicsDBConf.MPIHOSTFILE, hostfile);

    GenomicsDBContext gc = new GenomicsDBContext(hadoopConf, sc.sc());
    JavaRDD<VariantContext> myrdd = gc.getVariantContexts().toJavaRDD();
    System.out.println(myrdd.count());

    List<VariantContext> vlist = myrdd.collect();
    for (Object aVlist : vlist) {
      System.out.println(aVlist);
    }
  }

  public static void usingNewAPIHadoopRDD(String[] args) {
    String master = args[0];
    int port = Integer.parseInt(args[1]);
    String loaderJsonFile = args[2];
    String queryJsonFile = args[3];
    String hostfile = args[4];

    SparkConf conf = new SparkConf();
    conf.setMaster("spark://" + master + ":" + port);
    conf.setAppName("GenomicsDBTest using newAPIHadoopRDD");
    JavaSparkContext sc = new JavaSparkContext(conf);

    Configuration hadoopConf = sc.hadoopConfiguration();
    hadoopConf.set(GenomicsDBConf.LOADERJSON, loaderJsonFile);
    hadoopConf.set(GenomicsDBConf.QUERYJSON, queryJsonFile);
    hadoopConf.set(GenomicsDBConf.MPIHOSTFILE, hostfile);

    JavaPairRDD<String, VariantContext> pairedVariants;

    GenomicsDBInputFormat<VariantContext, PositionalBufferedStream> t = new
      GenomicsDBInputFormat<VariantContext, PositionalBufferedStream>();

    Class<GenomicsDBInputFormat<VariantContext, PositionalBufferedStream>> gformatClazz
      = ReflectionUtils.getClass(t);

    pairedVariants = sc.<String, VariantContext, GenomicsDBInputFormat<VariantContext, PositionalBufferedStream>>newAPIHadoopRDD(
      hadoopConf, gformatClazz, String.class, VariantContext.class);

    JavaRDD<VariantContext> variants = pairedVariants.map(GenomicsDBSparkFactory::tuple2variants);

    System.out.println(variants.count());
    List<VariantContext> variantList = variants.collect();
    for (Object variantObj : variantList) {
      System.out.println(variantObj);
    }
  }

  private static VariantContext tuple2variants(Tuple2<String, VariantContext> tuple) {
    return tuple._2();
  }

  public static void main(String[] args) throws Exception {
//    usingGenomicsRDD(args);
    usingNewAPIHadoopRDD(args);
  }
}
