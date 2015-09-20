package org.broadinstitute.hellbender.tools;

import com.google.common.io.CharStreams;
import org.apache.spark.launcher.SparkLauncher;

import java.io.File;
import java.io.InputStreamReader;

public class MySparkLauncher {
    // You must set SPARK_HOME in the environment.
    public static void main(String[] args) throws Exception {
        String path = new File("build/classes/main").getCanonicalPath();
        String jarLibPath = new File("build/lib/*").getCanonicalPath();
        System.out.println("starting... " + path);

        String inputBam = new File("src/test/resources/org/broadinstitute/hellbender/tools/BQSR/CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam").getCanonicalPath();
        String inputVcf = new File("src/test/resources/org/broadinstitute/hellbender/tools/BQSR/dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf").getCanonicalPath();
        String outputBam = "file:///" + new File("output.bam").getCanonicalPath();
        Process spark = new SparkLauncher()
                .setAppResource("")
                .setMainClass("org.broadinstitute.hellbender.Main")
                .setMaster("spark://davidada-macbookpro2.roam.corp.google.com:7077")
                .setConf(SparkLauncher.EXECUTOR_EXTRA_CLASSPATH, path + ":" + jarLibPath)
                .setConf(SparkLauncher.DRIVER_EXTRA_CLASSPATH, path + ":" + jarLibPath)
                .addAppArgs("ReadsPipelineSpark",
                        "-I", inputBam,
                        "-O", outputBam,
                        "-R", "gg://reference/EOSsjdnTicvzwAE",
                        "-BQSRKnownVariants", inputVcf,
                        "--apiKey", "AIzaSyDjKyy2KWcEVgtlOlJBdqy4O4oUVMxluE4",
                        "--sparkMaster", "spark://davidada-macbookpro2.roam.corp.google.com:7077")
                .launch();
        String stdout = CharStreams.toString(new InputStreamReader(spark.getInputStream()));// Drain the streams.
        String stderr = CharStreams.toString(new InputStreamReader(spark.getErrorStream()));// Drain the streams.
        int i = spark.waitFor();
        System.out.println(i);
        System.out.println(stdout);
        System.out.println(stderr);

    }
}
