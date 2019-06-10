package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredIntervalArgumentCollection;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
//TODO: This is not meant for prime time yet.  DO NOT MERGE INTO MASTER
public class MapTestReadWalkerIntegrationTest extends CommandLineProgramTest {
    @Test
    public void basicTest() throws IOException {

        final File outputFile = new File("/home/lichtens/broad_oncotator_configs/test_ct_test/test.txt"); //File.createTempFile("test_simple", ".txt");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        arguments.add("gs://fc-908984e1-32e4-4745-910f-1319506618ce/RCCBMS-00150-N.bam");
        arguments.add("--" + MapTestReadWalker.BWA_INDEX_PARAM_FULL_NAME);
        arguments.add("/home/lichtens/broad_oncotator_configs/hg19_bwa_idx/Homo_sapiens_assembly19.fasta.img");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add("/home/lichtens/broad_oncotator_configs/hg19_bwa_idx/Homo_sapiens_assembly19.fasta");
        arguments.add("-" + StandardArgumentDefinitions.INTERVALS_SHORT_NAME);
        arguments.add("/home/lichtens/broad_oncotator_configs/hg19_bwa_idx/ice_targets.tsv.preprocessed.interval_list");
//        arguments.add("-" + StandardArgumentDefinitions.INTERVALS_SHORT_NAME);
//        arguments.add("1:3000000-4500000");
        arguments.add("--" + RequiredIntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME);
        arguments.add("INTERSECTION");
        arguments.add("--" + StandardArgumentDefinitions.VERBOSITY_NAME);
        arguments.add(Log.LogLevel.INFO.name());
        runCommandLine(arguments);
    }
}
