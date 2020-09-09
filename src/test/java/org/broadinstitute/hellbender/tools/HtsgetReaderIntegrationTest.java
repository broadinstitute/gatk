package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Map;

import com.google.common.collect.ImmutableMap;


public class HtsgetReaderIntegrationTest extends CommandLineProgramTest {

    private static final String ENDPOINT = "https://htsget.ga4gh.org/reads";
    private static final String FILE_ID = "A1-B000168-3_57_F-1-1_R2.mus.Aligned.out.sorted.bam";
    private static final String LARGE_FILE_ID = "10X_P8_15_possorted_genome.bam";
    private static final String LARGE_FILE_MD5_FILE = "md5_checksum";

    @DataProvider(name = "successfulParameters")
    Object[][] successfulParameters() {
        return new Object[][] {
            { // basic download test
                ImmutableMap.of(
                    HtsgetReader.ID_LONG_NAME, FILE_ID),
                "A1-B000168-3_57_F-1-1_R2.mus.Aligned.out.sorted.bam"
            },
            { // parallel download test
                ImmutableMap.of(
                    HtsgetReader.ID_LONG_NAME, FILE_ID,
                    StandardArgumentDefinitions.INTERVALS_LONG_NAME, "chr1",
                    HtsgetReader.NUM_THREADS_LONG_NAME, "2"),
                "A1-B000168-3_57_F-1-1_R2.mus.Aligned.out.sorted.bam.refname"
            },
            { // header only
                ImmutableMap.of(
                    HtsgetReader.ID_LONG_NAME, FILE_ID,
                    HtsgetReader.CLASS_LONG_NAME, "header"),
                    "A1-B000168-3_57_F-1-1_R2.mus.Aligned.out.sorted.bam.header"
            },
            { // reference name without range
                ImmutableMap.of(
                    HtsgetReader.ID_LONG_NAME, FILE_ID,
                    StandardArgumentDefinitions.INTERVALS_LONG_NAME, "chr1"),
                "A1-B000168-3_57_F-1-1_R2.mus.Aligned.out.sorted.bam.refname"
            },
            { // reference name with range
                ImmutableMap.of(
                    HtsgetReader.ID_LONG_NAME, FILE_ID,
                    StandardArgumentDefinitions.INTERVALS_LONG_NAME, "chr1:24000000-25000000"),
                "A1-B000168-3_57_F-1-1_R2.mus.Aligned.out.sorted.bam.startend"
            },
//            TODO enable when server update goes in
//            { // filter by field
//                ImmutableMap.of(
//                    HtsgetReader.ID_LONG_NAME, FILE_ID,
//                    HtsgetReader.FIELDS_LONG_NAME, "QNAME"),
//                "A1-B000168-3_57_F-1-1_R2.mus.Aligned.out.sorted.bam.field"
//            }
        };
    }

    // Test successful combinations of query parameters
    // disabled until https://github.com/broadinstitute/gatk/pull/6799 is resolved
    @Test(enabled = false, dataProvider = "successfulParameters")
    public void testSuccessfulParameters(final Map<String, String> params, final String expectedFileName) throws IOException {
        final File expected = new File(getToolTestDataDir(), expectedFileName);
        final File output = createTempFile("output", ".bam");

        final ArgumentsBuilder args = new ArgumentsBuilder()
            .add(HtsgetReader.URL_LONG_NAME, ENDPOINT)
            .addOutput(output);
        params.forEach(args::add);
        
        runCommandLine(args);
        SamAssertionUtils.assertEqualBamFiles(output, expected, false, ValidationStringency.LENIENT);
    }

    //This test is disabled because it takes a long time.  It's included in order to run it manually.
    @Test(enabled = false)
    public void testLargeFileParallelDownload() throws IOException {
        final String expectedMd5 = new String(Files.readAllBytes(new File(getToolTestDataDir(), LARGE_FILE_MD5_FILE).toPath()));
        final File output = createTempFile("output", ".bam");

        final ArgumentsBuilder args = new ArgumentsBuilder()
            .add(HtsgetReader.URL_LONG_NAME, ENDPOINT)
            .add(HtsgetReader.ID_LONG_NAME, LARGE_FILE_ID)
            .add(HtsgetReader.NUM_THREADS_LONG_NAME, 4)
            .addOutput(output);

        runCommandLine(args);

        final String md5String = Utils.calculateFileMD5(output);
        Assert.assertEquals(md5String, expectedMd5);
    }

    @DataProvider(name = "failureParameters")
    Object[][] failureParameters() {
        return new Object[][] {
            { // Invalid sample id
                ImmutableMap.of(
                    HtsgetReader.URL_LONG_NAME, ENDPOINT,
                    HtsgetReader.ID_LONG_NAME, "nonexistentSample"
                )
            },
            { // Invalid endpoint
                ImmutableMap.of(
                    HtsgetReader.URL_LONG_NAME, "invalidEndpoint",
                    HtsgetReader.ID_LONG_NAME, FILE_ID
                )
            },
            { // Nonexistent contig
                ImmutableMap.of(
                    HtsgetReader.URL_LONG_NAME, ENDPOINT,
                    HtsgetReader.ID_LONG_NAME, FILE_ID,
                    StandardArgumentDefinitions.INTERVALS_LONG_NAME, "nonexistentInterval"
                )
            }
        };
    }

    // Test parameters that result in invalid request
    @Test(dataProvider = "failureParameters", expectedExceptions = UserException.class)
    public void testFailureParameters(final Map<String, String> params) {
        final File output = createTempFile("output", ".bam");

        final ArgumentsBuilder args = new ArgumentsBuilder()
            .addOutput(output);
        params.forEach(args::add);

        runCommandLine(args);
    }
}