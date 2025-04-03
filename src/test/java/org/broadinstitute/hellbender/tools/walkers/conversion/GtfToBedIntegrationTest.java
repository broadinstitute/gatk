package org.broadinstitute.hellbender.tools.walkers.conversion;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class GtfToBedIntegrationTest extends CommandLineProgramTest {
    private static final String mapk1 = ConversionTestUtils.getMapk1Gtf();
    private static final String decoys = ConversionTestUtils.getDecoySampleGtf();
    private static final String manyTranscripts = ConversionTestUtils.getManyTranscriptsGtf();
    private static final String chr14Gene = ConversionTestUtils.getChr14GeneGtf();
    private static final String mouse = ConversionTestUtils.getMouseGtf();
    private static final String mouseDictionary = ConversionTestUtils.getMouseDict();
    private static final String dictionary = ConversionTestUtils.getReferenceDict();
    private static final String chr14Fasta = ConversionTestUtils.getChr14Fasta();
    private static final String expectedDecoysGeneBed = ConversionTestUtils.getDecoySamplesGeneBed();
    private static final String expectedDecoysTranscriptBed = ConversionTestUtils.getDecoySamplesTranscriptBed();
    private static final String expectedMapk1GeneBed = ConversionTestUtils.getMapk1GeneBed();
    private static final String expectedMapk1TranscriptBed = ConversionTestUtils.getMapk1TranscriptBed();
    private static final String expectedManyTranscriptsBed = ConversionTestUtils.getManyTranscriptsBed();
    private static final String expectedNotBasicBed = ConversionTestUtils.getNotBasicBed();
    private static final String expectedMouseBed = ConversionTestUtils.getMouseBed();
    private static final String expectedChr14GeneBed = ConversionTestUtils.getChr14GeneBed();

    private static class GtfToBedTest {
        final String input;
        final String SD;
        final String transcript;
        final String basic;
        final String expected;

        private GtfToBedTest(String input, String SD, String transcript, String basic, String expected) {
            this.input = input;
            this.SD = SD;
            this.transcript = transcript;
            this.basic = basic;
            this.expected = expected;
        }
    }

    @DataProvider(name = "GtfToBedTestProvider")
    public Object[][] equalRangeData() {
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{new GtfToBedTest(mapk1, dictionary, String.valueOf(false), String.valueOf(true), expectedMapk1GeneBed)});
        tests.add(new Object[]{new GtfToBedTest(mapk1, dictionary, String.valueOf(true), String.valueOf(true), expectedMapk1TranscriptBed)});
        tests.add(new Object[]{new GtfToBedTest(decoys, dictionary, String.valueOf(false), String.valueOf(true), expectedDecoysGeneBed)});
        tests.add(new Object[]{new GtfToBedTest(decoys, dictionary, String.valueOf(true), String.valueOf(true), expectedDecoysTranscriptBed)});
        tests.add(new Object[]{new GtfToBedTest(manyTranscripts, dictionary, String.valueOf(true), String.valueOf(true), expectedManyTranscriptsBed)});
        tests.add(new Object[]{new GtfToBedTest(manyTranscripts, dictionary, String.valueOf(true), String.valueOf(false), expectedNotBasicBed)});
        tests.add(new Object[]{new GtfToBedTest(mouse, mouseDictionary, String.valueOf(true), String.valueOf(false), expectedMouseBed)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GtfToBedTestProvider")
    public void testGtfToBed(GtfToBedTest params) throws IOException {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArrayList<String> args = new ArrayList<>();

        args.add("--" + GtfToBed.INPUT_LONG_NAME);
        args.add(new File (params.input).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME);
        args.add(new File (params.SD).getAbsolutePath());
        args.add("--" + GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME);
        args.add(params.transcript);
        args.add("--" + GtfToBed.USE_BASIC_TRANSCRIPT_LONG_NAME);
        args.add(params.basic);
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outputFile.getAbsolutePath());

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(new File(params.expected), outputFile);
    }

    @Test(expectedExceptions = UserException.class)
    public void testGtfToBedNoSequenceDictionary() {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArrayList<String> args = new ArrayList<>();

        args.add("--" + GtfToBed.INPUT_LONG_NAME);
        args.add(new File (mapk1).getAbsolutePath());
        args.add("--" + GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME);
        args.add(String.valueOf(true));
        args.add("--" + GtfToBed.USE_BASIC_TRANSCRIPT_LONG_NAME);
        args.add(String.valueOf(true));
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outputFile.getAbsolutePath());

        runCommandLine(args);
    }

    @Test
    public void testGtfToBedSequenceDictionaryFromReference() throws IOException {
        final File outputFile = createTempFile("outputBed", ".bed");
        final ArrayList<String> args = new ArrayList<>();

        args.add("--" + GtfToBed.INPUT_LONG_NAME);
        args.add(new File (chr14Gene).getAbsolutePath());
        args.add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        args.add(chr14Fasta);
        args.add("--" + GtfToBed.SORT_BY_TRANSCRIPT_LONG_NAME);
        args.add(String.valueOf(false));
        args.add("--" + GtfToBed.USE_BASIC_TRANSCRIPT_LONG_NAME);
        args.add(String.valueOf(false));
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outputFile.getAbsolutePath());

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(new File(expectedChr14GeneBed), outputFile);
    }
}
