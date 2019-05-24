package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.codecs.AnnotatedIntervalCodec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SimpleAnnotatedIntervalWriterUnitTest extends GATKBaseTest {

    private static final File TEST_FILE_NO_SAMHEADER = new File(toolsTestDir,
            "copynumber/utils/combine-segment-breakpoints-no-samheader.seg");
    private static final File TEST_FILE = new File(toolsTestDir,
            "copynumber/utils/combine-segment-breakpoints.seg");
    private static final File TEST_FILE_OLD_HEADER = new File(toolsTestDir,
            "copynumber/utils/annotatedinterval/simple-annotated-interval-writer-replacement-header-comments.seg");
    private static final File TEST_FILE_OLD_HEADER_CONFIG = new File(toolsTestDir,
            "copynumber/utils/annotatedinterval/old-header.config");

    @Test
    public void testNoSamFileHeaderWithComments() throws IOException {
        final File outputFile = File.createTempFile("simpleannotatedintervalwriter_no_samheader", ".seg");
        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(TEST_FILE_NO_SAMHEADER.toPath(), null);

        final SimpleAnnotatedIntervalWriter writer = new SimpleAnnotatedIntervalWriter(outputFile);
        writer.writeHeader(AnnotatedIntervalCodec.createHeaderForWriter(collection.getAnnotations(), collection.getSamFileHeader()));
        collection.getRecords().forEach(r -> writer.add(r));
        writer.close();

        final AnnotatedIntervalCollection testCollection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        // Reminder: In this case, the output will have additional comments for the header
        Assert.assertEquals(testCollection.getComments().subList(0,2), collection.getComments());
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.CONTIG_COL_COMMENT + "CONTIG")));
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.START_COL_COMMENT + "START")));
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.END_COL_COMMENT + "END")));
        Assert.assertEquals(testCollection.getRecords(), collection.getRecords());
        Assert.assertEquals(testCollection.getSamFileHeader(), collection.getSamFileHeader());
        Assert.assertEquals(testCollection.getAnnotations(), collection.getAnnotations());
    }

    @Test
    public void testRoundTrip() throws IOException {
        final File outputFile = File.createTempFile("simpleannotatedintervalwriter_", ".seg");
        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(TEST_FILE.toPath(), null);

        final SimpleAnnotatedIntervalWriter writer = new SimpleAnnotatedIntervalWriter(outputFile);

        writer.writeHeader(AnnotatedIntervalCodec.createHeaderForWriter(collection.getAnnotations(), collection.getSamFileHeader()));
        collection.getRecords().forEach(r -> writer.add(r));
        writer.close();

        final AnnotatedIntervalCollection testCollection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        // Reminder: In this case, the output will have additional comments for the header.  The input had none.
        Assert.assertEquals(testCollection.getComments().size(), 3);
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.CONTIG_COL_COMMENT + "CONTIG")));
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.START_COL_COMMENT + "START")));
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.END_COL_COMMENT + "END")));
        Assert.assertEquals(testCollection.getRecords(), collection.getRecords());
        Assert.assertEquals(testCollection.getSamFileHeader(), collection.getSamFileHeader());
        Assert.assertEquals(testCollection.getAnnotations(), collection.getAnnotations());
    }

    @Test
    public void testRoundTripWithReplacementHeader() throws IOException {
        // Tests that we can read a file with a non-standard header (specified by a config file) and that we can output
        //  with the default headers (input headers would then not match output headers).
        //  Also, tests that the structured comments get replaced with the new value.
        final File outputFile = File.createTempFile("simpleannotatedintervalwriter_", ".seg");
        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(TEST_FILE_OLD_HEADER.toPath(),
                TEST_FILE_OLD_HEADER_CONFIG.toPath(),null);

        final String newContigColumn = "CONTIG";
        final String newStartColumn = "START";
        final String newEndColumn = "END";

        try (final SimpleAnnotatedIntervalWriter writer = new SimpleAnnotatedIntervalWriter(outputFile)) {
            writer.writeHeader(AnnotatedIntervalCodec.createHeaderForWriter(collection.getAnnotations(), collection.getSamFileHeader()));
            collection.getRecords().forEach(r -> writer.add(r));
        }

        final AnnotatedIntervalCollection testCollection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        // Reminder: In this case, the output will have additional comments for the header.  The input had none.
        Assert.assertEquals(testCollection.getComments().size(), 1+3);
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.CONTIG_COL_COMMENT + newContigColumn)));
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.START_COL_COMMENT + newStartColumn)));
        Assert.assertTrue(testCollection.getComments().stream().anyMatch(c -> c.equals(AnnotatedIntervalCodec.END_COL_COMMENT + newEndColumn)));
        Assert.assertEquals(testCollection.getRecords(), collection.getRecords());
        Assert.assertEquals(testCollection.getSamFileHeader(), collection.getSamFileHeader());
        Assert.assertEquals(testCollection.getAnnotations(), collection.getAnnotations());
    }
}
