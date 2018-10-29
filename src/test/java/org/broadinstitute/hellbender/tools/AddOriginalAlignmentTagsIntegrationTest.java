package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static org.testng.Assert.*;

public final class AddOriginalAlignmentTagsIntegrationTest extends CommandLineProgramTest {
    private File localTestData = new File(getTestDataDir(), "mitochondria/NA12878.bam");
    private static final ArrayList<String> expectedOATags = new ArrayList<>(Arrays.asList(
            "chrM,1,+,92S59M,60,0;",
            "chrM,1,+,19S132M,60,0;",
            "chrM,1,-,47S104M,60,0;",
            "*,0,*,*,0,0;",
            "chrM,2,-,116S35M,60,0;",
            "chrM,3,-,92S59M,60,0;",
            "*,0,*,*,0,0;",
            "*,0,*,*,0,0;"));

    private static final ArrayList<String> expectedXMTags = new ArrayList<>(Arrays.asList(
            "chrM",
            "chrM",
            "*",
            "chrM",
            "chrM",
            "chrM",
            "*",
            "*"));

    @Test()
    public void testAddingOATag() {
        final File oaAddedBam = createTempFile("oaAdded", ".bam");

        final List<String> args = Arrays.asList("-I", localTestData.getPath(),
                "-O", oaAddedBam.getPath());
        runCommandLine(args);

        SamReader out = SamReaderFactory.makeDefault().open(oaAddedBam);
        Iterator<String> expectedOAIter = expectedOATags.iterator();
        Iterator<String> expectedXMIter = expectedXMTags.iterator();
        int readsInBam = 0;

        for(SAMRecord r : out) {
            assertEquals(r.getAttribute(AddOriginalAlignmentTags.OA_TAG_NAME), expectedOAIter.next());
            assertEquals(r.getAttribute(AddOriginalAlignmentTags.MATE_CONTIG_TAG_NAME), expectedXMIter.next());
            readsInBam++;
        }
        assertEquals(readsInBam, 8);
    }

    @Test()
    public void testUsingIntervalList() {
        final File oaAddedBam = createTempFile("oaAdded", ".bam");

        final List<String> args = Arrays.asList("-I", localTestData.getPath(),
                "-O", oaAddedBam.getPath(),
                "-L", "chrM");
        runCommandLine(args);

        SamReader out = SamReaderFactory.makeDefault().open(oaAddedBam);
        Iterator<String> expectedOAIter = expectedOATags.iterator();
        Iterator<String> expectedXMIter = expectedXMTags.iterator();
        int readsInBam = 0;

        for(SAMRecord r : out) {
            assertEquals(r.getAttribute(AddOriginalAlignmentTags.OA_TAG_NAME), expectedOAIter.next());
            assertEquals(r.getAttribute(AddOriginalAlignmentTags.MATE_CONTIG_TAG_NAME), expectedXMIter.next());
            readsInBam++;
        }
        assertEquals(readsInBam, 6);
    }


}
