package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.api.client.util.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Iterator;
import java.util.List;

public class SubsettingRealignmentEngineTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return SubsettingRealignmentEngine.class.getSimpleName();
    }

    public static final String TEST_DATA_DIR = "src/test/resources/large/SubsettingRealignmentEngine/";
    public static final String BWA_IMAGE_PATH = TEST_DATA_DIR + "hg38_chr22_test.fasta.img";
    public static final String BAM_FILE = TEST_DATA_DIR + "test.realigned.bam";

    private static SubsettingRealignmentEngine getDefaultEngine(final SAMFileHeader header) {
        return new SubsettingRealignmentEngine(BWA_IMAGE_PATH, header, 10, 1, false);
    }

    private static boolean alwaysTrue(final GATKRead read) {
        return true;
    }

    @Test
    public void testSmall() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        try (final SubsettingRealignmentEngine engine = getDefaultEngine(inputHeader)) {
            final Iterator<SAMRecord> iter = reader.iterator();
            engine.addRead(new SAMRecordToGATKReadAdapter(iter.next()), SubsettingRealignmentEngineTest::alwaysTrue);
            engine.addRead(new SAMRecordToGATKReadAdapter(iter.next()), r -> false);
            final List<GATKRead> output = Lists.newArrayList(engine.alignAndMerge());

            Assert.assertEquals(engine.getSelectedReadsCount(), 1);
            Assert.assertEquals(engine.getNonselectedReadsCount(), 1);
            Assert.assertEquals(engine.getPairedAlignmentReadsCount(), 0);
            Assert.assertEquals(engine.getUnpairedAlignmentReadsCount(), 1);

            final long realignedOutput = output.stream().filter(SubsettingRealignmentEngine::readIsRealigned).count();
            Assert.assertEquals(realignedOutput, 1);

            final long notRealignedOutput = output.stream().filter(r -> !SubsettingRealignmentEngine.readIsRealigned(r)).count();
            Assert.assertEquals(notRealignedOutput, 1);

            Assert.assertFalse(output.get(0).isUnmapped());
            Assert.assertTrue(SubsettingRealignmentEngine.readIsRealigned(output.get(0)));
            Assert.assertFalse(output.get(1).isUnmapped());
            Assert.assertFalse(SubsettingRealignmentEngine.readIsRealigned(output.get(1)));
        }
    }

    @Test
    public void testFull() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        try (final SubsettingRealignmentEngine engine = getDefaultEngine(inputHeader)) {
            final Iterator<SAMRecord> iter = reader.iterator();
            while (iter.hasNext()) {
                final GATKRead read = new SAMRecordToGATKReadAdapter(iter.next());
                engine.addRead(read, r -> r.getMappingQuality() < 30);
            }
            final List<GATKRead> output = Lists.newArrayList(engine.alignAndMerge());

            Assert.assertEquals(engine.getSelectedReadsCount(), 103);
            Assert.assertEquals(engine.getNonselectedReadsCount(), 7713);
            Assert.assertEquals(engine.getPairedAlignmentReadsCount(), 34);
            Assert.assertEquals(engine.getUnpairedAlignmentReadsCount(), 69);

            final long realignedOutput = output.stream().filter(SubsettingRealignmentEngine::readIsRealigned).count();
            Assert.assertEquals(realignedOutput, 105);

            final long notRealignedOutput = output.stream().filter(r -> !SubsettingRealignmentEngine.readIsRealigned(r)).count();
            Assert.assertEquals(notRealignedOutput, 7713);

            final long properlyPairedOutput = output.stream().filter(GATKRead::isProperlyPaired).count();
            Assert.assertEquals(properlyPairedOutput, 7562);

            final long hasReadGroupOutput = output.stream().map(GATKRead::getReadGroup).filter(rg -> rg != null).count();
            Assert.assertEquals(hasReadGroupOutput, 7818);
        }
    }
    @Test
    public void testRetainDuplicateFlag() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        try (final SubsettingRealignmentEngine keepDuplicateEngine = new SubsettingRealignmentEngine(BWA_IMAGE_PATH, inputHeader, 10, 1, true);
             final SubsettingRealignmentEngine notKeepDuplicateEngine = new SubsettingRealignmentEngine(BWA_IMAGE_PATH, inputHeader, 10, 1, false)) {
            final Iterator<SAMRecord> iter = reader.iterator();
            final GATKRead duplicate = new SAMRecordToGATKReadAdapter(iter.next());
            duplicate.setIsDuplicate(true);
            final GATKRead nonDuplicate = new SAMRecordToGATKReadAdapter(iter.next());
            nonDuplicate.setIsDuplicate(false);

            keepDuplicateEngine.addRead(duplicate, SubsettingRealignmentEngineTest::alwaysTrue);
            keepDuplicateEngine.addRead(nonDuplicate, SubsettingRealignmentEngineTest::alwaysTrue);
            final List<GATKRead> keepDuplicateOutput = Lists.newArrayList(keepDuplicateEngine.alignAndMerge());
            Assert.assertEquals(keepDuplicateOutput.size(), 2);
            Assert.assertEquals(keepDuplicateOutput.stream().filter(GATKRead::isDuplicate).count(), 1);

            notKeepDuplicateEngine.addRead(duplicate, SubsettingRealignmentEngineTest::alwaysTrue);
            notKeepDuplicateEngine.addRead(nonDuplicate, SubsettingRealignmentEngineTest::alwaysTrue);
            final List<GATKRead> notKeepDuplicateOutput = Lists.newArrayList(notKeepDuplicateEngine.alignAndMerge());
            Assert.assertEquals(notKeepDuplicateOutput.size(), 2);
            Assert.assertEquals(notKeepDuplicateOutput.stream().filter(GATKRead::isDuplicate).count(), 0);
        }
    }

    @Test
    public void testEmptyInput() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        try (final SubsettingRealignmentEngine engine = getDefaultEngine(inputHeader)) {
            final List<GATKRead> output = Lists.newArrayList(engine.alignAndMerge());
            Assert.assertEquals(output.size(), 0);
            Assert.assertEquals(engine.getSelectedReadsCount(), 0);
            Assert.assertEquals(engine.getNonselectedReadsCount(), 0);
            Assert.assertEquals(engine.getPairedAlignmentReadsCount(), 0);
            Assert.assertEquals(engine.getUnpairedAlignmentReadsCount(), 0);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddingReadAfterAlignment() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        try (final SubsettingRealignmentEngine engine = getDefaultEngine(inputHeader)) {
            engine.alignAndMerge();
            engine.addRead(new SAMRecordToGATKReadAdapter(reader.iterator().next()), SubsettingRealignmentEngineTest::alwaysTrue);
        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testWrongSortType() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        // Must be coordinate sorted
        inputHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        getDefaultEngine(inputHeader);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testClosed1() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        final SubsettingRealignmentEngine engine = getDefaultEngine(inputHeader);
        engine.close();
        engine.addRead(new SAMRecordToGATKReadAdapter(reader.iterator().next()), SubsettingRealignmentEngineTest::alwaysTrue);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testClosed2() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        final SubsettingRealignmentEngine engine = getDefaultEngine(inputHeader);
        engine.addRead(new SAMRecordToGATKReadAdapter(reader.iterator().next()), SubsettingRealignmentEngineTest::alwaysTrue);
        engine.close();
        engine.alignAndMerge();
    }

}