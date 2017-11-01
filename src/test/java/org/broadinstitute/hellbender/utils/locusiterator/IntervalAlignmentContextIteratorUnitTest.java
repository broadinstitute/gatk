package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;


public class IntervalAlignmentContextIteratorUnitTest extends GATKBaseTest {

    private static final File TEST_LARGE_DATA_DIR = new File("src/test/resources/large/");
    private static final String BAM_FILE_NAME = TEST_LARGE_DATA_DIR.getAbsolutePath() + "/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";
    private static final File TEST_DATA_DIR = new File("src/test/resources/");
    private static final File MINI_BAM = new File(TEST_DATA_DIR.getAbsolutePath() + "/NA12878.chr17_69k_70k.dictFix.bam");

    @Test
    public void testSimple() {

        // Note: IGV was used to determine coverage in the input bam file. (chr20:9,999,875-9,999,945)
        //   Make sure that viewing is *not* downsampled and *includes* dupe reads

        // Totally uncovered
        final SimpleInterval record_20_9999600_9999600 = new SimpleInterval("20:9999600-9999600");

        // Partially covered
        final SimpleInterval record_20_9999900_9999910 = new SimpleInterval("20:9999900-9999910");

        // Totally covered
        final SimpleInterval record_20_9999910_9999913 = new SimpleInterval("20:9999910-9999913");


        final List<SimpleInterval> locusIntervals = new ArrayList<>(3);
        locusIntervals.add(record_20_9999600_9999600);
        locusIntervals.add(record_20_9999900_9999910);
        locusIntervals.add(record_20_9999910_9999913);

        final List<AlignmentContext> allAlignmentContexts = getAlignmentContexts(locusIntervals, BAM_FILE_NAME);

        Assert.assertEquals(allAlignmentContexts.size(), 16);
        Assert.assertTrue(allAlignmentContexts.stream().allMatch(ac -> ac != null));

        // No reads in the first three alignment contexts
        Assert.assertEquals(allAlignmentContexts.get(0).getBasePileup().getReads().size(), 0);
        Assert.assertEquals(allAlignmentContexts.get(1).getBasePileup().getReads().size(),  0);
        Assert.assertEquals(allAlignmentContexts.get(2).getBasePileup().getReads().size(),  0);

        // Make sure that at least one locus-interval has at least one read
        Assert.assertTrue(allAlignmentContexts.stream().anyMatch(ac -> ac.getBasePileup().getReads().size() > 0));
        Assert.assertTrue(allAlignmentContexts.stream().anyMatch(ac -> ac.getBasePileup().getReads().size() == 0));

        Assert.assertEquals(allAlignmentContexts.get(3).getBasePileup().getReads().size(), 2);
        Assert.assertEquals(allAlignmentContexts.get(4).getBasePileup().getReads().size(), 4);
        Assert.assertEquals(allAlignmentContexts.get(13).getBasePileup().getReads().size(), 17);
    }

    @Test
    public void testCoveredOnly() {

        // This test is good at finding places where the alignment contexts start behind the interval.
        // Totally covered
        final SimpleInterval record_20_9999910_9999913 = new SimpleInterval("20:9999910-9999913");
        final List<SimpleInterval> locusIntervals = new ArrayList<>(1);
        locusIntervals.add(record_20_9999910_9999913);

        final List<AlignmentContext> allAlignmentContexts = getAlignmentContexts(locusIntervals, BAM_FILE_NAME);

        Assert.assertEquals(allAlignmentContexts.size(), 4);
        Assert.assertTrue(allAlignmentContexts.stream().allMatch(ac -> ac != null));
    }

    @Test
    public void testUnCoveredOnly() {

        final SimpleInterval record_20_9999900_9999901 = new SimpleInterval("20:9999900-9999901");
        final List<SimpleInterval> locusIntervals = new ArrayList<>(1);
        locusIntervals.add(record_20_9999900_9999901);

        final List<AlignmentContext> allAlignmentContexts = getAlignmentContexts(locusIntervals, BAM_FILE_NAME);

        Assert.assertEquals(allAlignmentContexts.size(), 2);
        Assert.assertTrue(allAlignmentContexts.stream().allMatch(ac -> ac != null));
        Assert.assertEquals(allAlignmentContexts.get(0).getBasePileup().getReads().size(), 0);
        Assert.assertEquals(allAlignmentContexts.get(1).getBasePileup().getReads().size(),  0);
    }

    @Test
    public void testPileupToNextPileup(){
        // Test reading the end of one pileup then no coverage then beginning of one pileup in one interval
        //  IGV: chr20:10,098,446-10,098,587
        final SimpleInterval record_20_10098490_10098560 = new SimpleInterval("20:10098490-10098560");
        final List<SimpleInterval> locusIntervals = new ArrayList<>(1);
        locusIntervals.add(record_20_10098490_10098560);
        final List<AlignmentContext> allAlignmentContexts = getAlignmentContexts(locusIntervals, BAM_FILE_NAME);

        Assert.assertEquals(allAlignmentContexts.size(), 71);
        Assert.assertTrue(allAlignmentContexts.stream().allMatch(ac -> ac.getBasePileup().getReads().size() < 2));
        Assert.assertTrue(IntStream.range(0, 6).allMatch(i -> allAlignmentContexts.get(i).size() == 1));
        Assert.assertTrue(IntStream.range(6, 62).allMatch(i -> allAlignmentContexts.get(i).size() == 0));
        Assert.assertTrue(IntStream.range(62, 71).allMatch(i -> allAlignmentContexts.get(i).size() == 1));
        Assert.assertTrue(allAlignmentContexts.stream().allMatch(ac -> ac != null));
    }

    @Test
    public void testNoIntervalsSpecified() {
        // $ samtools view -H NA12878.chr17_69k_70k.dictFix.bam
        // @HD	VN:1.5	GO:none	SO:coordinate
        // @SQ	SN:17	LN:1000000	AS:GRCh37	UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	M5:351f64d4f4f9ddd45b35336ad97aa6de	SP:Homo Sapiens

        final List<AlignmentContext> allAlignmentContexts = getAlignmentContexts(null, MINI_BAM.getAbsolutePath());

        Assert.assertEquals(allAlignmentContexts.size(), 1000000);
    }

    @Test
    public void testIntervalIteratorIsEmpty() {
        // This test should not throw an exception.
        final List<AlignmentContext> allAlignmentContexts = getAlignmentContexts(Collections.<SimpleInterval>emptyList(),MINI_BAM.getAbsolutePath());
        Assert.assertEquals(allAlignmentContexts.size(), 0);
    }

    private List<AlignmentContext> getAlignmentContexts(final List<SimpleInterval> locusIntervals, final String bamPath) {
        final List<String> sampleNames = Collections.singletonList("NA12878");

        final ReadsDataSource gatkReads = new ReadsDataSource(IOUtils.getPath(bamPath));
        final SAMFileHeader header = gatkReads.getHeader();
        final Stream<GATKRead> filteredReads = Utils.stream(gatkReads).filter(new WellformedReadFilter(header).and(new ReadFilterLibrary.MappedReadFilter()));

        final SAMSequenceDictionary dictionary = header.getSequenceDictionary();

        final LocusIteratorByState locusIteratorByState = new LocusIteratorByState(filteredReads.iterator(), LocusIteratorByState.NO_DOWNSAMPLING, false, sampleNames, header, true);


        List<SimpleInterval> relevantIntervals = locusIntervals;
        if (relevantIntervals == null) {
            relevantIntervals = IntervalUtils.getAllIntervalsForReference(dictionary);
        }
        final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(relevantIntervals.iterator());

        final IntervalAlignmentContextIterator intervalAlignmentContextIterator = new IntervalAlignmentContextIterator(locusIteratorByState, intervalLocusIterator, dictionary);

        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(intervalAlignmentContextIterator, Spliterator.ORDERED), false).collect(Collectors.toList());
    }
}
