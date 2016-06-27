package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link SegmentMergeUtils}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentMergeUtilsUnitTest extends BaseTest {
    //segments for testing mergeSegments
    private static final SimpleInterval segment1 = new SimpleInterval("chr1", 1, 4);
    private static final SimpleInterval segment2 = new SimpleInterval("chr1", 5, 12);
    private static final SimpleInterval segment3 = new SimpleInterval("chr1", 11, 20);
    private static final SimpleInterval segment4 = new SimpleInterval("chr2", 1, 10);
    private static final SimpleInterval segment5 = new SimpleInterval("chr2", 2, 9);

    @Test
    public void testMergeSegments() {
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment1, segment2), new SimpleInterval("chr1", 1, 12));
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment2, segment1), new SimpleInterval("chr1", 1, 12));
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment2, segment3), new SimpleInterval("chr1", 5, 20));
        Assert.assertEquals(SegmentMergeUtils.mergeSegments(segment5, segment4), new SimpleInterval("chr2", 1, 10));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testJoinWithIllegalArgumentException() {
        SegmentMergeUtils.mergeSegments(segment1, segment4);
    }

    @DataProvider(name = "dataNeutralSegmentMerging")
    public Object[][] dataNeutralSegmentMerging() {
        return new Object[][]{
                //=================================================================================================
                {
                        //segments to merge
                        Arrays.asList(
                                //no calls, merge all
                                new ModeledSegment(new SimpleInterval("1", 1, 10), 1, 0),
                                new ModeledSegment(new SimpleInterval("1", 11, 20), 1, 0),
                                new ModeledSegment(new SimpleInterval("1", 21, 30), 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 1, 10), 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 11, 20), 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 21, 30), 1, 0)
                        ),
                        //expected merged segments
                        Arrays.asList(
                                new SimpleInterval("1", 1, 30),
                                new SimpleInterval("2", 1, 30)
                        )
                },
                //=================================================================================================
                {
                        //segments to merge
                        Arrays.asList(
                                new ModeledSegment(new SimpleInterval("1", 1, 10), ReCapSegCaller.AMPLIFICATION_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("1", 11, 20), ReCapSegCaller.NEUTRAL_CALL, 1, 0),  //merge right
                                new ModeledSegment(new SimpleInterval("1", 21, 30), ReCapSegCaller.NEUTRAL_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("1", 31, 40), ReCapSegCaller.DELETION_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("1", 41, 50), ReCapSegCaller.AMPLIFICATION_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("1", 51, 60), ReCapSegCaller.NEUTRAL_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 1, 10), ReCapSegCaller.AMPLIFICATION_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 11, 20), ReCapSegCaller.NEUTRAL_CALL, 1, 0),  //merge right twice
                                new ModeledSegment(new SimpleInterval("2", 21, 30), ReCapSegCaller.NEUTRAL_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 31, 40), ReCapSegCaller.NEUTRAL_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 41, 50), ReCapSegCaller.AMPLIFICATION_CALL, 1, 0),
                                new ModeledSegment(new SimpleInterval("2", 51, 60), ReCapSegCaller.NEUTRAL_CALL, 1, 0)
                        ),
                        //expected merged segments
                        Arrays.asList(
                                new SimpleInterval("1", 1, 10),
                                new SimpleInterval("1", 11, 30),
                                new SimpleInterval("1", 31, 40),
                                new SimpleInterval("1", 41, 50),
                                new SimpleInterval("1", 51, 60),
                                new SimpleInterval("2", 1, 10),
                                new SimpleInterval("2", 11, 40),
                                new SimpleInterval("2", 41, 50),
                                new SimpleInterval("2", 51, 60)
                        )
                }
        };
    }

    @Test(dataProvider = "dataNeutralSegmentMerging")
    public void testNeutralSegmentMerging(final List<ModeledSegment> segments,
                                          final List<SimpleInterval> expectedMergedSegments) {
        final List<SimpleInterval> resultMergedSegments =
                SegmentMergeUtils.mergeNeutralSegments(segments);
        Assert.assertEquals(resultMergedSegments, expectedMergedSegments);
    }

    /**
     * Tests for small-segment merging, along with methods for constructing test data.
     */
    public static final class SmallSegmentsTestHelper {
        //segment is small if number of targets it contains is strictly less than SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD
        private static final int SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD = 3;

        private static final SimpleInterval leftSegment = new SimpleInterval("chr2", 100, 200);
        private static final SimpleInterval leftSegmentCloser = new SimpleInterval("chr2", 100, 250);
        private static final SimpleInterval leftSegmentDifferentContig = new SimpleInterval("chr1", 100, 200);

        private static final SimpleInterval centerSegment = new SimpleInterval("chr2", 300, 400);

        private static final SimpleInterval rightSegment = new SimpleInterval("chr2", 500, 600);
        private static final SimpleInterval rightSegmentCloser = new SimpleInterval("chr2", 450, 600);
        private static final SimpleInterval rightSegmentDifferentContig = new SimpleInterval("chr3", 500, 600);

        private static List<SimpleInterval> makeSegments(final SimpleInterval leftSegment,
                                                         final SimpleInterval rightSegment) {
            return Arrays.asList(leftSegment, centerSegment, rightSegment);
        }

        private static Genome makeGenome(final List<List<ReadCountRecord.SingleSampleRecord>> targetCoverages,
                                         final List<List<AllelicCount>> snpCounts) {
            final List<ReadCountRecord.SingleSampleRecord> records = new ArrayList<>();
            final List<AllelicCount> snps = new ArrayList<>();
            targetCoverages.stream().forEach(records::addAll);
            snpCounts.stream().forEach(snps::addAll);

            final List<Target> targets = records.stream().map(ReadCountRecord::getTarget).collect(Collectors.toList());
            final RealMatrix counts = new Array2DRowRealMatrix(targets.size(), 1);
            counts.setColumn(0, records.stream().mapToDouble(r -> r.getCount()).toArray());
            return new Genome(new ReadCountCollection(targets, Arrays.asList("sample"), counts), snps);
        }

        private static final List<SimpleInterval> mergedSegmentsLeft =
                Arrays.asList(new SimpleInterval("chr2", 100, 400), rightSegment);
        private static final List<SimpleInterval> mergedSegmentsRight =
                Arrays.asList(leftSegment, new SimpleInterval("chr2", 300, 600));
        private static final List<SimpleInterval> mergedSegmentsAll =
                Arrays.asList(new SimpleInterval("chr2", 100, 600));

        @DataProvider(name = "dataSmallSegmentMerging")
        public Object[][] dataSmallSegmentMerging() {
            return new Object[][]{
                    //=================================================================================================
                    {
                            //case description
                            "missing targets and SNPs, different genomic distance = merge all on genomic distance",
                            //original segments
                            makeSegments(leftSegmentCloser, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 250)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            //no targets
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            //no SNPs
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 10, 10)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 10, 10)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsAll
                    },
                    //=================================================================================================
                    {
                            //case description
                            "missing targets and SNPs, different genomic distance = merge all on genomic distance",
                            //original segments
                            makeSegments(leftSegment, rightSegmentCloser),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            //no targets
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 450, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 10, 10)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 10, 10)
                                    ),
                                    Arrays.asList(
                                            //no SNPs
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsAll
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, missing SNPs, different genomic distance = merge left on genomic distance",
                            //original segments
                            makeSegments(leftSegmentCloser, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 250)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            //no SNPs
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 10, 10)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 10, 10)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsLeft
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, no SNPs, different genomic distance = merge right on genomic distance",
                            //original segments
                            makeSegments(leftSegment, rightSegmentCloser),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 450, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    //no SNPs
                            ),
                            //expected merged segments
                            mergedSegmentsRight
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, no SNPs, left on different contig = merge right on genomic distance",
                            //original segments
                            makeSegments(leftSegmentDifferentContig, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr1", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr1", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr1", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    //no SNPs
                            ),
                            //expected merged segments
                            Arrays.asList(leftSegmentDifferentContig, new SimpleInterval("chr2", 300, 600))
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, no SNPs, right on different contig = merge left on genomic distance",
                            //original segments
                            makeSegments(leftSegment, rightSegmentDifferentContig),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr3", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr3", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr3", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    //no SNPs
                            ),
                            //expected merged segments
                            Arrays.asList(new SimpleInterval("chr2", 100, 400), rightSegmentDifferentContig)
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, no SNPs, left and right on different contigs = no merge",
                            //original segments
                            makeSegments(leftSegmentDifferentContig, rightSegmentDifferentContig),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr1", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr1", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr1", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr3", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr3", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr3", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    //no SNPs
                            ),
                            //expected merged segments
                            Arrays.asList(leftSegmentDifferentContig, rightSegmentDifferentContig)
                    },
                    //=================================================================================================
                    {
                            //case description
                            "center coverage lowest, no SNPs, same genomic distance = merge left on coverage",
                            //original segments
                            makeSegments(leftSegment, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 2.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 2.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 200)), 2.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 3.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 3.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 3.0)
                                    )
                            ),
                            //SNP counts
                            Collections.emptyList(), //no SNPS
                            //expected merged segments
                            mergedSegmentsLeft
                    },
                    //=================================================================================================
                    {
                            //case description
                            "center coverage highest, no SNPs, same genomic distance = merge right on coverage",
                            //original segments
                            makeSegments(leftSegment, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 3.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 2.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 2.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 2.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    //no SNPs
                            ),
                            //expected merged segments
                            mergedSegmentsRight
                    },
                    //=================================================================================================
                    {
                            //case description
                            "center coverage highest, not enough SNPs, same genomic distance = merge right on coverage",
                            //original segments
                            makeSegments(leftSegment, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 3.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 2.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 2.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 2.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 10, 10)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 10, 10)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 10, 10)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsRight
                    },
                    //=================================================================================================
                    {
                            //case description
                            "right coverages closer, no SNPs, same genomic distance = merge right on coverage",
                            //original segments
                            makeSegments(leftSegment, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 0.7),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 0.8),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 129)), 0.9),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 130, 139)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 140, 200)), 1.1)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 310, 319)), 0.95),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 320, 329)), 1.05)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target8", new SimpleInterval("chr2", 500, 509)), 0.85),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target9", new SimpleInterval("chr2", 510, 519)), 0.95),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target10", new SimpleInterval("chr2", 520, 529)), 1.05),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target11", new SimpleInterval("chr2", 530, 539)), 1.15),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target12", new SimpleInterval("chr2", 540, 600)), 1.25)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    //no SNPs
                            ),
                            //expected merged segments
                            mergedSegmentsRight
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, non-overlapping AAFs, same genomic distance = merge right on inverse MAF",
                            //original segments
                            makeSegments(leftSegment, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 100 - 20, 20)

                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 100 - 30, 30),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 100 - 35, 35),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 100 - 40, 40)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 100 - 50, 50),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 100 - 55, 55),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 100 - 60, 60)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsRight
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages and AAFs, different genomic distance = merge left on genomic distance",
                            //original segments
                            makeSegments(leftSegmentCloser, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 250)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 103, 103), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 104, 104), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 105, 105), 100 - 20, 20)

                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 303, 303), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 304, 304), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 305, 305), 100 - 20, 20)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 503, 503), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 504, 504), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 505, 505), 100 - 20, 20)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsLeft
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, different number of same AAFs, different genomic distance = " +
                                    "merge left on genomic distance",
                            //original segments
                            makeSegments(leftSegmentCloser, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 250)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 100 - 11, 11),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 100 - 20, 20)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 303, 303), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 304, 304), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 305, 305), 100 - 20, 20)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 100 - 11, 11),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 503, 503), 100 - 11, 11),
                                            new AllelicCount(new SimpleInterval("chr2", 504, 504), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 505, 505), 100 - 20, 20)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsLeft
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, right closer in AAF but same KS distance, different genomic distance = " +
                                    "merge right on inverse MAF",
                            //original segments
                            makeSegments(leftSegmentCloser, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 250)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 100 - 7, 7),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 103, 103), 100 - 8, 8),
                                            new AllelicCount(new SimpleInterval("chr2", 104, 104), 100 - 13, 13),
                                            new AllelicCount(new SimpleInterval("chr2", 105, 105), 100 - 18, 18)

                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 303, 303), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 304, 304), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 305, 305), 100 - 20, 20)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 100 - 11, 11),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 100 - 21, 21),
                                            new AllelicCount(new SimpleInterval("chr2", 503, 503), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 504, 504), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 505, 505), 100 - 22, 22)
                                    )
                            ),
                            //expected merged segments
                            Arrays.asList(leftSegmentCloser, new SimpleInterval("chr2", 300, 600))
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, right closer in AAF variance, different genomic distance = " +
                                    "merge right on AAF",
                            //original segments
                            makeSegments(leftSegmentCloser, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 250)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            //15 alt counts drawn from normal with mean 15, variance 5, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 100 - 5, 5),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 100 - 6, 6),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 103, 103), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 104, 104), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 105, 105), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 106, 106), 100 - 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 107, 107), 100 - 11, 11),
                                            new AllelicCount(new SimpleInterval("chr2", 108, 108), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 109, 109), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 110, 110), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 111, 111), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 112, 112), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 113, 113), 100 - 18, 18),
                                            new AllelicCount(new SimpleInterval("chr2", 114, 114), 100 - 19, 19),
                                            //15 alt counts drawn from normal with mean 85, variance 5, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 115, 115), 100 - 77, 77),
                                            new AllelicCount(new SimpleInterval("chr2", 116, 116), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 117, 117), 100 - 81, 81),
                                            new AllelicCount(new SimpleInterval("chr2", 118, 118), 100 - 81, 81),
                                            new AllelicCount(new SimpleInterval("chr2", 119, 119), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 120, 120), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 121, 121), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 122, 122), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 123, 123), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 124, 124), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 125, 125), 100 - 87, 87),
                                            new AllelicCount(new SimpleInterval("chr2", 126, 126), 100 - 87, 87),
                                            new AllelicCount(new SimpleInterval("chr2", 127, 127), 100 - 91, 91),
                                            new AllelicCount(new SimpleInterval("chr2", 128, 128), 100 - 93, 93),
                                            new AllelicCount(new SimpleInterval("chr2", 129, 129), 100 - 93, 93)

                                    ),
                                    Arrays.asList(
                                            //10 alt counts drawn from normal with mean 15, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 303, 303), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 304, 304), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 305, 305), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 306, 306), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 307, 307), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 308, 308), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 309, 309), 100 - 18, 18),
                                            //10 alt counts drawn from normal with mean 85, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 310, 310), 100 - 82, 82),
                                            new AllelicCount(new SimpleInterval("chr2", 311, 311), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 312, 312), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 313, 313), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 314, 314), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 315, 315), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 316, 316), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 317, 317), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 318, 318), 100 - 87, 87),
                                            new AllelicCount(new SimpleInterval("chr2", 319, 319), 100 - 88, 88)
                                    ),
                                    Arrays.asList(
                                            //10 alt counts drawn from normal with mean 15, variance 3, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 100 - 11, 11),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 503, 503), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 504, 504), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 505, 505), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 506, 506), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 507, 507), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 508, 508), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 509, 509), 100 - 19, 19),
                                            //10 alt counts drawn from normal with mean 85, variance 3, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 510, 510), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 511, 511), 100 - 82, 82),
                                            new AllelicCount(new SimpleInterval("chr2", 512, 512), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 513, 513), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 514, 514), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 515, 515), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 516, 516), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 517, 517), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 518, 518), 100 - 86, 86),
                                            new AllelicCount(new SimpleInterval("chr2", 519, 519), 100 - 87, 87)
                                    )
                            ),
                            //expected merged segments
                            Arrays.asList(leftSegmentCloser, new SimpleInterval("chr2", 300, 600))
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, right closer in AAF mean, different genomic distance = merge right on AAF",
                            //original segments
                            makeSegments(leftSegmentCloser, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 250)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            //15 alt counts drawn from normal with mean 15, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 100 - 12, 12),
                                            new AllelicCount(new SimpleInterval("chr2", 103, 103), 100 - 14, 14),
                                            new AllelicCount(new SimpleInterval("chr2", 104, 104), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 105, 105), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 106, 106), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 107, 107), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 108, 108), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 109, 109), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 110, 110), 100 - 15, 15),
                                            new AllelicCount(new SimpleInterval("chr2", 111, 111), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 112, 112), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 113, 113), 100 - 16, 16),
                                            new AllelicCount(new SimpleInterval("chr2", 114, 114), 100 - 17, 17),
                                            //15 alt counts drawn from normal with mean 85, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 115, 115), 100 - 82, 82),
                                            new AllelicCount(new SimpleInterval("chr2", 116, 116), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 117, 117), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 118, 118), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 119, 119), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 120, 120), 100 - 83, 83),
                                            new AllelicCount(new SimpleInterval("chr2", 121, 121), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 122, 122), 100 - 84, 84),
                                            new AllelicCount(new SimpleInterval("chr2", 123, 123), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 124, 124), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 125, 125), 100 - 85, 85),
                                            new AllelicCount(new SimpleInterval("chr2", 126, 126), 100 - 86, 86),
                                            new AllelicCount(new SimpleInterval("chr2", 127, 127), 100 - 86, 86),
                                            new AllelicCount(new SimpleInterval("chr2", 128, 128), 100 - 87, 87),
                                            new AllelicCount(new SimpleInterval("chr2", 129, 129), 100 - 88, 88)

                                    ),
                                    Arrays.asList(
                                            //10 alt counts drawn from normal with mean 20, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 100 - 19, 19),
                                            new AllelicCount(new SimpleInterval("chr2", 303, 303), 100 - 19, 19),
                                            new AllelicCount(new SimpleInterval("chr2", 304, 304), 100 - 19, 19),
                                            new AllelicCount(new SimpleInterval("chr2", 305, 305), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 306, 306), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 307, 307), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 308, 308), 100 - 21, 21),
                                            new AllelicCount(new SimpleInterval("chr2", 309, 309), 100 - 24, 24),
                                            //10 alt counts drawn from normal with mean 80, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 310, 310), 100 - 74, 75),
                                            new AllelicCount(new SimpleInterval("chr2", 311, 311), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 312, 312), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 313, 313), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 314, 314), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 315, 315), 100 - 81, 81),
                                            new AllelicCount(new SimpleInterval("chr2", 316, 316), 100 - 81, 81),
                                            new AllelicCount(new SimpleInterval("chr2", 317, 317), 100 - 81, 81),
                                            new AllelicCount(new SimpleInterval("chr2", 318, 318), 100 - 82, 82),
                                            new AllelicCount(new SimpleInterval("chr2", 319, 319), 100 - 82, 82)
                                    ),
                                    Arrays.asList(
                                            //10 alt counts drawn from normal with mean 20, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 100 - 17, 17),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 100 - 18, 18),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 100 - 19, 19),
                                            new AllelicCount(new SimpleInterval("chr2", 503, 503), 100 - 19, 19),
                                            new AllelicCount(new SimpleInterval("chr2", 504, 504), 100 - 19, 19),
                                            new AllelicCount(new SimpleInterval("chr2", 505, 505), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 506, 506), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 507, 507), 100 - 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 508, 508), 100 - 21, 21),
                                            new AllelicCount(new SimpleInterval("chr2", 509, 509), 100 - 23, 23),
                                            //10 alt counts drawn from normal with mean 80, variance 2, rounded
                                            new AllelicCount(new SimpleInterval("chr2", 510, 510), 100 - 78, 78),
                                            new AllelicCount(new SimpleInterval("chr2", 511, 511), 100 - 78, 78),
                                            new AllelicCount(new SimpleInterval("chr2", 512, 512), 100 - 79, 79),
                                            new AllelicCount(new SimpleInterval("chr2", 513, 513), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 514, 514), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 515, 515), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 516, 516), 100 - 80, 80),
                                            new AllelicCount(new SimpleInterval("chr2", 517, 517), 100 - 81, 81),
                                            new AllelicCount(new SimpleInterval("chr2", 518, 518), 100 - 81, 81),
                                            new AllelicCount(new SimpleInterval("chr2", 519, 519), 100 - 82, 82)
                                    )
                            ),
                            //expected merged segments
                            Arrays.asList(leftSegmentCloser, new SimpleInterval("chr2", 300, 600))
                    },
                    //=================================================================================================
                    //tests with random results below; if any of these tests are removed or their order changed,
                    //or the random seed is changed from RANDOM_SEED = 42, recheck randomly generated expected segments
                    //=================================================================================================
                    {
                            //case description
                            "missing coverages, missing SNPs, same genomic distance = merge left",
                            //original segments
                            makeSegments(leftSegment, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 200)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 500, 509)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 510, 519)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 520, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            //no SNPs
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 10, 10)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 10, 10)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsLeft
                    },
                    //=================================================================================================
                    {
                            //case description
                            "same coverages, AAFs, and genomic distance = merge right",
                            //original segments
                            makeSegments(leftSegment, rightSegment),
                            //target coverages
                            Arrays.asList(
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target1", new SimpleInterval("chr2", 100, 109)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target2", new SimpleInterval("chr2", 110, 119)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target3", new SimpleInterval("chr2", 120, 299)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target4", new SimpleInterval("chr2", 300, 400)), 1.0)
                                    ),
                                    Arrays.asList(
                                            new ReadCountRecord.SingleSampleRecord(new Target("target5", new SimpleInterval("chr2", 401, 409)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target6", new SimpleInterval("chr2", 410, 419)), 1.0),
                                            new ReadCountRecord.SingleSampleRecord(new Target("target7", new SimpleInterval("chr2", 420, 600)), 1.0)
                                    )
                            ),
                            //SNP counts
                            Arrays.asList(
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 100, 100), 20, 20),
                                            new AllelicCount(new SimpleInterval("chr2", 101, 101), 30, 30),
                                            new AllelicCount(new SimpleInterval("chr2", 102, 102), 40, 40)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 300, 300), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 301, 301), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 302, 302), 10, 10)
                                    ),
                                    Arrays.asList(
                                            new AllelicCount(new SimpleInterval("chr2", 500, 500), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 501, 501), 10, 10),
                                            new AllelicCount(new SimpleInterval("chr2", 502, 502), 10, 10)
                                    )
                            ),
                            //expected merged segments
                            mergedSegmentsRight
                    }
            };
        }

        @Test(dataProvider = "dataSmallSegmentMerging")
        public void testSmallSegmentMerging(final String caseDescription,
                                            final List<SimpleInterval> segments,
                                            final List<List<ReadCountRecord.SingleSampleRecord>> targetCoverages,
                                            final List<List<AllelicCount>> snpCounts,
                                            final List<SimpleInterval> expectedMergedSegments) {
            logger.info("Testing case: " + caseDescription);
            final Genome genome = makeGenome(targetCoverages, snpCounts);
            final List<SimpleInterval> resultMergedSegments =
                    SegmentMergeUtils.mergeSmallSegments(segments, genome, SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD);
            Assert.assertEquals(resultMergedSegments, expectedMergedSegments);
        }
    }

    /**
     * Tests for similar-segment merging.
     */
    public static final class SimilarSegmentsTestHelper {
        private static final double INTERVAL_THRESHOLD = 1.;    //merge when parameter modes are within 95% HPD interval

        private static ACNVModeledSegment constructSegment(final String contig, final int start, final int end,
                                                          final double segmentMeanPosteriorCenter,
                                                          final double segmentMeanPosteriorHalfWidth,
                                                          final double minorAlleleFractionPosteriorCenter,
                                                          final double minorAlleleFractionPosteriorHalfWidth) {
            return new ACNVModeledSegment(new SimpleInterval(contig, start, end),
                    new PosteriorSummary(
                            segmentMeanPosteriorCenter,
                            segmentMeanPosteriorCenter - segmentMeanPosteriorHalfWidth,
                            segmentMeanPosteriorCenter + segmentMeanPosteriorHalfWidth),
                    new PosteriorSummary(
                            minorAlleleFractionPosteriorCenter,
                            minorAlleleFractionPosteriorCenter - minorAlleleFractionPosteriorHalfWidth,
                            minorAlleleFractionPosteriorCenter + minorAlleleFractionPosteriorHalfWidth));
        }

        @DataProvider(name = "dataSimilarSegmentMerging")
        public Object[][] dataSimilarSegmentMerging() {
            return new Object[][]{
                    //=================================================================================================
                    {
                            //segments to merge
                            Arrays.asList(
                                    constructSegment("1", 1, 10, 1., 0.05, 0.45, 0.1),      //merge right
                                    constructSegment("1", 11, 20, 1.05, 0.05, 0.25, 0.2),
                                    constructSegment("1", 21, 30, 5., 0.05, 0.25, 0.1),     //no merge
                                    constructSegment("2", 1, 10, 5., 0.05, 0.25, 0.05),     //no merge
                                    constructSegment("2", 11, 20, 1.05, 0.01, 0.25, 0.05),  //merge right
                                    constructSegment("2", 21, 30, 1., 0.05, 0.3, 0.01)
                            ),
                            //expected merged segments
                            Arrays.asList(
                                    constructSegment("1", 1, 20, 1.025, 0.035, 0.41, 0.089),
                                    constructSegment("1", 21, 30, 5., 0.05, 0.25, 0.1),
                                    constructSegment("2", 1, 10, 5., 0.05, 0.25, 0.05),
                                    constructSegment("2", 11, 30, 1.048, 0.0098, 0.298, 0.0098)
                            )
                    },
                    //=================================================================================================
                    {
                            //segments to merge
                            Arrays.asList(
                                    constructSegment("1", 1, 10, 1., 0.05, 0.45, 0.1),      //merge right twice
                                    constructSegment("1", 11, 20, 1.05, 0.05, 0.25, 0.2),
                                    constructSegment("1", 21, 30, 1.05, 0.01, 0.25, 0.05),
                                    constructSegment("1", 31, 40, 5., 0.05, 0.25, 0.1)      //no merge
                            ),
                            //expected merged segments
                            Arrays.asList(
                                    constructSegment("1", 1, 30, 1.048, 0.0096, 0.288, 0.044),
                                    constructSegment("1", 31, 40, 5., 0.05, 0.25, 0.1)
                            )
                    }
            };
        }

        @Test(dataProvider = "dataSimilarSegmentMerging")
        public void testSimilarSegmentMerging(final List<ACNVModeledSegment> segments,
                                              final List<ACNVModeledSegment> expectedMergedSegments) {
            final List<ACNVModeledSegment> resultMergedSegments =
                    SegmentMergeUtils.mergeSimilarSegments(segments, INTERVAL_THRESHOLD, INTERVAL_THRESHOLD);
            Assert.assertTrue(resultMergedSegments.size() == expectedMergedSegments.size());
            for (int i = 0; i < resultMergedSegments.size(); i++) {
                final ACNVModeledSegment resultSegment = resultMergedSegments.get(i);
                final ACNVModeledSegment expectedSegment = expectedMergedSegments.get(i);
                Assert.assertEquals(resultSegment.getInterval(), expectedSegment.getInterval());
                Assert.assertEquals(
                        resultSegment.getSegmentMeanPosteriorSummary().getCenter(),
                        expectedSegment.getSegmentMeanPosteriorSummary().getCenter(), 0.001);
                Assert.assertEquals(
                        resultSegment.getSegmentMeanPosteriorSummary().getLower(),
                        expectedSegment.getSegmentMeanPosteriorSummary().getLower(), 0.001);
                Assert.assertEquals(
                        resultSegment.getSegmentMeanPosteriorSummary().getUpper(),
                        expectedSegment.getSegmentMeanPosteriorSummary().getUpper(), 0.001);
                Assert.assertEquals(
                        resultSegment.getMinorAlleleFractionPosteriorSummary().getCenter(),
                        expectedSegment.getMinorAlleleFractionPosteriorSummary().getCenter(), 0.001);
                Assert.assertEquals(
                        resultSegment.getMinorAlleleFractionPosteriorSummary().getLower(),
                        expectedSegment.getMinorAlleleFractionPosteriorSummary().getLower(), 0.001);
                Assert.assertEquals(
                        resultSegment.getMinorAlleleFractionPosteriorSummary().getUpper(),
                        expectedSegment.getMinorAlleleFractionPosteriorSummary().getUpper(), 0.001);
            }
        }
    }
}