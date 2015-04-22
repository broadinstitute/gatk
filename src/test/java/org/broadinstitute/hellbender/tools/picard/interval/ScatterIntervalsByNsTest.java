package org.broadinstitute.hellbender.tools.picard.interval;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Created by farjoun on 3/26/14.
 */
public final class ScatterIntervalsByNsTest extends CommandLineProgramTest {

    @DataProvider(name = "testSegregate")
    public Object[][] testSegregate() {

        return new Object[][]{
                new Object[]{"NNNNNAAAAANNNNNN", 1, CollectionUtil.makeList(
                        new Interval("fake1", 1, 5),
                        new Interval("fake1", 6, 10),
                        new Interval("fake1", 11, 16))},

                new Object[]{"NNNNNAAAAANNNNNN", 5, CollectionUtil.makeList(
                        new Interval("fake1", 1, 5),
                        new Interval("fake1", 6, 10),
                        new Interval("fake1", 11, 16))},

                new Object[]{"NNNNNAAAnANNNNNN", 6, CollectionUtil.makeList(
                        new Interval("fake1", 1, 5),
                        new Interval("fake1", 6, 10),
                        new Interval("fake1", 11, 16))},

                new Object[]{"AnTGCNNNNNACGTCG", 1, CollectionUtil.makeList(
                        new Interval("fake1", 1, 5),
                        new Interval("fake1", 6, 10),
                        new Interval("fake1", 11, 16))},

                new Object[]{"ACTGCNNNNNACGTCG", 4, CollectionUtil.makeList(
                        new Interval("fake1", 1, 5),
                        new Interval("fake1", 6, 10),
                        new Interval("fake1", 11, 16))},

                new Object[]{"ACTGCNNNNNACgctG", 5, CollectionUtil.makeList(
                        new Interval("fake1", 1, 16))},

                new Object[]{"ACNGCNNNNNACnTNG", 5, CollectionUtil.makeList(
                        new Interval("fake1", 1, 16))},

                new Object[]{"ACNGCNNNNNACnTNG", 0, CollectionUtil.makeList(
                        new Interval("fake1", 1, 2),   //acgt
                        new Interval("fake1", 3, 3),   //N
                        new Interval("fake1", 4, 5),
                        new Interval("fake1", 6, 10),  //N
                        new Interval("fake1", 11, 12),
                        new Interval("fake1", 13, 13), //N
                        new Interval("fake1", 14, 14),
                        new Interval("fake1", 15, 15), //N
                        new Interval("fake1", 16, 16))},

                new Object[]{"AAAAAAAAAAACnTNG", 0, CollectionUtil.makeList(
                        new Interval("fake1", 1, 12),
                        new Interval("fake1", 13, 13), //N
                        new Interval("fake1", 14, 14),
                        new Interval("fake1", 15, 15), //N
                        new Interval("fake1", 16, 16))},

                new Object[]{"AAAAAAAAAAACATNG", 0, CollectionUtil.makeList(
                        new Interval("fake1", 1, 14),
                        new Interval("fake1", 15, 15), //N
                        new Interval("fake1", 16, 16))},

                new Object[]{"ANCNNGCNNNNNACNTNGN", 1, CollectionUtil.makeList(
                        new Interval("fake1", 1, 3),
                        new Interval("fake1", 4, 5),
                        new Interval("fake1", 6, 7),
                        new Interval("fake1", 8, 12),
                        new Interval("fake1", 13, 18),
                        new Interval("fake1", 19, 19)

                )},
                new Object[]{"ANCNNGCNNNNNACNTNGN", 0, CollectionUtil.makeList(
                        new Interval("fake1", 1, 1),
                        new Interval("fake1", 2, 2),
                        new Interval("fake1", 3, 3),
                        new Interval("fake1", 4, 5),
                        new Interval("fake1", 6, 7),
                        new Interval("fake1", 8, 12),
                        new Interval("fake1", 13, 14),
                        new Interval("fake1", 15, 15),
                        new Interval("fake1", 16, 16),
                        new Interval("fake1", 17, 17),
                        new Interval("fake1", 18, 18),
                        new Interval("fake1", 19, 19)
                )}


        };
    }

    @Test(dataProvider = "testSegregate")
    public void testSegregateReference(final String referenceString, final int maxNmerToMerge, final List<Interval> result) throws Exception {

        final SAMSequenceRecord record = new SAMSequenceRecord("fake1", referenceString.length());
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(record);

        final ReferenceSequenceFile reference = new ReferenceSequenceFile() {

            boolean done = false;

            @Override
            public SAMSequenceDictionary getSequenceDictionary() {
                return dictionary;
            }

            @Override
            public ReferenceSequence nextSequence() {
                if (!done) {
                    done = true;
                    return getSequence(record.getSequenceName());
                }
                return null;
            }

            @Override
            public void reset() {
                done = false;
            }

            @Override
            public boolean isIndexed() {
                return false;
            }

            @Override
            public ReferenceSequence getSequence(final String contig) {
                if (contig.equals(record.getSequenceName())) {
                    return new ReferenceSequence(record.getSequenceName(), 0, referenceString.getBytes());
                } else {
                    return null;
                }
            }

            @Override
            public ReferenceSequence getSubsequenceAt(final String contig, final long start, final long stop) {
                return null;
            }

            @Override
            public String toString() {
                return null;
            }

            @Override
            public void close() {}
        };
        Assert.assertEquals(ScatterIntervalsByNs.segregateReference(reference, maxNmerToMerge).getIntervals(), result);
    }
}
