package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.segment.SegmentExonOverlaps;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.segment.SegmentExonUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfGeneFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SegmentExonUtilsUnitTest extends GATKBaseTest {

    private static final ReferenceDataSource refDataSourceHg38;

    static {
        refDataSourceHg38 = ReferenceDataSource.of( IOUtils.getPath(hg38Reference) );
    }

        @DataProvider
    public Object[][] provideExonSegmentCalculation() {
        // 5'UTR 1000-1199
        // exon0 1200-1299
        // intron 1300-1499
        // exon1 1500-1599
        // intron 1600-1799
        // exon2 1800-1899
        // intron 1900-2099
        // exon3 2100-2199
        // 3'UTR 2200-2399
        final GencodeGtfGeneFeature geneFeaturePos = DataProviderForExampleGencodeGtfGene.dynamicallyCreateTestGencodeGtfGeneFeature("chr1",
                1000, "TESTING", Strand.POSITIVE, 4, 200, 200, 100, 200);
        final GencodeGtfGeneFeature geneFeatureNeg = DataProviderForExampleGencodeGtfGene.dynamicallyCreateTestGencodeGtfGeneFeature("chr1",
                1000, "TESTING", Strand.NEGATIVE, 4, 200, 200, 100, 200);

        return new Object[][]{

                //  Reminder UTRs still considered part of the exon for purposes of this test.
                // Segment start breakpoint is in the transcript (overlap)
                {new SimpleInterval("chr1", 1050, 5000), geneFeaturePos, "0+", ""},
                {new SimpleInterval("chr1", 1050, 5000), geneFeatureNeg, "3-", ""},
                {new SimpleInterval("chr1", 1550, 5000), geneFeaturePos, "1+", ""},
                {new SimpleInterval("chr1", 1550, 5000), geneFeatureNeg, "2-", ""},
                {new SimpleInterval("chr1", 1850, 5000), geneFeaturePos, "2+", ""},
                {new SimpleInterval("chr1", 1850, 5000), geneFeatureNeg, "1-", ""},
                {new SimpleInterval("chr1", 2150, 5000), geneFeaturePos, "3+", ""},
                {new SimpleInterval("chr1", 2150, 5000), geneFeatureNeg, "0-", ""},

                // Segment start breakpoint is in the transcript (no overlap/ in the intron)
                {new SimpleInterval("chr1", 1350, 5000), geneFeaturePos, "1+", ""},
                {new SimpleInterval("chr1", 1350, 5000), geneFeatureNeg, "2-", ""},
                {new SimpleInterval("chr1", 1650, 5000), geneFeaturePos, "2+", ""},
                {new SimpleInterval("chr1", 1650, 5000), geneFeatureNeg, "1-", ""},
                {new SimpleInterval("chr1", 1950, 5000), geneFeaturePos, "3+", ""},
                {new SimpleInterval("chr1", 1950, 5000), geneFeatureNeg, "0-", ""},

                // Segment start and end breakpoints are in the transcript (overlap)
                {new SimpleInterval("chr1", 1050, 1550), geneFeaturePos, "0+", "1-"},
                {new SimpleInterval("chr1", 1050, 1550), geneFeatureNeg, "3-", "2+"},
                {new SimpleInterval("chr1", 1550, 1850), geneFeaturePos, "1+", "2-"},
                {new SimpleInterval("chr1", 1550, 1850), geneFeatureNeg, "2-", "1+"},
                {new SimpleInterval("chr1", 1850, 2150), geneFeaturePos, "2+", "3-"},
                {new SimpleInterval("chr1", 1850, 2150), geneFeatureNeg, "1-", "0+"},
                {new SimpleInterval("chr1", 1050, 1850), geneFeaturePos, "0+", "2-"},
                {new SimpleInterval("chr1", 1050, 1850), geneFeatureNeg, "3-", "1+"},
                {new SimpleInterval("chr1", 1550, 2150), geneFeaturePos, "1+", "3-"},
                {new SimpleInterval("chr1", 1550, 2150), geneFeatureNeg, "2-", "0+"},


                // Segment start and end breakpoints are in the transcript (no overlap/in the intron)
                {new SimpleInterval("chr1", 1350, 1650), geneFeaturePos, "1+", "1-"},
                {new SimpleInterval("chr1", 1350, 1650), geneFeatureNeg, "2-", "2+"},
                {new SimpleInterval("chr1", 1650, 1950), geneFeaturePos, "2+", "2-"},
                {new SimpleInterval("chr1", 1650, 1950), geneFeatureNeg, "1-", "1+"},
                {new SimpleInterval("chr1", 1350, 1950), geneFeaturePos, "1+", "2-"},
                {new SimpleInterval("chr1", 1350, 1950), geneFeatureNeg, "2-", "1+"},


                // Segment end breakpoint is in the transcript (overlap)
                {new SimpleInterval("chr1", 50, 1050), geneFeaturePos, "", "0-"},
                {new SimpleInterval("chr1", 50, 1050), geneFeatureNeg, "", "3+"},
                {new SimpleInterval("chr1", 50, 1550), geneFeaturePos, "", "1-"},
                {new SimpleInterval("chr1", 50, 1550), geneFeatureNeg, "", "2+"},
                {new SimpleInterval("chr1", 50, 1850), geneFeaturePos, "", "2-"},
                {new SimpleInterval("chr1", 50, 1850), geneFeatureNeg, "", "1+"},
                {new SimpleInterval("chr1", 50, 2150), geneFeaturePos, "", "3-"},
                {new SimpleInterval("chr1", 50, 2150), geneFeatureNeg, "", "0+"},

                // Segment end breakpoint is in the transcript (no overlap/in the intron)
                {new SimpleInterval("chr1", 50, 1350), geneFeaturePos, "", "0-"},
                {new SimpleInterval("chr1", 50, 1350), geneFeatureNeg, "", "3+"},
                {new SimpleInterval("chr1", 50, 1650), geneFeaturePos, "", "1-"},
                {new SimpleInterval("chr1", 50, 1650), geneFeatureNeg, "", "2+"},
                {new SimpleInterval("chr1", 50, 1950), geneFeaturePos, "", "2-"},
                {new SimpleInterval("chr1", 50, 1950), geneFeatureNeg, "", "1+"},


                // Segment breakpoints do not overlap the transcript.
                {new SimpleInterval("chr1", 4000, 5000), geneFeaturePos, "", ""},
                {new SimpleInterval("chr1", 4000, 5000), geneFeatureNeg, "", ""},

                // Segment breakpoints do not overlap the transcript.
                {new SimpleInterval("chr1", 1, 500), geneFeaturePos, "", ""},
                {new SimpleInterval("chr1", 1, 500), geneFeatureNeg, "", ""},

                // Segment breakpoints do not overlap the transcript, though the transcript does.
                {new SimpleInterval("chr1", 1, 5000), geneFeaturePos, "", ""},
                {new SimpleInterval("chr1", 1, 5000), geneFeatureNeg, "", ""}
        };
    }

    @Test(dataProvider = "provideExonSegmentCalculation")
    public void testExonSegmentCalculation(final SimpleInterval segment, final GencodeGtfGeneFeature geneFeature, final String gtStart, final String gtEnd) {

        final SegmentExonOverlaps guess = SegmentExonUtils.determineSegmentExonPosition(geneFeature.getTranscripts().get(0), segment);

        Assert.assertEquals(guess.getSegmentStartExonOverlap(), gtStart, "Start segment exon overlap incorrect");
        Assert.assertEquals(guess.getSegmentEndExonOverlap(), gtEnd, "End segment exon overlap incorrect");
    }
}
