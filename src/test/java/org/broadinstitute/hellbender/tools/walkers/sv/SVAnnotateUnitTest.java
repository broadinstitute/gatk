package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.bed.FullBEDFeature;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfGeneFeature;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfTranscriptFeature;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class SVAnnotateUnitTest extends GATKBaseTest {
    protected final File TOY_GTF_FILE = new File(getToolTestDataDir() + "unittest.gtf");
    protected final File TINY_NONCODING_BED_FILE = new File(getToolTestDataDir() + "noncoding.unittest.bed");

    /**
     * Load toy GTF containing toy genes into FeatureDataSource
     * Visible for use in SVAnnotateEngineUnitTest
     * @return - toy GTF loaded into FeatureDataSource
     */
    protected static FeatureDataSource<GencodeGtfGeneFeature> loadToyGTFSource(File toyGtfFile) {
        return new FeatureDataSource<>(toyGtfFile);
    }

    // Load toy GTF and check for expected promoter and TSS intervals
    @Test
    public void testGetTranscriptionStartSiteAndPromoterInterval() {
        Map<String, Integer> tssByGene = new HashMap<>();
        tssByGene.put("EMMA1", 100);
        tssByGene.put("EMMA2", 3000);

        Map<String, SimpleInterval> promoterByGene = new HashMap<>();
        promoterByGene.put("EMMA1", new SimpleInterval("chr1", 1,99));
        promoterByGene.put("EMMA2", new SimpleInterval("chr1", 3001, 4000));

        int promoterWindow = 1000;

        FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource(TOY_GTF_FILE);
        for (final GencodeGtfGeneFeature gene : toyGTFSource) {
            GencodeGtfTranscriptFeature transcript = gene.getTranscripts().get(0); // each toy gene only has one transcript
            String geneName = transcript.getGeneName();
            if (tssByGene.containsKey(geneName)) {
                int expectedTSS = tssByGene.get(geneName);
                Assert.assertEquals(SVAnnotate.getTranscriptionStartSite(transcript), expectedTSS);

                Assert.assertEquals(SVAnnotate.getPromoterInterval(transcript, promoterWindow),
                        promoterByGene.get(geneName));
            }
        }
    }

    /**
     * Load tiny noncoding BED file into FeatureDataSource
     * Tiny noncoding BED file has a few noncoding elements to test noncoding annotation
     * Visible for use in SVAnnotateEngineUnitTest
     * @return - tiny noncoding BED file loaded into FeatureDataSource
     */
    protected static FeatureDataSource<FullBEDFeature> loadTinyNoncodingBEDSource(File tinyNonCodingBedFile) {
        return new FeatureDataSource<>(tinyNonCodingBedFile);
    }

    /**
     * Create synthetic SAMSequenceDictionary for provided contig list. Convenience function for testing.
     * Visible for use in tests for other classes.
     * @param contigs - list of contig names (Strings) in the desired order
     * @return - SAMSequenceDictionary containing the provided contigs, in the provided order,
     *          with a length of 100000 for each
     */
    public static SAMSequenceDictionary createSequenceDictionary(final List<String> contigs) {
        final int toyContigLen = 100000;
        final int n = contigs.size();
        final List<SAMSequenceRecord> sequenceRecords = new ArrayList<>(n);
        for (final String contig : contigs) {
            sequenceRecords.add(new SAMSequenceRecord(contig, toyContigLen));
        }
        return new SAMSequenceDictionary(sequenceRecords);
    }

    // Test building interval trees (GTF and BED) with different sequence dictionaries
    // and expected number of elements in BED interval tree & transcript interval tree
    @DataProvider(name = "buildIntervalTrees")
    public Object[][] getBuildIntervalTreesTestData() {
        return new Object[][] {
                { createSequenceDictionary(Arrays.asList("chr1")), 2, 2 },
                { createSequenceDictionary(Arrays.asList("chr1", "chr2")), 3, 3 }
        };
    }

    // Test building BED and GTF interval trees from toy files with different sequence dictionaries
    @Test(dataProvider = "buildIntervalTrees")
    public void testIgnoreUnknownContigsWhenBuildingIntervalTrees(
            final SAMSequenceDictionary sequenceDictionary,
            final int expectedBEDTreeSize,
            final int expectedTranscriptTreeSize
    ) {
        final FeatureDataSource<FullBEDFeature> tinyNoncodingBedSource = loadTinyNoncodingBEDSource(TINY_NONCODING_BED_FILE);
        final SVIntervalTree<String> nonCodingIntervalTree =
                SVAnnotate.buildIntervalTreeFromBED(tinyNoncodingBedSource, sequenceDictionary);
        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource(TOY_GTF_FILE);
        final SVAnnotateEngine.GTFIntervalTreesContainer gtfTrees =
                SVAnnotate.buildIntervalTreesFromGTF(toyGTFSource, sequenceDictionary, 100);
        // check size to ensure contigs not included in the map are excluded from the interval tree successfully
        Assert.assertEquals(nonCodingIntervalTree.size(), expectedBEDTreeSize);
        Assert.assertEquals(gtfTrees.getTranscriptIntervalTree().size(), expectedTranscriptTreeSize);
    }
}
