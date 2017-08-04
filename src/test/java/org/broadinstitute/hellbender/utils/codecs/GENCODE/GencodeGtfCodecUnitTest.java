package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Files;
import java.util.*;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;


/**
 * Test class for the GENCODE GTF Reader.
 * Modeled after the TableCodecUnitTest, with extras specific to this file format.
 * Created by jonn on 7/27/17.
 */
public class GencodeGtfCodecUnitTest extends BaseTest {

    private static final String testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/codecs/GENCODE/";
    private static final String xyzTestFile = largeFileTestDir + "gencode.v26.primary_assembly.annotation.XYZ.gtf";
    private static final String gencodeHg19TestFile = largeFileTestDir + "gencode.v19.LargeFile.gtf";

    /**
     * Checks the given feature and the given start and end positions for whether the regions overlap.
     * @param feature {@link GencodeGtfFeature} to use for overlap checking
     * @param contig Contig to check for overlap
     * @param start Interval start position to check for overlap
     * @param end Interval end position to check for overlap
     * @return {@code true} if the region in the {@link GencodeGtfFeature} and the given interval overlap, {@code false} otherwise.
     */
    private boolean checkForOverlap(final GencodeGtfFeature feature,
                                    final String contig,
                                    final int start,
                                    final int end) {

        boolean overlaps = feature.getChromosomeName().equals(contig);

        // Check for any overlap.
        // This includes overlapping on either end, as well as either region
        // being contained within the other.
        overlaps = overlaps &&
                ((start >= feature.getStart() ) && (start <= feature.getEnd())) ||
                ((end >= feature.getStart() ) && (end <= feature.getEnd())) ||
                ((start <= feature.getStart() ) && (end >= feature.getEnd()));

        return overlaps;
    }

    /**
     * Tests that a given query in a file returns the correct number of results
     * @param contig Contiguous region in the genome in which to search.
     * @param start Start position within {@code contig} in which to search.
     * @param end End position within {@code contig} in which to search.
     * @param numExpectedGenes The number of expected results from the query.
     * @param testFile A GENCODE GTF {@link File} to query against.
     */
    private void testIndexHelper(String contig, int start, int end, int numExpectedGenes, File testFile) {
        // Now we do our queries:
        try (FeatureDataSource<GencodeGtfFeature> featureDataSource = new FeatureDataSource<>(testFile) )
        {
            final Iterator<GencodeGtfFeature> it = featureDataSource.query( new SimpleInterval(contig, start, end) );

            int geneCount = 0;

            for ( ; it.hasNext() ; ) {

                GencodeGtfFeature feature = it.next();

                // Verify the bounds:
                Assert.assertTrue( checkForOverlap(feature, contig, start, end) );

                // Keep track of how many genes we've seen:
                ++geneCount;
            }

            Assert.assertEquals( geneCount, numExpectedGenes );
        }
    }

    /**
     * Creates a valid {@link GencodeGtfGeneFeature} that corresponds to the data in the file {@code gencode.valid1.gtf}
     * @return a {@link GencodeGtfGeneFeature} representing the data in the file {@code gencode.valid1.gtf}
     */
    private GencodeGtfGeneFeature createGencodeGtfGene_valid1() {

        // Create the Features as they exist in the test file:
        GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature)GencodeGtfFeature.create(6, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.gene,
                30366, 30503, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000284332.1", null, GencodeGtfFeature.GeneTranscriptType.miRNA,
                null, "MIR1302-2", null, null, null, -1, null, GencodeGtfFeature.LocusLevel.THREE, null, null);

        GencodeGtfTranscriptFeature transcript = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(7, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.transcript,
                30366, 30503, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000284332.1", "ENST00000607096.1", GencodeGtfFeature.GeneTranscriptType.miRNA,
                null, "MIR1302-2", GencodeGtfFeature.GeneTranscriptType.miRNA, null, "MIR1302-2-201", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.NA),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic)
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon = (GencodeGtfExonFeature) GencodeGtfFeature.create(8, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                30366, 30503, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000284332.1", "ENST00000607096.1", GencodeGtfFeature.GeneTranscriptType.miRNA,
                null, "MIR1302-2", GencodeGtfFeature.GeneTranscriptType.miRNA, null, "MIR1302-2-201", 1, "ENSE00003695741.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.NA),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic)
                        )
                ),
                null
        );

        // Aggregate the Features as they should be:
        transcript.addExon(exon);
        gene.addTranscript(transcript);

        return gene;
    }

    /**
     * Creates a valid {@link GencodeGtfGeneFeature} that corresponds to the data in the file {@code gencode.valid_gencode_file2.gtf}
     * @return a {@link GencodeGtfGeneFeature} representing the data in the file {@code gencode.valid_gencode_file2.gtf}
     */
    private GencodeGtfGeneFeature createGencodeGtfGene_file2() {

        // Let's define all our features up front:

        GencodeGtfGeneFeature gene1 = (GencodeGtfGeneFeature) GencodeGtfFeature.create(6, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.gene,
                50200979, 50217616, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", null, GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", null, null, null, -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript1 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(7, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.transcript,
                50200979, 50217615, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon1 = (GencodeGtfExonFeature) GencodeGtfFeature.create(8, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50200979, 50201590, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds1 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(9, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50201037, 50201590, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfStartCodonFeature start_codon1 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(10, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.start_codon,
                50201037, 50201039, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon2 = (GencodeGtfExonFeature) GencodeGtfFeature.create(11, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50206317, 50206520, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 2, "ENSE00001129529.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds2 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(12, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50206317, 50206520, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 2, "ENSE00001129529.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon3 = (GencodeGtfExonFeature) GencodeGtfFeature.create(13, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50208536, 50208716, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 3, "ENSE00001129524.2", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds3 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(14, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50208536, 50208716, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 3, "ENSE00001129524.2", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon4 = (GencodeGtfExonFeature) GencodeGtfFeature.create(15, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50210181, 50210311, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 4, "ENSE00003473644.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds4 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(16, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50210181, 50210311, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 4, "ENSE00003473644.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon5 = (GencodeGtfExonFeature) GencodeGtfFeature.create(17, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50210631, 50210911, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 5, "ENSE00003503715.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds5 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(18, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50210631, 50210911, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 5, "ENSE00003503715.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon6 = (GencodeGtfExonFeature) GencodeGtfFeature.create(19, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50215717, 50215867, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 6, "ENSE00003573348.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds6 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(20, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50215717, 50215867, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 6, "ENSE00003573348.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon7 = (GencodeGtfExonFeature) GencodeGtfFeature.create(21, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50216691, 50216876, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 7, "ENSE00003510005.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds7 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(22, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50216691, 50216876, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 7, "ENSE00003510005.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon8 = (GencodeGtfExonFeature) GencodeGtfFeature.create(23, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50216972, 50217128, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 8, "ENSE00003591346.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds8 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(24, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50216972, 50217128, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 8, "ENSE00003591346.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon9 = (GencodeGtfExonFeature) GencodeGtfFeature.create(25, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50217205, 50217357, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 9, "ENSE00003728455.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds9 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(26, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50217205, 50217357, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 9, "ENSE00003728455.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon10 = (GencodeGtfExonFeature) GencodeGtfFeature.create(27, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                50217361, 50217615, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 10, "ENSE00003739808.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds10 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(28, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                50217361, 50217366, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 10, "ENSE00003739808.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfStopCodonFeature stop_codon1 = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(29, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.stop_codon,
                50217367, 50217369, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 10, "ENSE00003739808.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr1 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(30, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                50200979, 50201036, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr2 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(31, "chr22", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                50217367, 50217615, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000611222.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-201", 10, "ENSE00003739808.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000483593.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.FIVE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_alternative_2),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript2 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(32, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.transcript,
                50200979, 50217616, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfSelenocysteineFeature selenocysteine1 = (GencodeGtfSelenocysteineFeature) GencodeGtfFeature.create(33, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.selenocysteine,
                50217358, 50217360, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon11 = (GencodeGtfExonFeature) GencodeGtfFeature.create(34, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50200979, 50201590, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds11 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(35, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50201037, 50201590, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfStartCodonFeature start_codon2 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(36, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.start_codon,
                50201037, 50201039, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon12 = (GencodeGtfExonFeature) GencodeGtfFeature.create(37, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50206317, 50206520, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 2, "ENSE00001129529.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds12 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(38, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50206317, 50206520, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 2, "ENSE00001129529.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon13 = (GencodeGtfExonFeature) GencodeGtfFeature.create(39, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50208536, 50208716, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 3, "ENSE00001129524.2", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds13 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(40, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50208536, 50208716, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 3, "ENSE00001129524.2", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon14 = (GencodeGtfExonFeature) GencodeGtfFeature.create(41, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50210181, 50210311, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 4, "ENSE00003473644.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds14 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(42, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50210181, 50210311, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 4, "ENSE00003473644.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon15 = (GencodeGtfExonFeature) GencodeGtfFeature.create(43, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50210631, 50210911, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 5, "ENSE00003503715.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds15 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(44, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50210631, 50210911, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 5, "ENSE00003503715.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon16 = (GencodeGtfExonFeature) GencodeGtfFeature.create(45, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50215717, 50215867, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 6, "ENSE00003573348.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds16 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(46, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50215717, 50215867, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 6, "ENSE00003573348.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon17 = (GencodeGtfExonFeature) GencodeGtfFeature.create(47, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50216691, 50216876, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 7, "ENSE00003510005.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds17 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(48, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50216691, 50216876, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 7, "ENSE00003510005.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon18 = (GencodeGtfExonFeature) GencodeGtfFeature.create(49, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50216972, 50217128, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 8, "ENSE00003591346.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds18 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(50, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50216972, 50217128, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 8, "ENSE00003591346.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon19 = (GencodeGtfExonFeature) GencodeGtfFeature.create(51, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50217205, 50217616, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 9, "ENSE00003512975.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds19 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(52, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50217205, 50217366, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 9, "ENSE00003512975.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfStopCodonFeature stop_codon2 = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(53, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.stop_codon,
                50217367, 50217369, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 9, "ENSE00003512975.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr3 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(54, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.utr,
                50200979, 50201036, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr4 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(55, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.utr,
                50217367, 50217616, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000380903.6", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "SELENOO-001", 9, "ENSE00003512975.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000370288.2"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal_2),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript3 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(56, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.transcript,
                50206442, 50217616, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon20 = (GencodeGtfExonFeature) GencodeGtfFeature.create(57, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50206442, 50206520, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 1, "ENSE00001890724.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon21 = (GencodeGtfExonFeature) GencodeGtfFeature.create(58, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50208488, 50208716, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 2, "ENSE00001952603.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon22 = (GencodeGtfExonFeature) GencodeGtfFeature.create(59, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50210181, 50210311, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 3, "ENSE00003583919.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon23 = (GencodeGtfExonFeature) GencodeGtfFeature.create(60, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50210631, 50210911, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 4, "ENSE00003620115.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon24 = (GencodeGtfExonFeature) GencodeGtfFeature.create(61, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50215717, 50215867, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 5, "ENSE00003636069.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon25 = (GencodeGtfExonFeature) GencodeGtfFeature.create(62, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50216691, 50216876, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 6, "ENSE00003579717.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon26 = (GencodeGtfExonFeature) GencodeGtfFeature.create(63, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50216972, 50217128, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 7, "ENSE00003650938.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon27 = (GencodeGtfExonFeature) GencodeGtfFeature.create(64, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50217205, 50217616, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.13", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "SELENOO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, null, "SELENOO-002", 8, "ENSE00003475904.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        // ======================

        // Now let's collapse these objects into their correct structure:
        exon1.setCds(cds1);
        exon1.setStartCodon(start_codon1);

        exon2.setCds(cds2);
        exon3.setCds(cds3);
        exon4.setCds(cds4);
        exon5.setCds(cds5);
        exon6.setCds(cds6);
        exon7.setCds(cds7);
        exon8.setCds(cds8);
        exon9.setCds(cds9);

        exon10.setCds(cds10);
        exon10.setStopCodon(stop_codon1);

        transcript1.addExon(exon1);
        transcript1.addExon(exon2);
        transcript1.addExon(exon3);
        transcript1.addExon(exon4);
        transcript1.addExon(exon5);
        transcript1.addExon(exon6);
        transcript1.addExon(exon7);
        transcript1.addExon(exon8);
        transcript1.addExon(exon9);
        transcript1.addExon(exon10);

        transcript1.addUtr(utr1);
        transcript1.addUtr(utr2);

        gene1.addTranscript(transcript1);

        // ======================

        transcript2.addSelenocysteine(selenocysteine1);

        exon11.setCds(cds11);
        exon11.setStartCodon(start_codon2);

        exon12.setCds(cds12);
        exon13.setCds(cds13);
        exon14.setCds(cds14);
        exon15.setCds(cds15);
        exon16.setCds(cds16);
        exon17.setCds(cds17);
        exon18.setCds(cds18);

        exon19.setCds(cds19);
        exon19.setStopCodon(stop_codon2);

        transcript2.addExon(exon11);
        transcript2.addExon(exon12);
        transcript2.addExon(exon13);
        transcript2.addExon(exon14);
        transcript2.addExon(exon15);
        transcript2.addExon(exon16);
        transcript2.addExon(exon17);
        transcript2.addExon(exon18);
        transcript2.addExon(exon19);

        transcript2.addUtr(utr3);
        transcript2.addUtr(utr4);

        gene1.addTranscript(transcript2);

        // ======================

        transcript3.addExon(exon20);
        transcript3.addExon(exon21);
        transcript3.addExon(exon22);
        transcript3.addExon(exon23);
        transcript3.addExon(exon24);
        transcript3.addExon(exon25);
        transcript3.addExon(exon26);
        transcript3.addExon(exon27);

        gene1.addTranscript(transcript3);

        // ======================

        return gene1;
    }

    /**
     * Creates a valid {@link GencodeGtfGeneFeature} that corresponds to the data in the file {@code gencode.valid_gencode_file2.gtf}
     * @return a {@link GencodeGtfGeneFeature} representing the data in the file {@code gencode.valid_gencode_file2.gtf}
     */
    private GencodeGtfGeneFeature createGencodeGtfGene_file3() {

        GencodeGtfGeneFeature gene1 = (GencodeGtfGeneFeature) GencodeGtfFeature.create(6, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.gene,
                138082, 161852, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", null, GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", null, null, null, -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(), null
        );


        GencodeGtfTranscriptFeature transcript1 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(7, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.transcript,
                138082, 161750, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon1 = (GencodeGtfExonFeature) GencodeGtfFeature.create(8, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                161689, 161750, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 1, "ENSE00003735197.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon2 = (GencodeGtfExonFeature) GencodeGtfFeature.create(9, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                156289, 156497, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 2, "ENSE00003737280.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds1 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(10, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                156289, 156446, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 2, "ENSE00003737280.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfStartCodonFeature start_codon1 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(11, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.start_codon,
                156444, 156446, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 2, "ENSE00003737280.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon3 = (GencodeGtfExonFeature) GencodeGtfFeature.create(12, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                150987, 151021, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 3, "ENSE00003731891.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds2 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(13, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                150987, 151021, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 3, "ENSE00003731891.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon4 = (GencodeGtfExonFeature) GencodeGtfFeature.create(14, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                150350, 150499, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 4, "ENSE00003724613.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds3 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(15, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                150350, 150499, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 4, "ENSE00003724613.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon5 = (GencodeGtfExonFeature) GencodeGtfFeature.create(16, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                148414, 148478, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 5, "ENSE00003732418.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds4 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(17, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                148414, 148478, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 5, "ENSE00003732418.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon6 = (GencodeGtfExonFeature) GencodeGtfFeature.create(18, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                148116, 148232, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 6, "ENSE00003733960.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds5 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(19, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                148116, 148232, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 6, "ENSE00003733960.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon7 = (GencodeGtfExonFeature) GencodeGtfFeature.create(20, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                147624, 147703, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 7, "ENSE00003727207.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds6 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(21, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                147624, 147703, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 7, "ENSE00003727207.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon8 = (GencodeGtfExonFeature) GencodeGtfFeature.create(22, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                146640, 146721, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 8, "ENSE00003728972.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds7 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(23, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                146640, 146721, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 8, "ENSE00003728972.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon9 = (GencodeGtfExonFeature) GencodeGtfFeature.create(24, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                145004, 145096, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 9, "ENSE00003733844.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds8 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(25, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                145004, 145096, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 9, "ENSE00003733844.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon10 = (GencodeGtfExonFeature) GencodeGtfFeature.create(26, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                144749, 144895, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 10, "ENSE00003752738.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds9 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(27, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                144749, 144895, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 10, "ENSE00003752738.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon11 = (GencodeGtfExonFeature) GencodeGtfFeature.create(28, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                143614, 143789, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 11, "ENSE00003720006.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds10 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(29, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                143614, 143789, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 11, "ENSE00003720006.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon12 = (GencodeGtfExonFeature) GencodeGtfFeature.create(30, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                142194, 142292, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 12, "ENSE00003719283.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds11 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(31, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                142194, 142292, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 12, "ENSE00003719283.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon13 = (GencodeGtfExonFeature) GencodeGtfFeature.create(32, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                138743, 138831, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 13, "ENSE00003751415.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds12 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(33, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                138743, 138831, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 13, "ENSE00003751415.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon14 = (GencodeGtfExonFeature) GencodeGtfFeature.create(34, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                138082, 138667, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 14, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds13 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(35, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                138483, 138667, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 14, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfStopCodonFeature stop_codon1 = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(36, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.stop_codon,
                138480, 138482, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 14, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr1 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(37, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                161689, 161750, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 1, "ENSE00003735197.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr2 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(38, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                156447, 156497, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 2, "ENSE00003737280.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr3 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(39, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                138082, 138482, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000615165.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-202", 14, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000482462.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript2 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(40, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.transcript,
                138082, 161852, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon15 = (GencodeGtfExonFeature) GencodeGtfFeature.create(41, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                161689, 161852, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 1, "ENSE00003746084.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon16 = (GencodeGtfExonFeature) GencodeGtfFeature.create(42, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                161314, 161626, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 2, "ENSE00003719550.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds14 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(43, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                161314, 161586, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 2, "ENSE00003719550.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfStartCodonFeature start_codon2 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(44, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.start_codon,
                161584, 161586, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 2, "ENSE00003719550.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon17 = (GencodeGtfExonFeature) GencodeGtfFeature.create(45, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                156289, 156497, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 3, "ENSE00003723757.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds15 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(46, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                156289, 156497, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 3, "ENSE00003723757.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon18 = (GencodeGtfExonFeature) GencodeGtfFeature.create(47, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                150987, 151021, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 4, "ENSE00003731891.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds16 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(48, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                150987, 151021, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 4, "ENSE00003731891.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon19 = (GencodeGtfExonFeature) GencodeGtfFeature.create(49, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                150350, 150499, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 5, "ENSE00003724613.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds17 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(50, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                150350, 150499, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 5, "ENSE00003724613.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon20 = (GencodeGtfExonFeature) GencodeGtfFeature.create(51, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                148414, 148478, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 6, "ENSE00003732418.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds18 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(52, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                148414, 148478, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 6, "ENSE00003732418.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon21 = (GencodeGtfExonFeature) GencodeGtfFeature.create(53, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                148116, 148232, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 7, "ENSE00003733960.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds19 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(54, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                148116, 148232, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 7, "ENSE00003733960.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon22 = (GencodeGtfExonFeature) GencodeGtfFeature.create(55, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                147624, 147703, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 8, "ENSE00003727207.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds20 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(56, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                147624, 147703, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 8, "ENSE00003727207.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon23 = (GencodeGtfExonFeature) GencodeGtfFeature.create(57, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                146640, 146721, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 9, "ENSE00003728972.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds21 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(58, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                146640, 146721, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 9, "ENSE00003728972.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon24 = (GencodeGtfExonFeature) GencodeGtfFeature.create(59, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                145004, 145096, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 10, "ENSE00003733844.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds22 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(60, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                145004, 145096, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 10, "ENSE00003733844.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon25 = (GencodeGtfExonFeature) GencodeGtfFeature.create(61, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                144749, 144895, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 11, "ENSE00003752738.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds23 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(62, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                144749, 144895, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 11, "ENSE00003752738.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon26 = (GencodeGtfExonFeature) GencodeGtfFeature.create(63, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                143614, 143789, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 12, "ENSE00003720006.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds24 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(64, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                143614, 143789, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 12, "ENSE00003720006.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon27 = (GencodeGtfExonFeature) GencodeGtfFeature.create(65, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                142194, 142292, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 13, "ENSE00003719283.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds25 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(66, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                142194, 142292, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 13, "ENSE00003719283.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon28 = (GencodeGtfExonFeature) GencodeGtfFeature.create(67, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                138743, 138831, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 14, "ENSE00003751415.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds26 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(68, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                138743, 138831, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 14, "ENSE00003751415.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon29 = (GencodeGtfExonFeature) GencodeGtfFeature.create(69, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                138082, 138667, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 15, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds27 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(70, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                138483, 138667, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 15, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfStopCodonFeature stop_codon2 = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(71, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.stop_codon,
                138480, 138482, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 15, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr4 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(72, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                161689, 161852, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 1, "ENSE00003746084.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr5 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(73, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                161587, 161626, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 2, "ENSE00003719550.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr6 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(74, "KI270734.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                138082, 138482, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000277196.4", "ENST00000621424.4", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                null, "AC007325.2", GencodeGtfFeature.GeneTranscriptType.protein_coding, null, "AC007325.2-201", 15, "ENSE00003753010.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("protein_id", "ENSP00000481127.1"),
                                new GencodeGtfFeature.OptionalField<>("transcript_support_level", GencodeGtfFeature.TranscriptSupportLevel.ONE),
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        // ======================
        // Now let's collapse these objects into their correct structure:

        exon2.setCds(cds1);
        exon2.setStartCodon(start_codon1);

        exon3.setCds(cds2);
        exon4.setCds(cds3);
        exon5.setCds(cds4);
        exon6.setCds(cds5);
        exon7.setCds(cds6);
        exon8.setCds(cds7);
        exon9.setCds(cds8);
        exon10.setCds(cds9);
        exon11.setCds(cds10);
        exon12.setCds(cds11);
        exon13.setCds(cds12);

        exon14.setCds(cds13);
        exon14.setStopCodon(stop_codon1);

        transcript1.addExon(exon1);
        transcript1.addExon(exon2);
        transcript1.addExon(exon3);
        transcript1.addExon(exon4);
        transcript1.addExon(exon5);
        transcript1.addExon(exon6);
        transcript1.addExon(exon7);
        transcript1.addExon(exon8);
        transcript1.addExon(exon9);
        transcript1.addExon(exon10);
        transcript1.addExon(exon11);
        transcript1.addExon(exon12);
        transcript1.addExon(exon13);
        transcript1.addExon(exon14);

        transcript1.addUtr(utr1);
        transcript1.addUtr(utr2);
        transcript1.addUtr(utr3);

        gene1.addTranscript(transcript1);

        // ======================

        exon16.setCds(cds14);
        exon16.setStartCodon(start_codon2);

        exon17.setCds(cds15);
        exon18.setCds(cds16);
        exon19.setCds(cds17);
        exon20.setCds(cds18);
        exon21.setCds(cds19);
        exon22.setCds(cds20);
        exon23.setCds(cds21);
        exon24.setCds(cds22);
        exon25.setCds(cds23);
        exon26.setCds(cds24);
        exon27.setCds(cds25);
        exon28.setCds(cds26);

        exon29.setCds(cds27);
        exon29.setStopCodon(stop_codon2);

        transcript2.addExon(exon15);
        transcript2.addExon(exon16);
        transcript2.addExon(exon17);
        transcript2.addExon(exon18);
        transcript2.addExon(exon19);
        transcript2.addExon(exon20);
        transcript2.addExon(exon21);
        transcript2.addExon(exon22);
        transcript2.addExon(exon23);
        transcript2.addExon(exon24);
        transcript2.addExon(exon25);
        transcript2.addExon(exon26);
        transcript2.addExon(exon27);
        transcript2.addExon(exon28);
        transcript2.addExon(exon29);

        transcript2.addUtr(utr4);
        transcript2.addUtr(utr5);
        transcript2.addUtr(utr6);

        gene1.addTranscript(transcript2);

        // ======================

        return gene1;
    }

    private GencodeGtfGeneFeature createGencodeGtfGene_v19_valid1() {
        // Create the Features as they exist in the test file:
        GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature)GencodeGtfFeature.create(6, "chr1", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.gene,
                11869, 14412, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000223972.4", "ENSG00000223972.4",
                GencodeGtfFeature.GeneTranscriptType.pseudogene, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "DDX11L1", GencodeGtfFeature.GeneTranscriptType.pseudogene,
                 GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "DDX11L1", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000000961.2")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(7, "chr1", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.transcript,
                11869, 14409, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000223972.4", "ENST00000456328.2",
                GencodeGtfFeature.GeneTranscriptType.pseudogene, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "DDX11L1", GencodeGtfFeature.GeneTranscriptType.processed_transcript,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "DDX11L1-002", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000000961.2"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000362751.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon = (GencodeGtfExonFeature) GencodeGtfFeature.create(8, "chr1", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                11869, 12227, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000223972.4", "ENST00000456328.2",
                GencodeGtfFeature.GeneTranscriptType.pseudogene, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "DDX11L1", GencodeGtfFeature.GeneTranscriptType.processed_transcript,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "DDX11L1-002", 1, "ENSE00002234944.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000000961.2"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000362751.1")
                        )
                ),
                null
        );

        // Aggregate the Features as they should be:
        transcript.addExon(exon);
        gene.addTranscript(transcript);

        return gene;
    }

    private ArrayList<GencodeGtfGeneFeature> createGencodeGtfGene_v19_file2() {

        GencodeGtfGeneFeature gene1 = (GencodeGtfGeneFeature) GencodeGtfFeature.create(6, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.gene,
                50637519, 50638976, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000273253.1", "ENSG00000273253.1", GencodeGtfFeature.GeneTranscriptType.antisense,
                GencodeGtfFeature.GeneTranscriptStatus.NOVEL, "RP3-402G11.26", GencodeGtfFeature.GeneTranscriptType.antisense, GencodeGtfFeature.GeneTranscriptStatus.NOVEL, "RP3-402G11.26", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000186123.2")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript1 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(7, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.transcript,
                50637519, 50638976, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000273253.1", "ENST00000608025.1", GencodeGtfFeature.GeneTranscriptType.antisense,
                GencodeGtfFeature.GeneTranscriptStatus.NOVEL, "RP3-402G11.26", GencodeGtfFeature.GeneTranscriptType.antisense, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "RP3-402G11.26-001", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000186123.2"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000472292.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon1 = (GencodeGtfExonFeature) GencodeGtfFeature.create(8, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50638505, 50638976, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000273253.1", "ENST00000608025.1", GencodeGtfFeature.GeneTranscriptType.antisense,
                GencodeGtfFeature.GeneTranscriptStatus.NOVEL, "RP3-402G11.26", GencodeGtfFeature.GeneTranscriptType.antisense, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "RP3-402G11.26-001", 1, "ENSE00003710600.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000186123.2"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000472292.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon2 = (GencodeGtfExonFeature) GencodeGtfFeature.create(9, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50637519, 50637757, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000273253.1", "ENST00000608025.1", GencodeGtfFeature.GeneTranscriptType.antisense,
                GencodeGtfFeature.GeneTranscriptStatus.NOVEL, "RP3-402G11.26", GencodeGtfFeature.GeneTranscriptType.antisense, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "RP3-402G11.26-001", 2, "ENSE00003710731.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000186123.2"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000472292.2")
                        )
                ),
                null
        );

        GencodeGtfGeneFeature gene2 = (GencodeGtfGeneFeature) GencodeGtfFeature.create(10, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.gene,
                50639408, 50656045, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENSG00000073169.9", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript2 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(11, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.transcript,
                50639408, 50656045, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfSelenocysteineFeature selenocysteine1 = (GencodeGtfSelenocysteineFeature) GencodeGtfFeature.create(12, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.selenocysteine,
                50655787, 50655789, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon3 = (GencodeGtfExonFeature) GencodeGtfFeature.create(13, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50639408, 50640019, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds1 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(14, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50639466, 50640019, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfStartCodonFeature start_codon1 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(15, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.start_codon,
                50639466, 50639468, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 1, "ENSE00001541223.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon4 = (GencodeGtfExonFeature) GencodeGtfFeature.create(16, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50644746, 50644949, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 2, "ENSE00001129529.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds2 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(17, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50644746, 50644949, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 2, "ENSE00001129529.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon5 = (GencodeGtfExonFeature) GencodeGtfFeature.create(18, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50646965, 50647145, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 3, "ENSE00001129524.2", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds3 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(19, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50646965, 50647145, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 3, "ENSE00001129524.2", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon6 = (GencodeGtfExonFeature) GencodeGtfFeature.create(20, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50648610, 50648740, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 4, "ENSE00003473644.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds4 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(21, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50648610, 50648740, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 4, "ENSE00003473644.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon7 = (GencodeGtfExonFeature) GencodeGtfFeature.create(22, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50649060, 50649340, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 5, "ENSE00003503715.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds5 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(23, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50649060, 50649340, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 5, "ENSE00003503715.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon8 = (GencodeGtfExonFeature) GencodeGtfFeature.create(24, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50654146, 50654296, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 6, "ENSE00003573348.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds6 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(25, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50654146, 50654296, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.TWO, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 6, "ENSE00003573348.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon9 = (GencodeGtfExonFeature) GencodeGtfFeature.create(26, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50655120, 50655305, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 7, "ENSE00003510005.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds7 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(27, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50655120, 50655305, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 7, "ENSE00003510005.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon10 = (GencodeGtfExonFeature) GencodeGtfFeature.create(28, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50655401, 50655557, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 8, "ENSE00003591346.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds8 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(29, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50655401, 50655557, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 8, "ENSE00003591346.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon11 = (GencodeGtfExonFeature) GencodeGtfFeature.create(30, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50655634, 50656045, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 9, "ENSE00003512975.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds9 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(31, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.cds,
                50655634, 50655795, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 9, "ENSE00003512975.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfStopCodonFeature stop_codon1 = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(32, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.stop_codon,
                50655796, 50655798, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", 9, "ENSE00003512975.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr1 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(33, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.utr,
                50639408, 50639465, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr2 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(34, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.utr,
                50655796, 50656045, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000380903.2", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-001", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.appris_principal),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.CCDS),
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.seleno),
                                new GencodeGtfFeature.OptionalField<>("ccdsid", "CCDS43034.1"),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000075003.2")
                        )
                ),
                null
        );

        GencodeGtfTranscriptFeature transcript3 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(35, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.transcript,
                50644871, 50656045, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", -1, null, GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon12 = (GencodeGtfExonFeature) GencodeGtfFeature.create(36, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50644871, 50644949, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 1, "ENSE00001890724.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon13 = (GencodeGtfExonFeature) GencodeGtfFeature.create(37, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50646917, 50647145, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 2, "ENSE00001952603.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon14 = (GencodeGtfExonFeature) GencodeGtfFeature.create(38, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50648610, 50648740, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 3, "ENSE00003583919.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon15 = (GencodeGtfExonFeature) GencodeGtfFeature.create(39, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50649060, 50649340, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 4, "ENSE00003620115.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon16 = (GencodeGtfExonFeature) GencodeGtfFeature.create(40, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50654146, 50654296, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 5, "ENSE00003636069.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon17 = (GencodeGtfExonFeature) GencodeGtfFeature.create(41, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50655120, 50655305, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 6, "ENSE00003579717.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon18 = (GencodeGtfExonFeature) GencodeGtfFeature.create(42, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50655401, 50655557, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 7, "ENSE00003650938.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon19 = (GencodeGtfExonFeature) GencodeGtfFeature.create(43, "chr22", GencodeGtfFeature.AnnotationSource.HAVANA, GencodeGtfFeature.FeatureType.exon,
                50655634, 50656045, GencodeGtfFeature.GenomicStrand.FORWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000073169.9", "ENST00000492092.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO", GencodeGtfFeature.GeneTranscriptType.processed_transcript, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "SELO-002", 8, "ENSE00003475904.1", GencodeGtfFeature.LocusLevel.TWO,
                new ArrayList<>(
                        Arrays.asList(
                                new GencodeGtfFeature.OptionalField<>("tag", GencodeGtfFeature.FeatureTag.basic),
                                new GencodeGtfFeature.OptionalField<>("havana_gene", "OTTHUMG00000044645.3"),
                                new GencodeGtfFeature.OptionalField<>("havana_transcript", "OTTHUMT00000316993.1")
                        )
                ),
                null
        );

        // ======================
        // Now let's collapse these objects into their correct structure:

        transcript1.addExon(exon1);
        transcript1.addExon(exon2);
        gene1.addTranscript(transcript1);

        // ======================
        // ======================

        transcript2.addSelenocysteine(selenocysteine1);

        exon3.setCds(cds1);
        exon3.setStartCodon(start_codon1);

        exon4.setCds(cds2);
        exon5.setCds(cds3);
        exon6.setCds(cds4);
        exon7.setCds(cds5);
        exon8.setCds(cds6);
        exon9.setCds(cds7);
        exon10.setCds(cds8);

        exon11.setCds(cds9);
        exon11.setStopCodon(stop_codon1);

        transcript2.addExon(exon3);
        transcript2.addExon(exon4);
        transcript2.addExon(exon5);
        transcript2.addExon(exon6);
        transcript2.addExon(exon7);
        transcript2.addExon(exon8);
        transcript2.addExon(exon9);
        transcript2.addExon(exon10);
        transcript2.addExon(exon11);

        transcript2.addUtr(utr1);
        transcript2.addUtr(utr2);

        // ======================

        transcript3.addExon(exon12);
        transcript3.addExon(exon13);
        transcript3.addExon(exon14);
        transcript3.addExon(exon15);
        transcript3.addExon(exon16);
        transcript3.addExon(exon17);
        transcript3.addExon(exon18);
        transcript3.addExon(exon19);

        // ======================

        gene2.addTranscript(transcript2);
        gene2.addTranscript(transcript3);

        return new ArrayList<>( Arrays.asList(gene1, gene2) );
    }

    private GencodeGtfGeneFeature createGencodeGtfGene_v19_theOtherFile() {

        GencodeGtfGeneFeature gene1 = (GencodeGtfGeneFeature) GencodeGtfFeature.create(6, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.gene,
                38792, 97421, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENSG00000215615.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(), null
        );


        GencodeGtfTranscriptFeature transcript1 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(7, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.transcript,
                38792, 97421, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon1 = (GencodeGtfExonFeature) GencodeGtfFeature.create(8, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                97368, 97421, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 1, "ENSE00001544212.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon2 = (GencodeGtfExonFeature) GencodeGtfFeature.create(9, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                95174, 95232, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 2, "ENSE00001849396.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds1 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(10, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                95174, 95230, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 2, "ENSE00001849396.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfStartCodonFeature start_codon1 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(11, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.start_codon,
                95228, 95230, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 2, "ENSE00001849396.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon3 = (GencodeGtfExonFeature) GencodeGtfFeature.create(12, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                40642, 40872, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 3, "ENSE00001900862.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds2 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(13, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                40642, 40872, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 3, "ENSE00001900862.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon4 = (GencodeGtfExonFeature) GencodeGtfFeature.create(14, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                39874, 40028, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 4, "ENSE00001544206.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds3 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(15, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                39874, 40028, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ZERO, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 4, "ENSE00001544206.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfExonFeature exon5 = (GencodeGtfExonFeature) GencodeGtfFeature.create(16, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.exon,
                38792, 39019, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 5, "ENSE00001544203.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfCDSFeature cds4 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(17, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.cds,
                38794, 39019, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.ONE, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", 5, "ENSE00001544203.1", GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr1 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(18, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                97368, 97421, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr2 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(19, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                95231, 95232, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        GencodeGtfUTRFeature utr3 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(20, "GL000218.1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.utr,
                38792, 38793, GencodeGtfFeature.GenomicStrand.BACKWARD, GencodeGtfFeature.GenomicPhase.DOT, "ENSG00000215615.1", "ENST00000400681.1", GencodeGtfFeature.GeneTranscriptType.protein_coding,
                GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1", GencodeGtfFeature.GeneTranscriptType.protein_coding, GencodeGtfFeature.GeneTranscriptStatus.KNOWN, "AL354822.1-201", -1, null, GencodeGtfFeature.LocusLevel.THREE,
                new ArrayList<>(
                        Collections.singletonList(
                                new GencodeGtfFeature.OptionalField<>("tag", "basic")
                        )
                ),
                null
        );

        // ======================
        // Now let's collapse these objects into their correct structure:

        exon2.setCds(cds1);
        exon2.setStartCodon(start_codon1);

        exon3.setCds(cds2);
        exon4.setCds(cds3);
        exon5.setCds(cds4);

        transcript1.addExon(exon1);
        transcript1.addExon(exon2);
        transcript1.addExon(exon3);
        transcript1.addExon(exon4);
        transcript1.addExon(exon5);

        transcript1.addUtr(utr1);
        transcript1.addUtr(utr2);
        transcript1.addUtr(utr3);

        gene1.addTranscript(transcript1);

        return gene1;
    }

    /**
     * Helper method to create data for the {@link DataProvider} {@link #decodeTestProvider()}
     * @return An {@link Object} array of size 2 containing a file name the object-representation of a file's data.
     */
    private Object[] createTestData_valid1() {
        return new Object[] {
            "gencode.valid1.gtf", new ArrayList<>( Collections.singletonList( createGencodeGtfGene_valid1()) )
        };
    }

    /**
     * Helper method to create data for the {@link DataProvider} {@link #decodeTestProvider()}
     * @return An {@link Object} array of size 2 containing a file name the object-representation of a file's data.
     */
    private Object[] createTestData_gencode_file2() {

        return new Object[] {
                "gencode.valid_gencode_file2.gtf",
                new ArrayList<>( Collections.singletonList( createGencodeGtfGene_file2()) )
        };
    }

    private Object[] createTestData_theOtherfile() {
        return new Object[] {
                "gencode.and.this.is.a.valid.one.too.table.gtf",
                new ArrayList<>( Collections.singletonList( createGencodeGtfGene_file3()) )
        };
    }

    private Object[] createTestData_v19_valid1() {
        return new Object[] {
                "gencode.v19.valid1.gtf",
                new ArrayList<>( Collections.singletonList( createGencodeGtfGene_v19_valid1() ) )
        };
    }

    private Object[] createTestData_v19_file2() {
        return new Object[] {
                "gencode.v19.valid_gencode_file2.gtf",
                createGencodeGtfGene_v19_file2()
        };
    }

    private Object[] createTestData_v19_theOtherFile() {
        return new Object[] {
                "gencode.v19.and.this.is.a.valid.one.too.gtf",
                new ArrayList<>( Collections.singletonList( createGencodeGtfGene_v19_theOtherFile() ) )
        };
    }

    // ============================================================================================================
    // ============================================================================================================
    // ============================================================================================================

    @DataProvider
    private Object[][] testIndexingProvider() {

        return new Object[][] {

                // HG38
                { xyzTestFile, "chrX", 105958620, 106070470, 1, },          // gene fits in region
                { xyzTestFile, "chrX", 106033442, 106034738, 1, },          // region fits in gene
                { xyzTestFile, "chrX", 106033442, 106070470, 1, },          // region start in gene
                { xyzTestFile, "chrX", 105958620, 106034738, 1, },          // region end in gene
                { xyzTestFile, "chrX", 1, 200000000, 2366, },               // Many genes in region
                { xyzTestFile, "chrX", 10000, 200000, 0, },                 // no genes in region

                // HG19
                { gencodeHg19TestFile, "chr1", 32198, 41014, 1, },          // gene fits in region
                { gencodeHg19TestFile, "chr1", 35001, 35500, 1, },          // region fits in gene
                { gencodeHg19TestFile, "chr1", 35500, 37000, 1, },          // region start in gene
                { gencodeHg19TestFile, "chr1", 32198, 35500, 1, },          // region end in gene
                { gencodeHg19TestFile, "chr1", 1, 200000000, 1982, },       // Many genes in region
                { gencodeHg19TestFile, "chr1", 33001, 34091, 0, },          // no genes in region

        };
    }

    @DataProvider
    private Object[][] toStringTestProvider() {

        // Hand-done results:
        GencodeGtfGeneFeature gene = createGencodeGtfGene_valid1();

        String expected = "chr1\tENSEMBL\tgene\t30366\t30503\t.\t+\t.\tgene_id \"ENSG00000284332.1\"; gene_type \"miRNA\"; gene_name \"MIR1302-2\"; level 3;\n" +
                "chr1\tENSEMBL\ttranscript\t30366\t30503\t.\t+\t.\tgene_id \"ENSG00000284332.1\"; transcript_id \"ENST00000607096.1\"; gene_type \"miRNA\"; gene_name \"MIR1302-2\"; transcript_type \"miRNA\"; transcript_name \"MIR1302-2-201\"; level 3; transcript_support_level \"NA\"; tag \"basic\";\n" +
                "chr1\tENSEMBL\texon\t30366\t30503\t.\t+\t.\tgene_id \"ENSG00000284332.1\"; transcript_id \"ENST00000607096.1\"; gene_type \"miRNA\"; gene_name \"MIR1302-2\"; transcript_type \"miRNA\"; transcript_name \"MIR1302-2-201\"; exon_number 1; exon_id \"ENSE00003695741.1\"; level 3; transcript_support_level \"NA\"; tag \"basic\";";

        return new Object[][] {
                {gene, expected},
        };
    }

    @DataProvider
    public Object[][] nameProvider() {

        return new Object[][] {
                { "a.tsv"     , false },                                    // Wrong File name / extension
                { "a.table.gz", false },                                    // Wrong File name / extension
                { "a.bed"     , false },                                    // Wrong File name / extension
                { "a.bcf"     , false },                                    // Wrong File name / extension
                { "a.hapmap"  , false },                                    // Wrong File name / extension
                { "a.refseq"  , false },                                    // Wrong File name / extension
                { "a.beagle"  , false },                                    // Wrong File name / extension
                { "a.table"   , false },                                    // Wrong File name / extension

                { "gencode.v26.annotation.gtf.tsv", false},                 // Wrong File name / extension
                { "gencode.v26.annotation.tgz"    , false},                 // Wrong File name / extension
                { "gencode.v26.annotation.tar.gz" , false},                 // Wrong File name / extension

                { "gencode.gtf"                                , false},    // File does not exist
                { "gencode.v26.primary_assembly.annotation.gtf", false},    // File does not exist
                { "gencode.v26.long_noncoding_RNAs.gtf"        , false},    // File does not exist

                { "gencode.invalid_short_header.gtf"           , false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header.gtf"       , false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_desc.gtf"  , false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_prov.gtf"  , false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_cont.gtf"  , false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_form.gtf"  , false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_date.gtf"  , false},    // File exists, has invalid header

                { "gencode.valid1.gtf"                           , true},   // Valid file
                { "gencode.valid_gencode_file2.gtf"              , true},   // Valid file
                { "gencode.and.this.is.a.valid.one.too.table.gtf", true},   // Valid file
        };
    }

    @DataProvider
    public Object[][] headerProvider() {
        return new Object[][] {

                { new String[] {}, false },                             // Wrong length header
                { new String[] { "",
                                 "",
                                 "",
                                 "",
                                 "" },
                                 false },                                // Bad content
                { new String[] { "##descr",
                                 "##provider: GENCODE",
                                 "##contact: gencode-help@sanger.ac.uk",
                                 "##format: gtf",
                                 "##date: 2017-04-08" },
                                 false },                                // Bad header - description
                { new String[] { "##description: THIS IS A SAMPLE",
                                 "##provider: GARBAGEDAY",
                                 "##contact: gencode-help@sanger.ac.uk",
                                 "##format: gtf",
                                 "##date: 2017-04-08" },
                                false },                                // Bad header - provider
                { new String[] { "##description: THIS IS A SAMPLE",
                                 "##provider: GENCODE",
                                 "##contact: gencode@NORTHPOLE.pl",
                                 "##format: gtf",
                                 "##date: 2017-04-08" },
                                 false },                                // Bad header - contact
                { new String[] { "##description: THIS IS A SAMPLE",
                                 "##provider: GENCODE",
                                 "##contact: SANTACLAUSE@sanger.ac.uk",
                                 "##format: gtf",
                                 "##date: 2017-04-08" },
                                 false },                                // Bad header - contact
                { new String[] { "##description: THIS IS A SAMPLE",
                                 "##provider: GENCODE",
                                 "##contact: gencode-help@sanger.ac.uk",
                                 "##format: dumpy",
                                 "##date: 2017-04-08" },
                                false },                                // Bad header - format
                { new String[] { "##description: THIS IS A SAMPLE",
                                 "##provider: GENCODE",
                                 "##contact: gencode-help@sanger.ac.uk",
                                 "##format: gtf",
                                 "##doom: ID Software" },
                                 false },                                // Bad header - date
                { new String[] { "##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)",
                                 "##provider: GENCODE",
                                 "##contact: gencode@sanger.ac.uk",
                                 "##format: gtf",
                                 "##date: 2014-07-25" },
                                 true },                                // Good Header!
                { new String[] { "##description: evidence-based annotation of the human genome (GRCh38), version 26 (Ensembl 88)",
                                 "##provider: GENCODE",
                                 "##contact: gencode-help@sanger.ac.uk",
                                 "##format: gtf",
                                 "##date: 2014-07-25" },
                                 true },                                 // Good Header!

        };
    }

    @DataProvider
    public Object[][] decodeTestProvider() {

        return new Object[][] {
                createTestData_valid1(),
                createTestData_gencode_file2(),
                createTestData_theOtherfile(),
                createTestData_v19_valid1(),
                createTestData_v19_file2(),
                createTestData_v19_theOtherFile()
        };
    }

    // =============================================================================================================

    @Test(dataProvider = "nameProvider")
    public void testCanDecode(final String name, final boolean expected) {
        GencodeGtfCodec gencodeGtfCodec = new GencodeGtfCodec();
        Assert.assertEquals(gencodeGtfCodec.canDecode(testResourceDir + name), expected, name);
    }

    @Test(dataProvider = "headerProvider")
    public void testValidateHeader( final String[] header, final boolean expected ) {
        Assert.assertEquals( GencodeGtfCodec.validateHeader(header), expected );
    }

    @Test(dataProvider = "decodeTestProvider")
    public void testDecode( final String filePath, final List<GencodeGtfFeature> expected) throws IOException {
        GencodeGtfCodec gencodeGtfCodec = new GencodeGtfCodec();

        try (BufferedInputStream bufferedInputStream =
                     new BufferedInputStream(
                             new FileInputStream(testResourceDir + filePath)
                     )
        ) {
            // Get the line iterator:
            LineIterator lineIterator = gencodeGtfCodec.makeSourceFromStream(bufferedInputStream);

            // Get the header (required for the read to work correctly):
            gencodeGtfCodec.readHeader(lineIterator);

            // Setup our expected data iterator:
            Iterator<GencodeGtfFeature> expectedIterator = expected.iterator();

            // Now read our features and make sure they're what we expect:
            while ( lineIterator.hasNext() ) {
                GencodeGtfFeature feature = gencodeGtfCodec.decode(lineIterator);
                Assert.assertEquals(feature, expectedIterator.next());
            }
        }
    }

    @Test(dataProvider = "toStringTestProvider")
    public void testToString( final GencodeGtfFeature feature, final String expected ) {
        Assert.assertEquals(feature.toString(), expected);
    }

    @Test(dataProvider = "testIndexingProvider")
    public void testIndexing( final String fileName, final String contig, final int start, final int end, final int numExpectedGenes ) {

        final File gencodeTestFile = new File(fileName);
        testIndexHelper(contig, start, end, numExpectedGenes, gencodeTestFile);
    }

    @Test(dataProvider = "testIndexingProvider")
    public void testIndexingAndIndexCreation( final String fileName,
                                              final String contig,
                                              final int start,
                                              final int end,
                                              final int numExpectedGenes ) throws IOException {

        GencodeGtfCodec codec = new GencodeGtfCodec();

        // Create a temp dir:
        final File tmpDir = createTempDir("testIndexingAndIndexCreation_" + start + "_" + end);
        tmpDir.deleteOnExit();

        // Create a copy of our index file:
        final File originalTestFile = new File(fileName);
        final File testFile = new File(tmpDir.getAbsolutePath() + File.separator + originalTestFile.getName());

        // Copy our file to the tmp dir:
        Files.copy(originalTestFile.toPath(), testFile.toPath(), REPLACE_EXISTING);

        // Create our Index:
        File indexFile = Tribble.indexFile(testFile);
        Index index = IndexFactory.createDynamicIndex(testFile, codec, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
        index.write(indexFile);

        // Make sure it works:
        testIndexHelper(contig, start, end, numExpectedGenes, testFile);
    }
}
