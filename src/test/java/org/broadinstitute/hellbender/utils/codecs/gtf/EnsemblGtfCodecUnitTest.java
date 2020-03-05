package org.broadinstitute.hellbender.utils.codecs.gtf;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

/**
 * Test class for the ENSEMBL GTF Reader.
 * Modeled after the TableCodecUnitTest, with extras specific to this file format.
 * Created by jonn on 2020 02 28.
 */
public class EnsemblGtfCodecUnitTest extends GATKBaseTest {

    private static final String ECOLI_UCSC_GENOME_VERSION = "ASM584v2";
    private static final String ECOLI_CONTIG_NAME         = "Chromosome";

    private static final String testResourceDir           = publicTestDir + "org/broadinstitute/hellbender/utils/codecs/gtf/";
    private static final String eColiTestDir              = largeFileTestDir + "funcotator/ecoli_ds/gencode/" + ECOLI_UCSC_GENOME_VERSION + "/";

    private GencodeGtfFeature createEcoliEnsemblGene(final int startingFeatureOrder, final String geneId, final String geneName, final int geneStart,
                                             final int geneEnd, final String transcriptId, final String transcriptName,
                                             final String exonId, final int cdsStart, final int cdsEnd) {

        // Placeholder for constant extra data:
        final String geneAnonymousOptionalFields = " gene_source \"ena\";";
        final String transcriptAnonymousOptionalFields = geneAnonymousOptionalFields + " transcript_source \"ena\";";

        int featureOrderNum = startingFeatureOrder;

        GencodeGtfFeatureBaseData data;

        data = new GencodeGtfFeatureBaseData(EnsemblGtfCodec.GTF_FILE_TYPE_STRING, featureOrderNum++, ECOLI_CONTIG_NAME, GencodeGtfFeature.AnnotationSource.ena, GencodeGtfFeature.FeatureType.GENE,
                geneStart, geneEnd, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, geneId, null, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, geneName, null, null, null, -1, null, null,
                null,
                geneAnonymousOptionalFields);
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature)GencodeGtfFeature.create(data);
        gene.setUcscGenomeVersion(ECOLI_UCSC_GENOME_VERSION);

        data = new GencodeGtfFeatureBaseData(EnsemblGtfCodec.GTF_FILE_TYPE_STRING, featureOrderNum++, ECOLI_CONTIG_NAME, GencodeGtfFeature.AnnotationSource.ena, GencodeGtfFeature.FeatureType.TRANSCRIPT,
                geneStart, geneEnd, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, geneId, transcriptId, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, geneName, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, transcriptName, -1, null, null,
                null,
                transcriptAnonymousOptionalFields
        );
        final GencodeGtfTranscriptFeature transcript = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(data);
        transcript.setUcscGenomeVersion(ECOLI_UCSC_GENOME_VERSION);

        data = new GencodeGtfFeatureBaseData(EnsemblGtfCodec.GTF_FILE_TYPE_STRING, featureOrderNum++, ECOLI_CONTIG_NAME, GencodeGtfFeature.AnnotationSource.ena, GencodeGtfFeature.FeatureType.EXON,
                geneStart, geneEnd, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, geneId, transcriptId, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, geneName, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, transcriptName, 1, exonId, null,
                null,
                transcriptAnonymousOptionalFields
        );
        final GencodeGtfExonFeature exon = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);
        exon.setUcscGenomeVersion(ECOLI_UCSC_GENOME_VERSION);

        data = new GencodeGtfFeatureBaseData(EnsemblGtfCodec.GTF_FILE_TYPE_STRING, featureOrderNum++, ECOLI_CONTIG_NAME, GencodeGtfFeature.AnnotationSource.ena, GencodeGtfFeature.FeatureType.CDS,
                cdsStart, cdsEnd, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.ZERO, geneId, transcriptId, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, geneName, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, transcriptName, 1, null, null,
                Collections.singletonList(new GencodeGtfFeature.OptionalField<>("protein_id", transcriptId)),
                transcriptAnonymousOptionalFields
        );
        final GencodeGtfCDSFeature cds = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);
        cds.setUcscGenomeVersion(ECOLI_UCSC_GENOME_VERSION);

        data = new GencodeGtfFeatureBaseData(EnsemblGtfCodec.GTF_FILE_TYPE_STRING, featureOrderNum++, ECOLI_CONTIG_NAME, GencodeGtfFeature.AnnotationSource.ena, GencodeGtfFeature.FeatureType.START_CODON,
                cdsStart, cdsStart+2, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.ZERO, geneId, transcriptId, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, geneName, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, transcriptName, 1, null, null,
                null,
                transcriptAnonymousOptionalFields
        );
        final GencodeGtfStartCodonFeature startCodon = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(data);
        startCodon.setUcscGenomeVersion(ECOLI_UCSC_GENOME_VERSION);

        data = new GencodeGtfFeatureBaseData(EnsemblGtfCodec.GTF_FILE_TYPE_STRING, featureOrderNum, ECOLI_CONTIG_NAME, GencodeGtfFeature.AnnotationSource.ena, GencodeGtfFeature.FeatureType.STOP_CODON,
                cdsEnd+1, cdsEnd+3, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.ZERO, geneId, transcriptId, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, geneName, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, transcriptName, 1, null, null,
                null,
                transcriptAnonymousOptionalFields
        );
        final GencodeGtfStopCodonFeature stopCodon = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(data);
        stopCodon.setUcscGenomeVersion(ECOLI_UCSC_GENOME_VERSION);

        // Aggregate the Features as they should be:
        exon.setStartCodon(startCodon);
        exon.setStopCodon(stopCodon);
        exon.setCds(cds);
        transcript.addExon(exon);
        gene.addTranscript(transcript);

        return gene;
    }

    private List<GencodeGtfFeature> createEnsemblGeneSet() {
        // Create the Features as they exist in the test file:
        return Arrays.asList(
                    createEcoliEnsemblGene(6,"b0001", "thrL", 190, 255, "AAC73112", "thrL-1", "AAC73112-1",190,252),
                    createEcoliEnsemblGene(12,"b0002", "thrA", 337, 2799, "AAC73113", "thrA-1", "AAC73113-1",337,2796),
                    createEcoliEnsemblGene(18,"b0003", "thrB", 2801, 3733, "AAC73114", "thrB-1", "AAC73114-1",2801,3730),
                    createEcoliEnsemblGene(24,"b0004", "thrC", 3734, 5020, "AAC73115", "thrC-1", "AAC73115-1",3734,5017),
                    createEcoliEnsemblGene(30,"b0005", "yaaX", 5234, 5530, "AAC73116", "yaaX-1", "AAC73116-1",5234,5527)
                );
    }

    @DataProvider
    public Object[][] canDecodeProvider() {

        return new Object[][] {
                { "a.tsv"     , testResourceDir, false },                                    // Wrong File name / extension
                { "a.table.gz", testResourceDir, false },                                    // Wrong File name / extension
                { "a.bed"     , testResourceDir, false },                                    // Wrong File name / extension
                { "a.bcf"     , testResourceDir, false },                                    // Wrong File name / extension
                { "a.hapmap"  , testResourceDir, false },                                    // Wrong File name / extension
                { "a.refseq"  , testResourceDir, false },                                    // Wrong File name / extension
                { "a.beagle"  , testResourceDir, false },                                    // Wrong File name / extension
                { "a.table"   , testResourceDir, false },                                    // Wrong File name / extension

                { "gencode.v26.annotation.gtf.tsv", testResourceDir, false},                 // Wrong File name / extension
                { "gencode.v26.annotation.tgz"    , testResourceDir, false},                 // Wrong File name / extension
                { "gencode.v26.annotation.tar.gz" , testResourceDir, false},                 // Wrong File name / extension

                { "gencode.gtf"                                , testResourceDir, false},    // File does not exist
                { "gencode.v26.primary_assembly.annotation.gtf", testResourceDir, false},    // File does not exist
                { "gencode.v26.long_noncoding_RNAs.gtf"        , testResourceDir, false},    // File does not exist

                { "gencode.invalid_short_header.gtf"           , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header.gtf"       , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_desc.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_prov.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_cont.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_form.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_date.gtf"  , testResourceDir, false},    // File exists, has invalid header

                { "gencode.valid1.gtf"                           , testResourceDir, false},   // Not an Ensembl GTF file
                { "gencode.valid_gencode_file2.gtf"              , testResourceDir, false},   // Not an Ensembl GTF file
                { "gencode.and.this.is.a.valid.one.too.table.gtf", testResourceDir, false},   // Not an Ensembl GTF file

                { "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.44.gtf", eColiTestDir, true},   // Name doesn't start with 'gencode'
        };
    }

    @Test(dataProvider = "canDecodeProvider")
    public void testCanDecode(final String fileName, final String containingFolder, final boolean expected) {
        final EnsemblGtfCodec ensemblGtfCodec = new EnsemblGtfCodec();
        Assert.assertEquals(ensemblGtfCodec.canDecode(containingFolder + fileName), expected, fileName);
    }

    @DataProvider
    public Object[][] headerProvider() {
        return new Object[][] {

                { new ArrayList<String>(), false },                             // Wrong length header
                { Arrays.asList( "",
                        "",
                        "",
                        "",
                        ""  ),
                        false },                                // Bad content
                { Arrays.asList( "##descr",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - description
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GARBAGEDAY",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - provider
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: gencode@NORTHPOLE.pl",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - contact
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: SANTACLAUSE@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - contact
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: dumpy",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - format
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##doom: ID Software" ),
                        false },                                // Bad header - date
                { Arrays.asList( "##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)",
                        "##provider: GENCODE",
                        "##contact: gencode@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2014-07-25" ),
                        false },                                // Good GENCODE Header, but BAD ENSEMBL header!
                { Arrays.asList( "##description: evidence-based annotation of the human genome (GRCh38), version 26 (Ensembl 88)",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2014-07-25" ),
                        false },                                 // Good GENCODE Header, but BAD ENSEMBL header!

                // -------------

                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        true },                                           // Good Ensembl GTF Header!

                { Arrays.asList( "ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-build
                { Arrays.asList( "#!genome-build ASM584v2",
                        "ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-version
                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-date
                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-build-accession
                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#2014-08" ),
                        false },                                           // Bad header - genebuild-last-updated

                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2"),
                        false },                                           // Bad Ensembl GTF Header - not enough lines!

        };
    }

    @Test(dataProvider = "headerProvider")
    public void testValidateHeader(final List<String> header, final boolean expected ) {
        final EnsemblGtfCodec ensemblGtfCodec = new EnsemblGtfCodec();
        Assert.assertEquals( ensemblGtfCodec.validateHeader(header), expected );
    }

    @DataProvider
    public Object[][] decodeTestProvider() {

        return new Object[][] {
                {
                        eColiTestDir + "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.44.gtf",
                        createEnsemblGeneSet(),
                        ECOLI_UCSC_GENOME_VERSION
                }
        };
    }

    @Test(dataProvider = "decodeTestProvider")
    public void testDecode( final String filePath, final List<GencodeGtfFeature> expected, final String expectedUcscVersion) throws IOException {
        final EnsemblGtfCodec ensemblGtfCodec = new EnsemblGtfCodec();

        try (final BufferedInputStream bufferedInputStream =
                     new BufferedInputStream(
                             new FileInputStream(filePath)
                     )
        ) {
            // Get the line iterator:
            final LineIterator lineIterator = ensemblGtfCodec.makeSourceFromStream(bufferedInputStream);

            // Get the header (required for the read to work correctly):
            ensemblGtfCodec.readHeader(lineIterator);

            // Setup our expected data iterator:
            final Iterator<GencodeGtfFeature> expectedIterator = expected.iterator();

            // Now read our features and make sure they're what we expect:
            // NOTE: We only decode the number of features expect to see.
            int numDecoded = 0;
            while ( lineIterator.hasNext() && (numDecoded < expected.size()) ) {
                final GencodeGtfFeature feature = ensemblGtfCodec.decode(lineIterator);

                Assert.assertTrue(expectedIterator.hasNext());

                for ( final GencodeGtfFeature subFeature : feature.getAllFeatures() ) {
                    Assert.assertEquals(subFeature.getUcscGenomeVersion(), expectedUcscVersion);
                }
                final GencodeGtfFeature expectedFeature = expectedIterator.next();

                // Big equals check:
                Assert.assertEquals(feature, expectedFeature);

                ++numDecoded;
            }
        }
    }

}
