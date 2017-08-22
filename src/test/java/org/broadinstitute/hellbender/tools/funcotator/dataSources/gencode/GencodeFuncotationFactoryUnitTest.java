package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.GENCODE.*;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Unit test class for the {@link GencodeFuncotationFactory} class.
 * Created by jonn on 9/1/17.
 */
public class GencodeFuncotationFactoryUnitTest extends BaseTest {

    //==================================================================================================================
    // Static Variables:

    private static final FeatureReader<GencodeGtfFeature> muc16FeatureReader;
    private static final FeatureReader<GencodeGtfFeature> pik3caFeatureReader;

    private static final ReferenceDataSource refDataSourceHg19Ch19;
    private static final ReferenceDataSource refDataSourceHg19Ch3;

    // Initialization of static variables:
    static {
        muc16FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec() );
        pik3caFeatureReader = AbstractFeatureReader.getFeatureReader( FuncotatorTestConstants.PIK3CA_GENCODE_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec() );
        refDataSourceHg19Ch19 = ReferenceDataSource.of( new File (FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME) );
        refDataSourceHg19Ch3 = ReferenceDataSource.of( new File (FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME) );
    }

    //==================================================================================================================
    // Helper Methods:

    private static GencodeGtfExonFeature getExonForVariant( final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator,
                                                            final SimpleInterval variantLocation ) {

        final Optional<GencodeGtfExonFeature> exonOption = gtfFeatureIterator.stream()
                                                    .map( f -> ((GencodeGtfGeneFeature) f))
                                                    .flatMap(g -> g.getTranscripts().stream())
                                                    .flatMap(t -> t.getExons().stream())
                                                    .filter( x -> x.overlaps(variantLocation) )
                                                    .findFirst();

        return exonOption.orElseThrow( () -> new GATKException("Could not get an exon associated with this variant's interval: " + variantLocation.toString()) );
    }

    private static GencodeGtfExonFeature getExonForVariant( final GencodeGtfGeneFeature gene,
                                                            final SimpleInterval variantLocation ) {

        final Optional<GencodeGtfExonFeature> exonOption = gene.getTranscripts().stream()
                                                                .flatMap(t -> t.getExons().stream())
                                                                .filter( x -> x.overlaps(variantLocation) )
                                                                .findFirst();

        return exonOption.orElseThrow( () -> new GATKException("Could not get an exon associated with this variant's interval: " + variantLocation.toString()) );
    }

    private static GencodeGtfTranscriptFeature getMuc16Transcript(final GencodeGtfGeneFeature gene) {
        final Optional<GencodeGtfTranscriptFeature> transcriptOption = gene.getTranscripts().stream()
                .filter( x -> x.getTranscriptId().equals(FuncotatorTestConstants.MUC_16_TRANSCRIPT) )
                .findFirst();

        return transcriptOption.orElseThrow( () -> new GATKException("Could not get the MUC16 transcript from the MUC16 gene!  The test has mutated!  No good!") );
    }

    private static List<Object[]> addReferenceDataToUnitTestData(final List<Object[]> unitTestData,
                                                                 final String referenceFileName,
                                                                 final FeatureReader<GencodeGtfFeature> featureReader,
                                                                 final ReferenceDataSource referenceDataSource,
                                                                 final String transcriptFastaFile) {

        final List<Object[]> outList = new ArrayList<>(unitTestData.size());

        for ( final Object[] rawData : unitTestData ) {
            final Object[] dataWithReference = new Object[rawData.length + 4];
            for ( int i = 0; i < rawData.length; ++i ) {
                dataWithReference[i] = rawData[i];
            }
            dataWithReference[dataWithReference.length-4] = referenceFileName;
            dataWithReference[dataWithReference.length-3] = featureReader;
            dataWithReference[dataWithReference.length-2] = referenceDataSource;
            dataWithReference[dataWithReference.length-1] = transcriptFastaFile;
            outList.add(dataWithReference);
        }

        return outList;
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideTranscriptForGetSortedExonAndStartStopPositions() {
        return new Object[][] {
                {
                        DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature().getTranscripts().get(0),
                        Arrays.asList(
                                new SimpleInterval("chr1", 99, 101),
                                new SimpleInterval("chr1", 201, 400),
                                new SimpleInterval("chr1", 401, 600),
                                new SimpleInterval("chr1", 601, 800),
                                new SimpleInterval("chr1", 900, 902)
                        )
                },
                {
                        DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature().getTranscripts().get(1),
                        Arrays.asList(
                                new SimpleInterval("chr1", 1099, 1200),
                                new SimpleInterval("chr1", 1201, 1400),
                                new SimpleInterval("chr1", 1401, 1600),
                                new SimpleInterval("chr1", 1601, 1800),
                                new SimpleInterval("chr1", 1801, 1899),
                                new SimpleInterval("chr1", 1900, 1902)
                        )
                },
                {
                        DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature().getTranscripts().get(2),
                        Collections.emptyList()
                },
        };
    }

    @DataProvider
    Object[][] provideMuc16SnpDataForGetVariantClassification() {
        final List<Object[]> l = new ArrayList<>();

        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_1() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_2() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_3() );
        
        return l.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideMuc16SnpDataForGetVariantClassificationWithOutOfCdsData() {
        final List<Object[]> l = new ArrayList<>();

        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_1() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_2() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_3() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_NotInCDS() );

        return l.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideDataForCreateFuncotations() {
        final List<Object[]> outList = new ArrayList<>();

        // MUC13 SNPs / DNPs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_1(), FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_2(), FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_3(), FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_4(), FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_5(), FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideEdgeCasesForMUC16Data_1(), FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME, muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        // PIK3CA SNPs / DNPs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caMnpData(), FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, pik3caFeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        // PIK3CA INDELs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData(), FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, pik3caFeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        // PIK3CA Other Indels:
//        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData2(), FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME, pik3caFeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        final Object[][] outArray = outList.toArray(new Object[][]{{}});
        return outArray;
    }

    @DataProvider
    Object[][] provideDataForTestIsVariantInCodingRegion() {
        return new Object[][] {

                // Trivially false cases:
                {GencodeFuncotation.VariantClassification.INTRON, null, false},
                {GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, null, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, null, false},
                {GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.IGR, null, false},
                {GencodeFuncotation.VariantClassification.IGR, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.IGR, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, null, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, null, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, null, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.RNA, null, false},
                {GencodeFuncotation.VariantClassification.RNA, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.RNA, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.LINCRNA, null, false},
                {GencodeFuncotation.VariantClassification.LINCRNA, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.LINCRNA, GencodeFuncotation.VariantClassification.MISSENSE, false},

                // Trivially true cases:
                {GencodeFuncotation.VariantClassification.MISSENSE, null, true},
                {GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.NONSENSE, null, true},
                {GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.NONSTOP, null, true},
                {GencodeFuncotation.VariantClassification.NONSTOP, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.NONSTOP, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.SILENT, null, true},
                {GencodeFuncotation.VariantClassification.SILENT, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.SILENT, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_DEL, null, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_DEL, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_DEL, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_INS, null, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_INS, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_INS, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, null, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, null, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.START_CODON_SNP, null, true},
                {GencodeFuncotation.VariantClassification.START_CODON_SNP, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.START_CODON_SNP, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.START_CODON_INS, null, true},
                {GencodeFuncotation.VariantClassification.START_CODON_INS, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.START_CODON_INS, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.START_CODON_DEL, null, true},
                {GencodeFuncotation.VariantClassification.START_CODON_DEL, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.START_CODON_DEL, GencodeFuncotation.VariantClassification.MISSENSE, true},

                // Splice site special cases:
                {GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantClassification.MISSENSE, true},
        };
    }

    //==================================================================================================================
    // Tests:

    @Test ( dataProvider = "provideTranscriptForGetSortedExonAndStartStopPositions")
    void testGetSortedExonAndStartStopPositions(final GencodeGtfTranscriptFeature transcript, final List<? extends Locatable> expected) {

        final List<? extends Locatable> exons = GencodeFuncotationFactory.getSortedExonAndStartStopPositions(transcript);

        Assert.assertEquals(exons.size(), expected.size());

        for( int i = 0; i < exons.size() ; ++i ) {
            final SimpleInterval interval = new SimpleInterval( exons.get(i).getContig(), exons.get(i).getStart(), exons.get(i).getEnd());
            Assert.assertEquals( interval, expected.get(i) );
        }
    }

    @Test (dataProvider = "provideMuc16SnpDataForGetVariantClassification")
    void testGetVariantClassificationForCodingRegions(final int chromosomeNumber,
                                      final int start,
                                      final int end,
                                      final GencodeFuncotation.VariantType variantType,
                                      final String ref,
                                      final String alt,
                                      final GencodeFuncotation.VariantClassification expectedVariantClassification) {

        // This test can only deal with variants in coding regions.
        // So we ignore any tests that are expected outside of coding regions.
        // i.e. expectedVariantClassification is one of:
        //     { INTRON, FIVE_PRIME_UTR, THREE_PRIME_UTR, IGR, FIVE_PRIME_FLANK, DE_NOVO_START_IN_FRAME, DE_NOVO_START_OUT_FRAME, RNA, LINCRNA }
        // We test these cases in another unit test.
        if ((expectedVariantClassification == GencodeFuncotation.VariantClassification.INTRON) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.THREE_PRIME_UTR) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.IGR) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.RNA) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.LINCRNA) )
        {
            return;
        }

        final String contig = "chr" + Integer.toString(chromosomeNumber);
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME,
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = muc16FeatureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature             gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final GencodeGtfTranscriptFeature transcript = getMuc16Transcript(gene);
        final GencodeGtfExonFeature             exon = getExonForVariant( gene, variantInterval );

        final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch19, variantInterval );

        final List<? extends Locatable> exonPositionList = GencodeFuncotationFactory.getSortedExonAndStartStopPositions(transcript);

        final ReferenceDataSource muc16TranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE));
        final Map<String, GencodeFuncotationFactory.MappedTranscriptIdInfo> muc16TranscriptIdMap = GencodeFuncotationFactory.createTranscriptIdMap(muc16TranscriptDataSource);

        final FuncotatorUtils.SequenceComparison seqComp =
                GencodeFuncotationFactory.createSequenceComparison(
                        variantContext,
                        altAllele,
                        referenceContext,
                        transcript,
                        exonPositionList,
                        muc16TranscriptIdMap,
                        muc16TranscriptDataSource);

        final GencodeFuncotation.VariantClassification varClass = GencodeFuncotationFactory.getVariantClassification(
                variantContext,
                altAllele,
                variantType,
                exon,
                seqComp
        );

        Assert.assertEquals( varClass, expectedVariantClassification );
    }

    @Test (dataProvider = "provideMuc16SnpDataForGetVariantClassificationWithOutOfCdsData")
    void testMuc16SnpCreateFuncotations(final int chromosomeNumber,
                                        final int start,
                                        final int end,
                                        final GencodeFuncotation.VariantType expectedVariantType,
                                        final String ref,
                                        final String alt,
                                        final GencodeFuncotation.VariantClassification expectedVariantClassification) {

        final String contig = "chr" + Integer.toString(chromosomeNumber);
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME,
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = muc16FeatureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch19, variantInterval );

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(new File(FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE))) {
            // Generate our funcotations:
            final List<GencodeFuncotation> funcotations = funcotationFactory.createFuncotations(variantContext, altAllele, gene, referenceContext);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 1);
            Assert.assertEquals(funcotations.get(0).getVariantClassification(), expectedVariantClassification);
            Assert.assertEquals(funcotations.get(0).getVariantType(), expectedVariantType);
        }
    }

    @Test (dataProvider = "provideDataForCreateFuncotations")
    void testCreateFuncotations(final String expectedGeneName,
                                final int chromosomeNumber,
                                final int start,
                                final int end,
                                final GencodeFuncotation.VariantClassification expectedVariantClassification,
                                final GencodeFuncotation.VariantType expectedVariantType,
                                final String ref,
                                final String alt,
                                final String expectedGenomeChange,
                                final String expectedStrand,
                                final String expectedCDnaChange,
                                final String expectedCodonChange,
                                final String expectedProteinChange,
                                final String referenceFileName,
                                final FeatureReader<GencodeGtfFeature> featureReader,
                                final ReferenceDataSource referenceDataSource,
                                final String transcriptFastaFile) {

        final String contig = "chr" + Integer.toString(chromosomeNumber);
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                referenceFileName,
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = featureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final ReferenceContext referenceContext = new ReferenceContext(referenceDataSource, variantInterval );

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(new File(transcriptFastaFile))) {

            // Generate our funcotations:
            final List<GencodeFuncotation> funcotations = funcotationFactory.createFuncotations(variantContext, altAllele, gene, referenceContext);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 1);

            final GencodeFuncotation funcotation = funcotations.get(0);

            final boolean geneNameCorrect              = Objects.equals( funcotation.getHugoSymbol(), expectedGeneName );
            final boolean variantClassificationCorrect = Objects.equals( funcotation.getVariantClassification(), expectedVariantClassification );
            final boolean variantTypeCorrect           = Objects.equals( funcotation.getVariantType(), expectedVariantType );
            final boolean genomeChangeCorrect          = Objects.equals( funcotation.getGenomeChange(), expectedGenomeChange );
            final boolean strandCorrect                = Objects.equals( funcotation.getTranscriptStrand(), expectedStrand );
            final boolean cDnaChangeCorrect            = Objects.equals( funcotation.getcDnaChange(), expectedCDnaChange );
            final boolean codonChangeCorrect           = Objects.equals( funcotation.getCodonChange(), expectedCodonChange );
            final boolean proteinChangeCorrect         = Objects.equals( funcotation.getProteinChange(), expectedProteinChange );

            final StringBuilder errorMessageStringBuilder = new StringBuilder();

            if (!geneNameCorrect) {
                errorMessageStringBuilder.append("\n\tGene Name is not correct!\n\t\tExpected:  [" + expectedGeneName + "]\n\t\tBut found: [" + funcotation.getHugoSymbol() + "]");
            }
            if (!variantClassificationCorrect) {
                errorMessageStringBuilder.append("\n\tVariant Classification is not correct!\n\t\tExpected:  [" + expectedVariantClassification + "]\n\t\tBut found: [" + funcotation.getVariantClassification() + "]");
            }
            if (!variantTypeCorrect) {
                errorMessageStringBuilder.append("\n\tVariant Type is not correct!\n\t\tExpected:  [" + expectedVariantType + "]\n\t\tBut found: [" + funcotation.getVariantType() + "]");
            }
            if (!genomeChangeCorrect) {
                errorMessageStringBuilder.append("\n\tGenome Change is not correct!\n\t\tExpected:  [" + expectedGenomeChange + "]\n\t\tBut found: [" + funcotation.getGenomeChange() + "]");
            }
            if (!strandCorrect) {
                errorMessageStringBuilder.append("\n\tStrand is not correct!\n\t\tExpected:  [" + expectedStrand + "]\n\t\tBut found: [" + funcotation.getTranscriptStrand() + "]");
            }
            if (!cDnaChangeCorrect) {
                errorMessageStringBuilder.append("\n\tCDna Change is not correct!\n\t\tExpected:  [" + expectedCDnaChange + "]\n\t\tBut found: [" + funcotation.getcDnaChange() + "]");
            }
            if (!codonChangeCorrect) {
                errorMessageStringBuilder.append("\n\tCodon Change is not correct!\n\t\tExpected:  [" + expectedCodonChange + "]\n\t\tBut found: [" + funcotation.getCodonChange() + "]");
            }
            if (!proteinChangeCorrect) {
                errorMessageStringBuilder.append("\n\tProtein Change is not correct!\n\t\tExpected:  [" + expectedProteinChange + "]\n\t\tBut found: [" + funcotation.getProteinChange() + "]");
            }

            Assert.assertTrue(
                    geneNameCorrect && variantClassificationCorrect && variantTypeCorrect &&
                    genomeChangeCorrect && strandCorrect && cDnaChangeCorrect && codonChangeCorrect && proteinChangeCorrect,
                    errorMessageStringBuilder.toString() + "\n"
            );
            
        }
    }

    @Test (dataProvider = "provideDataForTestIsVariantInCodingRegion")
    void testIsVariantInCodingRegion(final GencodeFuncotation.VariantClassification varClass, final GencodeFuncotation.VariantClassification secondaryVarClass, final boolean expected) {
        Assert.assertEquals( GencodeFuncotationFactory.isVariantInCodingRegion(varClass, secondaryVarClass), expected );
    }
}
