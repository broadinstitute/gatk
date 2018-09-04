package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gencode.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Unit test class for the {@link GencodeFuncotationFactory} class.
 * Created by jonn on 9/1/17.
 */
public class GencodeFuncotationFactoryUnitTest extends GATKBaseTest {
    // Gencode v19
    private static final String CNTN4_GENCODE_ANNOTATIONS_FILE_NAME = toolsTestDir + "funcotator/gencode.v19.CNTN4.annotation.gtf";
    private static final String CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE = toolsTestDir + "funcotator/gencode.v19.CNTN4.pc_transcripts.fasta";

    //==================================================================================================================
    // Multi-Test Static Variables:

    private static final double doubleEqualsEpsilon = 0.000001;

    private static final FeatureReader<GencodeGtfFeature> muc16FeatureReader;
    private static final FeatureReader<GencodeGtfFeature> muc16NonBasicFeatureReader;
    private static final FeatureReader<GencodeGtfFeature> pik3caFeatureReader;

    private static final ReferenceDataSource refDataSourceHg19Ch19;
    private static final ReferenceDataSource refDataSourceHg19Ch3;

    private static GencodeFuncotationFactory testMuc16SnpCreateFuncotationsFuncotationFactory;

    // Initialization of static variables:
    static {
        muc16NonBasicFeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.MUC16_GENCODE_NON_BASIC_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec() );
        muc16FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec() );
        pik3caFeatureReader = AbstractFeatureReader.getFeatureReader( FuncotatorTestConstants.PIK3CA_GENCODE_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec() );
        refDataSourceHg19Ch19 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref()) );
        refDataSourceHg19Ch3 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()) );

        // Gets cleaned up in `cleanupAfterTests()`
        // NOTE: This is initialized here to save time in testing.
        testMuc16SnpCreateFuncotationsFuncotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                new HashSet<>(),
                new LinkedHashMap<>());
    }

    //==================================================================================================================
    // Setup and Breakdown Methods with Hooks:

    @BeforeClass
    public static void setupBeforeTests() {
        System.out.println("Setting up before tests...");
    }

    @AfterClass
    public static void cleanupAfterTests() {
        System.out.println("Cleaning up after tests...");

        testMuc16SnpCreateFuncotationsFuncotationFactory.close();
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
                .filter( x -> x.getTranscriptId().equals(FuncotatorTestConstants.MUC16_TRANSCRIPT) )
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

    private Set<String> getValidTranscriptsForGene(final String expectedGeneName) {

        final Set<String> requestedTranscriptIds = new HashSet<>();
        if ( expectedGeneName.equals("PIK3CA") ) {
            requestedTranscriptIds.add( FuncotatorTestConstants.PIK3CA_TRANSCRIPT );
        }
        else if ( expectedGeneName.equals("MUC16") ) {
            requestedTranscriptIds.add( FuncotatorTestConstants.MUC16_TRANSCRIPT );
        }
        return requestedTranscriptIds;
    }

    private static GencodeFuncotation createFuncotationForTestGencodeFuncotationComparatorUnitTest(
            final String transcriptId,
            final GencodeGtfFeature.FeatureTag apprisLevel,
            final Integer locusLevel,
            final GencodeFuncotation.VariantClassification variantClassification,
            final Integer transcriptLength
            ) {

        final GencodeFuncotationBuilder builder = new GencodeFuncotationBuilder();

        builder.setAnnotationTranscript( transcriptId );
        builder.setApprisRank( apprisLevel );
        builder.setLocusLevel( locusLevel );
        builder.setVariantClassification( variantClassification );
        builder.setTranscriptLength( transcriptLength );

        return builder.build();
    }

    private static ReferenceContext referenceHelperForTestCalculateGcContent(final String sequence,
                                                                             final String contigName,
                                                                             final int refStartPos,
                                                                             final int refEndPos) {
        // Create an in-memory ReferenceContext:
        final SimpleInterval wholeReferenceInterval = new SimpleInterval( contigName, 1, sequence.length() );

        return new ReferenceContext(
                new ReferenceMemorySource(
                        new ReferenceBases(sequence.getBytes(), wholeReferenceInterval),
                        new SAMSequenceDictionary(
                                Collections.singletonList(
                                        new SAMSequenceRecord(contigName, sequence.length())
                                )
                        )
                ),
                new SimpleInterval( contigName, refStartPos, refEndPos )
        );
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideTranscriptForGetSortedCdsAndStartStopPositions() {
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

        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_0() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_1() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_2() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_3() );
        
        return l.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideMuc16SnpDataForGetVariantClassificationWithOutOfCdsData() {
        final List<Object[]> l = new ArrayList<>();

        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_0() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_1() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_2() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_3() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_NotInCDS() );

        return l.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideForCreateNonBasicFuncotations() {
        return new Object[][] {
                { 9002502, 9002502 },
                { 9002154, 9002154 },
                { 9001832, 9001832 },
                { 9000443, 9000443 },
                { 9000148, 9000148 },
                { 8999393, 8999393 },
                { 8999026, 8999026 },
                { 8998699, 8998699 },
                { 8997413, 8997413 },
                { 8997119, 8997119 },
                { 8996322, 8996322 },
                { 8995955, 8995955 },
                { 8995636, 8995636 },
                { 8994418, 8994418 },
                { 8994143, 8994143 },
                { 8993374, 8993374 },
                { 8993008, 8993008 },
                { 8987211, 8987211 },
                { 8987046, 8987046 },
                { 8982158, 8982158 },
                { 8973586, 8973586 },
                { 9006643, 9006643 },
                { 9006345, 9006345 },
                { 9005560, 9005560 },
                { 9005195, 9005195 },
                { 9004870, 9004870 },
                { 9003567, 9003567 },
                { 9003294, 9003294 },
                { 9002502, 9002502 },
                { 9002154, 9002154 },
                { 9001832, 9001832 },
                { 9000443, 9000443 },
                { 9000148, 9000148 },
                { 8999393, 8999393 },
                { 8999026, 8999026 },
                { 8998699, 8998699 },
                { 8997413, 8997413 },
                { 8997119, 8997119 },
                { 8996322, 8996322 },
                { 8995955, 8995955 },
                { 8994418, 8994418 },
                { 8994143, 8994143 },
                { 8993374, 8993374 },
                { 8993008, 8993008 },
                { 8987211, 8987211 },
                { 8987046, 8987046 },
                { 8982158, 8982158 },
                { 8979218, 8979218 },
                { 8977641, 8977641 },
                { 8976740, 8976740 },
                { 8976582, 8976582 },
                { 8976261, 8976261 },
                { 8973973, 8973973 },
                { 8973550, 8973550 },
                { 8971677, 8971677 },
                { 8969263, 8969263 },
                { 8968881, 8968881 },
                { 8966750, 8966750 },
                { 9006643, 9006643 },
                { 9006345, 9006345 },
                { 9005560, 9005560 },
                { 9005195, 9005195 },
                { 9004870, 9004870 },
                { 9003567, 9003567 },
                { 9003294, 9003294 },
                { 9002502, 9002502 },
                { 9002154, 9002154 },
                { 9001832, 9001832 },
                { 9000443, 9000443 },
                { 9000148, 9000148 },
                { 8999393, 8999393 },
                { 8999026, 8999026 },
                { 8998699, 8998699 },
                { 8997413, 8997413 },
                { 8997119, 8997119 },
                { 8996322, 8996322 },
                { 8995955, 8995955 },
                { 8995636, 8995636 },
                { 8994418, 8994418 },
                { 8994143, 8994143 },
                { 8993374, 8993374 },
                { 8993008, 8993008 },
                { 8987211, 8987211 },
                { 8987046, 8987046 },
                { 8982158, 8982158 },
                { 8979218, 8979218 },
                { 8976740, 8976740 },
                { 8976582, 8976582 },
                { 8976261, 8976261 },
                { 8973973, 8973973 },
                { 8973550, 8973550 },
                { 8971677, 8971677 },
                { 8969277, 8969277 },
                { 8968881, 8968881 },
                { 8966651, 8966651 },
                { 8962355, 8962355 },
                { 8961953, 8961953 },
                { 8959612, 8959612 },
        };
    }

    @DataProvider
    Object[][] provideDataForCreateFuncotations() {
        final List<Object[]> outList = new ArrayList<>();

        // MUC16 SNPs / DNPs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_1(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_2(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_3(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_4(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_5(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideEdgeCasesForMUC16Data_1(), FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        // MUC16 INDELs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16IndelData.provideIndelDataForMuc16(), FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), muc16FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        // PIK3CA SNPs / DNPs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caMnpData(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), pik3caFeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        // PIK3CA INDELs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), pik3caFeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        // PIK3CA Other Indels:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData2(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), pik3caFeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.PIK3CA_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        return outList.toArray(new Object[][]{{}});
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

    @DataProvider
    Object[][] provideDataForTestGencodeFuncotationComparatorUnitTest() {

        final String transcriptId1 = "ENST0000123456";
        final String transcriptId2 = "ENST0000987654";
        final Set<String> transcriptSet = new HashSet<>( Arrays.asList(transcriptId1, transcriptId2) );

        return new Object[][] {
                // ==================================================================================================
                // CanonicalGencodeFuncotationComparator
                // Transcript list:
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.CANONICAL.getComparator(new HashSet<>()),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.CANONICAL.getComparator(new HashSet<>(Collections.singletonList("TEST"))),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Locus Level:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL,  null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // IGR / non-IGR:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Variant Classification:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // Appris Annotation:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL.ordinal() - GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal()
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal() - GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL.ordinal()
                },
                // Transcript Length:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // ABC Order of Transcript ID:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId1.compareTo(transcriptId2)
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId2.compareTo(transcriptId1)
                },


                // ==================================================================================================
                // BestEffectGencodeFuncotationComparator
                // Transcript list:
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(new HashSet<>()),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(new HashSet<>(Collections.singletonList("TEST"))),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("NOT_USER_REQUESTED", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("NOT_USER_REQUESTED", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // IGR / non-IGR:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Variant Classification:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // Locus Level:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL,  null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Appris Annotation:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL.ordinal() - GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal()
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal() - GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL.ordinal()
                },
                // Transcript Length:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // ABC Order of Transcript ID:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId1.compareTo(transcriptId2)
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId2.compareTo(transcriptId1)
                }
        };
    }

    @DataProvider
    Object[][] provideDataForTestCalculateGcContent() {

                          // Base Position:
                          //
                          //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
                          //0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
                          //0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        final String seq = "AAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCC";
        final String contig = "test";

        return new Object[][] {
                // MNPs:
                { Allele.create("A", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig,100,100), 10, 10.0/11.0 },
                { Allele.create("T", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig, 50, 50),  1, 1.0/3.0 },
                { Allele.create("T", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig, 50, 50),  9, 9.0/19.0 },
                { Allele.create("T", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig, 50, 50), 10, 11.0/21.0 },
                { Allele.create("T", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig, 50, 50), 50, 50.0/100.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), referenceHelperForTestCalculateGcContent(seq, contig, 1,  5),   1, 0.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), referenceHelperForTestCalculateGcContent(seq, contig, 1,  5),   9, 4.0 / 14.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), referenceHelperForTestCalculateGcContent(seq, contig, 1,  5),  10, 5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), referenceHelperForTestCalculateGcContent(seq, contig,95,100),   1, 1 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), referenceHelperForTestCalculateGcContent(seq, contig,95,100),   9, 10.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), referenceHelperForTestCalculateGcContent(seq, contig,95,100),  10, 10.0/16.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), referenceHelperForTestCalculateGcContent(seq, contig,50, 55),   1, 6.0/8.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), referenceHelperForTestCalculateGcContent(seq, contig,50, 55),   9, 10.0/24.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), referenceHelperForTestCalculateGcContent(seq, contig,50, 55),  10, 11.0/26.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), referenceHelperForTestCalculateGcContent(seq, contig,50, 55),  50, 50.0/100.0 },

                // Insertions:
                { Allele.create("A", true), Allele.create("AG"),    referenceHelperForTestCalculateGcContent(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("AGG"),   referenceHelperForTestCalculateGcContent(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("AGGG"),  referenceHelperForTestCalculateGcContent(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("CG"),    referenceHelperForTestCalculateGcContent(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("CGG"),   referenceHelperForTestCalculateGcContent(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("CGGG"),  referenceHelperForTestCalculateGcContent(seq, contig,100,100), 10, 1 },
                { Allele.create("T", true), Allele.create("TG"),    referenceHelperForTestCalculateGcContent(seq, contig, 50, 50),  1, 1.0/2.0 },
                { Allele.create("T", true), Allele.create("TGG"),   referenceHelperForTestCalculateGcContent(seq, contig, 50, 50),  9, 9.0/18.0 },
                { Allele.create("T", true), Allele.create("TGGG"),  referenceHelperForTestCalculateGcContent(seq, contig, 50, 50), 10, 10.0/20.0 },
                { Allele.create("T", true), Allele.create("TGGGG"), referenceHelperForTestCalculateGcContent(seq, contig, 50, 50), 49, 49.0/98.0 },

                // Deletions:
                { Allele.create("AAAAA",  true), Allele.create("A"), referenceHelperForTestCalculateGcContent(seq, contig, 1,  5),  10,  5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("G"), referenceHelperForTestCalculateGcContent(seq, contig,95,100),  10, 10.0/15.0  },
                { Allele.create("TGGGGG", true), Allele.create("T"), referenceHelperForTestCalculateGcContent(seq, contig,50, 55),  10, 10.0/25.0 },
                { Allele.create("TG", true),     Allele.create("T"), referenceHelperForTestCalculateGcContent(seq, contig,50, 51),  10, 10.0/21.0 },
        };
    }

    @DataProvider
    private Object[][] provideDataForTestGetReferenceBases() {

        // NOTE: Genome positions start at 1, not 0.
        //                                                                                                                          1
        //                       0        1         2         3         4    *    5         6         7         8         9         0
        //                       1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        final String sequence = "CGTCGACGGAACAAAAGTAGACCATCCCTCTTGGTAAGTACGTCTTCATACTCTACAAATACCCATAGCACAATTCGGAGCCCAACGCCCGACGGGTCAT";
        //   REVERSE COMPLEMENT: ATGACCCGTCGGGCGTTGGGCTCCGAATTGTGCTATGGGTATTTGTAGAGTATGAAGACGTACTTACCAAGAGGGATGGTCTACTTTTGTTCCGTCGACG

        final String contig = "specialTestContig";

        final int startPos = 45;
        final int endPos = 45;

        // NOTE: Variants are always described in the + direction.
        //       The start/end positions are for the REFERENCE allele ONLY.

        return new Object[][] {
                // SNP in + direction
                {
                        Allele.create("T", true),
                        Allele.create("C"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTA"
                },
                // DNP in + direction
                {
                        Allele.create("TT", true),
                        Allele.create("CC"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos+1),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTAC"
                },
                // TNP in + direction
                {
                        Allele.create("TTC", true),
                        Allele.create("CGG"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos+2),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTACA"
                },
                // Insertion - 1 base in + direction
                {
                        Allele.create("T", true),
                        Allele.create("TC"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTAC"
                },
                // Insertion - 2 bases in + direction
                {
                        Allele.create("T", true),
                        Allele.create("TCC"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTACA"
                },
                // Insertion - 3 bases in + direction
                {
                        Allele.create("T", true),
                        Allele.create("TCCA"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.POSITIVE,
                        "TAAGTACGTCTTCATACTCTACAA"
                },
                // Deletion - 1 base in + direction
                {
                        Allele.create("CT", true),
                        Allele.create("T"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos-1,  endPos),
                        Strand.POSITIVE,
                        "GTAAGTACGTCTTCATACTCTA"
                },
                // Deletion - 2 bases in + direction
                {
                        Allele.create("CTT", true),
                        Allele.create("T"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos-1,  endPos+1),
                        Strand.POSITIVE,
                        "GTAAGTACGTCTTCATACTCTAC"
                },
                // Deletion - 3 bases in + direction
                {
                        Allele.create("CTTC", true),
                        Allele.create("T"),
                        referenceHelperForTestCalculateGcContent(sequence, contig, startPos - 1, endPos + 2),
                        Strand.POSITIVE,
                        "GTAAGTACGTCTTCATACTCTACA"
                },

                // ================================================================================

                // SNP in - direction
                {
                        Allele.create("T", true),
                        Allele.create("C"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTA"
                },
                // DNP in - direction
                {
                        Allele.create("CT", true),
                        Allele.create("TA"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos-1,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTAC"
                },
                // TNP in - direction
                {
                        Allele.create("TCT", true),
                        Allele.create("ATG"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos-2,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTACC"
                },
                // Insertion - 1 base in - direction
                {
                        Allele.create("T", true),
                        Allele.create("TC"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "GTAGAGTATGAAGACGTACTTA"
                },
                // Insertion - 2 bases in - direction
                {
                        Allele.create("T", true),
                        Allele.create("TCC"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "TGTAGAGTATGAAGACGTACTTA"
                },
                // Insertion - 3 bases in - direction
                {
                        Allele.create("T", true),
                        Allele.create("TCCA"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos,  endPos),
                        Strand.NEGATIVE,
                        "TTGTAGAGTATGAAGACGTACTTA"
                },
                // Deletion - 1 base in - direction
                {
                        Allele.create("TT", true),
                        Allele.create("T"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos-1,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTAC"
                },
                // Deletion - 2 bases in - direction
                {
                        Allele.create("CTT", true),
                        Allele.create("C"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos-2,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTACC"
                },
                // Deletion - 3 bases in - direction
                {
                        Allele.create("TCTT", true),
                        Allele.create("T"),
                        referenceHelperForTestCalculateGcContent(sequence, contig,  startPos-3,  endPos),
                        Strand.NEGATIVE,
                        "TAGAGTATGAAGACGTACTTACCA"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test ( dataProvider = "provideTranscriptForGetSortedCdsAndStartStopPositions")
    void testGetSortedExonAndStartStopPositions(final GencodeGtfTranscriptFeature transcript, final List<? extends Locatable> expected) {

        final List<? extends Locatable> exons = GencodeFuncotationFactory.getSortedCdsAndStartStopPositions(transcript);

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
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
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

        final ReferenceContext referenceContext = new ReferenceContext( refDataSourceHg19Ch19, variantInterval );

        final List<? extends Locatable> exonPositionList = GencodeFuncotationFactory.getSortedCdsAndStartStopPositions(transcript);

        final ReferenceDataSource muc16TranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE).toPath());
        final Map<String, GencodeFuncotationFactory.MappedTranscriptIdInfo> muc16TranscriptIdMap = GencodeFuncotationFactory. createTranscriptIdMap(muc16TranscriptDataSource);

        final SequenceComparison seqComp =
                GencodeFuncotationFactory.createSequenceComparison(
                        variantContext,
                        altAllele,
                        referenceContext,
                        transcript,
                        exonPositionList,
                        muc16TranscriptIdMap,
                        muc16TranscriptDataSource,
                        true);

        final GencodeFuncotation.VariantClassification varClass = GencodeFuncotationFactory.createVariantClassification(
                variantContext,
                altAllele,
                variantType,
                exon,
                transcript.getExons().size(),
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
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
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

        // TODO: Make this an input argument:
        final Set<String> requestedTranscriptIds = getValidTranscriptsForGene("MUC16");

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                requestedTranscriptIds,
                new LinkedHashMap<>())) {

            // Generate our funcotations:
            final List<Feature> featureList = new ArrayList<>();
            featureList.add( gene );
            final List<Funcotation> funcotations = funcotationFactory.createFuncotationsOnVariant(variantContext, referenceContext, featureList);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 1);

            final GencodeFuncotation funcotation = (GencodeFuncotation)funcotations.get(0);

            Assert.assertEquals(funcotation.getVariantClassification(), expectedVariantClassification);
            Assert.assertEquals(funcotation.getVariantType(), expectedVariantType);
        }
    }

    @Test(dataProvider = "provideForCreateNonBasicFuncotations")
    void createNonBasicFuncotations(final int start, final int end) {

        final String contig = "chr19";

        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create("C", true);
        final Allele altAllele = Allele.create("G", false);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = muc16NonBasicFeatureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch19, variantInterval );

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                                                                    IOUtils.getPath(FuncotatorTestConstants.MUC16_GENCODE_TRANSCRIPT_FASTA_FILE),
                                                                    "VERSION",
                                                                    GencodeFuncotationFactory.DEFAULT_NAME,
                                                                    FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                                                                    new HashSet<>(),
                                                                    new LinkedHashMap<>())) {

            // Generate our funcotations:
            final List<Feature> featureList = new ArrayList<>();
            featureList.add( gene );
            final List<Funcotation> funcotations = funcotationFactory.createFuncotationsOnVariant(variantContext, referenceContext, featureList);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 0);
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

        // TODO: Make this an input argument:
        final Set<String> requestedTranscriptIds = getValidTranscriptsForGene(expectedGeneName);

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(transcriptFastaFile),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                requestedTranscriptIds,
                new LinkedHashMap<>())) {

            final List<Feature> featureList = new ArrayList<>();
            featureList.add( gene );
            final List<Funcotation> funcotations = funcotationFactory.createFuncotationsOnVariant(variantContext, referenceContext, featureList);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 1);

            final GencodeFuncotation funcotation = (GencodeFuncotation)funcotations.get(0);

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
                errorMessageStringBuilder.append("\n\tGene Name is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedGeneName);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getHugoSymbol());
                errorMessageStringBuilder.append("]");
            }
            if (!variantClassificationCorrect) {
                errorMessageStringBuilder.append("\n\tVariant Classification is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedVariantClassification);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getVariantClassification());
                errorMessageStringBuilder.append("]");
            }
            if (!variantTypeCorrect) {
                errorMessageStringBuilder.append("\n\tVariant Type is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedVariantType);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getVariantType());
                errorMessageStringBuilder.append("]");
            }
            if (!genomeChangeCorrect) {
                errorMessageStringBuilder.append("\n\tGenome Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedGenomeChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getGenomeChange());
                errorMessageStringBuilder.append("]");
            }
            if (!strandCorrect) {
                errorMessageStringBuilder.append("\n\tStrand is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedStrand);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getTranscriptStrand());
                errorMessageStringBuilder.append("]");
            }
            if (!cDnaChangeCorrect) {
                errorMessageStringBuilder.append("\n\tCDna Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedCDnaChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getcDnaChange());
                errorMessageStringBuilder.append("]");
            }
            if (!codonChangeCorrect) {
                errorMessageStringBuilder.append("\n\tCodon Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedCodonChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getCodonChange());
                errorMessageStringBuilder.append("]");
            }
            if (!proteinChangeCorrect) {
                errorMessageStringBuilder.append("\n\tProtein Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedProteinChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getProteinChange());
                errorMessageStringBuilder.append("]");
            }
            errorMessageStringBuilder.append("\n");
            errorMessageStringBuilder.append(funcotation.getAnnotationTranscript());

            Assert.assertTrue(
                    geneNameCorrect && variantClassificationCorrect && variantTypeCorrect &&
                    genomeChangeCorrect && strandCorrect && cDnaChangeCorrect && codonChangeCorrect && proteinChangeCorrect,
                    errorMessageStringBuilder.toString() + "\n"
            );

            Assert.assertNotNull(funcotation.getMetadata(), "Metadata was null.");
            Assert.assertTrue(funcotation.getMetadata().retrieveAllHeaderInfo().size() > 0, "Metadata was empty.");

            final Set<String> metaDataFields =  funcotation.getMetadata().retrieveAllHeaderInfo().stream()
                    .map(f -> f.getID()).collect(Collectors.toSet());
            final Set<String> symmetricDifference = Sets.symmetricDifference(metaDataFields, funcotation.getFieldNames());

            Assert.assertEquals(symmetricDifference.size(), 0, "Metadata fields did not match exactly the funcotation field names: " +
                symmetricDifference.stream().collect(Collectors.joining(", ")));
        }
    }

    @Test (dataProvider = "provideDataForTestIsVariantInCodingRegion")
    void testIsVariantInCodingRegion(final GencodeFuncotation.VariantClassification varClass, final GencodeFuncotation.VariantClassification secondaryVarClass, final boolean expected) {
        Assert.assertEquals( GencodeFuncotationFactory.isVariantInCodingRegion(varClass, secondaryVarClass), expected );
    }

    @Test ( dataProvider = "provideDataForTestGencodeFuncotationComparatorUnitTest")
    void testGencodeFuncotationComparatorUnitTest( final Comparator<GencodeFuncotation> comparator,
                                                   final GencodeFuncotation a,
                                                   final GencodeFuncotation b,
                                                   final int expected) {
        Assert.assertEquals( comparator.compare(a,b), expected );
    }

    @Test (dataProvider = "provideDataForTestCalculateGcContent")
    void testCalculateGcContent(final Allele refAllele,
                                final Allele altAllele,
                                final ReferenceContext referenceContext,
                                final int windowSize,
                                final double expected) {
        Assert.assertEquals( GencodeFuncotationFactory.calculateGcContent( refAllele, altAllele, referenceContext, windowSize ), expected, doubleEqualsEpsilon);
    }

    @Test(dataProvider = "provideDataForTestGetReferenceBases")
    void testGetReferenceBases(final Allele refAllele, final Allele altAllele, final ReferenceContext reference, final Strand strand, final String expected ) {
        Assert.assertEquals( GencodeFuncotationFactory.getReferenceBases(refAllele, altAllele, reference, strand) , expected);
    }

    @DataProvider
    public Object[][] provideSortingOfUserRequestedTranscripts() {
        final String transcriptId1 = "ENST0000001.4";
        final String transcriptId2 = "ENST0000002.3";
        final String transcriptId3 = "ENST0000003.8";
        return new Object[][] {
                {       Arrays.asList(
                            createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                            createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                            createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        Collections.singleton(transcriptId2),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        Collections.singleton(transcriptId1),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 2000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 2000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGtfFeature.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),Collections.singleton(transcriptId3),
                        transcriptId2
                }
        };
    }

    /**
     * Also, tests some basic sorting independent of the user defined list.
     */
    @Test(dataProvider = "provideSortingOfUserRequestedTranscripts")
    public void testSortingOfUserRequestedTranscripts(final List<GencodeFuncotation> gencodeFuncotations, final Set<String> userRequestedTranscripts, final String gtFirstTranscript) {

        final List<TranscriptSelectionMode> transcriptSelectionModes = Arrays.asList(TranscriptSelectionMode.ALL, TranscriptSelectionMode.BEST_EFFECT, TranscriptSelectionMode.CANONICAL);
        for (final TranscriptSelectionMode transcriptSelectionMode : transcriptSelectionModes) {
            final Comparator<GencodeFuncotation> comparator = transcriptSelectionMode.getComparator(userRequestedTranscripts);
            gencodeFuncotations.sort(comparator);
            Assert.assertEquals(gencodeFuncotations.get(0).getAnnotationTranscript(), gtFirstTranscript, " Failed on " + transcriptSelectionMode.toString());
        }
    }


    /**
     * This test (of {@link GencodeFuncotationFactory#createFuncotationsOnVariant}) makes sure that if multiple gene features are detected, there is still only one transcript returned
     *  when in the BEST_EFFECT or CANONICAL selection mode.  And tests that ALL has multiple transcripts returned.
     *  The selected test data will ensure that ALL returns more than one transcript.
     * @throws IOException
     */
    @Test
    public void testMultipleGeneFeaturesOnlyProduceOneTranscript() throws IOException {
        final GencodeGtfCodec gencodeGtfCodec = new GencodeGtfCodec();
        Assert.assertTrue(gencodeGtfCodec.canDecode(CNTN4_GENCODE_ANNOTATIONS_FILE_NAME));

        final List<Feature> gencodeFeatures = new ArrayList<>();

        // Note the "chr" here to make this work.
        final SimpleInterval variantInterval = new SimpleInterval("chr3", 2944600, 2944600);
        final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch3, variantInterval);
        final VariantContext vc = new VariantContextBuilder()
                .alleles(Arrays.asList(Allele.create("T", true), Allele.create("AT", false)))
                .chr(variantInterval.getContig()).start(variantInterval.getStart()).stop(variantInterval.getEnd())
                .make();

        try (BufferedInputStream bufferedInputStream =
                     new BufferedInputStream(
                             new FileInputStream(CNTN4_GENCODE_ANNOTATIONS_FILE_NAME)
                     )
        ) {
            // Get the line iterator:
            final LineIterator lineIterator = gencodeGtfCodec.makeSourceFromStream(bufferedInputStream);

            // Get the header (required for the read to work correctly):
            gencodeGtfCodec.readHeader(lineIterator);

            while (lineIterator.hasNext()) {
                gencodeFeatures.add(gencodeGtfCodec.decode(lineIterator));
            }
            Assert.assertTrue(gencodeFeatures.size() > 1);
        }

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.CANONICAL,
                Collections.emptySet(),
                new LinkedHashMap<>())) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertEquals(gencodeFuncotations.size(), 1);
        }

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.BEST_EFFECT,
                Collections.emptySet(),
                new LinkedHashMap<>())) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertEquals(gencodeFuncotations.size(), 1);
        }

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.ALL,
                Collections.emptySet(),
                new LinkedHashMap<>())) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertTrue(gencodeFuncotations.size() > 1);
        }
    }
}
