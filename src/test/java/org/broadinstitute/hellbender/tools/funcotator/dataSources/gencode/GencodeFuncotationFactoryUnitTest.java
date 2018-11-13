package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import com.google.common.collect.Sets;
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
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gencode.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
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

    private static final FeatureReader<GencodeGtfFeature> gencodeHg19FeatureReader;
    private static final FeatureReader<GencodeGtfFeature> gencodeHg38FeatureReader;
    private static final FeatureReader<GencodeGtfFeature> muc16NonBasicFeatureReader;

    private static final ReferenceDataSource refDataSourceHg19Ch19;
    private static final ReferenceDataSource refDataSourceHg19Ch3;
    private static final ReferenceDataSource refDataSourceHg38;

    private static GencodeFuncotationFactory testMuc16SnpCreateFuncotationsFuncotationFactory;

    private static final List<AutoCloseable> autoCloseableList = new ArrayList<>();

    // Initialization of static variables:
    static {
        muc16NonBasicFeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.MUC16_GENCODE_NON_BASIC_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec());
        gencodeHg19FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19, new GencodeGtfCodec());
        gencodeHg38FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG38, new GencodeGtfCodec());
        refDataSourceHg19Ch19 = ReferenceDataSource.of(IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref()));
        refDataSourceHg19Ch3 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()) );
        refDataSourceHg38 = ReferenceDataSource.of( IOUtils.getPath(hg38Reference) );

        // Gets cleaned up in `cleanupAfterTests()`
        // NOTE: This is initialized here to save time in testing.
        testMuc16SnpCreateFuncotationsFuncotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                new HashSet<>(),
                new LinkedHashMap<>(),
                createFeatureInputForMuc16Ds(GencodeFuncotationFactory.DEFAULT_NAME));

        // Add all to the closeable list:
        autoCloseableList.add( muc16NonBasicFeatureReader );
        autoCloseableList.add( gencodeHg19FeatureReader );
        autoCloseableList.add( gencodeHg38FeatureReader );
        autoCloseableList.add( refDataSourceHg19Ch19 );
        autoCloseableList.add( refDataSourceHg19Ch3 );
        autoCloseableList.add( refDataSourceHg38 );
        autoCloseableList.add( testMuc16SnpCreateFuncotationsFuncotationFactory );
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

        for ( final AutoCloseable cl : autoCloseableList ) {

            try {
                cl.close();
            }
            catch ( final Exception ex ) {
                throw new GATKException("Could not close " + cl.toString(), ex);
            }
        }
    }

    //==================================================================================================================
    // Helper Methods:

    private static Object[] helpProvideForTestGetTranscriptEndPaddingBases( final String refAllele,
                                                                            final String altAllele,
                                                                            final int startPos,
                                                                            final List<? extends Locatable> exonPositionList,
                                                                            final String refSequence,
                                                                            final String expected) {
        final String contig = "TEST";

        final VariantContext vc = new VariantContextBuilder(
                "TEST_SOURCE",
                contig,
                startPos,
                startPos + refAllele.length() - 1,
                Arrays.asList(Allele.create(refAllele, true), Allele.create(altAllele, false))
        ).make();

        return new Object[] {
            vc,
            vc.getAlleles().get(1),
            exonPositionList,
            FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(refSequence, vc.getContig(),  vc.getEnd(),  vc.getEnd()),
            expected
        };
    }

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
                                                                 final String transcriptFastaFile,
                                                                 final String transcriptGtfFile) {

        final List<Object[]> outList = new ArrayList<>(unitTestData.size());

        for ( final Object[] rawData : unitTestData ) {
            final Object[] dataWithReference = new Object[rawData.length + 5];
            for ( int i = 0; i < rawData.length; ++i ) {
                dataWithReference[i] = rawData[i];
            }
            dataWithReference[dataWithReference.length-5] = referenceFileName;
            dataWithReference[dataWithReference.length-4] = featureReader;
            dataWithReference[dataWithReference.length-3] = referenceDataSource;
            dataWithReference[dataWithReference.length-2] = transcriptFastaFile;
            dataWithReference[dataWithReference.length-1] = transcriptGtfFile;
            outList.add(dataWithReference);
        }

        return outList;
    }

    private Set<String> getValidTranscriptsForGene(final String expectedGeneName) {

        final Set<String> requestedTranscriptIds = new HashSet<>();
        if ( expectedGeneName != null ){
            if ( expectedGeneName.equals("PIK3CA") ) {
                requestedTranscriptIds.add( FuncotatorTestConstants.PIK3CA_TRANSCRIPT );
            }
            else if ( expectedGeneName.equals("MUC16") ) {
                requestedTranscriptIds.add( FuncotatorTestConstants.MUC16_TRANSCRIPT );
            }
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

    private GencodeGtfFeatureBaseData createGtfBaseDataForTestIs5Prime(final SimpleInterval interval) {
        return new GencodeGtfFeatureBaseData(
                1,
                interval.getContig(),
                GencodeGtfFeature.AnnotationSource.ENSEMBL,
                GencodeGtfFeature.FeatureType.GENE,
                interval.getStart(),
                interval.getEnd(),
                Strand.POSITIVE,
                GencodeGtfFeature.GenomicPhase.DOT,
                "TEST-GENE-ID",
                "TEST-TX-ID",
                GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                GencodeGtfFeature.GeneTranscriptStatus.PUTATIVE,
                "TEST-GENE",
                GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                GencodeGtfFeature.GeneTranscriptStatus.PUTATIVE,
                "TEST-TX",
                1,
                "",
                GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                null,
                "");
    }

    private GencodeGtfExonFeature helpCreateExonFeature(final SimpleInterval interval,
                                                        final GencodeGtfStartCodonFeature startCodonFeature) {

        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.EXON;
        final GencodeGtfExonFeature exon = (GencodeGtfExonFeature)GencodeGtfExonFeature.create(baseData);

        if ( startCodonFeature != null ) {
            exon.setStartCodon(startCodonFeature);
        }

        return exon;
    }

    private GencodeGtfTranscriptFeature provideForTestIs5PrimeUtrTranscriptHelper(
                                            final SimpleInterval interval,
                                            final Strand strand,
                                            final GencodeGtfStartCodonFeature startCodonFeature) {
        //    All that is needed for the Transcript is:
        //      Strand
        //      Exons
        //          start codon within an exon

        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.TRANSCRIPT;
        baseData.genomicStrand = strand;

        final GencodeGtfTranscriptFeature transcriptFeature = (GencodeGtfTranscriptFeature) GencodeGtfTranscriptFeature.create(baseData);

        if ( startCodonFeature == null ) {
            // Always create a list of 3 exons:
            final int NUM_EXONS = 3;

            final int step = (interval.getEnd()-interval.getStart())/NUM_EXONS;
            for (int i = interval.getStart() ; i < interval.getEnd() ; i += step ) {
                transcriptFeature.addExon(
                        helpCreateExonFeature(new SimpleInterval(interval.getContig(), i, i+step-1), null)
                );
            }
        }
        else {
            // Create 2 exons - 1 to contain the start codon, the other to be a placeholder:
            final int NUM_EXONS = 2;

            boolean firstExon = true;
            final int step = (interval.getEnd()-interval.getStart())/NUM_EXONS;
            for (int i = interval.getStart() ; i < interval.getEnd() ; i += step ) {
                GencodeGtfStartCodonFeature startCodonLoopFeature = null;
                if (firstExon) {
                    startCodonLoopFeature = startCodonFeature;
                    firstExon = false;
                }

                transcriptFeature.addExon(
                        helpCreateExonFeature(new SimpleInterval(interval.getContig(), i, i+step-1), startCodonLoopFeature)
                );
            }
        }

        return transcriptFeature;
    }

    private GencodeGtfUTRFeature provideForTestIs5PrimeUtrUTRHelper(final SimpleInterval interval) {
        //    All that is needed for the UTR is:
        //      start
        //      end
        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.UTR;
        return (GencodeGtfUTRFeature) GencodeGtfUTRFeature.create(baseData);
    }

    private GencodeGtfStartCodonFeature provideForTestIs5PrimeUtrStartCodonHelper(final SimpleInterval interval) {
        //    All that is needed for the StartCodon is:
        //      start
        //      end
        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.START_CODON;
        return (GencodeGtfStartCodonFeature)GencodeGtfStartCodonFeature.create(baseData);
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
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_1(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_2(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_3(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_4(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_5(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideEdgeCasesForMUC16Data_1(), FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );

        // MUC16 INDELs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16IndelData.provideIndelDataForMuc16(), FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // PIK3CA SNPs / DNPs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caMnpData(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // PIK3CA INDELs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // PIK3CA Other Indels:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData2(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // Trouble variants:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForTroubleVariants.provideSymbolicAllelesAndMaskedBasesForHg38(), hg38Reference, gencodeHg38FeatureReader, refDataSourceHg38, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG38, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG38 ) );

        return outList.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideForTestGetTranscriptEndPaddingBases() {

                          // Base Position:
                          //
                          //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
                          //0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
                          //0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        final String seq = "AAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCC";
        final String contig = "test";

        final List<Locatable> exonPositionList = new ArrayList<>();
        exonPositionList.add( new SimpleInterval(contig, 5,15) );
        exonPositionList.add( new SimpleInterval(contig, 17,25) );
        exonPositionList.add( new SimpleInterval(contig, 35,55) );
        exonPositionList.add( new SimpleInterval(contig, 56,90) );

//        final VariantContext variant,
//        final Allele altAllele,
//        final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
//        final ReferenceContext reference,
//        final String expected

        return new Object[][] {
                // Not Indel:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "G", 5, exonPositionList, seq, ""),
                helpProvideForTestGetTranscriptEndPaddingBases(seq.substring(seq.length()-1, seq.length()), "A", seq.length(), exonPositionList, seq, ""),
                // Insertion 1 base, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AG", 5, exonPositionList, seq, ""),
                // Insertion 1 base, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CG", 80, exonPositionList, seq, ""),
                // Insertion 1 base, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CG", 82, exonPositionList, seq, "GGGGGC"),
                // Insertion 2 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AGC", 5, exonPositionList, seq, ""),
                // Insertion 2 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGC", 80, exonPositionList, seq, ""),
                // Insertion 2 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGA", 82, exonPositionList, seq, "GGGGGC"),
                // Insertion 3 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AGCG", 5, exonPositionList, seq, ""),
                // Insertion 3 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGCT", 80, exonPositionList, seq, ""),
                // Insertion 3 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGAT", 82, exonPositionList, seq, "GGGGGC"),
                // Insertion 4 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AGCGT", 5, exonPositionList, seq, ""),
                // Insertion 4 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGCTG", 80, exonPositionList, seq, ""),
                // Insertion 4 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGATT", 82, exonPositionList, seq, "GGGGGCCCC"),
                // Deletion 1 base, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("AT", "A", 5, exonPositionList, seq, ""),
                // Deletion 1 base, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CA", "C", 80, exonPositionList, seq, ""),
                // Deletion 1 base, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AA", "A", 82, exonPositionList, seq, "GGGGGC"),
                // Deletion 2 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("ATT", "A", 5, exonPositionList, seq, ""),
                // Deletion 2 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CAA", "C", 80, exonPositionList, seq, ""),
                // Deletion 2 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AAA", "A", 82, exonPositionList, seq, "GGGGGC"),
                // Deletion 3 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("ATTT", "A", 5, exonPositionList, seq, ""),
                // Deletion 3 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CAAA", "C", 80, exonPositionList, seq, ""),
                // Deletion 3 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AAAA", "A", 82, exonPositionList, seq, "GGGGGC"),
                // Deletion 3 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("ATTTT", "A", 5, exonPositionList, seq, ""),
                // Deletion 3 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CAAAA", "C", 80, exonPositionList, seq, ""),
                // Deletion 3 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AAAAT", "A", 82, exonPositionList, seq, "GGGGGCCCC"),
        };
    }

    @DataProvider
    Object[][] provideForTestCreateGencodeFuncotationBuilderWithTrivialFieldsPopulated() {

        //Line 1229 or so

        // Just simple field checks are OK.
        // Doesn't need to be fancy.

//        final VariantContext variant,
//        final Allele altAllele,
//        final GencodeGtfGeneFeature gtfFeature,
//        final GencodeGtfTranscriptFeature transcript

        final SimpleInterval variantInterval =  new SimpleInterval("chr3", 178921515, 178921517);
        final Allele refAllele = Allele.create("GCA", true);
        final Allele altAllele = Allele.create("TTG");
        final VariantContext variant = new VariantContextBuilder(
                    FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                    variantInterval.getContig(),
                    variantInterval.getStart(),
                    variantInterval.getEnd(),
                    Arrays.asList(refAllele, altAllele)
                ).make();


        final GencodeGtfGeneFeature gene = DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature();

        return new Object[][] {
                { variant, altAllele, gene, gene.getTranscripts().get(0) }
        };
    }

    @DataProvider
    Object[][] provideForCreateDefaultFuncotationsOnProblemVariant() {

//        variant, altAllele, gtfFeature, reference, transcript, version

        final SimpleInterval variantInterval =  new SimpleInterval("chr3", 178921515, 178921517);
        final Allele refAllele = Allele.create("GCA", true);
        final Allele altAllele = Allele.create("TTG");
        final VariantContext variant = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                variantInterval.getContig(),
                variantInterval.getStart(),
                variantInterval.getEnd(),
                Arrays.asList(refAllele, altAllele)
        ).make();

        final String versionString = "VERSION";

        // ======================
        // Create the GencodeGtfFeature:
        GencodeGtfFeatureBaseData data;

        data = new GencodeGtfFeatureBaseData(1, variantInterval.getContig(), GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.GENE,
                variantInterval.getStart()-2000, variantInterval.getEnd()+2000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", null, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", null, null, null, -1, null, GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED, null, null);
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature)GencodeGtfFeature.create(data);

        // ======================

        data = new GencodeGtfFeatureBaseData(2, variantInterval.getContig(), GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.TRANSCRIPT,
                variantInterval.getStart()-1000, variantInterval.getEnd()+1000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", -1, null, GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfTranscriptFeature transcript1 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(data);
        gene.addTranscript(transcript1);

        // ======================

        return new Object[][] {
                {
                    variant,
                    variant.getAlleles().get(1),
                    gene,
                    new ReferenceContext( refDataSourceHg19Ch3, variantInterval ),
                    gene.getTranscripts().get(0),
                    versionString,
                    new GencodeFuncotationBuilder()
                            .setHugoSymbol(gene.getGeneName())
                            .setChromosome(variant.getContig())
                            .setStart(variant.getStart())
                            .setEnd(variant.getEnd())
                            .setVariantClassification(GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE)
                            .setVariantType(GencodeFuncotation.VariantType.TNP)
                            .setRefAllele(variant.getReference())
                            .setTumorSeqAllele2(variant.getAlternateAllele(0).getBaseString())
                            .setGenomeChange("g.chr3:178921515_178921517GCA>TTG")
                            .setAnnotationTranscript(transcript1.getTranscriptName())
                            .setStrand(gene.getGenomicStrand())
                            .setReferenceContext("TATAAATAGTGCACTCAGAATAA")
                            .setGcContent(0.3399503722084367)
                            .setLocusLevel(Integer.valueOf(gene.getLocusLevel().toString()))
                            // This is OK because there are no exons:
                            .setTranscriptLength(0)
                            .setVersion(versionString)
                            .setGeneTranscriptType(transcript1.getTranscriptType())
                            .build()
                },
        };
    }

    @DataProvider
    Object[][] provideForTestIs5PrimeUtr() {

        // Need:
        //    Transcript:
        //      Strand
        //      Exons
        //          start codon within an exon
        //
        //    StartCodon:
        //      start
        //      end
        //
        //    UTR:
        //      start
        //      end

        final SimpleInterval transcriptInterval = new SimpleInterval("T", 1, 2000);

        final GencodeGtfStartCodonFeature startCodonFeatureFront = provideForTestIs5PrimeUtrStartCodonHelper(new SimpleInterval("T", 500, 502));
        final GencodeGtfStartCodonFeature startCodonFeatureBack = provideForTestIs5PrimeUtrStartCodonHelper(new SimpleInterval("T", 1498, 1500));

        final GencodeGtfUTRFeature utr = provideForTestIs5PrimeUtrUTRHelper(new SimpleInterval("T", 750, 1250));

        return new Object[][] {
                // Strand +
                    // Start Codon is Null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.POSITIVE, null), false },
                    // Start Codon non-null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.POSITIVE, startCodonFeatureFront), false },
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.POSITIVE, startCodonFeatureBack), true },
                // Strand -
                    // Start Codon is Null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.NEGATIVE, null), false },
                    // Start Codon non-null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.NEGATIVE, startCodonFeatureFront), true },
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.NEGATIVE, startCodonFeatureBack), false },
        };
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
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100), 10, 10.0/11.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  1, 1.0/3.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  9, 9.0/19.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 10, 11.0/21.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 50, 50.0/100.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),   1, 0.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),   9, 4.0 / 14.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),  10, 5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),   1, 1 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),   9, 10.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),  10, 10.0/16.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),   1, 6.0/8.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),   9, 10.0/24.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  10, 11.0/26.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  50, 50.0/100.0 },

                // Insertions:
                { Allele.create("A", true), Allele.create("AG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("AGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("AGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("CG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("CGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("CGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100), 10, 1 },
                { Allele.create("T", true), Allele.create("TG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  1, 1.0/2.0 },
                { Allele.create("T", true), Allele.create("TGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  9, 9.0/18.0 },
                { Allele.create("T", true), Allele.create("TGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 10, 10.0/20.0 },
                { Allele.create("T", true), Allele.create("TGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 49, 49.0/98.0 },

                // Deletions:
                { Allele.create("AAAAA",  true), Allele.create("A"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),  10,  5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),  10, 10.0/15.0  },
                { Allele.create("TGGGGG", true), Allele.create("T"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  10, 10.0/25.0 },
                { Allele.create("TG", true),     Allele.create("T"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 51),  10, 10.0/21.0 },
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
            gtfFeatureIterator = gencodeHg19FeatureReader.query(contig, start, end);
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

        final ReferenceDataSource muc16TranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19).toPath());
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
            gtfFeatureIterator = gencodeHg19FeatureReader.query(contig, start, end);
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
                IOUtils.getPath(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                requestedTranscriptIds,
                new LinkedHashMap<>(), createFeatureInputForMuc16Ds(GencodeFuncotationFactory.DEFAULT_NAME))) {

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
                                                                    IOUtils.getPath(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19),
                                                                    "VERSION",
                                                                    GencodeFuncotationFactory.DEFAULT_NAME,
                                                                    FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                                                                    new HashSet<>(),
                                                                    new LinkedHashMap<>(), createFeatureInputForMuc16Ds(GencodeFuncotationFactory.DEFAULT_NAME)
                                                                    )) {

            // Generate our funcotations:
            final List<Feature> featureList = new ArrayList<>();
            featureList.add( gene );
            final List<Funcotation> funcotations = funcotationFactory.createFuncotationsOnVariant(variantContext, referenceContext, featureList);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 0);
        }
    }

    @Test (dataProvider = "provideDataForCreateFuncotations")
    public void testCreateFuncotations(final String expectedGeneName,
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
                                final String transcriptFastaFile, final String transcriptGtfFile) {

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
                new LinkedHashMap<>(),
                new FeatureInput<>(transcriptGtfFile, GencodeFuncotationFactory.DEFAULT_NAME, Collections.emptyMap()))) {

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

    @Test(dataProvider = "provideForTestGetTranscriptEndPaddingBases")
    void testGetTranscriptEndPaddingBases(final VariantContext variant,
                                          final Allele altAllele,
                                          final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
                                          final ReferenceContext reference,
                                          final String expected) {
        Assert.assertEquals(GencodeFuncotationFactory.getTranscriptEndPaddingBases(variant, altAllele, exonPositionList, reference), expected);
    }

    @Test(dataProvider = "provideForTestCreateGencodeFuncotationBuilderWithTrivialFieldsPopulated")
    void testCreateGencodeFuncotationBuilderWithTrivialFieldsPopulated(final VariantContext variant,
                                                                       final Allele altAllele,
                                                                       final GencodeGtfGeneFeature gtfFeature,
                                                                       final GencodeGtfTranscriptFeature transcript ) {

        final GencodeFuncotationBuilder builder =
                GencodeFuncotationFactory.createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele,
                        gtfFeature, transcript);
        final GencodeFuncotation gf = builder.gencodeFuncotation;

        // Ultra-trivial checks:
        Assert.assertEquals( gf.getRefAllele(), variant.getReference().getBaseString() );
        Assert.assertEquals( gf.getTranscriptStrand(), transcript.getGenomicStrand().encode() );
        Assert.assertEquals( gf.getHugoSymbol(), gtfFeature.getGeneName() );
        Assert.assertEquals( gf.getNcbiBuild(), gtfFeature.getUcscGenomeVersion() );
        Assert.assertEquals( gf.getChromosome(), gtfFeature.getChromosomeName() );
        Assert.assertEquals( gf.getStart(), variant.getStart() );
        Assert.assertEquals( gf.getGeneTranscriptType(), transcript.getTranscriptType() );
        Assert.assertEquals( gf.getEnd(), variant.getStart() + altAllele.length() - 1 );
        Assert.assertEquals( gf.getTumorSeqAllele2(), altAllele.getBaseString() );
        Assert.assertEquals( gf.getAnnotationTranscript(), transcript.getTranscriptId() );
        Assert.assertEquals( gf.getLocusLevel(), Integer.valueOf(transcript.getLocusLevel().toString()) );
        Assert.assertEquals( gf.getEnd(), variant.getEnd() );

        // Still simple, but computational, checks:
        Assert.assertEquals( gf.getVariantType(), GencodeFuncotationFactory.getVariantType(variant.getReference(), altAllele));
        Assert.assertEquals( gf.getGenomeChange(), GencodeFuncotationFactory.getGenomeChangeString(variant, altAllele) );
        Assert.assertEquals( gf.getApprisRank(), GencodeFuncotationFactory.getApprisRank( transcript ) );
        Assert.assertTrue( gf.getTranscriptLength() == transcript.getExons().stream().mapToInt(Locatable::getLengthOnReference).sum() );
    }

    @Test(dataProvider = "provideForTestIs5PrimeUtr")
    void testIs5PrimeUtr(final GencodeGtfUTRFeature utr, final GencodeGtfTranscriptFeature transcript, final boolean expected) {
        Assert.assertEquals(GencodeFuncotationFactory.is5PrimeUtr( utr, transcript), expected);
    }

    @Test (dataProvider = "provideDataForTestIsVariantInCodingRegion")
    void testIsVariantInCodingRegion(final GencodeFuncotation.VariantClassification varClass, final GencodeFuncotation.VariantClassification secondaryVarClass, final boolean expected) {
        Assert.assertEquals( GencodeFuncotationFactory.isVariantInCodingRegion(varClass, secondaryVarClass), expected );
    }

    @Test ( dataProvider = "provideForCreateDefaultFuncotationsOnProblemVariant")
    void testCreateDefaultFuncotationsOnProblemVariant(final VariantContext variant,
                                                       final Allele altAllele,
                                                       final GencodeGtfGeneFeature gtfFeature,
                                                       final ReferenceContext reference,
                                                       final GencodeGtfTranscriptFeature transcript,
                                                       final String version,
                                                       final GencodeFuncotation expected) {

        final GencodeFuncotation funcotation = GencodeFuncotationFactory.createDefaultFuncotationsOnProblemVariant(variant, altAllele, gtfFeature, reference, transcript, version);
        Assert.assertEquals( funcotation, expected );
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
                new LinkedHashMap<>(), createFeatureInputForCntn4Ds(GencodeFuncotationFactory.DEFAULT_NAME))) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertEquals(gencodeFuncotations.size(), 1);
        }

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.BEST_EFFECT,
                Collections.emptySet(),
                new LinkedHashMap<>(), createFeatureInputForCntn4Ds(GencodeFuncotationFactory.DEFAULT_NAME))) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertEquals(gencodeFuncotations.size(), 1);
        }

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.ALL,
                Collections.emptySet(),
                new LinkedHashMap<>(), createFeatureInputForCntn4Ds(GencodeFuncotationFactory.DEFAULT_NAME))) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertTrue(gencodeFuncotations.size() > 1);
        }
    }

    private static FeatureInput<? extends Feature> createFeatureInputForMuc16Ds(final String dsName) {
        return new FeatureInput<>(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, dsName, Collections.emptyMap());
    }

    private static FeatureInput<? extends Feature> createFeatureInputForCntn4Ds(final String dsName) {
        return new FeatureInput<>(CNTN4_GENCODE_ANNOTATIONS_FILE_NAME, dsName, Collections.emptyMap());
    }
}
