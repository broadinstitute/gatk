package org.broadinstitute.hellbender.tools.funcotator.genelistoutput;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.AnnotatedIntervalToSegmentVariantContextConverter;
import org.broadinstitute.hellbender.tools.funcotator.FlankSettings;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class GeneListOutputRendererUnitTest extends GATKBaseTest {
    private static final String CNTN4_GENCODE_ANNOTATIONS_FILE_NAME = toolsTestDir + "funcotator/gencode.v19.CNTN4.annotation.gtf";
    private static final String CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE = toolsTestDir + "funcotator/gencode.v19.CNTN4.pc_transcripts.fasta";
    private static final VariantContext DUMMY_SEGMENT_VARIANT_CONTEXT = FuncotatorTestUtils.createDummySegmentVariantContext();

    @DataProvider
    public Object[][] provideGeneExonMapTest() {
        final VariantContext segmentVariantContext = FuncotatorTestUtils.createDummySegmentVariantContext();

        // Create a dummy funcotation factory to get real gencode funcotation field names for testing.
        final LinkedHashSet<String> funcotationFields = createDummyGencodeFuncotationFactory()
                .getSupportedFuncotationFieldsForSegments();

        final List<String>  testFieldValues1 = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1", "", "1-", "", "", "");
        final List<Pair<String,String>> gtKeys1  = Arrays.asList(Pair.of("GENE1", "1-"), Pair.of("GENE1-AS1", ""), Pair.of("GENE2", ""));
        final FuncotationMap funcotationMap1 = createTestFuncotationMap(funcotationFields, testFieldValues1);
        final Pair<VariantContext, FuncotationMap> pairValue1 = Pair.of(segmentVariantContext, funcotationMap1);

        final List<String> testFieldValues2 = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1-AS1", "GENE2", "5+", "2-", "", "");
        final List<Pair<String,String>> gtKeys2 = Arrays.asList(Pair.of("GENE1", ""), Pair.of("GENE1-AS1", "5+"), Pair.of("GENE2", "2-"));
        final FuncotationMap funcotationMap2 = createTestFuncotationMap(funcotationFields, testFieldValues2);
        final Pair<VariantContext, FuncotationMap> pairValue2 = Pair.of(segmentVariantContext, funcotationMap2);

        final List<String> testFieldValues3 = Arrays.asList("GENE1,GENE2,GENE1-AS1", "", "", "", "", "", "");
        final List<Pair<String,String>> gtKeys3 = Arrays.asList(Pair.of("GENE1", ""), Pair.of("GENE1-AS1", ""), Pair.of("GENE2", ""));
        final FuncotationMap funcotationMap3 = createTestFuncotationMap(funcotationFields, testFieldValues3);
        final Pair<VariantContext, FuncotationMap> pairValue3 = Pair.of(segmentVariantContext, funcotationMap3);

        final List<String> testFieldValues4 = Arrays.asList("", "", "", "", "", "", "");
        final List<Pair<String,String>> gtKeys4 = Collections.emptyList();
        final FuncotationMap funcotationMap4 = createTestFuncotationMap(funcotationFields, testFieldValues4);
        final Pair<VariantContext, FuncotationMap> pairValue4 = Pair.of(segmentVariantContext, funcotationMap4);

        return new Object[][] {
                // GT values for genes, start_gene, end_gene, start_exon, end_exon, ref_allele, alt_allele.  The last
                //   two are always blank.  Second entry is the GT for the internal map.
                {createGroundTruthGeneExonMap(gtKeys1, pairValue1), segmentVariantContext, funcotationMap1},
                {createGroundTruthGeneExonMap(gtKeys2, pairValue2), segmentVariantContext, funcotationMap2},
                {createGroundTruthGeneExonMap(gtKeys3, pairValue3), segmentVariantContext, funcotationMap3},
                {createGroundTruthGeneExonMap(gtKeys4, pairValue4), segmentVariantContext, funcotationMap4}
        };
    }

    private static SortedMap<Pair<String,String>, Pair<VariantContext, FuncotationMap>> createGroundTruthGeneExonMap(final List<Pair<String, String>> gtKeys, final Pair<VariantContext, FuncotationMap> pairValue3) {
        return gtKeys.stream().collect(Collectors.toMap(k -> k, k -> pairValue3, (x1, x2) -> {
                    throw new IllegalArgumentException("Should not be able to have duplicate field names.");
                }, TreeMap::new));
    }

    private FuncotationMap createTestFuncotationMap(final LinkedHashSet<String> funcotationFields, final List<String> testFuncotationFieldValues) {
        return FuncotationMap.createNoTranscriptInfo(
                Collections.singletonList(
                        TableFuncotation.create(funcotationFields, testFuncotationFieldValues,
                                AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                                "TEST",
                                FuncotationMetadataUtils.createWithUnknownAttributes(Lists.newArrayList(funcotationFields.iterator()))
                        )
                )
        );
    }

    private FuncotationMap createTestFuncotationMap(final List<String> funcotationFields, final List<String> testFuncotationFieldValues) {
        return createTestFuncotationMap(new LinkedHashSet<>(funcotationFields), testFuncotationFieldValues);
    }

    /**
     * Whitebox test making sure that the internal state of the geneExonToVCFuncotaitonMap is consistent.
     * Does not test default annotations, override annotations, or excluded fields.
     */
    @Test(dataProvider = "provideGeneExonMapTest")
    public void testWriteGeneExonToVariantFuncotationMap(SortedMap<Pair<String,String>, Pair<VariantContext, FuncotationMap>> gtMap,
        final VariantContext segmentVariantContext, final FuncotationMap funcotationMap) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");

        geneListOutputRenderer.write(segmentVariantContext, funcotationMap);
        final SortedMap<Pair<String,String>, Pair<VariantContext, FuncotationMap>> sortedMap = geneListOutputRenderer.getGeneExonToVariantFuncotationMap();
        Assert.assertEquals(sortedMap, gtMap);
    }

    //TODO: Add 2 write / validateAbleToWrite tests for segment data < 150 bases - .  Use data from issue.

    @DataProvider
    public Object[][] provideForSegmentLengthTests() {

        final LinkedHashSet<String> funcotationFields = createDummyGencodeFuncotationFactory()
                .getSupportedFuncotationFieldsForSegments();

        final List<String> testFieldValues = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1", "", "1-", "", "", "");
        final FuncotationMap funcotationMap = createTestFuncotationMap(funcotationFields, testFieldValues);

        final int start = 356000;

        return new Object[][] {
            {
                funcotationMap,
                    FuncotatorTestUtils.createDummySegmentVariantContext(
                            start,
                            start + FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT - 1,
                            "T"
                    )
            },
        };
    }

    @Test(dataProvider = "provideForSegmentLengthTests",
          expectedExceptions = UserException.BadInput.class)
    public void testValidateAbleToWriteFailureOnShortLength(final FuncotationMap funcotationMap,
                                                            final VariantContext segmentVariantContext) throws IOException {
        final File outputFile = File.createTempFile("testFileForSegmentLengthTesting", ".seg");
        final GeneListOutputRenderer geneListOutputRenderer =
                new GeneListOutputRenderer(
                        outputFile.toPath(),
                        new LinkedHashMap<>(),
                        new LinkedHashMap<>(),
                        new HashSet<>(),
                        "TEST_TOOL"
                );

        geneListOutputRenderer.validateAbleToWrite(segmentVariantContext, funcotationMap);
    }

    @Test(dataProvider = "provideForSegmentLengthTests")
    public void testValidateAbleToWriteSuccessOnShortLengthWithNonDefaultMinSegmentSize(
            final FuncotationMap funcotationMap, final VariantContext segmentVariantContext
    ) throws IOException {
        final File outputFile = File.createTempFile("testFileForSegmentLengthTesting", ".seg");
        final GeneListOutputRenderer geneListOutputRenderer =
                new GeneListOutputRenderer(
                        outputFile.toPath(),
                        new LinkedHashMap<>(),
                        new LinkedHashMap<>(),
                        new HashSet<>(),
                        "TEST_TOOL",
                        FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT/2
                );

        geneListOutputRenderer.validateAbleToWrite(segmentVariantContext, funcotationMap);
    }

    @DataProvider
    public Object[][] provideErroneousConfigurationsForGeneListMultipleGencodeFields() {

        final LinkedHashSet<String> gencodeFields = createDummyGencodeFuncotationFactory().getSupportedFuncotationFieldsForSegments();

        final List<String> multipleGencodeKeysFields = Stream.concat(gencodeFields.stream(), Stream.of("Gencode_000972_genes")).collect(Collectors.toList());
        final List<String> multipleGencodeKeysFieldValues = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1", "", "1-", "", "", "", "GO_BOOM,DENIED");


        return new Object[][] {
                // Multiple GencodeKeys
                {DUMMY_SEGMENT_VARIANT_CONTEXT, createTestFuncotationMap(multipleGencodeKeysFields, multipleGencodeKeysFieldValues)},
        };
    }

    @DataProvider
    public Object[][] provideErroneousWriteGeneExonToVariantFuncotationMapMissingGencodeFields() {
        final LinkedHashSet<String> gencodeFields = createDummyGencodeFuncotationFactory().getSupportedFuncotationFieldsForSegments();
        final List<String> gencodeTestFieldValues = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1-AS1", "GENE2", "5+", "2-", "", "");

        final List<String> missingGencodeFields = Arrays.asList("Gencode_28_genes", "Gencode_28_start_gene", "Gencode_28_end_gene", "Gencode_28_start_exon", "Gencode_28_end_exon");

        final List<FuncotationMap> missingGencodeFieldsFuncotationMaps = IntStream.range(0, missingGencodeFields.size()).boxed()
                .map(i -> createTestFuncotationMap(remove(gencodeFields, missingGencodeFields.get(i)), remove(gencodeTestFieldValues, gencodeTestFieldValues.get(i))))
                .collect(Collectors.toList());

        return new Object[][]{
                // Missing required fields, please note that this list must be updated manually if the list of required fields changes above.
                {DUMMY_SEGMENT_VARIANT_CONTEXT, missingGencodeFieldsFuncotationMaps.get(0)},
                {DUMMY_SEGMENT_VARIANT_CONTEXT, missingGencodeFieldsFuncotationMaps.get(1)},
                {DUMMY_SEGMENT_VARIANT_CONTEXT, missingGencodeFieldsFuncotationMaps.get(2)},
                {DUMMY_SEGMENT_VARIANT_CONTEXT, missingGencodeFieldsFuncotationMaps.get(3)},
                {DUMMY_SEGMENT_VARIANT_CONTEXT, missingGencodeFieldsFuncotationMaps.get(4)},

        };
    }

    @DataProvider
    public Object[][] provideErroneousWriteGeneExonToVariantFuncotationMapWrongAltAlleles() {
        final LinkedHashSet<String> gencodeFields = createDummyGencodeFuncotationFactory().getSupportedFuncotationFieldsForSegments();
        final List<String> gencodeTestFieldValues = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1-AS1", "GENE2", "5+", "2-", "", "");

        final FuncotationMap funcotationMapWithExtraAlleles = FuncotationMap.createNoTranscriptInfo(
                Collections.singletonList(TableFuncotation.create(gencodeFields,
                        gencodeTestFieldValues, AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                        "TEST_DS", null)));
        funcotationMapWithExtraAlleles.add(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, TableFuncotation.create(gencodeFields,
                gencodeTestFieldValues, Allele.create(SimpleSVType.createBracketedSymbAlleleString("INS"), false),
                "TEST_DS", null));


        return new Object[][]{
                // Too many alt alleles
                {DUMMY_SEGMENT_VARIANT_CONTEXT, funcotationMapWithExtraAlleles},
        };
    }
    @DataProvider
    public Object[][] provideErroneousWriteGeneExonToVariantFuncotationMapInvalidTranscriptID() {
        final LinkedHashSet<String> gencodeFields = createDummyGencodeFuncotationFactory().getSupportedFuncotationFieldsForSegments();
        final List<String> gencodeTestFieldValues = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1-AS1", "GENE2", "5+", "2-", "", "");

        final FuncotationMap funcotationMapWithOneRealTxId = FuncotationMap.createEmpty();
        funcotationMapWithOneRealTxId.add("TEST_TX_0000001", TableFuncotation.create(gencodeFields,
                gencodeTestFieldValues, AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                "TEST_DS", null));

        return new Object[][]{
                // No dummy txID, but only one txID
                {DUMMY_SEGMENT_VARIANT_CONTEXT, funcotationMapWithOneRealTxId},
        };
    }

    @DataProvider
    public Object[][] provideErroneousWriteGeneExonToVariantFuncotationMapInvalidVariantContext() {
        final LinkedHashSet<String> gencodeFields = createDummyGencodeFuncotationFactory().getSupportedFuncotationFieldsForSegments();
        final List<String> gencodeTestFieldValues = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1-AS1", "GENE2", "5+", "2-", "", "");

        return new Object[][]{
                // Variant context is not a segment (FuncotationMap is reasonable here)
                {FuncotatorTestUtils.createSimpleVariantContext(FuncotatorReferenceTestUtils.retrieveB37Chr3Ref(),"3", 150, 150, "A", "C"),
                        createTestFuncotationMap(gencodeFields, gencodeTestFieldValues)}
        };
    }

    @DataProvider
    public Object[][] provideErroneousWriteGeneExonToVariantFuncotationMapNotExactlyOneTranscript() {
        final LinkedHashSet<String> gencodeFields = createDummyGencodeFuncotationFactory().getSupportedFuncotationFieldsForSegments();
        final List<String> gencodeTestFieldValues = Arrays.asList("GENE1,GENE2,GENE1-AS1", "GENE1-AS1", "GENE2", "5+", "2-", "", "");

        final FuncotationMap funcotationMapWithExtraTxIds = FuncotationMap.createNoTranscriptInfo(
                Collections.singletonList(TableFuncotation.create(gencodeFields,
                        gencodeTestFieldValues, AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                        "TEST_DS", null)));
        funcotationMapWithExtraTxIds.add("TEST_TX_000001", TableFuncotation.create(gencodeFields,
                gencodeTestFieldValues, AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                "TEST_DS", null));

        final FuncotationMap funcotationMapWithManyRealTxId = FuncotationMap.createEmpty();
        funcotationMapWithManyRealTxId.add("TEST_TX_0000001", TableFuncotation.create(gencodeFields,
                gencodeTestFieldValues, AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                "TEST_DS", null));
        funcotationMapWithManyRealTxId.add("TEST_TX_0000002", TableFuncotation.create(gencodeFields,
                gencodeTestFieldValues, AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                "TEST_DS", null));

        return new Object[][]{
                // Too many txIDs, incl. dummy txID
                {DUMMY_SEGMENT_VARIANT_CONTEXT, funcotationMapWithExtraTxIds},

                // Too many txIDs, not incl. dummy txID
                {DUMMY_SEGMENT_VARIANT_CONTEXT, funcotationMapWithManyRealTxId},

                // No txIDs at all
                {DUMMY_SEGMENT_VARIANT_CONTEXT, FuncotationMap.createEmpty()},
        };
    }

    @Test(expectedExceptions = UserException.BadInput.class, dataProvider = "provideErroneousConfigurationsForGeneListMultipleGencodeFields")
    public void testErroneousWriteGeneExonToVariantFuncotationMapMultipleGencodeFields(final VariantContext segmentVariantContext, final FuncotationMap funcotationMap) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");
        geneListOutputRenderer.write(segmentVariantContext, funcotationMap);
    }

    @Test(expectedExceptions = UserException.BadInput.class, expectedExceptionsMessageRegExp = ".*" + GeneListOutputRenderer.MISSING_GENCODE_FIELD_ERROR_MSG + ".*",
            dataProvider = "provideErroneousWriteGeneExonToVariantFuncotationMapMissingGencodeFields")
    public void testErroneousWriteGeneExonToVariantFuncotationMapMissingGencodeField(final VariantContext segmentVariantContext, final FuncotationMap funcotationMap) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");
        geneListOutputRenderer.write(segmentVariantContext, funcotationMap);
    }

    @Test(expectedExceptions = UserException.BadInput.class, expectedExceptionsMessageRegExp = ".*" + GeneListOutputRenderer.ONE_ALTERNATE_ALLELE_ERR_MSG + ".*",
            dataProvider = "provideErroneousWriteGeneExonToVariantFuncotationMapWrongAltAlleles")
    public void testErroneousWriteGeneExonToVariantFuncotationMapOneAltAllele(final VariantContext segmentVariantContext, final FuncotationMap funcotationMap) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");
        geneListOutputRenderer.write(segmentVariantContext, funcotationMap);
    }

    @Test(expectedExceptions = UserException.BadInput.class, expectedExceptionsMessageRegExp = ".*" + GeneListOutputRenderer.INVALID_VARIANT_CONTEXT_ERR_MSG + ".*",
            dataProvider = "provideErroneousWriteGeneExonToVariantFuncotationMapInvalidVariantContext")
    public void testErroneousWriteGeneExonToVariantFuncotationMapInvalidVariantContext(final VariantContext segmentVariantContext, final FuncotationMap funcotationMap) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");
        geneListOutputRenderer.write(segmentVariantContext, funcotationMap);
    }

    @Test(expectedExceptions = UserException.BadInput.class, expectedExceptionsMessageRegExp = ".*" + GeneListOutputRenderer.INCORRECT_NUM_TRANSCRIPT_ERR_MSG + ".*",
            dataProvider = "provideErroneousWriteGeneExonToVariantFuncotationMapNotExactlyOneTranscript")
    public void testErroneousWriteGeneExonToVariantFuncotationMapNotExactlyOneTranscript(final VariantContext segmentVariantContext, final FuncotationMap funcotationMap) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");
        geneListOutputRenderer.write(segmentVariantContext, funcotationMap);
    }

    @Test(expectedExceptions = UserException.BadInput.class, expectedExceptionsMessageRegExp = ".*" + GeneListOutputRenderer.INVALID_TRANSCRIPT_ID_ERR_MSG + ".*",
            dataProvider = "provideErroneousWriteGeneExonToVariantFuncotationMapInvalidTranscriptID")
    public void testErroneousWriteGeneExonToVariantFuncotationMapNotValidTranscriptId(final VariantContext segmentVariantContext, final FuncotationMap funcotationMap) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");
        geneListOutputRenderer.write(segmentVariantContext, funcotationMap);
    }

    private GencodeFuncotationFactory createDummyGencodeFuncotationFactory() {
        final FeatureInput<GencodeGtfFeature> gencodeFeatureInput = new FeatureInput<>(CNTN4_GENCODE_ANNOTATIONS_FILE_NAME,
                GencodeFuncotationFactory.DEFAULT_NAME, Collections.emptyMap());
        return new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "28",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                Collections.emptySet(),
                new LinkedHashMap<>(),
                gencodeFeatureInput,
                new FlankSettings(0, 0),
                "TEST");
    }

    private static List<String> remove(final Collection<String> removeFrom, final String item) {
        return removeFrom.stream().filter(s -> !s.equals(item)).collect(Collectors.toList());
    }

    /**
     * @return Standard Object[][] for a data provider. Never {@code null}
     */
    @DataProvider
    public Object[][] provideGeneExonMapSortingTest() {
        final VariantContext segmentVariantContext = FuncotatorTestUtils.createDummySegmentVariantContext();

        // Create a dummy funcotation factory to get real gencode funcotation field names for testing.
        final LinkedHashSet<String> funcotationFields = createDummyGencodeFuncotationFactory()
                .getSupportedFuncotationFieldsForSegments();

        final List<String> testFieldValues1 = Arrays.asList("A223,A3,A3-AS1", "A223", "", "1-", "", "", "");
        final List<Pair<String,String>> gtSorting1 = Arrays.asList(
                Pair.of("A3", ""),
                Pair.of("A3-AS1", ""),
                Pair.of("A223", "1-"));
        final FuncotationMap funcotationMap1 = createTestFuncotationMap(funcotationFields, testFieldValues1);

        final List<String> testFieldValues2 = Arrays.asList("A223,A3,A3-AS1", "A223", "A223", "10-", "1+", "", "");
        final List<Pair<String,String>> gtSorting2 = Arrays.asList(
                Pair.of("A3", ""),
                Pair.of("A3-AS1", ""),
                Pair.of("A223", "1+"),
                Pair.of("A223", "10-"));
        final FuncotationMap funcotationMap2 = createTestFuncotationMap(funcotationFields, testFieldValues2);

        final List<String> testFieldValues3 = Arrays.asList("A223,A3,A3-AS1", "", "A223", "", "1-", "", "");
        final List<Pair<String,String>> gtSorting3 = Arrays.asList(
                Pair.of("A3", ""),
                Pair.of("A3-AS1", ""),
                Pair.of("A223", "1-"));
        final FuncotationMap funcotationMap3 = createTestFuncotationMap(funcotationFields, testFieldValues3);

        final List<String> testFieldValues4 = Arrays.asList("A223,A3,A3-AS1", "", "A223", "", "", "", "");
        final List<Pair<String,String>> gtSorting4 = Arrays.asList(
                Pair.of("A3", ""),
                Pair.of("A3-AS1", ""),
                Pair.of("A223", ""));
        final FuncotationMap funcotationMap4 = createTestFuncotationMap(funcotationFields, testFieldValues4);

        final List<String> testFieldValues5a = Arrays.asList("A223,A3,A3-AS1", "A223", "", "10+", "", "", "");
        final List<String> testFieldValues5b = Arrays.asList("A223,B4,B10", "A223", "A223", "5+", "9-", "", "");
        final List<String> testFieldValues5c = Arrays.asList("A223,B4,B10", "A223", "", "4-", "", "", "");
        final List<Pair<String,String>> gtSorting5 = Arrays.asList(
                Pair.of("A3", ""),
                Pair.of("A3-AS1", ""),
                Pair.of("A223", "4-"),
                Pair.of("A223", "5+"),
                Pair.of("A223", "9-"),
                Pair.of("A223", "10+"),
                Pair.of("B4", ""),
                Pair.of("B10", "")
                );
        final FuncotationMap funcotationMap5a = createTestFuncotationMap(funcotationFields, testFieldValues5a);
        final FuncotationMap funcotationMap5b = createTestFuncotationMap(funcotationFields, testFieldValues5b);
        final FuncotationMap funcotationMap5c = createTestFuncotationMap(funcotationFields, testFieldValues5c);

        return new Object[][] {
                // GT values for genes, start_gene, end_gene, start_exon, end_exon, ref_allele, alt_allele.  The last
                //   two are always blank.  Second entry is the GT for the internal map.
                {gtSorting1, Collections.singletonList(segmentVariantContext), Collections.singletonList(funcotationMap1)},
                {gtSorting2, Collections.singletonList(segmentVariantContext), Collections.singletonList(funcotationMap2)},
                {gtSorting3, Collections.singletonList(segmentVariantContext), Collections.singletonList(funcotationMap3)},
                {gtSorting4, Collections.singletonList(segmentVariantContext), Collections.singletonList(funcotationMap4)},
                {gtSorting5, createSegmentVariantContexts(segmentVariantContext, 2),
                        Arrays.asList(funcotationMap5a, funcotationMap5b, funcotationMap5c)},
        };
    }

    // numToCreate is the number to add to the starting varaint context.
    private List<VariantContext> createSegmentVariantContexts(final VariantContext startingVariantContext, int numToCreate) {
        ParamUtils.isPositiveOrZero(numToCreate, "Must create at least 0.");
        final List<VariantContext> result = Lists.newArrayList(startingVariantContext);
        IntStream.range(0,numToCreate).boxed().forEach(i -> {
            final int length = 100000;
            final int newStart = startingVariantContext.getEnd() + 10 + (length * i);
            result.add(
                new VariantContextBuilder(startingVariantContext)
                .start(newStart)
                .stop(newStart + length)
                .make()
            );
        });
        return result;
    }

    /**
     * More whitebox testing of the gene exon map.
     */
    @Test(dataProvider = "provideGeneExonMapSortingTest")
    public void testSortingOfGeneExonPair(final List<Pair<String,String>> gtSorting,
                                          final List<VariantContext> segmentVariantContexts, final List<FuncotationMap> funcotationMaps) throws IOException {
        final File outputFile = File.createTempFile("simpleSegFileWriting", ".seg");

        Assert.assertEquals(segmentVariantContexts.size(), funcotationMaps.size(), "The number of test VariantContexts and FuncotationMaps must be equal.  This needs to be fixed in the test code by a GATK developer.");

        final GeneListOutputRenderer geneListOutputRenderer = new GeneListOutputRenderer(outputFile.toPath(), new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), "TEST_TOOL");

        // Populate the gene-exon map
        IntStream.range(0, funcotationMaps.size()).boxed().forEach(i -> geneListOutputRenderer.write(segmentVariantContexts.get(i), funcotationMaps.get(i)));
        final SortedMap<Pair<String,String>, Pair<VariantContext, FuncotationMap>> sortedMap = geneListOutputRenderer.getGeneExonToVariantFuncotationMap();
        Assert.assertEquals(Lists.newArrayList(sortedMap.keySet()), gtSorting);
    }
}
