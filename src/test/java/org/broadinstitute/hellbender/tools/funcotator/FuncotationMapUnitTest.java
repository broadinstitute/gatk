package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.collect.ImmutableSortedMap;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfCodec;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class FuncotationMapUnitTest extends BaseTest{

    private static final String DS_PIK3CA_DIR = GATKBaseTest.largeFileTestDir + "funcotator/small_ds_pik3ca/";
    private static final String DS_PIK3CA_HG19_GENCODE_FASTA = DS_PIK3CA_DIR + "gencode_pik3ca/hg19/gencode.v19.PIK3CA_transcript.fasta";
    private static final String DS_MUC16_DIR = GATKBaseTest.largeFileTestDir + "funcotator/small_ds_muc16/";
    private static final String DS_MUC16_HG19_GENCODE_FASTA = DS_MUC16_DIR + "gencode_muc16/hg19/gencode.v19.MUC16_transcript.fasta";
    private static final FeatureReader<GencodeGtfFeature> pik3caFeatureReader = AbstractFeatureReader.getFeatureReader( FuncotatorTestConstants.PIK3CA_GENCODE_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec() );
    private static final FeatureReader<GencodeGtfFeature> muc16FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.MUC16_GENCODE_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec() );

    @DataProvider
    public Object[][] provideCreationFromFuncotationVcfHeaderString() {
        return new Object[][] {
                {
                    Allele.create("C"), "FOO",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|4]",
                        Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                        Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    // Test when the last field is blank
                    Allele.create("C"), "FOO",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", ""))
                }, {
                    // Test when the first field is blank
                    Allele.create("C"), "FOO",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[|hg19|4]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    // Test when we have two transcripts
                    Allele.create("C"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_annotationTranscript",
                    "[FOO|hg19|txID1]#[BAR|hg38|txID2]",
                    Arrays.asList("txID1", "txID2"),
                    Arrays.asList(
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_annotationTranscript", "txID1"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAR", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID2")
                        )
                }, {
                    // Test when we have three transcripts
                    Allele.create("C"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_annotationTranscript",
                    "[FOO|hg19|txID1]#[BAR|hg38|txID2]#[BAZ|hg38|txID3]",
                    Arrays.asList("txID1", "txID2", "txID3"),
                    Arrays.asList(
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_annotationTranscript", "txID1"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAR", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID2"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAZ", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID3")
                        )
                }, {
                    // The rest of the tests are testing same as above, but with different alleles and different datasource names
                    Allele.create("A"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|4]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    Allele.create("A"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", ""))
                }, {
                    Allele.create("A"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[|hg19|4]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    Allele.create("A"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_annotationTranscript",
                    "[FOO|hg19|txID1]#[BAR|hg38|txID2]",
                    Arrays.asList("txID1", "txID2"),
                    Arrays.asList(
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_annotationTranscript", "txID1"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAR", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID2")
                        )
                }, {
                    Allele.create("A"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_annotationTranscript",
                    "[FOO|hg19|txID1]#[BAR|hg38|txID2]#[BAZ|hg38|txID3]",
                    Arrays.asList("txID1", "txID2", "txID3"),
                    Arrays.asList(
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_annotationTranscript", "txID1"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAR", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID2"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAZ", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID3")
                        )
                }, {
                    Allele.create("AT"), "TEST",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|4]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    Allele.create("AT"), "BAZ",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", ""))
                }, {
                    Allele.create("AT"), "BAZ",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[|hg19|4]",
                    Collections.singletonList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Collections.singletonList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    Allele.create("AT"), "BAZ",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_annotationTranscript",
                    "[FOO|hg19|txID1]#[BAR|hg38|txID2]",
                    Arrays.asList("txID1", "txID2"),
                    Arrays.asList(
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_annotationTranscript", "txID1"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAR", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID2")
                        )
                }, {
                    Allele.create("AT"), "BAZ",
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_annotationTranscript",
                    "[FOO|hg19|txID1]#[BAR|hg38|txID2]#[BAZ|hg38|txID3]",
                    Arrays.asList("txID1", "txID2", "txID3"),
                    Arrays.asList(
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_annotationTranscript", "txID1"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAR", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID2"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAZ", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID3")
                        )
                }
        };
    }

    @Test(dataProvider = "provideCreationFromFuncotationVcfHeaderString")
    public void testCreationFromFuncotationVcfHeaderString(final Allele gtAllele, final String gtDatasourceName, final String headerDescription, final String funcotationValue, final List<String> gtTranscriptIDs, final List<SortedMap<String, String>> gtMaps) {

        final FuncotationMap testMap = FuncotationMap.createAsAllTableFuncotationsFromVcf("Gencode_19_annotationTranscript",
                FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(headerDescription),
                funcotationValue, gtAllele, gtDatasourceName);
        final List<String> transcriptIds = testMap.getTranscriptList();
        Assert.assertEquals(new HashSet<>(transcriptIds), new HashSet<>(gtTranscriptIDs));
        for (int i = 0; i < gtMaps.size(); i++){
            final SortedMap<String, String> gtMap = gtMaps.get(i);
            final String transcriptId = transcriptIds.get(i);
            Assert.assertEquals(testMap.get(transcriptId).get(0).getFieldNames().size(), gtMap.keySet().size());
            Assert.assertEquals(testMap.get(transcriptId).size(), 1); // We have one funcotation with X fields.
            Assert.assertTrue(gtMap.keySet().stream().allMatch(k -> gtMap.get(k).equals(testMap.getFieldValue(transcriptId, k, gtAllele))));
            Assert.assertTrue(testMap.get(transcriptId).stream().allMatch(f -> f.getAltAllele().equals(gtAllele)));
            Assert.assertTrue(testMap.get(transcriptId).stream().allMatch(f -> f.getDataSourceName().equals(gtDatasourceName)));
        }
    }

    @DataProvider
    public Object[][] provideGencodeFuncotationCreation() {
        // All of these were chosen not to be IGR.
        return new Object[][] {
                {"chr3", 178916538, 178916538, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref())),
                        pik3caFeatureReader, DS_PIK3CA_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Collections.singletonList("ENST00000263967.3")
                },{"chr3", 178916538, 178916538, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref())),
                        pik3caFeatureReader, DS_PIK3CA_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.CANONICAL, Collections.singletonList("ENST00000263967.3")
                },{"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Arrays.asList("ENST00000397910.4", "ENST00000380951.5")

                // Next one tests where we would be in a gene with more than one basic transcript, but variant only overlaps one.  And we still ask for all,
                //   but since one is IGR, it will never get added the the FuncotationMap.
                }, {"chr19", 9014550, 9014550, "T", "A", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of(IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Collections.singletonList("ENST00000397910.4")

                // Next one tests where we would be in a gene with more than one basic transcript, variant overlaps both, but we are in canonical mode.
                },{"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.CANONICAL, Collections.singletonList("ENST00000397910.4")

                // Next one tests where we would be in a gene with more than one basic transcript, variant overlaps both, but we are in effect mode.
                },{"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.BEST_EFFECT, Collections.singletonList("ENST00000397910.4")
                }
        };
    }
    @Test(dataProvider = "provideGencodeFuncotationCreation")
    public void testGencodeFuncotationCreation(final String contig,
                                               final int start,
                                               final int end,
                                               final String ref,
                                               final String alt,
                                               final String referenceFileName,
                                               final ReferenceDataSource referenceDataSource,
                                               final FeatureReader<GencodeGtfFeature> featureReader,
                                               final String transcriptFastaFile,
                                               final TranscriptSelectionMode transcriptSelectionMode,
                                               final List<String> gtTranscripts) {

        final List<GencodeFuncotation> gencodeFuncotations = createGencodeFuncotations(contig, start, end, ref, alt, referenceFileName, referenceDataSource, featureReader, transcriptFastaFile, transcriptSelectionMode);

        final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(gencodeFuncotations);

        Assert.assertEquals(funcotationMap.getTranscriptList(), gtTranscripts);
        Assert.assertTrue(funcotationMap.getTranscriptList().stream().allMatch(k -> funcotationMap.get(k).size() == 1));
        Assert.assertTrue(funcotationMap.getTranscriptList().stream()
                .noneMatch(k -> ((GencodeFuncotation) funcotationMap.get(k).get(0)).getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR) ));
        Assert.assertTrue(funcotationMap.getTranscriptList().stream()
                .noneMatch(k -> ((GencodeFuncotation) funcotationMap.get(k).get(0)).getVariantClassification().equals(GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE) ));
    }

    private static List<GencodeFuncotation> createGencodeFuncotations(final String contig, final int start, final int end, final String ref, final String alt, final String referenceFileName, final ReferenceDataSource referenceDataSource, final FeatureReader<GencodeGtfFeature> featureReader, final String transcriptFastaFile, final TranscriptSelectionMode transcriptSelectionMode) {
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );
        final VariantContext variantContext = createVariantContext(contig, start, end, ref, alt, referenceFileName);

        final ReferenceContext referenceContext = new ReferenceContext(referenceDataSource, variantInterval );

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = featureReader.query(variantContext.getContig(), variantContext.getStart(), variantContext.getEnd());
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }
        final List<Feature> featureList = Collections.singletonList(gtfFeatureIterator.next());

        final String gencode_test = "GENCODE_TEST";
        final GencodeFuncotationFactory gencodeFactory = new GencodeFuncotationFactory(Paths.get(transcriptFastaFile),
        "TEST", gencode_test, transcriptSelectionMode, new HashSet<>(), new LinkedHashMap<>());

        return gencodeFactory.createFuncotations(variantContext, referenceContext, Collections.singletonMap(gencode_test, featureList)).stream()
        .map(f -> (GencodeFuncotation) f).collect(Collectors.toList());
    }

    private List<String> createFieldValuesFromNameList(final String prefix, final List<String> baseFieldList) {
        return baseFieldList.stream().map(f -> prefix + f + "value").collect(Collectors.toList());
    }

    @DataProvider
    public Object[][] provideTableFuncotations() {
        final List<String> baseFieldNameList = Arrays.asList("FOO", "BAR");
        return new Object[][]{
        // NOTE: The data field names must match data sources that are checked in for this to work in an expected way:
                {
                        Collections.singletonList(
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("A", baseFieldNameList),
                                        Allele.create("T"),
                                        GencodeFuncotationFactory.DEFAULT_NAME, null
                                )
                        )},{
                        Collections.singletonList(
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("B", baseFieldNameList),
                                        Allele.create("C"),
                                        GencodeFuncotationFactory.DEFAULT_NAME, null
                                )
                        )},{
                        Collections.singletonList(
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("C", baseFieldNameList),
                                        Allele.create("GG"),
                                        GencodeFuncotationFactory.DEFAULT_NAME, null
                                )
                        )},{
                        Collections.singletonList(
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("D", baseFieldNameList),
                                        Allele.create("T"),
                                        "TestDataSource4", null
                                )
                        )},{
                        Arrays.asList(
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("E", baseFieldNameList),
                                        Allele.create("A"),
                                        "TestDataSource5", null
                                ),
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("F", baseFieldNameList),
                                        Allele.create("AG"),
                                        "TestDataSource5", null
                                ),
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("G", baseFieldNameList),
                                        Allele.create("AT"),
                                        "TestDataSource5", null
                                )
                        )
                }

        };
    }


    @DataProvider
    public Object[][] provideTestAdd() {
        final List<String> baseFieldNameList = Arrays.asList("TESTFIELD1", "TESTADD");
        return new Object[][]{
                {"chr3", 178916538, 178916538, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref())),
                        pik3caFeatureReader, DS_PIK3CA_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Collections.singletonList("ENST00000263967.3"),
                        Arrays.asList(
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("E", baseFieldNameList),
                                        Allele.create("A"),
                                        "TestDataSource5", null
                                ),
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("F", baseFieldNameList),
                                        Allele.create("AG"),
                                        "TestDataSource5", null
                                ),
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("G", baseFieldNameList),
                                        Allele.create("AT"),
                                        "TestDataSource5", null
                                )
                        )
                },{"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.BEST_EFFECT, Collections.singletonList("ENST00000397910.4"),
                        Arrays.asList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("E", baseFieldNameList),
                                Allele.create("A"),
                                "TestDataSource5", null
                        ),
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("F", baseFieldNameList),
                                Allele.create("AG"),
                                "TestDataSource5", null
                        ),
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("G", baseFieldNameList),
                                Allele.create("AT"),
                                "TestDataSource5", null
                        )
                )
                }, {"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Arrays.asList("ENST00000397910.4", "ENST00000380951.5"),
                        Arrays.asList(
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("E", baseFieldNameList),
                                        Allele.create("A"),
                                        "TestDataSource5", null
                                ),
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("F", baseFieldNameList),
                                        Allele.create("AG"),
                                        "TestDataSource5", null
                                ),
                                TableFuncotation.create(
                                        baseFieldNameList,
                                        createFieldValuesFromNameList("G", baseFieldNameList),
                                        Allele.create("AT"),
                                        "TestDataSource5", null
                                )
                        )
                }
        };
    }

    /**
     * Also tests that {@link FuncotationMap#getGencodeFuncotations(String)} will return a correct list.
     */
    @Test(dataProvider = "provideTestAdd")
    public void testAddAndGet(final String contig,
                              final int start,
                              final int end,
                              final String ref,
                              final String alt,
                              final String referenceFileName,
                              final ReferenceDataSource referenceDataSource,
                              final FeatureReader<GencodeGtfFeature> featureReader,
                              final String transcriptFastaFile,
                              final TranscriptSelectionMode transcriptSelectionMode,
                              final List<String> gtTranscripts, final List<Funcotation> funcotationsToAdd){
        final List<GencodeFuncotation> gencodeFuncotations = createGencodeFuncotations(contig, start, end, ref, alt, referenceFileName, referenceDataSource, featureReader, transcriptFastaFile, transcriptSelectionMode);
        final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(gencodeFuncotations);

        // Let's make sure that the gtTranscripts match what is in the map, even if this is tested elsewhere
        Assert.assertEquals(funcotationMap.getTranscriptList(), gtTranscripts);

        for (final String transcriptId : funcotationMap.getTranscriptList()) {
            funcotationMap.add(transcriptId, funcotationsToAdd);
        }

        for (final Funcotation funcotation : funcotationsToAdd) {
            for (final String transcriptId : funcotationMap.getTranscriptList()) {
                Assert.assertTrue(funcotationMap.get(transcriptId).contains(funcotation), "Missing funcotation for " + transcriptId + ":" + funcotation);
            }
        }

        for (final String transcriptId : funcotationMap.getTranscriptList()) {
            final List<GencodeFuncotation> gencodeFuncotationList = funcotationMap.getGencodeFuncotations(transcriptId);
            Assert.assertEquals(gencodeFuncotationList.size(), 1);
            Assert.assertEquals(gencodeFuncotationList.get(0).getAnnotationTranscript(), transcriptId);
        }

    }

    private static VariantContext createVariantContext(final String contig, final int start, final int end, final String ref, final String alt, final String referenceFileName) {
        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                referenceFileName,
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        return variantContextBuilder.make();
    }

    @DataProvider
    public Object[][] provideFakeGencodeData() {

        return new Object[][] {
                {Collections.singletonList(new GencodeFuncotationBuilder().setAnnotationTranscript("TXID").build())},
                {Arrays.asList(new GencodeFuncotationBuilder().setAnnotationTranscript("TXID1").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID2").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID3").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID4").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID5").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID6").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID7").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID8").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID9").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID10").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID11").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID12").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID13").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID14").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID15").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID16").build(),
                        new GencodeFuncotationBuilder().setAnnotationTranscript("TXID17").build())
                }
        };
    }

    @Test(dataProvider = "provideFakeGencodeData")
    public void testTranscriptListAndSetAreOrdered(final List<GencodeFuncotation> gencodeFuncotations) {
        // Tests that both the list and the set are ordered the same way and produce equals values.
        final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(gencodeFuncotations);
        Assert.assertEquals(funcotationMap.getTranscriptList(), new ArrayList<>(funcotationMap.getTranscriptSet()));
    }

    @Test(expectedExceptions = GATKException.ShouldNeverReachHereException.class, expectedExceptionsMessageRegExp = ".*a Gencode Funcotation cannot be added.*")
    public void testAddingGencodeFuncotationToFuncotationMap() {
        // Create some gencode funcotations.  The content does not really matter here.
        final List<Funcotation> gencodeFuncotations = createGencodeFuncotations("chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                TranscriptSelectionMode.ALL).stream().map(gf -> (Funcotation) gf).collect(Collectors.toList());

        // Create a funcotationMap with some pre-made funcotations.  Content does not really matter.
        final FuncotationMap funcotationMap = FuncotationMap.createNoTranscriptInfo(Arrays.asList(TableFuncotation.create(
                        Arrays.asList("TESTFIELD1", "TESTADD1"),
                        createFieldValuesFromNameList("E", Arrays.asList("TESTFIELD1", "TESTADD1")),
                        Allele.create("A"),
                        "TestDataSource5", null
                ),
                TableFuncotation.create(
                        Arrays.asList("TESTFIELD2", "TESTADD2"),
                        createFieldValuesFromNameList("F", Arrays.asList("TESTFIELD2", "TESTADD2")),
                        Allele.create("AG"),
                        "TestDataSource5", null
                ),
                TableFuncotation.create(
                        Arrays.asList("TESTFIELD3", "TESTADD3"),
                        createFieldValuesFromNameList("G", Arrays.asList("TESTFIELD3", "TESTADD3")),
                        Allele.create("AT"),
                        "TestDataSource5", null
                )));

        // Attempt to add the Gencode funcotations to the funcotation map.  This should cause an exception.
        funcotationMap.add(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, gencodeFuncotations);
    }

    /**
     * Also tests that {@link FuncotationMap#getGencodeFuncotations(String)} will return an empty list.
     * @param funcotations funcotations to add to the FuncotationMap.  None are GencodeFuncotations.
     */
    @Test(dataProvider = "provideTableFuncotations")
    public void testCreateNoTranscriptInfo(final List<Funcotation> funcotations) {
        final FuncotationMap funcotationMap = FuncotationMap.createNoTranscriptInfo(funcotations);
        Assert.assertEquals(funcotationMap.getTranscriptList().size(), 1);
        final String transcriptId = funcotationMap.getTranscriptList().get(0);
        Assert.assertEquals(transcriptId, FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY);
        Assert.assertEquals(funcotationMap.get(transcriptId), funcotations);
        Assert.assertEquals(funcotationMap.getGencodeFuncotations(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).size(), 0);
    }
}
