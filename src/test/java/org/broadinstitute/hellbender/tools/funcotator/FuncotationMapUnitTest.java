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
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfCodec;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.FuncotatorReferenceTestUtils;
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
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|4]",
                        Arrays.asList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Arrays.asList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[FOO|hg19|]",
                    Arrays.asList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Arrays.asList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", ""))
                }, {
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_chromosome",
                    "[|hg19|4]",
                    Arrays.asList(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY),
                    Arrays.asList(ImmutableSortedMap.of("Gencode_19_hugoSymbol", "", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_chromosome", "4"))
                }, {
                    "Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_19_hugoSymbol|Gencode_19_ncbiBuild|Gencode_19_annotationTranscript",
                    "[FOO|hg19|txID1]#[BAR|hg38|txID2]",
                    Arrays.asList("txID1", "txID2"),
                    Arrays.asList(
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "FOO", "Gencode_19_ncbiBuild", "hg19", "Gencode_19_annotationTranscript", "txID1"),
                        ImmutableSortedMap.of("Gencode_19_hugoSymbol", "BAR", "Gencode_19_ncbiBuild", "hg38", "Gencode_19_annotationTranscript", "txID2")
                        )
                }, {
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
    public void testCreationFromFuncotationVcfHeaderString(final String headerDescription, final String funcotationValue, final List<String> gtTranscriptIDs, final List<SortedMap<String, String>> gtMaps) {
        final Allele dummyAllele = Allele.create("A");
        final FuncotationMap testMap = FuncotationMap.createAsAllTableFuncotationsFromVcf("Gencode_19_annotationTranscript",
                FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(headerDescription),
                funcotationValue, dummyAllele, "TEST");
        final List<String> transcriptIds = testMap.keyList();
        Assert.assertEquals(new HashSet<>(transcriptIds), new HashSet<>(gtTranscriptIDs));
        for (int i = 0; i < gtMaps.size(); i++){
            final SortedMap<String, String> gtMap = gtMaps.get(i);
            final String transcriptId = transcriptIds.get(i);
            Assert.assertEquals(testMap.get(transcriptId).get(0).getFieldNames().size(), gtMap.keySet().size());
            Assert.assertEquals(testMap.get(transcriptId).size(), 1); // We have one funcotation with X fields.
            Assert.assertTrue(gtMap.keySet().stream().allMatch(k -> gtMap.get(k).equals(testMap.getFieldValue(transcriptId, k, dummyAllele))));
        }
    }

    @DataProvider
    public Object[][] provideGencodeFuncotationCreation() {
        // All of these were chosen not to be IGR.
        return new Object[][] {
                {"chr3", 178916538, 178916538, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref())),
                        pik3caFeatureReader, DS_PIK3CA_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Arrays.asList("ENST00000263967.3")
                },{"chr3", 178916538, 178916538, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref())),
                        pik3caFeatureReader, DS_PIK3CA_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.CANONICAL, Arrays.asList("ENST00000263967.3")
                },{"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Arrays.asList("ENST00000397910.4", "ENST00000380951.5")

                // Next one tests where we would be in a gene with more than one basic transcript, but variant only overlaps one.  And we still ask for all,
                //   but since one is IGR, it will never get added the the FuncotationMap.
                }, {"chr19", 9014550, 9014550, "T", "A", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of(IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.ALL, Arrays.asList("ENST00000397910.4")

                // Next one tests where we would be in a gene with more than one basic transcript, variant overlaps both, but we are in canonical mode.
                },{"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.CANONICAL, Arrays.asList("ENST00000397910.4")

                // Next one tests where we would be in a gene with more than one basic transcript, variant overlaps both, but we are in effect mode.
                },{"chr19", 8994200, 8994200, "G", "C", FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                        ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref())),
                        muc16FeatureReader, DS_MUC16_HG19_GENCODE_FASTA,
                        TranscriptSelectionMode.BEST_EFFECT, Arrays.asList("ENST00000397910.4")
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

        final ReferenceContext referenceContext = new ReferenceContext(referenceDataSource, variantInterval );

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = featureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }
        final List<Feature> featureList = Collections.singletonList(gtfFeatureIterator.next());

        final String gencode_test = "GENCODE_TEST";
        final GencodeFuncotationFactory gencodeFactory = new GencodeFuncotationFactory(Paths.get(transcriptFastaFile),
        "TEST", gencode_test, transcriptSelectionMode, new HashSet<>(), new LinkedHashMap<>(),
                true);

        final List<GencodeFuncotation> gencodeFuncotations =
                gencodeFactory.createFuncotations(variantContext, referenceContext, Collections.singletonMap(gencode_test, featureList)).stream()
                .map(f -> (GencodeFuncotation) f).collect(Collectors.toList());

        final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(gencodeFuncotations);

        Assert.assertEquals(funcotationMap.keyList(), gtTranscripts);
        Assert.assertTrue(funcotationMap.keyList().stream().allMatch(k -> funcotationMap.get(k).size() == 1));
        Assert.assertTrue(funcotationMap.keyList().stream()
                .noneMatch(k -> ((GencodeFuncotation) funcotationMap.get(k).get(0)).getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR) ));
        Assert.assertTrue(funcotationMap.keyList().stream()
                .noneMatch(k -> ((GencodeFuncotation) funcotationMap.get(k).get(0)).getVariantClassification().equals(GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE) ));
    }
}
