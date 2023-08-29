package org.broadinstitute.hellbender.tools.funcotator.compositeoutput;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.AnnotatedIntervalToSegmentVariantContextConverter;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput.SimpleTsvOutputRenderer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CompositeOutputRendererUnitTest extends GATKBaseTest {

    private static final String TEST_SUB_DIR = toolsTestDir + "/funcotator/";
    private static final String SEG_CONFIG_FILE = TEST_SUB_DIR + "test_tsv_output_composite_fields.config";

    @Test
    public void testWriteTwoFiles() throws IOException {
        final File outputFile1 = IOUtils.createTempFile("compositeWriteTwoFiles1", ".seg");
        final File outputFile2 = IOUtils.createTempFile("compositeWriteTwoFiles2", ".seg");

        final SimpleInterval interval = new SimpleInterval("3", 1000000, 2000000);
        final SimpleTsvOutputRenderer renderer1 = SimpleTsvOutputRenderer.createFromFile(outputFile1.toPath(),
                new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), Paths.get(SEG_CONFIG_FILE), "TEST", true);

        final SimpleTsvOutputRenderer renderer2 = SimpleTsvOutputRenderer.createFromFile(outputFile2.toPath(),
                new LinkedHashMap<>(),
                FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("FIELD1", "FIELD2"), Arrays.asList("VAL1", "VAL2")),
                new HashSet<>(Collections.singletonList("foo2")), Paths.get(SEG_CONFIG_FILE), "TEST", true);

        final OutputRenderer compositeOutputRenderer = new CompositeOutputRenderer(Arrays.asList(renderer1, renderer2), "TEST_TOOL");

        final int numLinesInOutput = 2;

        // Create a list of several funcotation maps
        final List<FuncotationMap> funcotationMaps = IntStream.range(0, numLinesInOutput).boxed()
                .map(i -> createSimpleFuncotationMap(Arrays.asList("foo", "foo2"), Arrays.asList("FOO_VALUE" + i, "FOO2_VALUE" + i)))
                .collect(Collectors.toList());

        // Create a corresponding list of segment-like variant contexts.
        final List<VariantContext> segmentVariantContexts = IntStream.range(0, numLinesInOutput).boxed()
                .map(i -> FuncotatorTestUtils.createSimpleVariantContext(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        interval.getContig(), interval.getStart(), interval.getEnd(), "C",
                        AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE_STRING))
                .collect(Collectors.toList());

        // Write'em
        IntStream.range(0, numLinesInOutput).boxed().forEach(i -> compositeOutputRenderer.write(segmentVariantContexts.get(i), funcotationMaps.get(i)));

        // Close the renderer
        compositeOutputRenderer.close();

        // Load up the output files and check contents
        final List<String> fields1 = Arrays.asList("foo1", "foo2", "CONTIG", "END", "START");
        final List<LinkedHashMap<String, String>> gt1 = Arrays.asList(
            FuncotatorUtils.createLinkedHashMapFromLists(fields1,
                    Arrays.asList("FOO_VALUE0", "FOO2_VALUE0", interval.getContig(), String.valueOf(interval.getEnd()),
                            String.valueOf(interval.getStart()))),
            FuncotatorUtils.createLinkedHashMapFromLists(fields1,
                    Arrays.asList("FOO_VALUE1", "FOO2_VALUE1", interval.getContig(), String.valueOf(interval.getEnd()),
                            String.valueOf(interval.getStart())))
        );
        final List<String> fields2 = Arrays.asList("foo1", "CONTIG", "END", "FIELD1", "FIELD2", "START");
        final List<LinkedHashMap<String, String>> gt2 = Arrays.asList(
            FuncotatorUtils.createLinkedHashMapFromLists(fields2,
                    Arrays.asList("FOO_VALUE0", interval.getContig(), String.valueOf(interval.getEnd()), "VAL1", "VAL2",
                            String.valueOf(interval.getStart()))),
            FuncotatorUtils.createLinkedHashMapFromLists(fields2,
                    Arrays.asList("FOO_VALUE1", interval.getContig(), String.valueOf(interval.getEnd()), "VAL1", "VAL2",
                            String.valueOf(interval.getStart())))
        );
        FuncotatorTestUtils.assertTsvFile(outputFile1, gt1);
        FuncotatorTestUtils.assertTsvFile(outputFile2, gt2);
    }

    // Creates a no-transcript funcotation map with one funcotation.
    private FuncotationMap createSimpleFuncotationMap(final List<String> fieldNames, final List<String> fieldValues) {
        return FuncotationMap.createNoTranscriptInfo(
                Collections.singletonList(
                        TableFuncotation.create(fieldNames, fieldValues,
                                AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                                "TEST",
                                null
                        )));
    }
}
