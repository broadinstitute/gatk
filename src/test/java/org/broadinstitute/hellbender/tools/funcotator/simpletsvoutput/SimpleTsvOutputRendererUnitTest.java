package org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.funcotator.AnnotatedIntervalToSegmentVariantContextConverter;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorConstants;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class SimpleTsvOutputRendererUnitTest extends GATKBaseTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "/funcotator/";
    private static final String SEG_CONFIG_FILE = TEST_SUB_DIR + "test_tsv_output.config";
    private static final String SEG_RESOURCE_FILE = "org/broadinstitute/hellbender/tools/funcotator/simple_funcotator_seg_file.config";
    private static final String FUNCOTATION_FIELD_1 = "TO_BE_ALIASED1";

    @DataProvider
    public Object[][] provideForSimpleSegFileWriting() {
        return new Object[][]{
                {
                        FuncotatorTestUtils.createDummySegmentVariantContext()
                }
        };
    }

    @Test(dataProvider = "provideForSimpleSegFileWriting")
    public void testSimpleSegFileWriting(final VariantContext segVC) throws IOException {
        final File outputFile = IOUtils.createTempFile("simpleSegFileWriting", ".seg");

        final SimpleTsvOutputRenderer renderer = SimpleTsvOutputRenderer.createFromFile(outputFile.toPath(),
                new LinkedHashMap<>(),
                new LinkedHashMap<>(), new HashSet<>(), Paths.get(SEG_CONFIG_FILE), "TEST", true);

        final List<String> fieldNames = Arrays.asList("Gencode_19_genes", "foo1", "foobar2", "TEST3", "foo3!!");
        final List<String> fieldValues = Arrays.asList("GENE1,GENE2", "stuff1", "stuff2", "stuff3", "stuff4");
        final FuncotationMap funcotationMap = FuncotationMap.createNoTranscriptInfo(
                Collections.singletonList(
                        TableFuncotation.create(fieldNames, fieldValues, segVC.getAlternateAllele(0),
                                "TEST", createDummySegmentFuncotationMetadata())

                ));
        renderer.write(segVC, funcotationMap);
        renderer.close();

        // Create the entire list of OUTPUT field names in order.
        //   Note that the locatable fields are always emitted in alphabetical order if no aliases are found.
        final List<String> outputFieldNames = Arrays.asList("genes", "foo1", "foo2", "foo3", "foo3!!", "CONTIG", "END", "START");
        final List<String> outputFieldValues = new ArrayList<>();
        outputFieldValues.addAll(fieldValues);
        outputFieldValues.addAll(Arrays.asList(segVC.getContig(), String.valueOf(segVC.getEnd()), String.valueOf(segVC.getStart())));
        final List<LinkedHashMap<String, String>> gtOutputRecords = Collections.singletonList(
                FuncotatorUtils.createLinkedHashMapFromLists(outputFieldNames, outputFieldValues));

        FuncotatorTestUtils.assertTsvFile(outputFile, gtOutputRecords);
    }

    private static FuncotationMetadata createDummySegmentFuncotationMetadata() {
        return VcfFuncotationMetadata.create(
                Arrays.asList(
                        new VCFInfoHeaderLine("Gencode_19_genes", 1, VCFHeaderLineType.String, "The genes overlapping the segment."),
                        new VCFInfoHeaderLine("foo1", 1, VCFHeaderLineType.String, "foo1"),
                        new VCFInfoHeaderLine("foobar2", 1, VCFHeaderLineType.String, "foobar2 (an alias relative to the config file)"),
                        new VCFInfoHeaderLine("TEST3", 1, VCFHeaderLineType.String, "Note that this has no spaces"),
                        new VCFInfoHeaderLine("foo3!!", 1, VCFHeaderLineType.String, "special character....")
                )
        );
    }

    @DataProvider()
    public Object[][] provideConfigFiles() throws IOException {

        // test_tsv_output.seg
        LinkedHashMap<String, List<String>> gt1 = new LinkedHashMap<>();
        gt1.put("genes", Collections.singletonList("Gencode_19_genes"));
        gt1.put("foo1", Collections.emptyList());
        gt1.put("foo2", Arrays.asList("TEST2", "foobar2"));
        gt1.put("foo3", Arrays.asList("TEST3", "foobar2"));
        gt1.put("foo3!!", Collections.emptyList());

        // Actual resource for funcotator seg files.
        LinkedHashMap<String, List<String>> gt2 = new LinkedHashMap<>();
        gt2.put("alt_allele", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("Gencode_19_alt_allele,Gencode_27_alt_allele,Gencode_28_alt_allele", ",")));
        gt2.put("end_gene", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("Gencode_19_end_gene,Gencode_27_end_gene,Gencode_28_end_gene", ",")));
        gt2.put("end", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("END,End,End_Position,end_position,chromEnd,segment_end,End_position,target_end,stop,Stop,Position,position,pos,POS,segment_end", ",")));
        gt2.put("start_gene", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("Gencode_19_start_gene,Gencode_27_start_gene,Gencode_28_start_gene", ",")));
        gt2.put("Segment_Mean", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("MEAN_LOG2_COPY_RATIO", ",")));
        gt2.put("genes", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("Gencode_19_genes,Gencode_27_genes,Gencode_28_genes", ",")));
        gt2.put("Sample", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("sample,sample_id", ",")));
        gt2.put("start", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("START,Start,Start_Position,start_position,chromStart,segment_start,Start_position,target_start,Position,position,pos,POS,segment_start", ",")));
        gt2.put("chr", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("CONTIG,contig,Chromosome,chrom,chromosome,Chrom,seqname,seqnames,CHROM,target_contig,segment_contig", ",")));
        gt2.put("build", Collections.emptyList());
        gt2.put("Num_Probes", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("NUM_POINTS_COPY_RATIO", ",")));
        gt2.put("start_exon", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("Gencode_19_start_exon,Gencode_27_start_exon,Gencode_28_start_exon", ",")));
        gt2.put("end_exon", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("Gencode_19_end_exon,Gencode_27_end_exon,Gencode_28_end_exon", ",")));
        gt2.put("ref_allele", Arrays.asList(SimpleTsvOutputRenderer.splitAndTrim("Gencode_19_ref_allele,Gencode_27_ref_allele,Gencode_28_ref_allele", ",")));
        gt2.put("Segment_Call", Collections.emptyList());

        return new Object[][]{
                {Paths.get(SEG_CONFIG_FILE), gt1},
                {Resource.getResourceContentsAsFile(SEG_RESOURCE_FILE).toPath(), gt2}
        };
    }

    @Test(dataProvider = "provideConfigFiles")
    public void testLoadingConfigFile(Path configFile, LinkedHashMap<String, List<String>> gt) {
        final LinkedHashMap<String, List<String>> guess = SimpleTsvOutputRenderer.createColumnNameToAliasesMap(configFile);
        FuncotatorTestUtils.assertLinkedHashMapsEqual(guess, gt);
    }

    @DataProvider
    public Object[][] provideForAliasing() {
        return new Object[][]{
                {
                        // No aliases, so remaining fields will be in alphabetical order
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), true, new LinkedHashMap<>(),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("CONTIG", "END", "START", FUNCOTATION_FIELD_1),
                        Arrays.asList("CONTIG", "END", "START", FUNCOTATION_FIELD_1),
                        Arrays.asList("3", "3200000", "2100000", "value1")
                },
                {
                        // Simple alias of foo:FUNCOTATION_FIELD_1
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Collections.singletonList(FUNCOTATION_FIELD_1))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "CONTIG", "END", "START"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "CONTIG", "END", "START"),
                        Arrays.asList("value1", "3", "3200000", "2100000")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "CONTIG", "END", "START"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "CONTIG", "END", "START"),
                        Arrays.asList("value1", "3", "3200000", "2100000")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1, but we exclude the "CONTIG" and "START" fields
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "START")), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "END"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "END"),
                        Arrays.asList("value1", "3200000")
                },
                {
                        //Test #5 (index 4)  In this case "foo" should not appear in the output.
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1, but we exclude the "CONTIG", "foo", and "START" fields
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("foo", "CONTIG", "START")), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Collections.singletonList("END"),
                        Collections.singletonList("END"),
                        Collections.singletonList("3200000")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1 with default.
                        FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "foo2", "foo3"), Arrays.asList("BAD_DEFAULT", "val2", "val3")), new LinkedHashMap<>(), new HashSet<>(), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "CONTIG", "END", "START", "foo2", "foo3"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "CONTIG", "END", "START", "foo2", "foo3"),
                        Arrays.asList("value1", "3", "3200000", "2100000", "val2", "val3")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1 with overrides.
                        new LinkedHashMap<>(), FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "foo2", "foo3"), Arrays.asList("OVERRIDE", "val2", "val3")), new HashSet<>(), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "CONTIG", "END", "START", "foo2", "foo3"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "CONTIG", "END", "START", "foo2", "foo3"),
                        Arrays.asList("OVERRIDE", "3", "3200000", "2100000", "val2", "val3")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1 with overrides and excluding some of the overrides.
                        new LinkedHashMap<>(), FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "foo2", "foo3"), Arrays.asList("OVERRIDE", "val2", "val3")), new HashSet<>(Arrays.asList("foo2", "foo3")), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "CONTIG", "END", "START"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "CONTIG", "END", "START"),
                        Arrays.asList("OVERRIDE", "3", "3200000", "2100000")
                },

                // Same tests as above, but this time disallowing funcotation fields to percolate into the output
                //  Some examples from above (rightfully) generate an error, so those have been moved to a different test.

                {
                        // Simple alias of foo:FUNCOTATION_FIELD_1
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), false,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Collections.singletonList(FUNCOTATION_FIELD_1))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo"),
                        Arrays.asList(FUNCOTATION_FIELD_1),
                        Arrays.asList("value1")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), false,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo"),
                        Arrays.asList(FUNCOTATION_FIELD_1),
                        Arrays.asList("value1")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1, but we exclude the "CONTIG" and "START" fields
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "START")), false,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo"),
                        Arrays.asList(FUNCOTATION_FIELD_1),
                        Arrays.asList("value1")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1 with default.
                        FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "foo2", "foo3"), Arrays.asList("BAD_DEFAULT", "val2", "val3")), new LinkedHashMap<>(), new HashSet<>(), false,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "foo2", "foo3"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "foo2", "foo3"),
                        Arrays.asList("value1", "val2", "val3")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1 with overrides.
                        new LinkedHashMap<>(), FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "foo2", "foo3"), Arrays.asList("OVERRIDE", "val2", "val3")), new HashSet<>(), false,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "foo2", "foo3"),
                        Arrays.asList(FUNCOTATION_FIELD_1, "foo2", "foo3"),
                        Arrays.asList("OVERRIDE", "val2", "val3")
                },
                {
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1 with overrides and excluding some of the overrides.
                        new LinkedHashMap<>(), FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "foo2", "foo3"), Arrays.asList("OVERRIDE", "val2", "val3")), new HashSet<>(Arrays.asList("foo2", "foo3")), false,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo"),
                        Arrays.asList(FUNCOTATION_FIELD_1),
                        Arrays.asList("OVERRIDE")
                },
                {
                        // In this case "foo" should appear as __UNKNOWN__in the output.
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1, but we exclude the "CONTIG", "foo", and "START" fields
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "START")), true,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Arrays.asList("foo", "END", "TO_BE_ALIASED1"),
                        Arrays.asList("", "END", "TO_BE_ALIASED1"),
                        Arrays.asList(FuncotatorConstants.UNKNOWN_VALUE_STRING, "3200000", "value1")
                },
                {
                        // In this case "foo" should appear as __UNKNOWN__in the output, but since we do not allow funcotations to pass, it should be the only field
                        // Simple alias of foo:...,FUNCOTATION_FIELD_1, but we exclude the "CONTIG", "foo", and "START" fields
                        new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "START")), false,
                        FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"))),
                        FuncotatorTestUtils.createDummySegmentVariantContext(),
                        createSimpleFuncotationMap(),
                        Collections.singletonList("foo"),
                        Collections.singletonList(""),
                        Collections.singletonList(FuncotatorConstants.UNKNOWN_VALUE_STRING)
                },{
                    // In this case "foo" should not appear (excluded).  But "Bar" should appear as __UNKNOWN__
                    new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "START", "foo")), true,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList("DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar", "END", "TO_BE_ALIASED1"),
                    Arrays.asList("", "END", "TO_BE_ALIASED1"),
                    Arrays.asList(FuncotatorConstants.UNKNOWN_VALUE_STRING, "3200000", "value1")
                },{
                    // In this case "foo" should not appear (excluded).  But "Bar" should appear as __UNKNOWN__
                    new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "START", "foo")), false,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList("DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar"),
                    Arrays.asList(""),
                    Arrays.asList(FuncotatorConstants.UNKNOWN_VALUE_STRING)
                },{
                    // In this case "foo" should not appear (excluded).  But "Bar" should appear with the alias along with the remaining funcotation field.
                    new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "START", "foo")), true,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList(FUNCOTATION_FIELD_1, "DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar", "END"),
                    Arrays.asList(FUNCOTATION_FIELD_1, "END"),
                    Arrays.asList("value1", "3200000")
                },{
                    // In this case "foo" should not appear (excluded).  But "Bar" should appear with the alias but no other funcotation fields.
                    new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG","START", "foo")), false,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList(FUNCOTATION_FIELD_1, "DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar"),
                    Arrays.asList(FUNCOTATION_FIELD_1),
                    Arrays.asList("value1")
                },{
                    // In this case "foo" should not appear (excluded).  But "Bar" should appear with the alias with two other funcotation fields sorted properly.
                    new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "foo")), true,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList(FUNCOTATION_FIELD_1, "DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar", "END", "START"),
                    Arrays.asList(FUNCOTATION_FIELD_1, "END", "START"),
                    Arrays.asList("value1", "3200000", "2100000")
                },{
                    // Test sorting when we have a mixture of default, funcotation fields, and excludes.
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar", "foo3"), Arrays.asList("DEFAULT", "val2", "val3")),
                    new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "foo")), true,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList(FUNCOTATION_FIELD_1, "DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar", "END", "START", "foo3"),
                    Arrays.asList(FUNCOTATION_FIELD_1, "END", "START", "foo3"),
                    Arrays.asList("value1", "3200000", "2100000", "val3")
                },{
                    // Test sorting of f1, f10, and f2 when we have a mixture of default, funcotation fields, and excludes.
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar", "foo3", "foo1", "foo10"), Arrays.asList("DEFAULT", "val2", "val3", "val1", "val10")),
                    new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "foo")), true,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList(FUNCOTATION_FIELD_1, "DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar", "END", "START", "foo1", "foo3", "foo10"),
                    Arrays.asList(FUNCOTATION_FIELD_1, "END", "START", "foo1","foo3", "foo10"),
                    Arrays.asList("value1", "3200000", "2100000", "val1", "val3", "val10")
                },{
                    // Test sorting of f1, f10, and f2 when we have a mixture of default, override, funcotation fields, and excludes.
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar",  "foo1", "foo10"), Arrays.asList("DEFAULT", "val2", "val1", "val10")),
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo3"), Arrays.asList("val3")), new HashSet<>(Arrays.asList("CONTIG", "foo")), true,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList(FUNCOTATION_FIELD_1, "DUMMY3", "DUMMY4", "DUMMY5"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("bar", "END", "START", "foo1", "foo3", "foo10"),
                    Arrays.asList(FUNCOTATION_FIELD_1, "END", "START", "foo1","foo3", "foo10"),
                    Arrays.asList("value1", "3200000", "2100000", "val1", "val3", "val10")
                },{
                    // Test that the order is preserved on explicit columns
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "bar", "foo3", "foo1", "foo10"), Arrays.asList("DEFAULT", "val2", "val3", "val1", "val10")),
                    new LinkedHashMap<>(), new HashSet<>(Arrays.asList("CONTIG", "foo")), true,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("foo", "foo3", "foo1", "foo10", "bar"), Arrays.asList(Arrays.asList("DUMMY1", "DUMMY2", "DUMMY3"), Arrays.asList("DUMMY03", "DUMMY04", "DUMMY05"), Arrays.asList("DUMMY13", "DUMMY14", "DUMMY15"), Arrays.asList("DUMMY23", "DUMMY24", "DUMMY25"), Arrays.asList(FUNCOTATION_FIELD_1, "DUMMY33", "DUMMY34", "DUMMY35"))),
                    FuncotatorTestUtils.createDummySegmentVariantContext(),
                    createSimpleFuncotationMap(),
                    Arrays.asList("foo3", "foo1", "foo10", "bar", "END", "START"),
                    Arrays.asList("", "", "", FUNCOTATION_FIELD_1, "END", "START"),
                    Arrays.asList("val3", "val1", "val10", "value1", "3200000", "2100000")
                }

        };
    }

    private FuncotationMap createSimpleFuncotationMap() {
        return FuncotationMap.createNoTranscriptInfo(
                Collections.singletonList(
                        TableFuncotation.create(Collections.singletonList(FUNCOTATION_FIELD_1), Collections.singletonList("value1"),
                                AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE,
                                "TEST",
                                VcfFuncotationMetadata.create(
                                        Collections.singletonList(
                                                new VCFInfoHeaderLine(FUNCOTATION_FIELD_1, 1, VCFHeaderLineType.String, "Unknown")))
                        )));
    }

    // This tests both getColumnNameToFuncotationFieldMap, createColumnNameToValueMap, and the output file itself.
    //  So, this is a whitebox test.
    @Test(dataProvider = "provideForAliasing", description = "Test aliasing both when funcotation fields get to the output and when that is disallowed.")
    public void testAliasing(final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                             final Set<String> excludedOutputFields, final boolean isWriteFuncotationFieldsNotInConfig,
                             final LinkedHashMap<String, List<String>> columnNameToAliasesMap,
                             final VariantContext segVC, final FuncotationMap funcotationMap,
                             final List<String> gtAliasKeys, final List<String> gtAliasFuncotationFields, final List<String> gtFinalValues) throws IOException {
        final File outputFile = IOUtils.createTempFile("testAliasing", ".seg");
        final SimpleTsvOutputRenderer simpleTsvOutputRenderer = new SimpleTsvOutputRenderer(outputFile.toPath(),
                unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, excludedOutputFields,
                columnNameToAliasesMap, "TESTING_VERSION", isWriteFuncotationFieldsNotInConfig);

        // You must write one record since SimpleTsvOutputRenderer lazy loads the writer.
        simpleTsvOutputRenderer.write(segVC, funcotationMap);

        // Test that all columns are in the output and the proper aliases are set up.  A whitebox test.
        final LinkedHashMap<String, String> columnNameToFuncotationFieldMap = simpleTsvOutputRenderer.getColumnNameToFuncotationFieldMap();
        Assert.assertEquals(columnNameToFuncotationFieldMap.keySet(), new LinkedHashSet<>(gtAliasKeys), "Alias key did not correspond to ground truth.");
        Assert.assertEquals(Lists.newArrayList(columnNameToFuncotationFieldMap.values()), gtAliasFuncotationFields, "Alias value did not correspond to ground truth.");

        // Check the aliasing to the values.  A whitebox test
        // This next line assumes that all test funcotation maps have the same tx ID (None) and alternate allele (copy neutral)
        final LinkedHashMap<String, String> columnNameToFuncotationValueMap = SimpleTsvOutputRenderer.createColumnNameToValueMap(columnNameToFuncotationFieldMap,
                funcotationMap, FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE, unaccountedForDefaultAnnotations,
                unaccountedForOverrideAnnotations, excludedOutputFields);
        Assert.assertEquals(columnNameToFuncotationValueMap.values(), gtFinalValues);

        // Get the actual values in the columns
        simpleTsvOutputRenderer.close();

        // Check the output file contents
        final List<LinkedHashMap<String, String>> gtOutputRecords = Collections.singletonList(
                FuncotatorUtils.createLinkedHashMapFromLists(gtAliasKeys, gtFinalValues));
        FuncotatorTestUtils.assertTsvFile(outputFile, gtOutputRecords);
    }


    @DataProvider
    public Object[][] provideForAliasingErrors() {
        return new Object[][]{
                // The next three entries test cases where the output file would be empty, which triggers an error.
            {
                // No aliases, so remaining fields will be in alphabetical order.  Test that VC attributes are
                //   dropped, so we should still get an error.  We are not allowing columns not in the alias list
                // (or default or override).
                new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), false, new LinkedHashMap<>(),
                FuncotatorTestUtils.createDummySegmentVariantContextWithAttributes(ImmutableMap.of("ATTR1", "VAL_ATTR1")),
                createSimpleFuncotationMap()
            },
            {
                    // Excluded the only column in the alias map and we are not rendering columns not in the alias map.
                new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(Arrays.asList("foo", "CONTIG", "START")), false,
                FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("foo"), Collections.singletonList(Arrays.asList("DUMMY1", "DUMMY2", FUNCOTATION_FIELD_1, "DUMMY3"))),
                FuncotatorTestUtils.createDummySegmentVariantContext(),
                createSimpleFuncotationMap()
            },{
                // No aliases, so remaining fields will be in alphabetical order.  We are not allowing columns
                // not in the alias list (or default or override).
                new LinkedHashMap<>(), new LinkedHashMap<>(), new HashSet<>(), false, new LinkedHashMap<>(),
                FuncotatorTestUtils.createDummySegmentVariantContext(),
                createSimpleFuncotationMap()
            }
        };
    }

    @Test(dataProvider = "provideForAliasingErrors", description = "Test aliasing configurations that are illegal.", expectedExceptions = IllegalArgumentException.class)
    public void testAliasingErrors(final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                                   final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                                   final Set<String> excludedOutputFields, final boolean isWriteFuncotationFieldsNotInConfig,
                                   final LinkedHashMap<String, List<String>> columnNameToAliasesMap,
                                   final VariantContext segVC, final FuncotationMap funcotationMap) {
        final File outputFile = IOUtils.createTempFile("testAliasing", ".seg");
        final SimpleTsvOutputRenderer simpleTsvOutputRenderer = new SimpleTsvOutputRenderer(outputFile.toPath(),
                unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, excludedOutputFields,
                columnNameToAliasesMap, "TESTING_VERSION", isWriteFuncotationFieldsNotInConfig);

        // You must write one record since SimpleTsvOutputRenderer lazy loads the writer.
        simpleTsvOutputRenderer.write(segVC, funcotationMap);

    }
}
