package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.collect.Lists;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.builder.fluent.Configurations;
import org.apache.commons.configuration2.ex.ConfigurationException;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.funcotator.genelistoutput.GeneListOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput.SimpleTsvOutputRenderer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.collections.Sets;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils.createLinkedHashMapListTableReader;

public class FuncotateSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "/funcotator/";
    private static final String SIMPLE_TEST_FILE = TEST_SUB_DIR + "simple.seg";
    private static final String SIMPLE_TEST_CNTN4_FILE = TEST_SUB_DIR + "simple_cntn4_overlap.seg";
    private static final String TEST_GATK_FILE_B37 = TEST_SUB_DIR + "SM-74NF5.called.seg";
    private static final String TEST_GATK_EMPTY_FILE_B37 = TEST_SUB_DIR + "empty_b37.seg";
    private static final String DS_PIK3CA_DIR  = largeFileTestDir + "funcotator" + File.separator + "small_ds_pik3ca" + File.separator;
    // This has transcripts with multiple gene names...
    private static final String DS_CNTN4_DIR  = toolsTestDir + "funcotator" + File.separator + "small_cntn4_ds" + File.separator;
    private static final String SEG_RESOURCE_FILE = "org/broadinstitute/hellbender/tools/funcotator/simple_funcotator_seg_file.config";
    private static final String GENE_LIST_RESOURCE_FILE = "org/broadinstitute/hellbender/tools/funcotator/gene_list_output.config";
    private static final String funcotatorLargeFileTestDir = largeFileTestDir + "funcotator" + File.separator;
    private static final String CALLED_CR_GATK_HG38_FILE = funcotatorLargeFileTestDir + "calledcopyratio_gatk_output_hg38_integration_test.seg";
    private static final String MODEL_SEGMENTS_FINAL_GATK_HG38_FILE = funcotatorLargeFileTestDir + "modelFinal_gatk_output_hg38_integration_test.seg";

    @Test(description = "Test simple tsv and gene list when no segment overlaps a gene.")
    public void testSimpleNoOverlap() throws IOException {
        final File outputFile = IOUtils.createTempFile("funcotatesegs_simple", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addRaw("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.addRaw(SIMPLE_TEST_FILE);
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.addRaw(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.addRaw("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.addRaw(b37Reference);
        arguments.addRaw("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.addRaw(outputFile.getAbsolutePath());
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.addRaw(BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19);
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);

        runCommandLine(arguments);

        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(collection.getRecords().size(), 3);
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.hasAnnotation("genes")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("genes").equals("")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("start_gene").equals("")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("end_gene").equals("")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("start_exon").equals("")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("end_exon").equals("")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("ref_allele").equals("")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("alt_allele").equals("")));

        // Check the gene list file.  It should only have headers.
        FuncotatorTestUtils.assertTsvFile(new File(outputFile.getAbsolutePath() + FuncotatorEngine.GENE_LIST_FILE_SUFFIX),
                Collections.emptyList(),
                getGeneListOutputFields()
        );
    }

    private static ArrayList<String> getGeneListOutputFields() throws IOException {
        return Lists.newArrayList(
                SimpleTsvOutputRenderer.createColumnNameToAliasesMap(Resource.getResourceContentsAsFile(GeneListOutputRenderer.CONFIG_RESOURCE).toPath()).keySet());
    }

    /**
     * Very dependent on the data in "simple_cntn4_overlap.seg"
     */
    @DataProvider
    public Object[][] cntn4GroundTruth() {
        return new Object[][] {
                {
                            // Locatable info
                            Arrays.asList(new SimpleInterval("3", 2000000, 2500000),
                                    new SimpleInterval("3", 3000000,3500000),
                                    new SimpleInterval("3",3500001,3900000)),

                            // genes
                            Arrays.asList("CNTN4,CNTN4-AS2", "CNTN4,CNTN4-AS1", ""),
                            //start_gene
                            Arrays.asList("", "CNTN4", ""),
                            //end_gene
                            Arrays.asList("CNTN4", "", ""),
                            // ref_allele (always blank)
                            Arrays.asList("", "", ""),
                            // alt_allele (always blank)
                            Arrays.asList("", "", ""),
                            // Call
                            Arrays.asList("0", "-", "+"),
                            // Segment_Mean
                            Arrays.asList("0.037099", "0.001748", "0.501748"),
                            // Num_Probes
                            Arrays.asList("2000", "3000", "4000"),
                            // Sample
                            Arrays.asList("SAMPLE1", "SAMPLE1", "SAMPLE1")
                }
        };
    }

    @Test(dataProvider = "cntn4GroundTruth", description = "Test simple tsv and gene list when multiple segments overlap a gene.")
    public void testSimpleMultipleGenesOverlap(List<Locatable> gtInterval, List<String> gtGenesValues, List<String> gtStartGeneValues, List<String> gtEndGeneValues,
                                               List<String> gtRefAlleles, List<String> gtAltAlleles, List<String> gtCalls,
                                               List<String> gtSegmentMeans, List<String> gtNumProbes, List<String> gtSamples)
            throws IOException {
        final File outputFile = IOUtils.createTempFile("funcotatesegs_simple_cntn4", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addRaw("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.addRaw(SIMPLE_TEST_CNTN4_FILE);
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.addRaw(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.addRaw("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.addRaw(b37Reference);
        arguments.addRaw("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.addRaw(outputFile.getAbsolutePath());
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.addRaw(BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19);
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_CNTN4_DIR);

        runCommandLine(arguments);

        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        final List<AnnotatedInterval> segments = collection.getRecords();
        Assert.assertEquals(segments.size(), 3);

        final List<String> testGenesValues = segments.stream().map(r -> r.getAnnotationValue("genes")).collect(Collectors.toList());
        Assert.assertEquals(testGenesValues, gtGenesValues);
        final List<String> testStartGeneValues = segments.stream().map(r -> r.getAnnotationValue("start_gene")).collect(Collectors.toList());
        Assert.assertEquals(testStartGeneValues, gtStartGeneValues);
        final List<String> testEndGeneValues = segments.stream().map(r -> r.getAnnotationValue("end_gene")).collect(Collectors.toList());
        Assert.assertEquals(testEndGeneValues, gtEndGeneValues);
        final List<String> testRefAlleleValues = segments.stream().map(r -> r.getAnnotationValue("ref_allele")).collect(Collectors.toList());
        Assert.assertEquals(testRefAlleleValues, gtRefAlleles);
        final List<String> testAltAlleleValues = segments.stream().map(r -> r.getAnnotationValue("alt_allele")).collect(Collectors.toList());
        Assert.assertEquals(testAltAlleleValues, gtAltAlleles);

        Assert.assertEquals(segments.stream().map(r -> r.getAnnotationValue("Num_Probes")).collect(Collectors.toList()), gtNumProbes);
        Assert.assertEquals(segments.stream().map(r -> r.getAnnotationValue("Segment_Mean")).collect(Collectors.toList()), gtSegmentMeans);
        Assert.assertEquals(segments.stream().map(r -> r.getAnnotationValue("Segment_Call")).collect(Collectors.toList()), gtCalls);
        Assert.assertEquals(segments.stream().map(r -> r.getAnnotationValue("Sample")).collect(Collectors.toList()), gtSamples);
        Assert.assertEquals(segments.stream().map(AnnotatedInterval::getInterval).collect(Collectors.toList()), gtInterval);

        // Check gene list.  All values except the first two, since these are per-segment and tend to get repeated in a gene list.
        final List<String> orderedSegmentValues0 = Arrays.asList(segments.get(0).getContig(), String.valueOf(segments.get(0).getStart()), "", "", String.valueOf(segments.get(0).getEnd()),
                "CNTN4", "1-", segments.get(0).getAnnotationValue("Num_Probes"), segments.get(0).getAnnotationValue("Segment_Mean"), segments.get(0).getAnnotationValue("Segment_Call"),
                segments.get(0).getAnnotationValue("build"), segments.get(0).getAnnotationValue("Sample"));
        final List<String> orderedSegmentValues1 = Arrays.asList( segments.get(1).getContig(), String.valueOf(segments.get(1).getStart()), "CNTN4", "12+", String.valueOf(segments.get(1).getEnd()),
                "", "", segments.get(1).getAnnotationValue("Num_Probes"), segments.get(1).getAnnotationValue("Segment_Mean"), segments.get(1).getAnnotationValue("Segment_Call"),
                segments.get(1).getAnnotationValue("build"), segments.get(1).getAnnotationValue("Sample"));
        FuncotatorTestUtils.assertTsvFile(new File(outputFile.getAbsolutePath() + FuncotatorEngine.GENE_LIST_FILE_SUFFIX),
                /* See gene_list_output.config, but as of this writing:
                gene
                exon
                segment_contig
                segment_start
                segment_start_gene
                segment_start_exon
                segment_end
                segment_end_gene
                segment_end_exon
                segment_num_probes
                segment_mean
                segment_call
                build
                sample
                 */

                // Accotding to UCSC genome browser hg19 CNTN4 end of 2500000 is 1- and the 3000000 start is 12+
                //  Note  that the third segment (index 2) did not overlap any of the 3 genes in the truncated datasource
                Arrays.asList(

                        FuncotatorUtils.createLinkedHashMapFromLists(getGeneListOutputFields(), Stream.concat(Stream.of("CNTN4", "1-"), orderedSegmentValues0.stream()).collect(Collectors.toList())),
                        FuncotatorUtils.createLinkedHashMapFromLists(getGeneListOutputFields(), Stream.concat(Stream.of("CNTN4", "12+"), orderedSegmentValues1.stream()).collect(Collectors.toList())),
                        FuncotatorUtils.createLinkedHashMapFromLists(getGeneListOutputFields(), Stream.concat(Stream.of("CNTN4-AS1", ""), orderedSegmentValues1.stream()).collect(Collectors.toList())),
                        FuncotatorUtils.createLinkedHashMapFromLists(getGeneListOutputFields(), Stream.concat(Stream.of("CNTN4-AS2", ""), orderedSegmentValues0.stream()).collect(Collectors.toList()))
                ), getGeneListOutputFields()
        );
    }

    @Test(description = "Test simple tsv and gene list for GATK seg file where all genes completely overlap one segment.")
    public void testGatkCalledSegmentFile() throws IOException {
        final File outputFile = IOUtils.createTempFile("funcotatesegs_gatk_called", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addRaw("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.addRaw(TEST_GATK_FILE_B37);
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.addRaw(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.addRaw("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.addRaw(b37Reference);
        arguments.addRaw("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.addRaw(outputFile.getAbsolutePath());
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.addRaw("hg19");
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_CNTN4_DIR);

        runCommandLine(arguments);

        final AnnotatedIntervalCollection outputSegmentCollection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        final List<AnnotatedInterval> segments = outputSegmentCollection.getRecords();
        Assert.assertEquals(segments.size(), 404);

        // This tests that all should-be-populated fields have a value
        final List<String> allSegmentFieldsThatShouldBePopulated = Arrays.asList("Num_Probes", "Segment_Mean", "Segment_Call");
        final List<String> allSegmentFieldsThatShouldBeEmpty = Arrays.asList("ref_allele", "alt_allele");
        allSegmentFieldsThatShouldBePopulated.forEach(f ->
            Assert.assertTrue(segments.stream().noneMatch(r -> StringUtils.isEmpty(r.getAnnotationValue(f))), f + " was not populated and it should be.")
        );
        allSegmentFieldsThatShouldBeEmpty.forEach(f ->
            Assert.assertTrue(segments.stream().allMatch(r -> StringUtils.isEmpty(r.getAnnotationValue(f))), f + " was populated and it should not be.")
        );

        // Left is the output name.  Right is the input name.  These fields should have the exact same values in the input and output
        final List<Pair<String,String>> testPairedFieldNames = Arrays.asList(
          Pair.of("Num_Probes", "NUM_POINTS_COPY_RATIO"),
                Pair.of("Segment_Mean", "MEAN_LOG2_COPY_RATIO"),
                Pair.of("Segment_Call", "CALL"),
                Pair.of("chr", "CONTIG"),
                Pair.of("start", "START"),
                Pair.of("end", "END")
        );

        for (final Pair<String,String> pairedFields : testPairedFieldNames) {
            final AnnotatedIntervalCollection inputSegmentCollection = AnnotatedIntervalCollection.create(Paths.get(TEST_GATK_FILE_B37), null);
            Assert.assertEquals(segments.stream().map(r -> r.getAnnotationValue(pairedFields.getLeft())).collect(Collectors.toList()),
                    inputSegmentCollection.getRecords().stream().map(r -> r.getAnnotationValue(pairedFields.getRight())).collect(Collectors.toList()));
        }
        // We should have one segment that has all CNTN4 genes (index 67 segment)
        final AnnotatedInterval cntn4Segment = segments.stream().filter(r -> !StringUtils.isEmpty(r.getAnnotationValue("genes"))).findFirst().get();
        Assert.assertEquals(cntn4Segment.getAnnotationValue("genes"), "CNTN4,CNTN4-AS1,CNTN4-AS2");

        final AnnotatedInterval segmentWithAllGenes = segments.get(67);
        final List<String> orderedSegmentValues0 = Arrays.asList(segmentWithAllGenes.getContig(), String.valueOf(segmentWithAllGenes.getStart()), "", "", String.valueOf(segmentWithAllGenes.getEnd()),
                "", "", segmentWithAllGenes.getAnnotationValue("Num_Probes"), segmentWithAllGenes.getAnnotationValue("Segment_Mean"), segmentWithAllGenes.getAnnotationValue("Segment_Call"),
                segmentWithAllGenes.getAnnotationValue("build"), segmentWithAllGenes.getAnnotationValue("Sample"));

        // Check the gene list.
        FuncotatorTestUtils.assertTsvFile(new File(outputFile.getAbsolutePath() + FuncotatorEngine.GENE_LIST_FILE_SUFFIX),
                Arrays.asList(
                        FuncotatorUtils.createLinkedHashMapFromLists(getGeneListOutputFields(), Stream.concat(Stream.of("CNTN4", ""), orderedSegmentValues0.stream()).collect(Collectors.toList())),
                        FuncotatorUtils.createLinkedHashMapFromLists(getGeneListOutputFields(), Stream.concat(Stream.of("CNTN4-AS1", ""), orderedSegmentValues0.stream()).collect(Collectors.toList())),
                        FuncotatorUtils.createLinkedHashMapFromLists(getGeneListOutputFields(), Stream.concat(Stream.of("CNTN4-AS2", ""), orderedSegmentValues0.stream()).collect(Collectors.toList()))
                ), getGeneListOutputFields()
        );
    }

    @Test(description = "Test simple tsv and gene list when input file has no segments.  Output files should just be headers.")
    public void testEmptyGatkCalledSegmentFile() throws IOException {
        final File outputFile = IOUtils.createTempFile("funcotatesegs_gatk_called", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addRaw("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.addRaw(TEST_GATK_EMPTY_FILE_B37);
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.addRaw(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.addRaw("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.addRaw(b37Reference);
        arguments.addRaw("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.addRaw(outputFile.getAbsolutePath());
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.addRaw("hg19");
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_CNTN4_DIR);

        runCommandLine(arguments);

        // Just read the resource config file for segments and gene file to get the column lists
        assertEmptyTsvFileHasCorrectHeaders(outputFile, SEG_RESOURCE_FILE);
        assertEmptyTsvFileHasCorrectHeaders(new File(outputFile.getAbsoluteFile() + FuncotatorEngine.GENE_LIST_FILE_SUFFIX), GENE_LIST_RESOURCE_FILE);
    }

    // "Empty" means that it has no records, but does have a header.
    private static void assertEmptyTsvFileHasCorrectHeaders(final File outputFile, final String segResourceFile) throws IOException {
        final Path configFile = Resource.getResourceContentsAsFile(segResourceFile).toPath();
        try {
            final Configuration configFileContents = new Configurations().properties(configFile.toFile());
            final List<String> expectedColumns = Lists.newArrayList(configFileContents.getKeys());
            FuncotatorTestUtils.assertTsvFieldNames(outputFile, expectedColumns);
        } catch (final ConfigurationException ce) {
            throw new UserException.BadInput("Unable to read from XSV config file: " + configFile.toUri().toString(), ce);
        }
    }

    /**
     * $ egrep -o "gene_name \"[A-Z0-9a-z]+\"" ~/IdeaProjects/gatk/src/test/resources/large/funcotator/funcotator_dataSources/gencode/hg38/gencode.v28.regressionTestVariantSet.gtf | sort | uniq
     * gene_name "CCDC40"
     * gene_name "ENAH"
     * gene_name "NPB"
     * gene_name "PCSK5"
     * gene_name "PCYT2"
     */
    @Test(description = "Simple smoke test for a hg38 GATK called segments file.")
    public void testHg38CalledCopyRatioSegmentsGatkFile() throws IOException {
        final File outputFile = IOUtils.createTempFile("funcotatehg38segs_gatk_called", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addRaw("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.addRaw(CALLED_CR_GATK_HG38_FILE);
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.addRaw(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.addRaw("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.addRaw(hg38Reference);
        arguments.addRaw("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.addRaw(outputFile.getAbsolutePath());
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.addRaw("hg38");
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER);

        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final File outputGeneList = new File(outputFile.getAbsolutePath() + FuncotatorEngine.GENE_LIST_FILE_SUFFIX);
        Assert.assertTrue(outputGeneList.exists());

        final AnnotatedIntervalCollection outputSegmentCollection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        final List<AnnotatedInterval> segments = outputSegmentCollection.getRecords();
        Assert.assertEquals(segments.size(), 5919 - 3369 );

        assertHCC1143TestTrimmedGeneListFile(outputGeneList);
    }

    @Test(description = "Simple smoke test for a hg38 GATK modeled segments file.")
    public void testModelSegmentsGatkFile() throws IOException {
        final File outputFile = IOUtils.createTempFile("funcotatehg38segs_gatk_modeled_segs", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.addRaw("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.addRaw(MODEL_SEGMENTS_FINAL_GATK_HG38_FILE);
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.addRaw(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.addRaw("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.addRaw(hg38Reference);
        arguments.addRaw("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.addRaw(outputFile.getAbsolutePath());
        arguments.addRaw("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.addRaw("hg38");
        arguments.add(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER);

        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final File outputGeneList = new File(outputFile.getAbsolutePath() + FuncotatorEngine.GENE_LIST_FILE_SUFFIX);
        assertHCC1143TestTrimmedGeneListFile(outputGeneList);
    }

    // Asserts the gene list when using the CALLED_CR_GATK_HG38_FILE and MODEL_SEGMENTS_FINAL_GATK_HG38_FILE with the
    //  funcotator_datasources hg38.  This works because the two input files were from the same run (HCC1143T).
    private void assertHCC1143TestTrimmedGeneListFile(final File outputGeneList) throws IOException {
        Assert.assertTrue(outputGeneList.exists());

        // Test gene list
        final Set<String> gtGenesToSee = Sets.newHashSet(Arrays.asList("CCDC40", "ENAH", "NPB", "PCSK5", "PCYT2"));

        try (final TableReader<LinkedHashMap<String, String>> outputReader = createLinkedHashMapListTableReader(outputGeneList)) {

            // Check that the ordering of the column is correct and that the values are all correct.
            final List<LinkedHashMap<String, String>> geneListRecords = outputReader.toList();
            // 5 genes, but one breakpoint in the middle of one of the genes, so a total of 6 entries.
            Assert.assertEquals(geneListRecords.size(), 6);

            final Set<String> genesSeen = geneListRecords.stream().map(r -> r.get("gene")).collect(Collectors.toSet());
            Assert.assertEquals(genesSeen, gtGenesToSee);
        }
    }
}
