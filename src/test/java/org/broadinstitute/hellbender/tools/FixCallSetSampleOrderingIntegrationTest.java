package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class FixCallSetSampleOrderingIntegrationTest extends CommandLineProgramTest{

    private static final String SAMPLE_NAME_KEY = "SN";
    private static final String REVERSE_ORDERED_SAMPLE_MAP = "reverseOrdered.sample_map";
    private static final String BADLY_SORTED1000_BATCH_SIZE50_VCF = "badlySorted1000-batch-size50.vcf";

    private static File writeManyVCFs(int howMany, boolean samplesInHeaderCanMismatch, boolean addDuplicateHeaderSamples) throws IOException {
        final Map<String,File> nameToFileMap = new LinkedHashMap<>();
        final String contig = "chr20";
        final SAMSequenceDictionary dict = new SAMSequenceDictionary(
                Collections.singletonList(new SAMSequenceRecord(contig, 64444167)));

        final VCFFormatHeaderLine formatField = new VCFFormatHeaderLine(SAMPLE_NAME_KEY, 1, VCFHeaderLineType.String,
                                                                        "the name of the sample this genotype came from");
        final Set<VCFHeaderLine> headerLines = new HashSet<>();
        headerLines.add(formatField);
        headerLines.add(VCFStandardHeaderLines.getFormatLine("GT"));

        for( int i = 0; i < howMany ; i++) {
            final File out = File.createTempFile("tiny_", ".vcf");

            try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(out, dict, false, Options.INDEX_ON_THE_FLY)) {
                final String sampleName = "Sample_" + String.valueOf(i);
                String sampleNameInHeader;
                if ( samplesInHeaderCanMismatch && addDuplicateHeaderSamples && (i == 100 || i == 700) ) {
                    sampleNameInHeader = "FromHeader_DUPLICATE_NAME";
                } else if ( samplesInHeaderCanMismatch && i % 10 == 0 ) {
                    sampleNameInHeader = "FromHeader_" + sampleName;
                } else {
                    sampleNameInHeader = sampleName;
                }

                final VCFHeader vcfHeader = new VCFHeader(headerLines, Collections.singleton(sampleNameInHeader));
                vcfHeader.setSequenceDictionary(dict);
                writer.writeHeader(vcfHeader);
                final Allele Aref = Allele.create("A", true);
                final Allele C = Allele.create("C");
                final List<Allele> alleles = Arrays.asList(Aref, C);
                final VariantContext variant = new VariantContextBuilder("invented", contig, 100, 100, alleles)
                        .genotypes(new GenotypeBuilder(sampleNameInHeader, alleles).attribute(SAMPLE_NAME_KEY, sampleName).make())
                        .make();
                writer.add(variant);
                nameToFileMap.put(sampleName, out);
            }
        }

        return IOUtils.writeTempFile(
                nameToFileMap.entrySet().stream()
                        .map(pair -> pair.getKey() + "\t" + pair.getValue())
                        .sorted(Comparator.reverseOrder())
                        .collect(Collectors.joining("\n")), "mapping", ".sample_map");
    }

    /**
     * This was used to generate some of the vcf files used in {@link #testUnshufflingRestoresCorrectSamples(File, File, int, int , File)}
     * it's kept around in case we need to generate new test cases or examine how they were made.
     *
     * This method generates VCFs in which the sample names in the VCF headers all match the sample names in the
     * generated sample name map file.
     */
    @Test(enabled = false)
    public void generateVCFs() throws IOException {
        final File file = writeManyVCFs(1000, false, false);
        System.out.println(file);
    }

    /**
     * this was used to generate some of the vcf files used in {@link #testUnshufflingRestoresCorrectSamples(File, File, int, int , File)}
     * it's kept around in case we need to generate new test cases or examine how they were made
     *
     * This was run two ways: once with addDuplicateHeaderSamples == false, and once with addDuplicateHeaderSamples == true.
     *
     * This method generates VCFs in which some of the sample names in the VCF headers don't match the sample names in the
     * generated sample name map file.
     */
    @Test(enabled = false)
    public void generateVCFsWithMismatchingSampleNamesInHeaders() throws IOException {
        final File file = writeManyVCFs(1000, true, true);
        System.out.println(file);
    }

    /**
     *   Test files were created by generating mangled files using the {@link #writeManyVCFs(int, boolean, boolean)} methods.
     *   These were then imported to GenomicsDB using GATK 4.beta.5-66-g9609cb3-SNAPSHOT which is before the fix for
     *   GenomicsDBImport sample name ordering was applied.
     *
     *   For the VCF set in which the sample names in the VCF headers do not all match the sample names in the
     *   sample name map provided to GenomicsDBImport, GenomicsDBImport was run with --readerThreads 4 and
     *   --batchSize 50, to trigger the buggy multithreaded vcf loading code path.
     *
     *   The combined vcfs were then extracted using SelectVariants, these are the badlySorted1000 files.
     *   Each vcf contains a single record with 1000 genotypes that have a single format field "SN: sample name"
     *   these were created with SN's that matched their assigned sample names, but the badly sorted versions are shuffled.
     */
    @DataProvider
    public Object[][] getBadlySortedFiles(){
        return new Object[][]{
                {getTestFile(REVERSE_ORDERED_SAMPLE_MAP), getTestFile(BADLY_SORTED1000_BATCH_SIZE50_VCF), 50, 1, null},
                {getTestFile(REVERSE_ORDERED_SAMPLE_MAP), getTestFile("badlySorted1000-batch-size50.vcf.gz"), 50, 1, null},
                {getTestFile(REVERSE_ORDERED_SAMPLE_MAP), getTestFile("badlySorted1000-batch-size13.vcf"), 13, 1, null},
                {getTestFile("reverseOrderedWithHeaderMismatches.sample_map"), getTestFile("badlysorted1000-batchSize50-readerThreads5-samplesInHeaderAreDifferent.vcf"), 50, 5, getTestFile("reverseOrderedWithHeaderMismatches.gvcfToHeaderSampleMapFile")},
                {getTestFile("reverseOrderedWithHeaderMismatchesAndDuplicateSample.sample_map"), getTestFile("badlysorted1000-batchSize50-readerThreads5-samplesInHeaderAreDifferent-withDuplicateSample.vcf"), 50, 5, getTestFile("reverseOrderedWithHeaderMismatchesAndDuplicateSample.gvcfToHeaderSampleMapFile")}
        };
    }

    @Test(dataProvider = "getBadlySortedFiles")
    public void testUnshufflingRestoresCorrectSamples(final File sampleMap, final File shuffled, final int batchSize, final int readerThreads, final File gvcfToHeaderSampleMapFile) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File fixed = createTempFile("fixed_", ".vcf");
        args.addVCF(shuffled)
                .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, sampleMap)
                .addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, String.valueOf(readerThreads))
                .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                .addBooleanArgument(FixCallSetSampleOrdering.SKIP_PROMPT_LONG_NAME, true)
                .addOutput(fixed);
        if ( gvcfToHeaderSampleMapFile != null ) {
            args.addFileArgument("gvcf-to-header-sample-map-file", gvcfToHeaderSampleMapFile);
        }
        
        runCommandLine(args);

        try(final FeatureReader<VariantContext> reader =
                AbstractFeatureReader.getFeatureReader(fixed.getAbsolutePath(), new VCFCodec(), true)){
            final VariantContext vc = reader.iterator().next();
            Assert.assertEquals(vc.getGenotypes().size(), 1000);
            vc.getGenotypes().iterator().forEachRemaining( g ->
                    Assert.assertEquals(g.getSampleName(), g.getAnyAttribute(SAMPLE_NAME_KEY))
            );
        }
    }

    @DataProvider
    public Object[][] getBadInputs() {
        return new Object[][]{
                // smaller number of samples than batch size
                {getTestFile(REVERSE_ORDERED_SAMPLE_MAP), getTestFile(BADLY_SORTED1000_BATCH_SIZE50_VCF), 10000, 1, null},
                // batch size 0
                {getTestFile(REVERSE_ORDERED_SAMPLE_MAP), getTestFile(BADLY_SORTED1000_BATCH_SIZE50_VCF), 0, 1, null},
                // samples in sample map don't match badly sorted vcf
                {getTestFile("mismatched.sample_map"), getTestFile(BADLY_SORTED1000_BATCH_SIZE50_VCF), 1, 1, null},
                // readerThreads > 1, but no gvcfToHeaderSampleMapFile provided
                {getTestFile("reverseOrderedWithHeaderMismatches.sample_map"), getTestFile("badlysorted1000-batchSize50-readerThreads5-samplesInHeaderAreDifferent.vcf"), 50, 5, null},
                // readerThreads == 1, but a gvcfToHeaderSampleMapFile provided
                {getTestFile("reverseOrderedWithHeaderMismatches.sample_map"), getTestFile("badlysorted1000-batchSize50-readerThreads5-samplesInHeaderAreDifferent.vcf"), 50, 1, getTestFile("reverseOrderedWithHeaderMismatches.gvcfToHeaderSampleMapFile")}
        };
    }

    @Test(expectedExceptions = FixCallSetSampleOrdering.SampleNameFixingCannotProceedException.class, dataProvider = "getBadInputs")
    public void testErrorConditions(final File sampleMap, final File shuffled, final int batchSize, final int readerThreads, final File gvcfToHeaderSampleMapFile){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File fixed = createTempFile("fixed_", ".vcf");
        args.addVCF(shuffled)
                .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, sampleMap)
                .addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, String.valueOf(readerThreads))
                .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                .addBooleanArgument(FixCallSetSampleOrdering.SKIP_PROMPT_LONG_NAME, true)
                .addOutput(fixed);
        if ( gvcfToHeaderSampleMapFile != null ) {
            args.addFileArgument("gvcf-to-header-sample-map-file", gvcfToHeaderSampleMapFile);
        }

        runCommandLine(args);
    }
}