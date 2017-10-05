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
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
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

    private static File writeManyVCFs(int howMany) throws IOException {
        Map<String,File> nameToFileMap = new LinkedHashMap<>();
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

                final VCFHeader vcfHeader = new VCFHeader(headerLines, Collections.singleton(sampleName));
                vcfHeader.setSequenceDictionary(dict);
                writer.writeHeader(vcfHeader);
                final Allele Aref = Allele.create("A", true);
                final Allele C = Allele.create("C");
                final List<Allele> alleles = Arrays.asList(Aref, C);
                final VariantContext variant = new VariantContextBuilder("invented", contig, 100, 100, alleles)
                        .genotypes(new GenotypeBuilder(sampleName, alleles).attribute(SAMPLE_NAME_KEY, sampleName).make())
                        .make();
                writer.add(variant);
                nameToFileMap.put(sampleName, out);
            }
        }

        final File mapping = IOUtils.writeTempFile(
                nameToFileMap.entrySet().stream()
                        .map(pair -> pair.getKey() + "\t" + pair.getValue())
                        .sorted(Comparator.reverseOrder())
                        .collect(Collectors.joining("\n")), "mapping", ".sample_map");
     //   mapping.deleteOnExit();
        return mapping;
    }

    @Test(enabled = false)
    public void generateVCFs() throws IOException {
        final File file = writeManyVCFs(1000);
        System.out.println(file);
    }

    @DataProvider
    public Object[][] getBadlySortedFiles(){
        return new Object[][]{
                /**
                 *   test files were created by generating mangled files using the {@link #writeManyVCFs(1000)
                 *   these were then imported to GenomicsDB using GATK 4.beta.5-66-g9609cb3-SNAPSHOT which is before the fix for GenomicsDBImport
                 *   sample name ordering was applied
                 *
                 *   the combined vcfs were then extracted using SelectVariants, these are the badlySorted1000 files
                 *   each vcf contains a single record with 1000 genotypes that have a single format field "SN: sample name"
                 *   these were created with SN's that matched their assigned sample names, but the badly sorted versions are shuffled
                */

                {getTestFile("reverseOrdered.sample_map"), getTestFile("badlySorted1000-batch-size-50.vcf"), 50},
                {getTestFile("reverseOrdered.sample_map"), getTestFile("badlySorted1000-batch-size13.vcf"), 13},
        };
    }

    @Test(dataProvider = "getBadlySortedFiles")
    public void testUnshufflingRestoresCorrectSamples(final File sampleMap, final File shuffled, final int batchSize) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File fixed = createTempFile("fixed_", ".vcf");
        args.addVCF(shuffled)
                .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, sampleMap)
                .addArgument(GenomicsDBImport.BATCHSIZE_ARG_NAME, String.valueOf(batchSize))
                .addOutput(fixed);
        runCommandLine(args);

        try(final FeatureReader<VariantContext> reader =
                AbstractFeatureReader.getFeatureReader(fixed.getAbsolutePath(), new VCFCodec(), true)){
            final VariantContext vc = reader.iterator().next();
            vc.getGenotypes().iterator().forEachRemaining( g ->
                    Assert.assertEquals(g.getSampleName(), g.getAnyAttribute(SAMPLE_NAME_KEY))
            );
        }

    }
}