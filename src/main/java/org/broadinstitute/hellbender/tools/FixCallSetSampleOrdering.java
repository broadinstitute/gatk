package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Lists;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.collections.CollectionUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.io.IOUtils;


import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

@BetaFeature
@CommandLineProgramProperties(summary = "fix the sample names in a vcf that ran through the GenomicsDBImport when " +
        "the batch ordering bug was present, this will restore the correct sample names provided it is given the exact sample name mapping / vcf ordering" +
        "and batch sized that was used in the initial import" +
        "see https://github.com/broadinstitute/gatk/issues/3682 for more information",
        oneLineSummary = "fix sample names in a shuffled callset",
        programGroup = VariantProgramGroup.class,
        omitFromCommandLine = true)
public final class FixCallSetSampleOrdering extends VariantWalker {
    public static final String SKIP_PROMPT_LONG_NAME = "skipPrompt";

    @Argument(fullName = GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME,
            doc="the same sampleNameMap file which was used to import the callset using GenomicsDBImport",
            optional = false)
    public String sampleNameMapPath;

    @Argument(fullName = GenomicsDBImport.BATCHSIZE_ARG_NAME,
            doc="the exact batch size that was used to import the callset using GenomicsDBImport",
            minValue = 0,
            optional = false)
    public Integer batchSize;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="where to write a reheadered version of the input VCF with the sample names in the correct order",
            optional = false)
    public File output;

    @Argument(fullName = SKIP_PROMPT_LONG_NAME,
            shortName = "Y",
            doc = "skip the commandline prompt as if you had entered Y",
            optional = true)
    public boolean skipCommandLinePrompt = false;

    private VariantContextWriter writer;

    @Override
    public void onTraversalStart() {
        assertThatTheyReallyWantToProceed();

        if (batchSize == 0) {
            throw new UserException("your callset is not affected by the bug if you ran with "+ GenomicsDBImport.BATCHSIZE_ARG_NAME +" 0");
        }

        final VCFHeader originalHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> originalHeaderLines = originalHeader.getMetaDataInInputOrder();
        final Set<VCFHeaderLine> newHeaderLines = new LinkedHashSet<>(originalHeaderLines);
        newHeaderLines.addAll(getDefaultToolVCFHeaderLines());

        final List<String> sampleNamesOriginalOrdering = getSampleNamesInInputOrder(sampleNameMapPath);
        if( sampleNamesOriginalOrdering.size() <= batchSize ){
            throw new UserException("you are not affected by the sample name ordering bug if your batch size is >= the number of samples in your callset. \n"
                                            + "batch size: " + batchSize + "\n"
                                            + "number of samples: " + sampleNamesOriginalOrdering.size());
        }

        assertSampleNamesMatchInputVCF(originalHeader.getSampleNamesInOrder(), sampleNamesOriginalOrdering);
        final List<String> batchSortedSampleNames = getBatchSortedList(sampleNamesOriginalOrdering, batchSize);

        final VCFHeader remappedHeader = new VCFHeader(newHeaderLines, batchSortedSampleNames);
        logger.info("Writing the new header with corrected sample names");
        writer = createVCFWriter(output);
        writer.writeHeader(remappedHeader);
        logger.info("Copying the rest of the VCF");
    }

    private static void assertSampleNamesMatchInputVCF(final ArrayList<String> sampleNamesFromVCF, final List<String> sampleNameFromSampleMapFile) {
        if(!CollectionUtils.isEqualCollection(sampleNamesFromVCF, sampleNameFromSampleMapFile)){
            throw new UserException("The sample names in the provided sample name map file do not match the sample names in the provided vcf.\n" +
                    "It is important this tool is run with exactly the same sample name map file that was used to create the vcf.");
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        writer.add(variant);
    }

    @Override
    public void closeTool() {
        logger.info("Finished writing the new vcf with corrected sampled names.");
        if( writer != null) {
            writer.close();
        }
    }
    private void assertThatTheyReallyWantToProceed(){
        if( !skipCommandLinePrompt ){
            final Scanner scanner = new Scanner(System.in);
            System.out.println(
                    "\n\nYou are about to create a new VCF that has it's header corrected for the sample swapping bug described in https://github.com/broadinstitute/gatk/issues/3682.\n"+ "" +
                            "You should be certain you want to do this before proceeding.\n" +
                            "If the following description does not apply to your VCF then the newly generated vcf will be \n\n \t\tHORRIBLY CORRUPTED\n\n" +
                            "1: your vcf was generated using a GenomicsDBImport released before gatk version 4.beta.6\n" +
                            "2: you set --batchSize != 0 when running GenomicsDBImport\n" +
                            "3: your callset was imported in multiple batches, i.e. your number of samples > batchSize\n" +
                            "4: you supplied the exact same sampleNameMap file and batch size you used in the initial GenomicsDBImport\n\n" +

                            "If you don't know if this applies to you, please contact GATK support at https://gatkforums.broadinstitute.org/gatk/ for assistance\n" +
                            "Would you like to proceed? You can rerun with --"+SKIP_PROMPT_LONG_NAME+" to skip this check (y/N)");
            final String next = scanner.next();
            if (!next.equalsIgnoreCase("y")){
                System.out.println("Aborting");
                System.exit(0);
            }
        }
    }

    private static List<String> getSampleNamesInInputOrder(final String sampleNameMapPath) {
        final Map<String, Path> stringPathMap = GenomicsDBImport.loadSampleNameMapFile(IOUtils.getPath(sampleNameMapPath));
        return new ArrayList<>(stringPathMap.keySet());
    }

    /**
     * Recreate the sort order that the buggy GenomicsDBImport used.  It takes the sample name in input order, partitions
     * them into batches, and then sorts within each batch.
     */
    private static List<String> getBatchSortedList(final List<String> sampleNames, final int batchSize){
        final List<List<String>> partition = Lists.partition(sampleNames, batchSize);
        return partition.stream()
                .flatMap( (List<String> batch) -> batch.stream().sorted() )
                .collect(Collectors.toList());
    }

}
