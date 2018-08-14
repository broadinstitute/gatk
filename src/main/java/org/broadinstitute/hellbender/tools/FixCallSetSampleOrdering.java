package org.broadinstitute.hellbender.tools;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.collections.CollectionUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

@ExperimentalFeature
@DocumentedFeature
@CommandLineProgramProperties(summary = "Fix the sample names in a vcf that ran through the GenomicsDBImport when " +
        "the batch ordering bug was present, this will restore the correct sample names provided it is given the exact sample name mapping / vcf ordering" +
        "and batch sized that was used in the initial import. " +
        "See https://github.com/broadinstitute/gatk/issues/3682 for more information",
        oneLineSummary = "fix sample names in a shuffled callset",
        programGroup = OtherProgramGroup.class,
        omitFromCommandLine = true)
public final class FixCallSetSampleOrdering extends VariantWalker {
    public static final String SKIP_PROMPT_LONG_NAME = "skipPrompt";

    @Argument(fullName = GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME,
            doc="the same sampleNameMap file which was used to import the callset using GenomicsDBImport",
            optional = false)
    public String sampleNameMapPath;

    @Argument(fullName = GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME,
            doc="the exact batch size that was used to import the callset using GenomicsDBImport",
            minValue = 0,
            optional = false)
    public Integer batchSize;

    @Argument(fullName = GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME,
            doc="the value of the --" + GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME + " argument used when importing the callset via GenomicsDBImport",
            optional = false)
    public Integer readerThreads;

    @Argument(fullName = "gvcf-to-header-sample-map-file",
            doc="A mapping of GVCF to actual sample name in the GVCF header. You must provide this file if and only if " +
                    "GenomicsDBImport was run with --" + GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME + " > 1. " +
                    "Each line in this file should contain a GVCF path (matching the one from the sample name " +
                    "map file), followed by whitespace, followed by the actual sample name declared in the header " +
                    "for that GVCF.",
            optional = true)
    public String gvcfToHeaderSampleMapFile;

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

    private LinkedHashMap<String, Path> sampleNameMapFromGenomicsDBImport;

    private Map<Path, String> gvcfToHeaderSampleMap;

    @Override
    public void onTraversalStart() {
        assertThatTheyReallyWantToProceed();

        if (batchSize == 0) {
            throw new SampleNameFixingCannotProceedException("your callset is not affected by the bug if you ran with --"+ GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME +" 0");
        }

        if ( readerThreads > 1 ) {
            if ( gvcfToHeaderSampleMapFile == null ) {
                throw new SampleNameFixingCannotProceedException("You must provide a --gvcfToHeaderSampleMapFile if GenomicsDBImport was run with --" + GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME + " > 1");
            }
        } else if ( gvcfToHeaderSampleMapFile != null ) {
            throw new SampleNameFixingCannotProceedException("You must NOT provide a --gvcfToHeaderSampleMapFile if GenomicsDBImport was run with --" + GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME + " 1");
        }

        final VCFHeader originalHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> originalHeaderLines = originalHeader.getMetaDataInInputOrder();
        final Set<VCFHeaderLine> newHeaderLines = new LinkedHashSet<>(originalHeaderLines);
        newHeaderLines.addAll(getDefaultToolVCFHeaderLines());

        loadSampleNameMappings();
        final List<String> sampleNamesOriginalOrdering = new ArrayList<>(sampleNameMapFromGenomicsDBImport.keySet());
        if( sampleNamesOriginalOrdering.size() <= batchSize ){
            throw new SampleNameFixingCannotProceedException("you are not affected by the sample name ordering bug if your batch size is >= the number of samples in your callset. \n"
                                            + "batch size: " + batchSize + "\n"
                                            + "number of samples: " + sampleNamesOriginalOrdering.size());
        }
        assertSampleNamesMatchInputVCF(originalHeader.getSampleNamesInOrder(), sampleNamesOriginalOrdering);

        final List<String> batchSortedSampleNames = getBatchSortedList();

        final VCFHeader remappedHeader = new VCFHeader(newHeaderLines, batchSortedSampleNames);
        logger.info("Writing the new header with corrected sample names");
        writer = createVCFWriter(output);
        writer.writeHeader(remappedHeader);
        logger.info("Copying the rest of the VCF");
    }

    private void loadSampleNameMappings() {
        sampleNameMapFromGenomicsDBImport = GenomicsDBImport.loadSampleNameMapFile(IOUtils.getPath(sampleNameMapPath));
        gvcfToHeaderSampleMap = loadGvcfToHeaderSampleMap();
    }

    private Map<Path, String> loadGvcfToHeaderSampleMap() {
        if ( gvcfToHeaderSampleMapFile == null ) {
            return null;
        }

        Map<Path, String> mapping = new HashMap<>();
        final Set<Path> gvcfPathsFromSampleNameMap = new HashSet<>(sampleNameMapFromGenomicsDBImport.values());

        try {
            final List<String> lines = Files.readAllLines(IOUtils.getPath(gvcfToHeaderSampleMapFile));
            if ( lines.size() != sampleNameMapFromGenomicsDBImport.size() ) {
                throw new SampleNameFixingCannotProceedException("Number of lines in the provided --gvcfToHeaderSampleMapFile (" + lines.size() +
                        ") does not match the number of lines in the original --" + GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME +
                        " file (" + sampleNameMapFromGenomicsDBImport.size() + ")");
            }

            for ( final String line : lines ) {
                final String[] tokens = line.split("\\s+", -1);

                if ( tokens.length != 2 ) {
                    throw new SampleNameFixingCannotProceedException("Malformed line in " + gvcfToHeaderSampleMapFile + " does not have exactly two fields: " + line);
                }

                if ( tokens[0].isEmpty() || tokens[1].isEmpty() ) {
                    throw new SampleNameFixingCannotProceedException("Malformed line in " + gvcfToHeaderSampleMapFile + " contains an empty key or value: " + line);
                }

                final Path gvcfPath = IOUtils.getPath(tokens[0]);
                final String sampleName = tokens[1];

                if ( ! gvcfPathsFromSampleNameMap.contains(gvcfPath) ) {
                    throw new SampleNameFixingCannotProceedException("GVCF path " + tokens[0] + " from provided --gvcfToHeaderSampleMapFile file " +
                            "is not present in the original --" + GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME + " file");
                }
                else {
                    gvcfPathsFromSampleNameMap.remove(gvcfPath);
                }

                if ( mapping.containsKey(gvcfPath) ) {
                    throw new SampleNameFixingCannotProceedException("Duplicate GVCF path specified in the provided --gvcfToHeaderSampleMapFile: " + tokens[0]);
                }

                mapping.put(gvcfPath, sampleName);
            }

            if ( ! gvcfPathsFromSampleNameMap.isEmpty() ) {
                throw new SampleNameFixingCannotProceedException("Not all GVCF paths from the --" + GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME +
                        " were found in the provided --gvcf-to-header-sample-map-file");
            }

            return mapping;
        }
        catch ( final IOException e ) {
            throw new UserException.CouldNotReadInputFile("Error loading " + gvcfToHeaderSampleMapFile, e);
        }

    }

    private static void assertSampleNamesMatchInputVCF(final ArrayList<String> sampleNamesFromVCF, final List<String> sampleNameFromSampleMapFile) {
        if(!CollectionUtils.isEqualCollection(sampleNamesFromVCF, sampleNameFromSampleMapFile)){
            throw new SampleNameFixingCannotProceedException("The sample names in the provided sample name map file do not match the sample names in the provided vcf.\n" +
                    "It is important this tool is run with exactly the same sample name map file that was used to create the vcf.");
        }
    }

    /**
     * Recreate the sort order that the buggy GenomicsDBImport used.  It takes the sample name in input order, partitions
     * them into batches, and then sorts within each batch.
     */
    private List<String> getBatchSortedList() {
        final List<String> sampleNamesOriginalOrdering = new ArrayList<>(sampleNameMapFromGenomicsDBImport.keySet());
        final List<List<String>> batches = Lists.partition(sampleNamesOriginalOrdering, batchSize);

        if ( readerThreads == 1 ) {
            return batches.stream()
                    .flatMap(( List<String> batch ) -> batch.stream().sorted())
                    .collect(Collectors.toList());
        }
        else {
            final List<String> batchSortedSamples = new ArrayList<>();

            for ( final List<String> batch : batches ) {
                // We have to compute these lookup tables per-batch, to allow for the possibility
                // of duplicate sample names across different batches (but not within the same batch!)
                Map<String, String> sampleNameInMapToSampleNameInHeader = new HashMap<>();
                Map<String, String> sampleNameInHeaderToSampleNameInMap = new HashMap<>();

                for ( final String batchSample : batch ) {
                    final Path vcfPath = sampleNameMapFromGenomicsDBImport.get(batchSample);
                    if ( vcfPath == null ) {
                        throw new GATKException("Hash lookup failed for sample " + batchSample);
                    }

                    final String sampleNameInHeader = gvcfToHeaderSampleMap.get(vcfPath);
                    if ( sampleNameInHeader == null ) {
                        throw new GATKException("Hash lookup failed for path " + vcfPath);
                    }

                    if ( sampleNameInHeaderToSampleNameInMap.containsKey(sampleNameInHeader) ) {
                        throw new SampleNameFixingCannotProceedException("Duplicate sample name from the VCF headers detected within the same batch! " +
                                "This tool is currently unable to repair callsets in this scenario. " +
                                "Sample name was: " + sampleNameInHeader + " (as declared in the vcf header)");
                    }
                    if ( sampleNameInMapToSampleNameInHeader.containsKey(batchSample) ) {
                        throw new GATKException("Duplicate sample name from the sample name map file (" +
                                batchSample + ") detected within a batch. This should never happen!");
                    }

                    sampleNameInMapToSampleNameInHeader.put(batchSample, sampleNameInHeader);
                    sampleNameInHeaderToSampleNameInMap.put(sampleNameInHeader, batchSample);
                }

                final List<String> currentBatchSortedSamples =
                        batch.stream()
                        .map(batchSample -> sampleNameInMapToSampleNameInHeader.get(batchSample))
                        .sorted()
                        .map(headerSample -> sampleNameInHeaderToSampleNameInMap.get(headerSample))
                        .collect(Collectors.toList());
                batchSortedSamples.addAll(currentBatchSortedSamples);
            }

            return batchSortedSamples;
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        writer.add(variant);
    }

    @Override
    public void closeTool() {
        logger.info("Finished writing the new vcf with corrected sample names.");
        if( writer != null) {
            writer.close();
        }
    }
    private void assertThatTheyReallyWantToProceed(){
        if( !skipCommandLinePrompt ) {
            try (final Scanner scanner = new Scanner(System.in)) {
                System.out.println(
                        "\n\nYou are about to create a new VCF that has its header corrected for the sample swapping bug described in https://github.com/broadinstitute/gatk/issues/3682.\n" + "" +
                                "You should be certain you want to do this before proceeding.\n" +
                                "If the following description does not apply to your VCF then the newly generated vcf will be \n\n \t\tHORRIBLY CORRUPTED: by having its sample names shuffled so that the genotypes don't correspond to the correct samples\n\n" +
                                "1: your vcf was generated using a GenomicsDBImport released before gatk version 4.beta.6\n" +
                                "2: you set --batch-size != 0 when running GenomicsDBImport\n" +
                                "3: your callset was imported in multiple batches, i.e. your number of samples > batchSize\n" +
                                "4: you supplied the exact same sampleNameMap file and batch size you used in the initial GenomicsDBImport\n" +
                                "or:\n" +
                                "1. you ran GenomicsDBImport with --reader-threads > 1, and at least one sample name as declared\n" +
                                "   in a GVCF header did not match the sample name specified for that file in the sample name map file\n" +
                                "   provided to GenomicsDBImport\n\n" +

                                "If you don't know if this applies to you, please contact GATK support at https://gatkforums.broadinstitute.org/gatk/ for assistance\n" +
                                "Would you like to proceed? You can rerun with --" + SKIP_PROMPT_LONG_NAME + " to skip this check (y/N)");
                final String next = scanner.nextLine();
                if (!next.equalsIgnoreCase("y")) {
                    System.out.println("Aborting");
                    System.exit(1);
                }
            }
        }
    }

    @VisibleForTesting
    static class SampleNameFixingCannotProceedException extends UserException {
        private static final long serialVersionUID = 1L;
        SampleNameFixingCannotProceedException(final String message){
            super(message);
        }
    }

}
