package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.pgen.PgenWriter;
import org.broadinstitute.pgen.PgenEmptyPgenException;

import java.io.IOException;
import java.io.OutputStream;
import java.util.EnumSet;


@SuppressWarnings("unused")
@CommandLineProgramProperties(
        summary = "(\"ExtractCohortToPgen\") - Filter and extract variants from BigQuery to a PLINK 2.0 (3 files: PGEN, PSAM, and PVAR) output.",
        oneLineSummary = "Tool to extract variants from BigQuery to a PLINK 2.0 output for a subset of samples.",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractCohortToPgen extends ExtractCohort {
    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output PGEN file to which annotated variants should be written."
    )
    private GATKPath outputPgenPath = null;

    @Argument(
            fullName = "pgen-chromosome-code",
            shortName = "pcc",
            doc = "Plink defines a set of chromosome codes that correspond to different sets of contig names for" +
                    " chromosomes. This tool supports codes 'chrM' and 'MT'.  Chromosome code reference for plink2" +
                    " can be found here: https://www.cog-genomics.org/plink/2.0/data#irreg_output"
    )
    private ChromosomeCode pgenChromosomeCode;

    @Argument(
            shortName = "wm",
            fullName = "write-mode",
            doc = "Write mode for the PGEN writer.",
            optional = true
    )
    private WriteMode writeMode = WriteMode.WRITE_AND_COPY;

    @Argument(
            fullName = "max-alt-alleles",
            shortName = "maa",
            doc = "Maximum alt alleles per site.  Cannot exceed 254.  Sites with more than the max will be excluded.",
            maxValue = PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
            optional = true
    )
    private int maxAltAlleles = PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES;

    @Argument(
            fullName = "lenient-ploidy-validation",
            shortName = "lpv",
            doc = "PGEN requires individual sample to be diploid (except for sex chromosomes, which may be haploid - " +
                    "these are accepted and recoded for pgen as heterozygous/diploid). By default, any ploidy failure " +
                    "will result in an exception to be thrown. Use this flag to tolerate ploidy failures " +
                    "(samples will be recoded as missing, and logged if a log file is provided using --writerLogFile)",
            optional = true
    )
    private boolean lenientPloidyValidation = false;

    @Argument(
            fullName = "writer-log-file",
            shortName = "wlf",
            doc = "An optional second log file that the PGEN Writer will use to log sites dropped for exceeding " +
                    "maxAltAlleles, and also, if --lenientPloidyValidation is used, logs sites and samples where " +
                    "samples have to be marked as missing for non-conforming ploidy.",
            optional = true
    )
    private String writerLogFile = null;

    @Argument(
            fullName = "allow-empty-pgen",
            doc = "In the case that no variants are written to the PGEN file, this flag allows the tool to write an " +
                    "empty PGEN file.  By default, the tool will throw an exception if no variants are written, because" +
                    " an empty PGEN file is not technically a valid PGEN file.",
            optional = true
    )
    private boolean allowEmptyPgen = false;

    private PgenWriter pgenWriter = null;

    /**
     * Enum for user input options for pgen chromosome code.  Converts to an equivalent PgenWriter.PgenChromosomeCode
     */
    private enum ChromosomeCode implements CommandLineParser.ClpEnum {
        chrM("Corresponds to hg38 contig names (e.g. 'chrX', 'chrY', etc.)", PgenWriter.PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM),
        MT("Corresponds to b37 contig names (e.g. 'X', 'Y', etc.)", PgenWriter.PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT);

        final String description;
        final PgenWriter.PgenChromosomeCode equivalentCode;

        ChromosomeCode(String description, PgenWriter.PgenChromosomeCode equivalentCode) {
            this.description = description;
            this.equivalentCode = equivalentCode;
        }

        @Override
        public String getHelpDoc() {
            return description;
        }

        public PgenWriter.PgenChromosomeCode getEquivalentCode() {
            return this.equivalentCode;
        }
    }

    /**
     * Enum for user input options for pgen write mode.  Converts to an equivalent PgenWriter.PgenWriteMode
     * Excluding backward seek because we won't know the variant count beforehand
     */
    private enum WriteMode implements CommandLineParser.ClpEnum {
        WRITE_SEPARATE_INDEX(
                "Write a separate .pgi index file",
                PgenWriter.PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX
        ),
        WRITE_AND_COPY(
                "The final real .pgen is only created at the end, by writing the index and then appending " +
                        "the body of the first temporary .pgen (which is then deleted)",
                PgenWriter.PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY
        );

        final String description;
        final PgenWriter.PgenWriteMode equivalentMode;

        WriteMode(String description, PgenWriter.PgenWriteMode equivalentMode) {
            this.description = description;
            this.equivalentMode = equivalentMode;
        }

        @Override
        public String getHelpDoc() {
            return description;
        }

        public PgenWriter.PgenWriteMode getEquivalentMode() {
            return this.equivalentMode;
        }
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        pgenWriter = new PgenWriter(
                outputPgenPath,
                header,
                writeMode.getEquivalentMode(),
                EnumSet.noneOf(PgenWriter.PgenWriteFlag.class),
                pgenChromosomeCode.getEquivalentCode(),
                lenientPloidyValidation,
                PgenWriter.VARIANT_COUNT_UNKNOWN,
                maxAltAlleles,
                writerLogFile
        );

    }

    @Override
    protected void apply(VariantContext variantContext) {
        if (variantContext != null) {
            // Add the variant contexts that aren't filtered or add everything if we aren't excluding anything
            if (variantContext.isNotFiltered() || !excludeFilteredSites) {
                try {
                    pgenWriter.add(variantContext);
                }
                catch(IllegalStateException e) {
                    logger.debug("Encountered an error.  Here's some debug info:\n" +
                            "ID: " + variantContext.getID() + "\n" +
                            "NAlleles: " + variantContext.getNAlleles() + "\n" +
                            "Contig: " + variantContext.getContig() + "\n" +
                            "Start: " + variantContext.getStart() + "\n" +
                            "End: " + variantContext.getEnd() + "\n" +
                            "NoCallCount: " + variantContext.getNoCallCount() + "\n" +
                            "HomRefCount: " + variantContext.getHomRefCount() + "\n" +
                            "HetCount: " + variantContext.getHetCount() + "\n" +
                            "HomVarCount: " + variantContext.getHomVarCount() + "\n" +
                            "MixedCount: " + variantContext.getMixedCount() + "\n" +
                            "NSamples: " + variantContext.getNSamples() + "\n" +
                            "Genotypes.size(): " + variantContext.getGenotypes().size());
                    throw e;
                }
            }
            progressMeter.update(variantContext);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if(pgenWriter.getDroppedVariantCount() > 0l) {
            logger.info(pgenWriter.getDroppedVariantCount() + " variants dropped by writer for exceeding the" +
                    " maximum number of allowed alt alleles.");
        }

        return null;
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        // Close up our writer if we have to:
        if (pgenWriter != null) {
            try {
                pgenWriter.close();
            }
            catch(PgenEmptyPgenException e) {
                if(allowEmptyPgen) {
                    logger.warn("No variants were written to the PGEN file.  This is not technically a valid PGEN file, " +
                            "but the --allow-empty-pgen flag was used, so an empty PGEN file was written.");
                    writeEmptyPgenFiles();
                }
                else {
                    throw e;
                }
            }
            catch(Exception e) {
                throw e;
            }
        }
    }

    /**
     * Writes empty .pgen, .pvar.zst, and .psam files to the output path
     */
    private void writeEmptyPgenFiles() {
        // Make output stream for the pgen file
        final OutputStream emptyPgenOutputStream = outputPgenPath.getOutputStream();
        // Build paths for the pvar and psam files and make output streams for them
        final String pgenAbsolutePath = outputPgenPath.toPath().toAbsolutePath().toString();
        final GATKPath pvarPath = new GATKPath(pgenAbsolutePath.substring(0, pgenAbsolutePath.lastIndexOf(".pgen")) + ".pvar.zst");
        final OutputStream emptyPvarOutputStream = pvarPath.getOutputStream();

        final GATKPath psamPath = new GATKPath(pgenAbsolutePath.substring(0, pgenAbsolutePath.lastIndexOf(".pgen")) + ".psam");
        final OutputStream emptyPsamOutputStream = psamPath.getOutputStream();

        // Write empty files
        try {
            emptyPgenOutputStream.close();
            emptyPvarOutputStream.close();
            emptyPsamOutputStream.close();
        }
        catch (IOException e) {
            throw new UserException("Failed to create empty pgen/pvar.zst/psam files from pgen path " + outputPgenPath.getRawInputString() + "\n Caused by:" + e.getMessage(), e);
        }
        
    }
}
