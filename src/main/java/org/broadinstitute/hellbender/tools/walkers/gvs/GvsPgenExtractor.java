package org.broadinstitute.hellbender.tools.walkers.gvs;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.pgen.PgenWriter;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.LinkedHashSet;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Reads from a VCF and writes the variants to a PGEN file",
        oneLineSummary = "Reads from a VCF and writes the variants to a PGEN file",
        programGroup = VariantFilteringProgramGroup.class
)
public final class GvsPgenExtractor extends TwoPassVariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output pgen file", optional = true)
    private GATKPath outputFile = null;

    @Argument(fullName = "maxAltAlleles", shortName = "maa", doc = "Maximum alt alleles to write", optional = true)
    private int maxAltAlleles = 254;

    @Argument(fullName = "writeMode", shortName = "wm", doc = "Write mode for the PGEN writer", optional = true)
    private PgenWriter.PgenWriteMode writeMode = PgenWriter.PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY;

    /**
     * This argument can be specified multiple times in order to provide multiple sample names, or to specify
     * the name of one or more files containing sample names. File names must use the extension ".args", and the
     * expected file format is simply plain text with one sample name per line. Note that sample exclusion takes
     * precedence over inclusion, so that if a sample is in both lists it will be excluded.
     */
    @Argument(fullName=StandardArgumentDefinitions.SAMPLE_NAME_LONG_NAME, shortName=StandardArgumentDefinitions.SAMPLE_NAME_SHORT_NAME, doc="Include genotypes from this sample", optional=true)
    private Set<String> sampleNames = new LinkedHashSet<>(0);

    private PgenWriter pgenWriter;

    // Count of variants in the input file.  We need to count this up to supply to the PgenWriter
    private int numberOfVariants = 0;

    @Override
    protected void firstPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        numberOfVariants += 1;
    }

    @Override
    protected void afterFirstPass() {
        final VCFHeader inputHeader = getHeaderForVariants();
        pgenWriter = new PgenWriter(outputFile, inputHeader, writeMode, numberOfVariants, maxAltAlleles);
    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        pgenWriter.add(variant);
    }

    @Override
    public Object onTraversalSuccess() {
        // Log the number of dropped variants
        logger.info("Number of variants dropped as a result of exceeding the maximum number of alt alleles: " + pgenWriter.getDroppedVariantCount());
        return null;
    }

    @Override
    public void closeTool() {
        if ( pgenWriter != null ) {
            pgenWriter.close();
        }
    }
}
