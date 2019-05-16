package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.ChromosomeCounts;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Example/toy program that shows how to implement the VariantWalker interface. Prints supplied variants
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints variants supplied to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Example tool that prints variants with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleVariantWalker extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    private PrintStream outputStream = null;

    private VariantContextWriter vcfWriter = null;

    @Override
    public void onTraversalStart() {
//        try {
//            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
//
//
//        }
//        catch ( final FileNotFoundException e ) {
//            throw new UserException.CouldNotReadInputFile(outputFile, e);
//        }


        final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());

        // Initialize VCF header lines
        final Set<VCFHeaderLine> headerLines = createVCFHeaderLineList(vcfHeaders);
        VCFInfoHeaderLine v = new VCFInfoHeaderLine("source", 1, VCFHeaderLineType.String, "what caller the variant came from");
        headerLines.add(v);

        SortedSet<String> samples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);


        vcfWriter = createVCFWriter(outputFile);
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        Set<String> sample_names = variant.getSampleNames();
        Boolean write = true;
        int total_het_samples = 0;
        int ad_failures = 0;

        for(String sample : sample_names) {
            if (write == false){
                break;
            }

            Genotype genotype = variant.getGenotype(sample);
            if (genotype.getGQ() < 25) {
                write = false;
                break;
            }
//            if (genotype.isHomVar()){
//                logger.warn("FOUND HOM VAR VARIANT???????????");
//            }
            if (genotype.isHet()) {
                total_het_samples += 1;
                BinomialDistribution b = new BinomialDistribution(genotype.getDP(), .5);
                int cutoff = b.inverseCumulativeProbability(.01);
                for (int ad : genotype.getAD()) {
                    if (ad < cutoff) {
                        ad_failures += 1;
                    }
                }
            }
        }

        if ((ad_failures * 1.0) / total_het_samples > .5) {
            write = false;
        }

        if (variant.getNoCallCount() == 0 && write) {
            vcfWriter.add(variant);
        }

    }

    /**
     * Prepare the VCF header lines
     */
    private Set<VCFHeaderLine> createVCFHeaderLineList(Map<String, VCFHeader> vcfHeaders) {

        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfHeaders.values(), true);
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AC_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AF_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AN_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_DP_KEY));

        headerLines.addAll(Arrays.asList(ChromosomeCounts.descriptions));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));

        return headerLines;
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
        if (vcfWriter != null){
            vcfWriter.close();
        }
    }
}
