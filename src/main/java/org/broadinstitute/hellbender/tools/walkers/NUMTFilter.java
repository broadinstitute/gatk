package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


@CommandLineProgramProperties(
        summary = "Creates a BAM with reads that correspond to known NUMTs",
        oneLineSummary = "Creates a BAM with reads that correspond to known NUMTs",
        programGroup = CoverageAnalysisProgramGroup.class
)
public final class NUMTFilter extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false)
    private File outputFile;

    @Argument(fullName = "groupVariantsBySource", shortName = "groupVariantsBySource", doc = "If true, group overlapping variants by their source when outputting them", optional = true)
    private boolean groupVariantsBySource = false;

    private PrintStream outputStream = null;

    static private Map<Integer, Allele> NUMTList = new HashMap<>();

    private SAMFileGATKReadWriter outputWriter;

    static {
        NUMTList.put(16301, Allele.ALT_T);
    }


    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(outputFile.getAbsolutePath()), true);
    }

    private String getNUMTAlternateBase(int position) {
        return NUMTList.get(position).getBaseString();
    }

    private byte getNUMTAlternateBaseAsByte(int position) {
        return NUMTList.get(position).getBases()[0];
    }

    @Override
    public void apply( AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext ) {

        // Passing in all FeatureInputs at once to featureContext.getValues() lets us get
        // all overlapping variants without regard to source
        for ( VariantContext variant : featureContext.getValues(variants) ) {
            // if variant is NUMT and base matches NUMT allele
            if (NUMTList.containsKey(variant.getStart())) {
                for (Allele allele : variant.getAlleles()) {
                    if (allele.isNonReference() ) {
                        String altBase = getNUMTAlternateBase(variant.getStart());
                        // TODO check if variant is homoplasmic for the NUMT allele - if so, don't add to new bam
                        // i think we will want a configurable threshold
                        if (altBase != null && allele.getBaseString().equals(altBase)) {
                             ReadPileup readsWeWant = alignmentContext.getBasePileup().makeFilteredPileup(
                                     element -> element.getBase() == getNUMTAlternateBaseAsByte(variant.getStart()));
                             readsWeWant.getReads().stream().forEach(read -> outputWriter.addRead(read));
                        }
                    }
                }
            }
        }

    }


    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}