package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Stratifies the variants by whether they overlap an interval in the set provided on the command line.
 *
 * The primary use of this stratification is to provide a mechanism to divide asssessment of a call set up
 * by whether a variant overlaps an interval or not.  I use this to differentiate between variants occurring
 * in CCDS exons vs. those in non-coding regions, in the 1000G call set, using a command line that looks like:
 *
 * -T VariantEval -R human_g1k_v37.fasta -eval 1000G.vcf -stratIntervals:BED ccds.bed -ST IntervalStratification
 *
 * Note that the overlap algorithm properly handles symbolic alleles with an INFO field END value.  In order to
 * safely use this module you should provide entire contigs worth of variants, and let the interval strat decide
 * overlap, as opposed to using -L which will not properly work with symbolic variants.
 */
public class IntervalStratification extends VariantStratifier {
    final protected static Logger logger = LogManager.getLogger(IntervalStratification.class);

    final List<Object> OVERLAPPING = Arrays.asList((Object)"all", (Object)"overlaps.intervals");
    final List<Object> NOT_OVERLAPPING = Arrays.asList((Object)"all", (Object)"outside.intervals");

    
    @Override
    public void initialize() {
        if ( getVariantEvalWalker().intervalsFile == null )
            throw new CommandLineException.MissingArgument("stratIntervals", "Must be provided when IntervalStratification is enabled");

        states.addAll(Arrays.asList("all", "overlaps.intervals", "outside.intervals"));
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval != null) {
            List<Feature> overlapping = featureContext.getValues(getVariantEvalWalker().intervalsFile);
            if ( !overlapping.isEmpty() )
                return OVERLAPPING;
            else
                return NOT_OVERLAPPING;
        }

        return Collections.emptyList();
    }
}
