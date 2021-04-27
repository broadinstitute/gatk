package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantSummary;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Stratifies the eval RODs by the allele count of the alternate allele
 *
 * Looks first at the MLEAC value in the INFO field, and uses that value if present.
 * If not present, it then looks for the AC value in the INFO field.  If both are absent,
 * it computes the AC from the genotypes themselves.  If no AC can be computed, 0 is used.
 */
public class AlleleCount extends VariantStratifier {
    private static final Logger log = LogManager.getLogger(AlleleCount.class);
    private int nchrom;

    public AlleleCount(VariantEvalEngine engine) {
        super(engine);

        // There are ploidy x n sample chromosomes
        // TODO -- generalize to handle multiple ploidy
        nchrom = getEngine().getNumberOfSamplesForEvaluation() * getEngine().getVariantEvalArgs().getPloidy();

        // create an array containing each of the allele counts
        for( int ac = 0; ac <= nchrom; ac++ ) {
            states.add(ac);
        }

        log.info("AlleleCount using " + nchrom + " chromosomes");
    }

    @Override
    public void validateArgs() {
        // we can only work with a single eval VCF, and it must have genotypes
        if ( getEngine().getVariantEvalArgs().getEvals().size() != 1 && !getEngine().getVariantEvalArgs().isMergeEvals() )
            throw new CommandLineException.BadArgumentValue("AlleleCount", "AlleleCount stratification only works with a single eval vcf");

        if ( nchrom < 2 )
            throw new CommandLineException.BadArgumentValue("AlleleCount", "AlleleCount stratification requires an eval vcf with at least one sample");
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null) {
            int AC = 0; // by default, the site is considered monomorphic

            try {
                if ( eval.isBiallelic() ) {
                    if ( eval.hasAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) ) {
                        // the MLEAC is allowed to be larger than the AN (e.g. in the case of all PLs being 0, the GT is ./. but the exact model may arbitrarily choose an AC>1)
                        AC = Math.min(eval.getAttributeAsInt(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, 0), nchrom);
                    } else if ( eval.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
                        AC = eval.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0);
                    }
                }
            } catch ( ClassCastException e ) {
                // protect ourselves from bad inputs
                // TODO -- fully decode VC
            }

            if ( AC == 0 && eval.isVariant() ) {
                // fall back to the direct calculation
                for (Allele allele : eval.getAlternateAlleles())
                    AC = Math.max(AC, eval.getCalledChrCount(allele));
            }

            // make sure that the AC isn't invalid
            if ( AC > nchrom )
                throw new UserException(String.format("The AC value (%d) at position %s:%d " +
                        "is larger than the number of chromosomes over all samples (%d)", AC,
                        eval.getContig(), eval.getStart(), nchrom));

            return Collections.singletonList((Object) AC);
        } else {
            return Collections.emptyList();
        }
    }

    @Override
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return new HashSet<Class<? extends VariantEvaluator>>(Arrays.asList(VariantSummary.class));
    }

    @Override
    public String getFormat() {
        return "%d";
    }
}
