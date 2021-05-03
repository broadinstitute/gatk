package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;

/**
 * Stratifies the evaluation by each contig in the reference sequence. Note: if the user supplies custom intervals, it will defer to these rather than the full sequence dictionary
 */
public class Contig extends VariantStratifier {
    public Contig(VariantEvalEngine engine) {
        super(engine);

        states.addAll(getContigNames());
        states.add("all");
    }

    /**
     * @return The list of contig names to be traversed, preferentially taking user supplied intervals, but otherwise defaulting to driving variants
     */
    private List<String> getContigNames() {
        final TreeSet<String> contigs = new TreeSet<>();
        if (getEngine().getTraversalIntervals() == null) {
            getEngine().getSequenceDictionaryForDrivingVariants().getSequences().stream().map(SAMSequenceRecord::getSequenceName).forEach(contigs::add);
        }
        else {
            getEngine().getTraversalIntervals().stream().map(SimpleInterval::getContig).forEach(contigs::add);
        }

        return new ArrayList<>(contigs);
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null) {
            return Arrays.asList("all", eval.getContig());
        } else {
            return Collections.emptyList();
        }
    }
}
