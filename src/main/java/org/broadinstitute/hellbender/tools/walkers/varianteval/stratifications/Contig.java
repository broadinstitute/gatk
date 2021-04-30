package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.*;

/**
 * Stratifies the evaluation by each contig in the reference sequence
 */
public class Contig extends VariantStratifier {
    public Contig(VariantEvalEngine engine) {
        super(engine);

        states.addAll(getContigNames(getEngine().getSequenceDictionaryForDrivingVariants()));
        states.add("all");
    }

    private Set<String> getContigNames(SAMSequenceDictionary dict) {
        final TreeSet<String> contigs = new TreeSet<>();
        for( final SAMSequenceRecord r :  dict.getSequences()) {
            contigs.add(r.getSequenceName());
        }
        return contigs;
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
