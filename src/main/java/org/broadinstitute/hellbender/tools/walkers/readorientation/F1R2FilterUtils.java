package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

public class F1R2FilterUtils {
    static Histogram<Integer> createAltHistogram(final String refContext, final Nucleotide altAllele, final ReadOrientation type,
                                                 final int maxDepth){
        final Histogram<Integer> h = new Histogram<>(F1R2FilterConstants.binName, tripletToLabel(refContext, altAllele, type));
        h.prefillBins(F1R2FilterConstants.getEmptyBins(maxDepth));
        return h;
    }

    static Histogram<Integer> createRefHistogram(final String refContext, final int maxDepth){
        final Histogram<Integer> h = new Histogram<>(F1R2FilterConstants.binName, refContext);
        h.prefillBins(F1R2FilterConstants.getEmptyBins(maxDepth));
        return h;
    }

    // Separates an alt histogram label into components
    // e.g. "ATG_C_F1R2" becomes {Ref Context = ATG, Alt Allele = C, Read Orientaiton = F1R2}
    static Triple<String, Nucleotide, ReadOrientation> labelToTriplet(final String label){
        final String[] parts = label.split(F1R2FilterConstants.FIELD_SEPARATOR);
        Utils.validate(parts.length == 3, "Invalid label: " + label);
        return new ImmutableTriple<>(parts[0], Nucleotide.valueOf(parts[1]), ReadOrientation.valueOf(parts[2]));
    }

    // Inverse of {@link labelToTriplet}
    // e.g. {Ref Context = ATG, Alt Allele = C, Read Orientaiton = F1R2} becomes "ATG_C_F1R2"
    static String tripletToLabel(final String context, final Nucleotide altAllele, final ReadOrientation type){
        return String.join(F1R2FilterConstants.FIELD_SEPARATOR, context, altAllele.toString(), type.toString());
    }

    public static Nucleotide getMiddleBase(final String refContext){
        return Nucleotide.valueOf(refContext.substring(F1R2FilterConstants.MIDDLE_INDEX, F1R2FilterConstants.MIDDLE_INDEX+1));
    }
}
