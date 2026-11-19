package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * The reference contexts of loci, read from a window of reference bases that is fetched again only when a locus's
 * context falls outside it.  Loci come in order, so this replaces the intervals, sequence dictionary lookups, base copies
 * and strings of {@link ReferenceContext#getKmerAround} at every locus with an array lookup.
 */
final class ReferenceContextWindow {
    static final int DEFAULT_WINDOW_SIZE = 10_000;

    // each context of F1R2FilterConstants.ALL_KMERS, indexed by its bases as a number in base 4, with A, C, G and T as
    // 0, 1, 2 and 3, which BaseUtils.simpleBaseToBaseIndex gives them
    private static final String[] CONTEXTS_BY_CODE = new String[F1R2FilterConstants.NUM_KMERS];
    static {
        for (final String context : F1R2FilterConstants.ALL_KMERS) {
            CONTEXTS_BY_CODE[code(context.getBytes(), 0)] = context;
        }
    }

    private final int windowSize;

    // the reference bases of contig:[start, start + bases.length)
    private String contig = null;
    private int start = 0;
    private byte[] bases = new byte[0];

    ReferenceContextWindow(final int windowSize) {
        Utils.validateArg(windowSize >= F1R2FilterConstants.REFERENCE_CONTEXT_SIZE, "the window must hold a reference context");
        this.windowSize = windowSize;
    }

    /**
     * The reference context centered on a position of a contig, as {@link ReferenceContext#getKmerAround} would give it,
     * or null if the context would run off either end of the contig or contains an N
     *
     * @param referenceContext a reference context on the contig, whose data source supplies the window's bases
     * @return one of {@link F1R2FilterConstants#ALL_KMERS}, or null
     */
    String contextAt(final ReferenceContext referenceContext, final String contig, final int position) {
        final int contextStart = position - F1R2FilterConstants.REF_CONTEXT_PADDING;
        final int contextEnd = position + F1R2FilterConstants.REF_CONTEXT_PADDING;
        if (contextStart < 1) {
            return null;
        }

        if (!contig.equals(this.contig) || contextStart < start || contextEnd >= start + bases.length) {
            // the window starts at the context, and getBases trims it at the end of the contig
            bases = referenceContext.getBases(new SimpleInterval(contig, contextStart, contextStart + windowSize - 1));
            this.contig = contig;
            start = contextStart;
            if (contextEnd >= start + bases.length) {
                return null;
            }
        }

        final int code = code(bases, contextStart - start);
        return code < 0 ? null : CONTEXTS_BY_CODE[code];
    }

    /**
     * The reference context of bases starting at an offset as a number in base 4, or -1 if a base isn't A, C, G or T.
     * GATK's reference has only upper-case A, C, G, T and N.
     */
    private static int code(final byte[] bases, final int offset) {
        int code = 0;
        for (int n = offset; n < offset + F1R2FilterConstants.REFERENCE_CONTEXT_SIZE; n++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[n]);
            if (baseIndex == -1) {
                return -1;
            }
            code = 4 * code + baseIndex;
        }
        return code;
    }
}
