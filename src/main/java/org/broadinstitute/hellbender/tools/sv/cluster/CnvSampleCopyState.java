package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * Memory-compact per-sample copy-number state for a pass-1 CNV clustering item, used by
 * {@link SVClusterLinkage#computeSampleOverlap} instead of retaining one htsjdk {@link Genotype} per
 * sample. At very large sample counts (e.g. ~240k) a dense region can hold ~1000+ CNV items in the
 * clustering active window simultaneously; storing each as ~240k {@code Genotype} objects exhausts the
 * heap. This holds two {@code int[]} (per sample, indexed by a shared {@link Dictionary}) plus a presence
 * {@link BitSet}, ~48x smaller, and lets the overlap computation run as primitive array scans.
 *
 * <p>The two stored values mirror exactly what {@link SVClusterLinkage}'s {@code getCopyState} reads:
 * <ul>
 *   <li>{@link #copyState}{@code [i]} = {@code CN} else {@code RD_CN} else {@code -1} &mdash; the value
 *       used when this record is the <i>present</i> side of a sample comparison.</li>
 *   <li>{@link #expectedCopyState}{@code [i]} = {@code ECN} else {@code -1} &mdash; the value supplied to
 *       the other record when this record's sample is the <i>missing</i> side.</li>
 * </ul>
 * Because the arrays are built from the same reduced genotypes the non-compact path would see, and the
 * extraction uses the same {@link VariantContextGetters#getAttributeAsInt} calls, the overlap result is
 * byte-identical to the genotype-based computation (locked by a differential unit test).</p>
 */
public final class CnvSampleCopyState {

    /** Shared, stable sample-name &harr; index mapping. One instance per run, referenced by every item. */
    public static final class Dictionary {
        private final String[] byIndex;
        private final Map<String, Integer> toIndex;

        public Dictionary(final Iterable<String> samples) {
            final java.util.List<String> names = new java.util.ArrayList<>();
            for (final String s : samples) {
                names.add(s);
            }
            this.byIndex = names.toArray(new String[0]);
            this.toIndex = new HashMap<>(SVUtils.hashMapCapacity(byIndex.length));
            for (int i = 0; i < byIndex.length; i++) {
                this.toIndex.put(byIndex[i], i);
            }
        }

        public int size() {
            return byIndex.length;
        }

        public String sampleAt(final int index) {
            return byIndex[index];
        }

        /** Global index of a sample, or -1 if not in the dictionary. */
        public int indexOf(final String sample) {
            final Integer i = toIndex.get(sample);
            return i == null ? -1 : i;
        }
    }

    private final Dictionary dictionary;
    private final int[] copyState;
    private final int[] expectedCopyState;
    private final BitSet present;

    private CnvSampleCopyState(final Dictionary dictionary, final int[] copyState,
                               final int[] expectedCopyState, final BitSet present) {
        this.dictionary = dictionary;
        this.copyState = copyState;
        this.expectedCopyState = expectedCopyState;
        this.present = present;
    }

    /**
     * Builds the compact copy-state from a record's genotypes (typically the CN-only reduced genotypes a
     * pass-1 CNV item would otherwise retain). Each genotype's sample must be in {@code dictionary}.
     */
    public static CnvSampleCopyState fromGenotypes(final GenotypesContext genotypes, final Dictionary dictionary) {
        final int n = dictionary.size();
        final int[] copyState = new int[n];
        final int[] expectedCopyState = new int[n];
        final BitSet present = new BitSet(n);
        for (final Genotype g : genotypes) {
            final int idx = dictionary.indexOf(g.getSampleName());
            Utils.validate(idx >= 0, () -> "Genotype sample not in dictionary: " + g.getSampleName());
            present.set(idx);
            copyState[idx] = copyStateOf(g);
            expectedCopyState[idx] = VariantContextGetters.getAttributeAsInt(
                    g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1);
        }
        return new CnvSampleCopyState(dictionary, copyState, expectedCopyState, present);
    }

    /**
     * Present-side copy state for a genotype: {@code CN}, else {@code RD_CN}, else {@code -1}. Must match
     * {@link SVClusterLinkage}'s {@code getCopyState} present branch exactly.
     */
    static int copyStateOf(final Genotype g) {
        return VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT,
                VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, -1));
    }

    Dictionary getDictionary() {
        return dictionary;
    }

    boolean isPresent(final int index) {
        return present.get(index);
    }

    int copyStateAt(final int index) {
        return copyState[index];
    }

    int expectedCopyStateAt(final int index) {
        return expectedCopyState[index];
    }

    /** Lowest set/clear bookkeeping for union iteration: number of samples genotyped in this record. */
    int presentCount() {
        return present.cardinality();
    }

    /** Index of the next present sample at or after {@code from}, or -1 if none. */
    int nextPresent(final int from) {
        return present.nextSetBit(from);
    }
}
