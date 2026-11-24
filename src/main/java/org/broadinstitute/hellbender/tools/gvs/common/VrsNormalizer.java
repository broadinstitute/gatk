package org.broadinstitute.hellbender.tools.gvs.common;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * VRS allele normalization implementation.
 *
 * This implementation handles substitutions and full indel normalization
 * (rolling, ambiguity-bound expansion, and state classification).
 *
 * <p>Reference sequence is accessed via a {@link ReferenceContext} that is
 * passed per {@link #normalize} call, matching the GATK walker pattern where
 * reference data is provided through {@code apply()} rather than stored at
 * construction time.</p>
 */
public class VrsNormalizer {

    // No stored reference — reference is passed per normalize() call.
    public VrsNormalizer() {}

    /**
     * Normalize an allele according to VRS specification.
     *
     * @param contig     Chromosome/contig name (e.g., "19", "chr19")
     * @param start      0-based interbase start position
     * @param end        0-based interbase end position (half-open interval)
     * @param alt        Alternate allele sequence
     * @param refContext ReferenceContext for the current variant (from apply())
     * @return Normalized allele with adjusted coordinates and state
     */
    public NormalizedAllele normalize(String contig, long start, long end, String alt, ReferenceContext refContext) {
        // Step 0: Get reference sequence
        String ref = getRefSequence(contig, start, end, refContext);

        // Step 1: Trim common prefix and suffix
        TrimResult trimmed = trimCommonFlanks(ref, alt, start, end);

        // Step 2: Classify variant type after trimming
        int refLen = trimmed.ref.length();
        int altLen = trimmed.alt.length();

        // Reference allele (no variation)
        if (trimmed.ref.isEmpty() && trimmed.alt.isEmpty()) {
            int referenceLength = safeToInt(end - start, "Reference allele length exceeds int range");
            return new NormalizedAllele(
                    start,
                    end,
                    NormalizedAllele.SequenceState.referenceLength(referenceLength, referenceLength)
            );
        }

        if (refLen > 0 && altLen > 0) {
            // SNP, multi-nucleotide substitution (MNP), or delins: both trimmed sides non-empty.
            // A delins (e.g. ACA -> T) is not ambiguous, so VRS 2.0 does not roll it
            // (vrs-python normalize step 2.b). All are a LiteralSequenceExpression of the
            // trimmed ALT over the trimmed interval.
            return new NormalizedAllele(
                trimmed.start,
                trimmed.end,
                NormalizedAllele.SequenceState.literalSequence(trimmed.alt)
            );
        }

        // Indel path (Steps 3-5): exactly one trimmed side is empty
        boolean isInsertion = refLen == 0;
        String seed = isInsertion ? trimmed.alt : trimmed.ref;
        int seedLength = seed.length();

        RollBounds bounds = determineRollBounds(contig, trimmed.start, seed, refContext);
        ExpansionResult expansion = expandAllele(contig, trimmed, bounds, refContext);

        // Unambiguous insertion
        if (isInsertion && bounds.left == bounds.right) {
            return new NormalizedAllele(
                    expansion.start,
                    expansion.end,
                    NormalizedAllele.SequenceState.literalSequence(expansion.alt)
            );
        }

        // Deletion (reference-derived)
        if (!isInsertion) {
            return new NormalizedAllele(
                    expansion.start,
                    expansion.end,
                    NormalizedAllele.SequenceState.referenceLength(
                            safeToInt(expansion.alt.length(), "Expanded deletion alt length exceeds int range"),
                            seedLength
                    )
            );
        }

        // Ambiguous insertion: reference-derived detection
        Integer repeatSubunitLength = findReferenceDerivedRepeatSubunitLength(
                expansion.ref,
                expansion.alt,
                seedLength
        );

        if (repeatSubunitLength != null) {
            return new NormalizedAllele(
                    expansion.start,
                    expansion.end,
                    NormalizedAllele.SequenceState.referenceLength(
                            safeToInt(expansion.alt.length(), "Expanded insertion alt length exceeds int range"),
                            repeatSubunitLength
                    )
            );
        }

        return new NormalizedAllele(
                expansion.start,
                expansion.end,
                NormalizedAllele.SequenceState.literalSequence(expansion.alt)
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Step 1: Trim common prefix and suffix
    // ─────────────────────────────────────────────────────────────────────────

    private static class TrimResult {
        final String ref;
        final String alt;
        final long start;
        final long end;

        TrimResult(String ref, String alt, long start, long end) {
            this.ref = ref;
            this.alt = alt;
            this.start = start;
            this.end = end;
        }
    }

    private TrimResult trimCommonFlanks(String ref, String alt, long start, long end) {
        // Trim common prefix
        int prefixLen = 0;
        int minLen = Math.min(ref.length(), alt.length());

        while (prefixLen < minLen && ref.charAt(prefixLen) == alt.charAt(prefixLen)) {
            prefixLen++;
        }

        // Trim common suffix (only from the remaining sequences)
        int suffixLen = 0;
        int refRemaining = ref.length() - prefixLen;
        int altRemaining = alt.length() - prefixLen;
        int minRemaining = Math.min(refRemaining, altRemaining);

        while (suffixLen < minRemaining &&
               ref.charAt(ref.length() - 1 - suffixLen) == alt.charAt(alt.length() - 1 - suffixLen)) {
            suffixLen++;
        }

        // Extract trimmed sequences
        String trimmedRef = ref.substring(prefixLen, ref.length() - suffixLen);
        String trimmedAlt = alt.substring(prefixLen, alt.length() - suffixLen);

        // Adjust coordinates
        long newStart = start + prefixLen;
        long newEnd = end - suffixLen;

        return new TrimResult(trimmedRef, trimmedAlt, newStart, newEnd);
    }

    private static class RollBounds {
        final long left;
        final long right;

        RollBounds(long left, long right) {
            this.left = left;
            this.right = right;
        }
    }

    private RollBounds determineRollBounds(String contig, long start, String seed, ReferenceContext refContext) {
        long left = start;
        long right = start;

        String rotated = seed;
        while (left > 0 && getReferenceBase(contig, left - 1, refContext) == rotated.charAt(rotated.length() - 1)) {
            left--;
            rotated = rotateRight(rotated);
        }

        rotated = seed;
        long contigLength = getContigLength(contig, refContext);
        while (right < contigLength && getReferenceBase(contig, right, refContext) == rotated.charAt(0)) {
            right++;
            rotated = rotateLeft(rotated);
        }

        return new RollBounds(left, right);
    }

    private static class ExpansionResult {
        final long start;
        final long end;
        final String ref;
        final String alt;

        ExpansionResult(long start, long end, String ref, String alt) {
            this.start = start;
            this.end = end;
            this.ref = ref;
            this.alt = alt;
        }
    }

    private ExpansionResult expandAllele(String contig, TrimResult trimmed, RollBounds bounds, ReferenceContext refContext) {
        String prefix = getRefSequence(contig, bounds.left, trimmed.start, refContext);
        // Start after the deleted bases (a no-op for insertions, where start == end)
        String suffix = getRefSequence(contig, trimmed.end, bounds.right, refContext);

        String expandedRef = prefix + trimmed.ref + suffix;
        String expandedAlt = prefix + trimmed.alt + suffix;

        return new ExpansionResult(bounds.left, bounds.right, expandedRef, expandedAlt);
    }

    private Integer findReferenceDerivedRepeatSubunitLength(String expandedRef, String expandedAlt, int seedLength) {
        for (int divisor : getDivisorsDescending(seedLength)) {
            if (divisor > expandedRef.length()) {
                continue;
            }

            for (int i = 0; i + divisor <= expandedRef.length(); i++) {
                String unit = expandedRef.substring(i, i + divisor);
                if (canBuildByCircularRepeat(unit, expandedAlt)) {
                    return divisor;
                }
            }
        }

        return null;
    }

    private static List<Integer> getDivisorsDescending(int n) {
        List<Integer> divisors = new ArrayList<>();
        for (int d = 1; d <= n; d++) {
            if (n % d == 0) {
                divisors.add(d);
            }
        }
        divisors.sort(Collections.reverseOrder());
        return divisors;
    }

    private static boolean canBuildByCircularRepeat(String unit, String target) {
        int d = unit.length();
        if (d == 0) {
            return target.isEmpty();
        }

        for (int offset = 0; offset < d; offset++) {
            boolean matches = true;
            for (int i = 0; i < target.length(); i++) {
                if (target.charAt(i) != unit.charAt((offset + i) % d)) {
                    matches = false;
                    break;
                }
            }
            if (matches) {
                return true;
            }
        }

        return false;
    }

    private static String rotateLeft(String value) {
        if (value.length() <= 1) {
            return value;
        }
        return value.substring(1) + value.charAt(0);
    }

    private static String rotateRight(String value) {
        if (value.length() <= 1) {
            return value;
        }
        return value.charAt(value.length() - 1) + value.substring(0, value.length() - 1);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Reference sequence access
    // ─────────────────────────────────────────────────────────────────────────

    /**
     * Get reference sequence via the per-variant ReferenceContext.
     *
     * @param contig     Chromosome/contig name
     * @param start      0-based interbase start (VRS coordinates)
     * @param end        0-based interbase end (VRS coordinates, half-open)
     * @param refContext ReferenceContext from the current apply() call
     * @return Reference sequence in uppercase
     */
    private String getRefSequence(String contig, long start, long end, ReferenceContext refContext) {
        if (start >= end) {
            return "";  // Empty interval (insertion site)
        }
        // ReferenceContext uses 1-based inclusive coordinates.
        // VRS 0-based [start, end) → 1-based [start+1, end]
        final SimpleInterval interval = new SimpleInterval(contig, (int)(start + 1), (int) end);
        final byte[] bases = refContext.getBases(interval);
        return new String(bases).toUpperCase();
    }

    /**
     * One reference base, fetched on its own.
     *
     * <p>Rolling ({@link #determineRollBounds}) calls this once per base it moves, so an indel in a
     * repeat makes one {@code getBases} call per base of the repeat. That is deliberately left simple.
     * The alternative, fetching a window around the variant and extending it when a roll reaches its
     * edge, adds window-boundary logic to the part of normalization that decides the VRS ID, where a
     * mistake would change IDs rather than speed.
     *
     * <p>Each call is an in-memory copy, not file I/O. With {@code --reference}, the walker's
     * {@code ReferenceContext} reads from a {@code ReferenceFileSource}, whose
     * {@link org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile} holds the last
     * 1 Mb of the contig in memory ({@code DEFAULT_CACHE_SIZE}). Variants arrive in position order, so
     * the FASTA is read about once per megabase, and a per-base call costs a few small allocations: two
     * intervals, a one-byte array copy and a String.
     *
     * <p>What we know of the cost (VS-2058): on the 25 integration samples, {@code CreateVariantIngestFiles}
     * took 167 minutes summed with VRS IDs on and 129 with them off, about 30% more. That is one run
     * of each, and the 30% also covers trimming, SHA-512 digests and writing the VRS rows, so how much
     * of it is reference fetching has not been measured. Profile before replacing this with a windowed
     * fetch, and check any replacement against the vrs-python concordance test.
     */
    private char getReferenceBase(String contig, long position, ReferenceContext refContext) {
        String base = getRefSequence(contig, position, position + 1, refContext);
        if (base.isEmpty()) {
            throw new GATKException("Reference base lookup returned empty sequence at " + contig + ":" + position);
        }
        return base.charAt(0);
    }

    private long getContigLength(String contig, ReferenceContext refContext) {
        final htsjdk.samtools.SAMSequenceRecord seq = refContext.getSequenceRecord();
        if (seq == null || !seq.getSequenceName().equals(contig)) {
            // Fall back: look up by name in the sequence dictionary via a temporary context call
            // getSequenceRecord() only returns the record for the context's own contig;
            // for cross-contig lookups during rolling we use the same contig, so this is fine.
            throw new GATKException("Cannot determine length for contig: '" + contig +
                "' from ReferenceContext on '" + refContext.getContig() + "'");
        }
        return seq.getSequenceLength();
    }

    private static int safeToInt(long value, String message) {
        if (value > Integer.MAX_VALUE || value < Integer.MIN_VALUE) {
            throw new GATKException(message + ": " + value);
        }
        return (int) value;
    }
}
