package org.broadinstitute.hellbender.tools.copynumber.datacollection;

import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Collects reference/alternate allele counts at specified sites.  The alt count is defined as the total count minus the ref count,
 * and the alt nucleotide is defined as the non-ref base with the highest count, with ties broken by the order of the
 * bases in {@link AllelicCountCollector#BASES}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCountCollector {

    private static final Logger logger = LogManager.getLogger(AllelicCountCollector.class);

    public static final List<Nucleotide> BASES = Collections.unmodifiableList(Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T));

    private final SampleLocatableMetadata metadata;
    private final List<AllelicCount> allelicCounts = new ArrayList<>();

    public AllelicCountCollector(final SampleLocatableMetadata metadata) {
        this.metadata = Utils.nonNull(metadata);
    }

    /**
     * Add counts to this class for a specific locus.
     *
     * @param refBase single nucleotide of the reference.  Not {@code null}
     * @param pileup associated pileup at the locus.  Not {@code null}
     * @param locus position in genome to collect alellic counts.  Not {@code null}
     * @param minBaseQuality minimum base quality in the read for that read to count at that position.  Must be greater than or equal to 0.
     */
    public void collectAtLocus(final Nucleotide refBase, final ReadPileup pileup, final Locatable locus, final int minBaseQuality) {
        Utils.nonNull(refBase);
        Utils.nonNull(pileup);
        Utils.nonNull(locus);
        ParamUtils.isPositiveOrZero(minBaseQuality, "Minimum base quality must be zero or higher.");

        if (!BASES.contains(refBase)) {
            logger.warn(String.format("The reference position at %s has an unknown base call (value: %s). Skipping...",
                    locus, refBase.toString()));
            return;
        }

        final Nucleotide.Counter nucleotideCounter = new Nucleotide.Counter();

        Utils.stream(pileup.iterator())
                .filter(r -> !r.isDeletion())
                .filter(r -> r.getQual() >= minBaseQuality)
                .forEach(r -> nucleotideCounter.add(r.getBase()));

        final int totalBaseCount = BASES.stream().mapToInt(b -> (int) nucleotideCounter.get(b)).sum();  //only include total ACGT counts (exclude N, etc.)
        final int refReadCount = (int) nucleotideCounter.get(refBase);
        final int altReadCount = totalBaseCount - refReadCount;                                         //we take alt = total - ref instead of the actual alt count
        final Nucleotide altBase = altReadCount == 0 ? Nucleotide.N : inferAltFromPileupBaseCounts(nucleotideCounter, refBase);

        allelicCounts.add(new AllelicCount(
                new SimpleInterval(locus.getContig(), locus.getStart(), locus.getEnd()),
                refReadCount, altReadCount, refBase, altBase));
    }

    /**
     * Get the allelic counts gathered so far.
     *
     * @return a <em>reference</em> to the AllelicCountCollection
     */
    public AllelicCountCollection getAllelicCounts() {
        return new AllelicCountCollection(metadata, allelicCounts);
    }

    /**
     * Returns the non-ref base with highest count (if there is a tie, the first base in the order given in
     * {@link AllelicCountCollector#BASES} will be returned).
     */
    private static Nucleotide inferAltFromPileupBaseCounts(final Nucleotide.Counter baseCounts,
                                                           final Nucleotide refNucleotide) {
        return BASES.stream()
                .filter(b -> b != refNucleotide)
                .sorted((b1, b2) -> Long.compare(baseCounts.get(b2), baseCounts.get(b1)))
                .findFirst().get();
    }
}
