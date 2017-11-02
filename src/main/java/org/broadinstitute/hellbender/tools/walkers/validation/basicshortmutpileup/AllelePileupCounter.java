package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Useful when you know the interval and the alleles of interest ahead of the counting.
 *
 * Supports simple Indels, SNP, and MNP.
 *
 */
public class AllelePileupCounter {

    private final Allele referenceAllele;
    private final List<Allele> alternateAlleles;
    private int minBaseQualityCutoff;
    private Map<Allele, MutableInt> countMap = new HashMap<>();

    /**
     *
     * @param referenceAllele allele to treat as reference.  Create with {@link Allele#create(String, boolean)}, where
     *                  second parameter is {@code true}.  Never {@code null}.   If the reference is symbolic, exception will be thrown.
     * @param alternateAlleles List of alleles to treat as the alternates.  Easy to create with {@link Allele#create(String, boolean)}, where
     *                  second parameter is {@code false}.  Never {@code null}
     * @param minBaseQualityCutoff minimum base quality for the bases that match the allele in order to be counted.
     *                             Must be positive or zero.
     */
    public AllelePileupCounter(final Allele referenceAllele, final List<Allele> alternateAlleles, int minBaseQualityCutoff) {

        this.referenceAllele = Utils.nonNull(referenceAllele);
        this.alternateAlleles = Utils.nonNull(alternateAlleles);

        // Additional checks
        if (referenceAllele.isSymbolic()) {
            throw new UserException.BadInput("A symbolic reference allele was specified.");
        }

        Utils.validateArg(!referenceAllele.isNonReference(), "Reference allele was non-reference: " + referenceAllele);
        Utils.validateArg(alternateAlleles.stream().allMatch(a -> a.isNonReference()),
                "One or more alternate alleles were reference: " + alternateAlleles.stream().map(a-> a.toString()).collect(Collectors.joining(", ")));

        this.minBaseQualityCutoff = ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Minimum base quality must be positive or zero.");;

        alternateAlleles.forEach(a -> countMap.put(a, new MutableInt(0)));
        countMap.put(referenceAllele, new MutableInt(0));
    }

    /**
     * See {@link AllelePileupCounter#AllelePileupCounter(Allele, List, int)}, except this will immediately add the given pileup as well.
     *
     * @param referenceAllele See {@link AllelePileupCounter#AllelePileupCounter(Allele, List, int)}
     * @param alternateAlleles See {@link AllelePileupCounter#AllelePileupCounter(Allele, List, int)}
     * @param minBaseQualityCutoff See {@link AllelePileupCounter#AllelePileupCounter(Allele, List, int)}
     * @param readPileup Initial pileup to add to this allele pileup counter
     */
    public AllelePileupCounter(final Allele referenceAllele, final List<Allele> alternateAlleles, int minBaseQualityCutoff, final ReadPileup readPileup) {
        this(referenceAllele, alternateAlleles, minBaseQualityCutoff);
        this.addPileup(readPileup);
    }

    /**
     * @param pileup read pileup to base the counts of each allele.  If null, it will not be considered.
     */
    public void addPileup(final ReadPileup pileup) {

        if ((pileup != null) && !referenceAllele.isSymbolic()) {
            Utils.stream(pileup)
                    .filter(pe -> isUsableRead(pe.getRead()))
                    .forEach(pe -> incrementAlleleCountMap(pe, referenceAllele, alternateAlleles, minBaseQualityCutoff, countMap));
        }
    }

    private static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    /**
     * Increment the given count map with one of the alleles passed in (reference or alt).  Note that this method
     * will not add keys to the count map.
     *
     * @param pileupElement Read to consider.  Never {@code null}
     * @param referenceAllele Never {@code null}
     * @param altAlleles Never {@code null}
     * @param minBaseQualityCutoff minimum base quality for the bases that match the allele in order to be counted.
     *                             Must be positive or zero.
     * @param countMap modified in place with an increment (if applicable) for the allele that is seen in the read.
     */
    private static void incrementAlleleCountMap(final PileupElement pileupElement, final Allele referenceAllele,
                                                final List<Allele> altAlleles, int minBaseQualityCutoff,
                                                final Map<Allele, MutableInt> countMap) {
        ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Minimum base quality must be positive or zero.");
        Utils.nonNull(countMap);
        Utils.nonNull(referenceAllele);
        Utils.nonNull(altAlleles);
        Utils.nonNull(pileupElement);

        final Allele pileupAllele = GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele, altAlleles, minBaseQualityCutoff);

        if ((pileupAllele != null) && (countMap.containsKey(pileupAllele))) {
            countMap.get(pileupAllele).increment();
        }
    }

    public Map<Allele, MutableInt> getCountMap() {
        return countMap;
    }
}
