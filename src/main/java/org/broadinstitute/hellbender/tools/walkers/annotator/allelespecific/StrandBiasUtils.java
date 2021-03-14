package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Common strand bias utilities used by allele specific strand bias annotators
 */
public class StrandBiasUtils {
    public static final int FORWARD = 0;
    public static final int REVERSE = 1;
    public static final int MIN_COUNT = 2;
    private static final List<Integer> ZERO_LIST = new ArrayList<>(Arrays.asList(0,0));

    public static Map<String, Object> computeSBAnnotation(VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods, String key) {
        // calculate the annotation from the likelihoods
        // likelihoods can come from HaplotypeCaller or Mutect2 call to VariantAnnotatorEngine
        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
        getStrandCountsFromLikelihoodMap(vc, likelihoods, myData, MIN_COUNT);
        final Map<Allele, List<Integer>> perAlleleValues = myData.getAttributeMap();
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), perAlleleValues);
        annotations.put(key, annotationString);
        return annotations;
    }

    /**
     * Helper method to output raw allele-specific strand counts as a string
     * @param vcAlleles relevant alleles
     * @param perAlleleValues forward and reverse read counts for each allele
     * @return a String appropriate to use for annotating a GVCF
     */
    public static String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Integer>> perAlleleValues) {
        String annotationString = "";
        for (final Allele a : vcAlleles) {
            if (!annotationString.isEmpty()) {
                annotationString += AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM;
            }
            List<Integer> alleleValues = perAlleleValues.get(a);
            if (alleleValues == null) {
                alleleValues = ZERO_LIST;
            }
            annotationString += encode(alleleValues);
        }
        return annotationString;
    }

    /**
     * Helper method to convert StrandBias values to string format
     * @param alleleValues values to  format
     * @return formatted string
     */
    public static String encode(List<Integer> alleleValues) {
        String annotationString = "";
        for (int j =0; j < alleleValues.size(); j++) {
            annotationString += alleleValues.get(j);
            if (j < alleleValues.size()-1) {
                annotationString += ",";
            }
        }
        return annotationString;
    }

    /**
     Allocate and fill a Nx2 strand contingency table where N is the number of alleles.  In the end, it'll look something like this:
     *             fwd      rev
     *   allele1   #       #
     *   allele2   #       #
     *
     *   NOTE:Only use informative reads
     *
     * @param vc    VariantContext from which to get alleles
     * @param likelihoods per-read allele likelihoods to determine if each read is informative
     * @param perAlleleValues modified to store the output counts
     * @param minCount minimum threshold of counts to use
     */
    public static void getStrandCountsFromLikelihoodMap( final VariantContext vc,
                                                  final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                  final ReducibleAnnotationData<List<Integer>> perAlleleValues,
                                                  final int minCount) {
        if( likelihoods == null || vc == null ) {
            return;
        }

        final Allele ref = vc.getReference();
        final List<Allele> allAlts = vc.getAlternateAlleles();

        for (final String sample : likelihoods.samples()) {
            final ReducibleAnnotationData<List<Integer>> sampleTable = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
            likelihoods.bestAllelesBreakingTies(sample).stream()
                    .filter(ba -> ba.isInformative())
                    .forEach(ba -> updateTable(ba.allele, ba.evidence, ref, allAlts, sampleTable));
            if (passesMinimumThreshold(sampleTable, minCount)) {
                combineAttributeMap(sampleTable, perAlleleValues);
            }
        }
    }

    /**
     * Combine allele-specific data from two ReducibleAnnotationData data structures
     * @param toAdd input values
     * @param combined  modified to return the combined values
     */
    protected static void combineAttributeMap(final ReducibleAnnotationData<List<Integer>> toAdd, final ReducibleAnnotationData<List<Integer>> combined) {
        for (final Allele a : combined.getAlleles()) {
            if (toAdd.hasAttribute(a) && toAdd.getAttribute(a) != null) {
                if (combined.getAttribute(a) != null) {
                    combined.getAttribute(a).set(FORWARD, (int) combined.getAttribute(a).get(FORWARD) + (int) toAdd.getAttribute(a).get(FORWARD));
                    combined.getAttribute(a).set(REVERSE, (int) combined.getAttribute(a).get(REVERSE) + (int) toAdd.getAttribute(a).get(REVERSE));
                }
                else {
                    List<Integer> alleleData = new ArrayList<>();
                    alleleData.add(FORWARD, toAdd.getAttribute(a).get(FORWARD));
                    alleleData.add(REVERSE, toAdd.getAttribute(a).get(REVERSE));
                    combined.putAttribute(a,alleleData);
                }
            }
        }
    }

    /**
     * Add another read to the strand count table
     * @param bestAllele    the Allele best supported by {@param read}
     * @param read  read to add
     * @param ref   reference Allele
     * @param allAlts   the (subset of) alternate alleles for the associated variant
     * @param perAlleleValues updated to return the values including @{value read}
     */
    private static void updateTable(final Allele bestAllele, final GATKRead read, final Allele ref, final List<Allele> allAlts, final ReducibleAnnotationData<List<Integer>> perAlleleValues) {

        final boolean matchesRef = bestAllele.equals(ref, true);
        final boolean matchesAnyAlt = allAlts.contains(bestAllele);

        //can happen if a read's most likely allele has been removed when --max_alternate_alleles is exceeded
        if (!( matchesRef || matchesAnyAlt )) {
            return;
        }

        final List<Integer> alleleStrandCounts;
        if (perAlleleValues.hasAttribute(bestAllele) && perAlleleValues.getAttribute(bestAllele) != null) {
            alleleStrandCounts = perAlleleValues.getAttribute(bestAllele);
        } else {
            alleleStrandCounts = new ArrayList<>();
            alleleStrandCounts.add(0,0);
            alleleStrandCounts.add(1,0);
        }

        final int strand = read.isReverseStrand() ? REVERSE : FORWARD;
        alleleStrandCounts.set(strand, alleleStrandCounts.get(strand) + 1);
        perAlleleValues.putAttribute(bestAllele, alleleStrandCounts);
    }

    /**
     * Does this strand data array pass the minimum threshold for inclusion?
     *
     * @param sampleTable  the per-allele fwd/rev read counts for a single sample
     * @param minCount The minimum threshold of counts in the array
     * @return true if it passes the minimum threshold, false otherwise
     */
    protected static boolean passesMinimumThreshold(final ReducibleAnnotationData<List<Integer>> sampleTable, final int minCount) {
        final int readCount = sampleTable.getAttributeMap().values().stream()
                .filter(alleleValues -> alleleValues != null)
                .mapToInt(alleleValues -> alleleValues.get(FORWARD) + alleleValues.get(REVERSE))
                .sum();
        return readCount > minCount;
    }


    /**
     * Get an allele-specific strand bias contingency table from a VariantContext
     * @param vc
     * @return a contingency table with forward and reverse counts for each allele
     */
    public static List<List<Integer>> getSBsForAlleles(VariantContext vc) {
        String sbStr = vc.getCommonInfo().getAttributeAsString(GATKVCFConstants.AS_SB_TABLE_KEY, null);
        if (sbStr == null || sbStr.isEmpty()) {
            return Collections.emptyList();
        }
        List<String> asb = AnnotationUtils.decodeAnyASListWithRawDelim(sbStr);
        return asb.stream()
                .map(fwdrev -> AnnotationUtils.decodeAnyASList(fwdrev).stream().map(String::trim)
                .mapToInt(Integer::parseInt).boxed().collect(Collectors.toList())).collect(Collectors.toList());

    }
}
