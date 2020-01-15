package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

public class StrandBiasUtils {
    public static final int FORWARD = 0;
    public static final int REVERSE = 1;
    public static final int MIN_COUNT = 2;
    public static final String PRINT_DELIM = "|";
    private static final List<Integer> ZERO_LIST = new ArrayList<>(Arrays.asList(0,0));

    public static Map<String, Object> computeSBAnnotation(VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods, String key) {
        // calculate the annotation from the likelihoods
        // likelihoods can come from HaplotypeCaller call to VariantAnnotatorEngine
        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
        getStrandCountsFromLikelihoodMap(vc, likelihoods, myData, MIN_COUNT);
        final Map<Allele, List<Integer>> perAlleleValues = myData.getAttributeMap();
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), perAlleleValues);
        annotations.put(key, annotationString);
        return annotations;
    }

    protected static String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Integer>> perAlleleValues) {
        String annotationString = "";
        for (final Allele a : vcAlleles) {
            if (!annotationString.isEmpty()) {
                annotationString += PRINT_DELIM;
            }
            List<Integer> alleleValues = perAlleleValues.get(a);
            if (alleleValues == null) {
                alleleValues = ZERO_LIST;
            }
            annotationString += encode(alleleValues);
        }
        return annotationString;
    }

    protected static String encode(List<Integer> alleleValues) {
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
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
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
        final boolean isForward = !read.isReverseStrand();
        if (isForward) {
            alleleStrandCounts.set(FORWARD, alleleStrandCounts.get(FORWARD) + 1);
        } else {
            alleleStrandCounts.set(REVERSE, alleleStrandCounts.get(REVERSE) + 1);
        }
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



}
