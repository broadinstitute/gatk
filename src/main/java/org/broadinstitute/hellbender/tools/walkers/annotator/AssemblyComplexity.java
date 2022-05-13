package org.broadinstitute.hellbender.tools.walkers.annotator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.mutable.MutableInt;
import org.broadinstitute.barclay.argparser.Argument;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY,
        summary="Describe the complexity of an assembly region")
public class AssemblyComplexity implements JumboInfoAnnotation {

    @Argument(fullName = "assembly-complexity-reference-mode",
            doc="If enabled will treat the reference as the basis for assembly complexity as opposed to estimated germline haplotypes",
            optional=true)
    public boolean germlineMode = false;

    public AssemblyComplexity() { }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final FeatureContext features,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                        final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods,
                                        final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods) {
        final Triple<int[], int[], double[]> annotations = annotate(vc, haplotypeLikelihoods);
        final Map<String, Object> result = new HashMap<>();

        result.put(GATKVCFConstants.HAPLOTYPE_EQUIVALENCE_COUNTS_KEY , annotations.getLeft());
        result.put(GATKVCFConstants.HAPLOTYPE_COMPLEXITY_KEY , annotations.getMiddle());
        result.put(GATKVCFConstants.HAPLOTYPE_DOMINANCE_KEY , annotations.getRight());
        return result;
    }

    public Triple<int[], int[], double[]> annotate(final VariantContext vc, final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods) {

        // count best-read support for each haplotype
        final Map<Haplotype, MutableInt> haplotypeSupportCounts = haplotypeLikelihoods.alleles().stream()
                .collect(Collectors.toMap(hap -> hap, label -> new MutableInt(0)));
        haplotypeLikelihoods.bestAllelesBreakingTies()
                .forEach(bestHaplotype -> haplotypeSupportCounts.get(bestHaplotype.allele).increment());

        // encode each haplotype as a string of variant starts and alt allele strings, excluding the locus of vc (to avoid reference/germline bias)
        // note that VariantContexts in an EventMap are always biallelic, so var.getAlternateAllele(0) is valid
        final Map<String, List<Haplotype>> haplotypeGroups = haplotypeLikelihoods.alleles().stream()
                .collect(Collectors.groupingBy(hap -> hap.getEventMap().getVariantContexts().stream()
                        .filter(var -> var.getStart() != vc.getStart())
                        .map(var -> var.getStart() + var.getAlternateAllele(0).getBaseString())
                        .collect(Collectors.joining())));

        // sum the read support counts for all haplotypes within each group
        final int[] equivalenceCounts = haplotypeGroups.values().stream()
                .map(haps -> haps.stream().map(haplotypeSupportCounts::get).mapToInt(MutableInt::intValue).sum())
                .sorted(Comparator.reverseOrder())
                .mapToInt(n->n)
                .toArray();

        // we're going to calculate the complexity of this variant's haplotype (that is, the variant-supporting haplotype
        // with the most reads) versus the closest (in terms of edit distance) germline haplotype.  The haplotype
        // with the greatest read support is considered germline, and as a heuristic we consider the second-most-supported
        // haplotype to be germline as well if it is at least half as supported (in terms of best read count) as the most-supported.
        // as above, we exclude differences at the variant site itself to avoid reference and germline bias
        final List<Haplotype> haplotypesByDescendingSupport = haplotypeSupportCounts.entrySet().stream()
                .sorted(Comparator.comparingInt(entry -> -entry.getValue().intValue()))
                .map(entry -> entry.getKey())
                .collect(Collectors.toList());

        final List<Haplotype> germlineHaplotypes;
        if (germlineMode) {
            germlineHaplotypes = Collections.singletonList(haplotypeLikelihoods.getAllele(haplotypeLikelihoods.indexOfReference()));
        } else {
            germlineHaplotypes = new ArrayList<>();
            germlineHaplotypes.add(haplotypesByDescendingSupport.get(0));
            if (haplotypesByDescendingSupport.size() > 1 && haplotypeSupportCounts.get(haplotypesByDescendingSupport.get(1)).intValue() >= haplotypeSupportCounts.get(haplotypesByDescendingSupport.get(0)).intValue() / 2) {
                germlineHaplotypes.add(haplotypesByDescendingSupport.get(1));
            }
        }

        final int[] editDistances = IntStream.range(0, vc.getNAlleles() - 1).map(altAlleleIndex -> {
            if (vc.getAlternateAllele(altAlleleIndex).isSymbolic() || vc.getAlternateAllele(altAlleleIndex).getBases()[0] == '*') {
                return 0;
            }
            final Haplotype mostSupportedHaplotypeWithAllele = haplotypesByDescendingSupport.stream()
                    .filter(hap -> containsAltAllele(hap.getEventMap(), vc, altAlleleIndex))
                    .findFirst().get();

            return germlineHaplotypes.stream().mapToInt(gh -> editDistance(gh, mostSupportedHaplotypeWithAllele, vc.getStart())).min().getAsInt();
        }).toArray();

        // measure which proportion of reads supporting each alt allele fit the most-supported haplotype for that allele
        final double[] haplotypeDominance = IntStream.range(0, vc.getNAlleles() - 1).mapToDouble(altAlleleIndex -> {
            if (vc.getAlternateAllele(altAlleleIndex).isSymbolic() || vc.getAlternateAllele(altAlleleIndex).getBases()[0] == '*') {
                return 0;
            }
            final int[] counts = haplotypesByDescendingSupport.stream()
                    .filter(hap -> containsAltAllele(hap.getEventMap(), vc, altAlleleIndex))
                    .mapToInt(hap -> haplotypeSupportCounts.get(hap).intValue())
                    .toArray();
            return MathUtils.arrayMax(counts) / (double) MathUtils.sum(counts);
        }).toArray();

        return Triple.of(equivalenceCounts, editDistances, haplotypeDominance);
    }


    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.HAPLOTYPE_EQUIVALENCE_COUNTS_KEY, GATKVCFConstants.HAPLOTYPE_COMPLEXITY_KEY, GATKVCFConstants.HAPLOTYPE_DOMINANCE_KEY);
    }

    // does an EventMap contain a variant allele
    // we assume that everything is derived from a common assembly and thus variant start positions
    // are consistent.  However, the variant might not have the same reference allele; for example if one
    // haplotype had a CACAC -> C and another had a CAC -> C, then in the merged variant context the latter would
    // be represented as CACAC -> CAC
    // since the EventMap is biallelic its alleles are as reduced as possible; hence if the variant context's reference is
    // k bases longer we simply check that the event map alt allele matches the variant context alt allele excluding the
    // latter's last k bases
    private static boolean containsAltAllele(final EventMap eventMap, final VariantContext vc, final int altAlleleIndex) {
        final List<VariantContext> overlapping =  eventMap.getOverlappingEvents(vc.getStart());
        if (overlapping.isEmpty()) {
            return false;
        } else if (overlapping.get(0).getStart() != vc.getStart()) {
            return false;
        } else {
            final VariantContext eventMapVC = overlapping.get(0);
            final int excessBases = vc.getReference().length() - eventMapVC.getReference().length();

            return equalBasesExcludingSuffix(eventMapVC.getAlternateAllele(0).getBases(),
                    vc.getAlternateAllele(altAlleleIndex).getBases(), excessBases);
        }
    }

    // does a minimally-represented EventMap allele match a VariantContextAllele excluding a certain number of suffix bases
    // in the latter
    private static boolean equalBasesExcludingSuffix(final byte[] eventMapBases, final byte[] variantContextBases, final int suffixSize) {
        if (eventMapBases.length + suffixSize != variantContextBases.length) {
            return false;
        } else if (eventMapBases.length > variantContextBases.length) {
            return false;   // edge case -- event map is longer, though minimal, because it is a MNP
                            // even if the leading bases match, let's call this not a match
        } else {
            for (int n = 0; n < eventMapBases.length; n++) {
                if (eventMapBases[n] != variantContextBases[n]) {
                    return false;
                }
            }
            return true;
        }

    }

    // count variants in one haplotype but not the other
    // note that we use the fact that EventMap VariantContexts are biallelic
    private static int uniqueVariants(final Haplotype hap1, final Haplotype hap2, final int excludedPosition) {
        final EventMap eventMap2 = hap2.getEventMap();
        return (int) hap1.getEventMap().getVariantContexts().stream()
                .filter(vc -> vc.getStart() != excludedPosition)
                .filter(vc -> !containsAltAllele(eventMap2, vc, 0))
                .count();
    }

    private static int editDistance(final Haplotype hap1, final Haplotype hap2, final int excludedPosition) {
        return uniqueVariants(hap1, hap2, excludedPosition) + uniqueVariants(hap2, hap1, excludedPosition);
    }

}
