package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Class for collapsing a collection of similar {@link SVCallRecord} objects, such as clusters produced by
 * {@link CanonicalSVLinkage}, into a single representative call.
 */
public class CanonicalSVCollapser implements SVCollapser<SVCallRecord, BasicOutputCluster<SVCallRecord>> {

    /**
     * Define strategies for collapsing variant intervals.
     */
    public enum BreakpointSummaryStrategy {
        /**
         * Use the (first) middle value to summarize cluster starts and ends, such that the start and end were seen in the data
         */
        MEDIAN_START_MEDIAN_END,

        /**
         * A conservative strategy to summarize a cluster by its smallest extent
         */
        MIN_START_MAX_END,

        /**
         * A permissive strategy to summarize a cluster by it largest extent
         */
        MAX_START_MIN_END,

        /**
         * Summarize a cluster using the mean value for each end, even if that value was not represented in any sample
         */
        MEAN_START_MEAN_END

    }

    /**
     * Define strategies for collapsing alt alleles with different subtypes.
     */
    public enum AltAlleleSummaryStrategy {
        /**
         * Use the most specific subtype that doesn't conflict with any of the other alleles.
         * For example, (&lt;INS&gt, &lt;INS:MEI:SVA&gt, &lt;INS:MEI:LINE&gt) results in &lt;INS:MEI&gt.
         */
        MOST_SPECIFIC_SUBTYPE,

        /**
         * Use subtypes in common among all alleles.
         * For example, (&lt;INS&gt, &lt;INS:MEI:SVA&gt, &lt;INS:MEI&gt) results in &lt;INS&gt.
         */
        COMMON_SUBTYPE

    }

    /**
     * Strategies for collapsing the insertion LENGTH attribute. Only applied when there are multiple defined values.
     */
    public enum InsertionLengthSummaryStrategy {
        MEDIAN,
        MEAN,
        MIN,
        MAX,
        UNDEFINED
    }

    // TODO Add support for complex (CPX)
    private static final Set<StructuralVariantType> SUPPORTED_SV_TYPES = Sets.newHashSet(
            StructuralVariantType.DEL,
            StructuralVariantType.DUP,
            StructuralVariantType.CNV,
            StructuralVariantType.INS,
            StructuralVariantType.INV,
            StructuralVariantType.BND
    );

    private final AltAlleleSummaryStrategy altAlleleSummaryStrategy;
    private final BreakpointSummaryStrategy breakpointSummaryStrategy;
    private final InsertionLengthSummaryStrategy insertionLengthSummaryStrategy;
    private final ReferenceSequenceFile reference;
    private final SAMSequenceDictionary dictionary;
    private final Comparator<SVLocatable> comparator;

    private static final AlleleCollectionCollapserComparator ALLELE_COMPARATOR = new AlleleCollectionCollapserComparator();

    public CanonicalSVCollapser(final ReferenceSequenceFile reference,
                                final AltAlleleSummaryStrategy altAlleleSummaryStrategy,
                                final BreakpointSummaryStrategy breakpointSummaryStrategy,
                                final InsertionLengthSummaryStrategy insertionLengthSummaryStrategy) {
        this.reference = Utils.nonNull(reference);
        this.dictionary = reference.getSequenceDictionary();
        this.altAlleleSummaryStrategy = altAlleleSummaryStrategy;
        this.breakpointSummaryStrategy = breakpointSummaryStrategy;
        this.insertionLengthSummaryStrategy = insertionLengthSummaryStrategy;
        this.comparator = SVCallRecordUtils.getSVLocatableComparator(dictionary);
    }

    @Override
    public SVCallRecord collapse(final BasicOutputCluster<SVCallRecord> cluster) {
        final Collection<SVCallRecord> items = cluster.getMembers();
        validateRecords(items);
        final String id = collapseIds(items);
        final List<String> algorithms = collapseAlgorithms(items);
        final StructuralVariantType type = collapseTypes(items);
        final Map<String, Object> attributes = collapseVariantAttributes(items);

        // Prefer using variants generated with PESR callers, which tend to generate more precise breakpoints
        final Collection<SVCallRecord> mostPreciseCalls = getRecordsWithMostPreciseBreakpoints(items);
        final SVCallRecord exampleCall = mostPreciseCalls.iterator().next();

        final Pair<Integer, Integer> coordinates = collapseInterval(mostPreciseCalls);
        final int start = coordinates.getKey();
        final int end = coordinates.getValue();
        final Integer length = collapseLength(mostPreciseCalls, type);

        final Allele refAllele = collapseRefAlleles(exampleCall.getContigA(), start);

        final List<Allele> altAlleles = collapseAltAlleles(items);
        final int numAlleles = 1 + altAlleles.size();
        final List<Allele> alleles = new ArrayList<>(numAlleles);
        alleles.add(refAllele);
        alleles.addAll(altAlleles);

        final List<Genotype> genotypes = collapseAllGenotypes(items, refAllele, altAlleles);
        final List<Genotype> harmonizedGenotypes = harmonizeAltAlleles(altAlleles, genotypes);

        return new SVCallRecord(id, exampleCall.getContigA(), start, exampleCall.getStrandA(), exampleCall.getContigB(),
                end, exampleCall.getStrandB(), type, length, algorithms, alleles, harmonizedGenotypes, attributes, dictionary);
    }

    /**
     * Asserts that the given records are valid for collapsing.
     */
    protected void validateRecords(final Collection<SVCallRecord> records) {
        for (final SVCallRecord r : records) {
        Utils.validateArg(SUPPORTED_SV_TYPES.contains(r.getType()),
                "Unsupported SV type: " + r.getType());
        }
    }

    protected List<Genotype> harmonizeAltAlleles(final List<Allele> sortedAltAlleles, final List<Genotype> collapsedGenotypes) {
        Utils.nonNull(sortedAltAlleles);
        Utils.nonNull(collapsedGenotypes);
        final Set<Allele> genotypeAltAlleles = collapsedGenotypes.stream()
                .map(Genotype::getAlleles)
                .flatMap(List::stream)
                .filter(SVCallRecordUtils::isAltAllele)
                .collect(Collectors.toSet());
        // Alt alleles match already, or some alts vanished in genotype collapsing (and we keep them)
        if (sortedAltAlleles.containsAll(genotypeAltAlleles)) {
            return collapsedGenotypes;
        }
        // One or more subtypes were collapsed - need to replace genotype alt alleles
        Utils.validate(sortedAltAlleles.size() == 1, "Multi-allelic variants with subtyped alleles are " +
                "not supported.");
        final Allele newAlt = sortedAltAlleles.get(0);
        return collapsedGenotypes.stream()
                .map(g -> new GenotypeBuilder(g).alleles(replaceAltAlleles(g.getAlleles(), newAlt)).make())
                .collect(Collectors.toList());
    }

    private List<Allele> replaceAltAlleles(final List<Allele> alleles, final Allele replacement) {
        return alleles.stream().map(a -> SVCallRecordUtils.isAltAllele(a) ? replacement : a).collect(Collectors.toList());
    }

    /**
     * Finds representative expected copy number for a sample. This method simply checks whether this is defined in all
     * genotypes and verifies that there is only one unique value and returns it.
     * @param sampleGenotypes genotypes belonging to a single sample
     * @return representative copy number for the reference state
     */
    protected int collapseExpectedCopyNumber(final Collection<Genotype> sampleGenotypes) {
        Utils.nonNull(sampleGenotypes);
        Utils.nonEmpty(sampleGenotypes);
        final List<Integer> expectedCopyNumberValues = sampleGenotypes.stream()
                .map(SVCallRecord::getExpectedCopyNumber)
                .distinct()
                .collect(Collectors.toList());
        Utils.validate(expectedCopyNumberValues.size() == 1,
                "Expected 1 unique expected copy number but found " + expectedCopyNumberValues.size());
        return expectedCopyNumberValues.get(0);
    }

    /**
     * Returns the ref allele at the given locus.
     */
    protected Allele collapseRefAlleles(final String contig, final int pos) {
        final byte[] bases = ReferenceUtils.getRefBaseAtPosition(reference, contig, pos);
        Utils.validate(bases != null && bases.length == 1, "Invalid reference locus " + contig + ":" + pos);
        return Allele.create(bases[0], true);
    }

    /**
     * Collapses alternate alleles into a list of representative alleles. Note this supports sub-typed alleles such as
     * &lt;INS:MEI&gt;. If multiple alt alleles are found, the variant must either be a CNV or sub-typed alleles with the
     * same base symbol (e.g. &lt;INS:MEI&gt; and &lt;INS:MEI:SVA&gt; would result in &lt;INS&gt;).
     * @param items records whose alt alleles should be collapsed
     * @return collapsed alt alleles
     */
    protected List<Allele> collapseAltAlleles(final Collection<SVCallRecord> items) {
        Utils.nonNull(items);
        final List<Allele> altAlleles = items.stream().map(SVCallRecord::getAltAlleles)
                .flatMap(List::stream)
                .distinct()
                .sorted()
                .collect(Collectors.toList());
        if (altAlleles.isEmpty()) {
            return Collections.emptyList();
        } else if (altAlleles.size() == 1) {
            return Collections.singletonList(altAlleles.get(0));
        } else {
            // Multiple non-ref alleles need collapsing
            // TODO does not search for subtypes e.g. <DUP:TANDEM>
            if (altAlleles.size() == 2 && altAlleles.contains(Allele.SV_SIMPLE_DEL) && altAlleles.contains(Allele.SV_SIMPLE_DUP)) {
                // CNVs
                return altAlleles;
            }

            final String[] collapsedAlleleTokens;
            if (altAlleleSummaryStrategy == AltAlleleSummaryStrategy.COMMON_SUBTYPE) {
                collapsedAlleleTokens = collapseAltAllelesCommon(altAlleles);
            } else if (altAlleleSummaryStrategy == AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE) {
                collapsedAlleleTokens = collapseAltAllelesMostSpecific(altAlleles);
            } else {
                throw new UnsupportedOperationException("Unimplemented alt allele summary strategy: " + altAlleleSummaryStrategy.name());
            }
            Utils.validate(collapsedAlleleTokens.length > 0, "Encountered multiple symbolic allele base symbols for non-CNV");
            return Collections.singletonList(Allele.create("<" + String.join(":", collapsedAlleleTokens) + ">", false));
        }
    }

    private String[] collapseAltAllelesCommon(final List<Allele> alleles) {
        final List<String[]> alleleTokens = alleles.stream()
                .map(GATKSVVariantContextUtils::getSymbolicAlleleSymbols)
                .collect(Collectors.toList());
        final String[] firstTokens = alleleTokens.get(0);
        int alleleSize = 0;
        for (int i = 0; i < firstTokens.length; i++) {
            final String subtype = firstTokens[i];
            for (int j = 1; j < alleleTokens.size(); j++) {
                final String[] tokens = alleleTokens.get(j);
                if (i < tokens.length && subtype.equals(tokens[i])) {
                    alleleSize = i + 1;
                } else {
                    return Arrays.copyOf(firstTokens, alleleSize);
                }
            }
        }
        return Arrays.copyOf(firstTokens, alleleSize);
    }

    private String[] collapseAltAllelesMostSpecific(final List<Allele> alleles) {
        final List<String[]> alleleTokens = alleles.stream()
                .map(GATKSVVariantContextUtils::getSymbolicAlleleSymbols)
                .collect(Collectors.toList());
        final int maxSize = alleleTokens.stream().mapToInt(arr -> arr.length).max().getAsInt();
        int alleleIndex = 0;
        int alleleSize = 0;
        outerloop:
        for (int i = 0; i < maxSize; i++) {
            String subtype = null;
            for (int j = 0; j < alleleTokens.size(); j++) {
                final String[] tokens = alleleTokens.get(j);
                if (i < tokens.length) {
                    if (subtype == null) {
                        subtype = tokens[i];
                        alleleIndex = j;
                        alleleSize = i + 1;
                    } else if (!subtype.equals(tokens[i])) {
                        alleleIndex = j;
                        alleleSize = i;
                        break outerloop;
                    }
                }
            }
        }
        return Arrays.copyOf(alleleTokens.get(alleleIndex), alleleSize);
    }

    private List<Genotype> collapseAllGenotypes(final Collection<SVCallRecord> items,
                                                final Allele refAllele,
                                                final List<Allele> altAlleles) {
        return items.stream()
                .map(SVCallRecord::getGenotypes)
                .flatMap(GenotypesContext::stream)
                .collect(Collectors.groupingBy(Genotype::getSampleName))
                .values()
                .stream()
                .map(g -> collapseSampleGenotypes(g, refAllele, altAlleles))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    protected List<Allele> collapseSampleGenotypeAlleles(final Collection<Genotype> genotypes,
                                                         final int expectedCopyNumber,
                                                         final Integer copyNumber,
                                                         final List<Allele> altAlleles,
                                                         final Allele refAllele) {
        final List<Allele> alleles;
        if (altAlleles.contains(Allele.SV_SIMPLE_DEL) || altAlleles.contains(Allele.SV_SIMPLE_DUP)) {
            // For CNVs, collapse genotype using copy number
            alleles = getCNVGenotypeAllelesFromCopyNumber(altAlleles, refAllele, expectedCopyNumber, copyNumber);
        } else {
            alleles = collapseNonCNVGenotypeAlleles(genotypes, expectedCopyNumber, refAllele);
        }
        return alleles;
    }

    /**
     * Collapses a set of genotypes from a single sample into a list of alleles for a representative genotype. The
     * result is the most common non-ref genotype alleles list, potentially augmented with additional ref alleles
     * to match the ploidy.
     * @param genotypes genotypes from a single variant and single sample
     * @param refAllele ref allele for the genotype
     * @return list of alleles for a representative genotype (may be empty)
     */
    private List<Allele> collapseNonCNVGenotypeAlleles(final Collection<Genotype> genotypes,
                                                       final int expectedCopyNumber,
                                                       final Allele refAllele) {
        Utils.nonNull(genotypes);

        // Trivial case
        if (genotypes.size() == 1) {
            return genotypes.iterator().next().getAlleles();
        }

        // Empty case
        if (genotypes.isEmpty()) {
            return Collections.emptyList();
        }

        // Ploidy 0
        // TODO : we may not want to always throw away alleles here. For example, in the case of segmental
        //  duplications on chromosome Y.
        if (expectedCopyNumber == 0) {
            return Collections.emptyList();
        }

        // Determine frequency of each genotype
        final Map<List<Allele>, Integer> genotypeCounts = genotypes.stream().map(Genotype::getAlleles)
                .collect(Collectors.groupingBy(l -> l, Collectors.collectingAndThen(Collectors.toList(), List::size)));
        // Sort keys somehow, for stability
        final List<List<Allele>> sortedAlleleLists = genotypeCounts.keySet().stream()
                .sorted(ALLELE_COMPARATOR).collect(Collectors.toList());
        // Find most common non-ref genotype
        // Ties go to the first result in the list, as determined by the comparator
        // TODO : implement an allele comparator that can take into account GQ if available
        List<Allele> bestGenotypeAltAlleles = null;
        int bestGenotypeFrequency = 0;
        for (final List<Allele> alleles : sortedAlleleLists) {
            final List<Allele> genotypeAltAlleles = alleles.stream().filter(SVCallRecordUtils::isAltAllele).collect(Collectors.toList());
            final int genotypeFrequency = genotypeCounts.get(alleles);
            if (!genotypeAltAlleles.isEmpty() && (bestGenotypeAltAlleles == null || genotypeFrequency > bestGenotypeFrequency)) {
                bestGenotypeAltAlleles = genotypeAltAlleles;
                bestGenotypeFrequency = genotypeFrequency;
            }
        }
        // No non-ref genotypes
        if (bestGenotypeAltAlleles == null) {
            return sortedAlleleLists.get(0);
        }

        // Create alleles list with the selected alt alleles, and fill in remaining with ref
        final List<Allele> alleles = new ArrayList<>(expectedCopyNumber);
        final int numCollapsedRefAlleles = expectedCopyNumber - bestGenotypeAltAlleles.size();
        Utils.validate(numCollapsedRefAlleles >= 0, "Expected copy number is less than number of alt alleles");
        for (int i = 0; i < numCollapsedRefAlleles; i++) {
            alleles.add(refAllele);
        }
        alleles.addAll(bestGenotypeAltAlleles);
        return alleles;
    }

    private Integer collapseSampleCopyNumber(final Collection<Genotype> genotypes,
                                             final int expectedCopyNumber) {
        final Map<Integer, Integer> nonRefCopyNumberCounts = genotypes.stream()
                .map(g -> VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, expectedCopyNumber))
                .filter(c -> c != expectedCopyNumber)
                .collect(Collectors.groupingBy(i -> i, Collectors.collectingAndThen(Collectors.toList(), Collection::size)));
        if (nonRefCopyNumberCounts.isEmpty()) {
            return expectedCopyNumber;
        } else if (nonRefCopyNumberCounts.size() == 1) {
            return nonRefCopyNumberCounts.keySet().iterator().next();
        }

        // If we have conflicting sv types, let them cancel out to no-call
        // TODO expose option to control this behavior
        final int numLoss = (int) nonRefCopyNumberCounts.keySet().stream().filter(c -> c < expectedCopyNumber).count();
        if (numLoss != 0 && numLoss != nonRefCopyNumberCounts.size()) {
            return null;
        }

        // Sort such that ties go to state closest to ref
        final List<Integer> sortedCopyStates = nonRefCopyNumberCounts.keySet().stream()
                .sorted(getCopyNumberCollapsingComparator(expectedCopyNumber)).collect(Collectors.toList());
        int maxFreqCopyState = expectedCopyNumber;
        int maxFreq = 0;
        for (final Integer c : sortedCopyStates) {
            final int freq = nonRefCopyNumberCounts.get(c);
            if (freq > maxFreq) {
                maxFreqCopyState = c;
                maxFreq = freq;
            }
        }
        return maxFreqCopyState;
    }


    private static Comparator<Integer> getCopyNumberCollapsingComparator(final int expectedCopyNumber) {
        return (o1, o2) -> compareCopyNumberForCollapsing(o1, o2, expectedCopyNumber);
    }

    private static int compareCopyNumberForCollapsing(final Integer first, final Integer second, final Integer expectedCopyNumber) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        Utils.nonNull(expectedCopyNumber);
        final int distanceResult = Integer.compare(Math.abs(first - expectedCopyNumber), Math.abs(second - expectedCopyNumber));
        if (distanceResult != 0) {
            return distanceResult;
        }
        // Shouldn't happen if copy states are either all > or all < expected state
        return Integer.compare(first, second);
    }

    /***
     * Collapses collection of genotypes belonging to a single sample.
     */
    protected Genotype collapseSampleGenotypes(final Collection<Genotype> genotypes,
                                               final Allele refAllele,
                                               final List<Allele> altAlleles) {
        final GenotypeBuilder builder = new GenotypeBuilder(genotypes.iterator().next());

        // Reset attributes and collapse extended attributes
        builder.noAttributes();
        final int expectedCopyNumber = collapseExpectedCopyNumber(genotypes);
        builder.attributes(collapseGenotypeAttributes(genotypes, expectedCopyNumber));

        final Integer copyNumber;
        if (altAlleles.contains(Allele.SV_SIMPLE_DEL) || altAlleles.contains(Allele.SV_SIMPLE_DUP)) {
            // For CNVs, collapse genotype using copy number
            copyNumber = collapseSampleCopyNumber(genotypes, expectedCopyNumber);
            builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
        } else {
            copyNumber = null;
        }
        final List<Allele> collapsedAlleles = collapseSampleGenotypeAlleles(genotypes, expectedCopyNumber, copyNumber, altAlleles, refAllele);
        // Replace ref alleles with the possibly-new one
        final List<Allele> collapsedAllelesNewRef = collapsedAlleles.stream()
                .map(a -> (a != null && a.isReference()) ? refAllele : a)
                .collect(Collectors.toList());
        builder.alleles(collapsedAllelesNewRef);
        return builder.make();
    }

    /**
     * Generates genotype alleles, i.e. for the GT field, for CNVs (DEL and/or DUP). Genotypes that cannot be determined
     * unambiguously (e.g. at diploid DUP sites) result in {@link Allele#NO_CALL} alleles.
     * @param siteAltAlleles  unique alt alleles for the sample at the site
     * @param refAllele  reference allele
     * @param expectedCopyNumber  expected copy number (i.e. ploidy)
     * @param copyNumber  copy number state
     * @return  alleles for the sample at the site
     * @throws {@link IllegalArgumentException} if the alt allele(s) are not CNV(s)
     */
    public static List<Allele> getCNVGenotypeAllelesFromCopyNumber(final List<Allele> siteAltAlleles,
                                                                   final Allele refAllele,
                                                                   final int expectedCopyNumber,
                                                                   final Integer copyNumber) {
        Utils.nonNull(siteAltAlleles);
        Utils.nonNull(refAllele);
        Utils.validateArg(siteAltAlleles.size() <= 2, "No support for variants with over 2 alleles");
        if (copyNumber == null) {
            // Had conflicting genotypes
            return Collections.nCopies(expectedCopyNumber, Allele.NO_CALL);
        }
        Utils.validateArg(copyNumber >= 0, "Invalid negative copy number: " + copyNumber);
        Utils.validateArg(expectedCopyNumber >= 0, "Invalid expected copy number: " + expectedCopyNumber);
        if (siteAltAlleles.isEmpty()) {
            return Collections.nCopies(expectedCopyNumber, refAllele);
        } else if (siteAltAlleles.size() == 2) {
            if (siteAltAlleles.contains(Allele.SV_SIMPLE_DUP) && siteAltAlleles.contains(Allele.SV_SIMPLE_DEL)) {
                return getMultiallelicCNVAllelesFromCopyNumber(refAllele, expectedCopyNumber, copyNumber);
            } else {
                final String messageAlleles = String.join(", ",
                        siteAltAlleles.stream().map(Allele::getDisplayString).collect(Collectors.toList()));
                throw new IllegalArgumentException("Unsupported CNV alt alleles: " + messageAlleles);
            }
        }
        final Allele altAllele = siteAltAlleles.get(0);
        if (altAllele.equals(Allele.SV_SIMPLE_DEL)) {
            return getDeletionAllelesFromCopyNumber(refAllele, expectedCopyNumber, copyNumber);
        } else if (altAllele.equals(Allele.SV_SIMPLE_DUP)) {
            return getDuplicationAllelesFromCopyNumber(refAllele, expectedCopyNumber, copyNumber);
        } else {
            throw new IllegalArgumentException("Unsupported CNV alt allele: " + altAllele.getDisplayString());
        }
    }

    /**
     * Generates genotype alleles for deletion genotypes from the given copy number.
     * @param refAllele  reference allele for the site
     * @param expectedCopyNumber  expected copy number for the genotype
     * @param copyNumber  copy number for the genotype
     * @return  genotype alleles
     */
    public static List<Allele> getDeletionAllelesFromCopyNumber(final Allele refAllele, final int expectedCopyNumber,
                                                                final int copyNumber) {
        Utils.nonNull(refAllele);
        Utils.validateArg(copyNumber <= expectedCopyNumber, "Invalid deletion copy number " +
                copyNumber + " when ploidy is " + expectedCopyNumber);
        Utils.validateArg(expectedCopyNumber >= 0, "Invalid expected copy number: " + expectedCopyNumber);
        if (expectedCopyNumber == 0) {
            return Collections.emptyList();
        } else if (expectedCopyNumber == copyNumber) {
            // Most common in practice - use faster method
            return Collections.nCopies(expectedCopyNumber, refAllele);
        } else {
            final List<Allele> alleles = new ArrayList<>(expectedCopyNumber);
            for (int i = 0; i < copyNumber; i++) {
                alleles.add(refAllele);
            }
            for (int i = copyNumber; i < expectedCopyNumber; i++) {
                alleles.add(Allele.SV_SIMPLE_DEL);
            }
            return alleles;
        }
    }

    /**
     * Generates genotype alleles for duplication genotypes from the given copy number. Genotypes that cannot be
     * determined unambiguously (e.g. diploid sites) result in {@link Allele#NO_CALL} alleles.
     * @param refAllele  reference allele for the site
     * @param expectedCopyNumber  expected copy number for the genotype
     * @param copyNumber  copy number for the genotype
     * @return  genotype alleles
     */
    public static List<Allele> getDuplicationAllelesFromCopyNumber(final Allele refAllele, final int expectedCopyNumber,
                                                                   final int copyNumber) {
        Utils.nonNull(refAllele);
        Utils.validateArg(copyNumber >= expectedCopyNumber, "Invalid duplication copy number " +
                copyNumber + " when ploidy is " + expectedCopyNumber);
        Utils.validateArg(expectedCopyNumber >= 0, "Invalid expected copy number: " + expectedCopyNumber);
        if (expectedCopyNumber == copyNumber) {
            // Most common in practice - use faster method
            return Collections.nCopies(expectedCopyNumber, refAllele);
        } else if (expectedCopyNumber == 0) {
            return Collections.emptyList();
        } else if (expectedCopyNumber == 1) {
            // At this point know it's a DUP, and since it's haploid the phasing is unambiguous
            // TODO : May not be a simple DUP, need a more general Allele for possible multi-copy DUP
            return Collections.singletonList(Allele.SV_SIMPLE_DUP);
        } else if (copyNumber == expectedCopyNumber + 1) {
            // Case where we can resolve alleles
            return makeBiallelicList(Allele.SV_SIMPLE_DUP, refAllele, 1, expectedCopyNumber);
        } else {
            // DUP alleles cannot be resolved in other cases
            return Collections.nCopies(expectedCopyNumber, Allele.NO_CALL);
        }
    }

    /**
     * Generates genotype alleles for multi-allelic CNV genotypes from the given copy number. Genotypes that cannot be
     * determined unambiguously (e.g. diploid sites) result in {@link Allele#NO_CALL} alleles.
     * @param refAllele  reference allele for the site
     * @param expectedCopyNumber  expected copy number for the genotype
     * @param copyNumber  copy number for the genotype
     * @return  genotype alleles
     */
    public static List<Allele> getMultiallelicCNVAllelesFromCopyNumber(final Allele refAllele,
                                                                       final int expectedCopyNumber,
                                                                       final int copyNumber) {
        Utils.nonNull(refAllele);
        Utils.validateArg(copyNumber >= 0, "Invalid multi-allelic CNV copy number: " + copyNumber);
        Utils.validateArg(expectedCopyNumber >= 0, "Invalid expected copy number: " + expectedCopyNumber);
        if (expectedCopyNumber == 0) {
            return Collections.emptyList();
        } else if (expectedCopyNumber == 1) {
            if (copyNumber == 0) {
                return Collections.singletonList(Allele.SV_SIMPLE_DEL);
            } else if (copyNumber == 1) {
                return Collections.singletonList(refAllele);
            } else {
                // copyNumber >= 2 - haploid dups have non-ambiguous phasing
                return Collections.singletonList(Allele.SV_SIMPLE_DUP);
            }
        } else if (copyNumber == 0) {
            return Collections.nCopies(expectedCopyNumber, Allele.SV_SIMPLE_DEL);
        } else if (copyNumber == 1) {
            return makeBiallelicList(Allele.SV_SIMPLE_DEL, refAllele, expectedCopyNumber - 1, expectedCopyNumber);
        } else {
            // Alleles cannot be resolved in other cases
            return Collections.nCopies(expectedCopyNumber, Allele.NO_CALL);
        }
    }

    /**
     * Creates list for biallelic sites with the specified number of each allele.
     * @param alt  alt allele
     * @param ref  ref allele
     * @param numAlt  number of alt alleles
     * @param ploidy  number of alleles (total)
     * @return  resulting list
     */
    public static List<Allele> makeBiallelicList(final Allele alt, final Allele ref, final int numAlt, final int ploidy) {
        Utils.nonNull(alt);
        Utils.nonNull(ref);
        Utils.validateArg(numAlt >= 0, "Negative number of alt alleles");
        Utils.validateArg(ploidy >= 0, "Negative number of ref alleles");
        if (ploidy == 0) {
            return Collections.emptyList();
        } else if (numAlt == 0) {
            return Collections.nCopies(ploidy, ref);
        } else if (numAlt == ploidy) {
            return Collections.nCopies(ploidy, alt);
        }
        final int numRef = ploidy - numAlt;
        final List<Allele> alleles = new ArrayList<>(numAlt + numRef);
        for (int i = 0; i < numRef; i++) {
            alleles.add(ref);
        }
        for (int i = 0; i < numAlt; i++) {
            alleles.add(alt);
        }
        return alleles;
    }

    protected Map<String, Object> collapseGenotypeAttributes(final Collection<Genotype> genotypes,
                                                             final int expectedCopyNumber) {
        Utils.nonNull(genotypes);
        Utils.nonEmpty(genotypes);
        final Map<String, Object> collapsedAttributes = new HashMap<>();
        final Map<String, Set<Object>> genotypeFields = genotypes.stream().map(Genotype::getExtendedAttributes)
                .map(Map::entrySet)
                .flatMap(Set::stream)
                .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toSet())));
        for (final Map.Entry<String, Set<Object>> entry : genotypeFields.entrySet()) {
            if (entry.getKey().equals(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                // Copy number attribute set later using special logic in collapseSampleCopyNumber
                // Unit testing is easier with this design
                continue;
            } else if (entry.getKey().equals(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT)) {
                // Set after loop
                continue;
            } else {
                collapsedAttributes.put(entry.getKey(), collapseSampleGenotypeAttribute(entry.getKey(), entry.getValue()));
            }
        }
        collapsedAttributes.put(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
        return collapsedAttributes;
    }

    protected Object collapseSampleGenotypeAttribute(final String key, final Set<Object> values) {
        if (values.size() == 1) {
            return values.iterator().next();
        } else {
            return null;
        }
    }

    protected Map<String, Object> collapseVariantAttributes(final Collection<SVCallRecord> items) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        final Map<String, Object> collapsedAttributes = new HashMap<>();
        final Map<String, Set<Object>> attributes = items.stream().map(SVCallRecord::getAttributes)
                .map(Map::entrySet)
                .flatMap(Set::stream)
                .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toSet())));
        for (final Map.Entry<String, Set<Object>> entry : attributes.entrySet()) {
            collapsedAttributes.put(entry.getKey(), collapseSingleVariantAttribute(entry.getValue()));
        }
        collapsedAttributes.put(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, items.stream().map(SVCallRecord::getId).sorted().collect(Collectors.toList()));
        return collapsedAttributes;
    }

    protected Object collapseSingleVariantAttribute(final Set<Object> values) {
        if (values.size() == 1) {
            return values.iterator().next();
        } else {
            return null;
        }
    }

    /***
     * Calculates new SVLEN value.
     * @param items  records to collapse
     * @param newType   collapsed sv type
     * @return
     */
    protected final Integer collapseLength(final Collection<SVCallRecord> items, final StructuralVariantType newType) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        if (newType.equals(StructuralVariantType.INS)) {
            return collapseInsertionLength(items);
        } else {
            return null;
        }
    }

    protected Integer collapseInsertionLength(final Collection<SVCallRecord> items) {
        // No need for strategy when the answer is obvious
        if (items.size() == 1) {
            return items.iterator().next().getLength();
        }

        // We only want to use strategies when there are multiple defined lengths
        final int[] definedLengths = items.stream()
                .map(SVCallRecord::getLength)
                .filter(len -> len != null)
                .mapToInt(Integer::intValue)
                .toArray();

        // All lengths were undefined
        if (definedLengths.length == 0) {
            return null;
        }

        // Only one length is defined - no ambiguity here
        if (definedLengths.length == 1) {
            return definedLengths[0];
        }

        // Employ strategies when there are two or more defined lengths
        if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MEDIAN) {
            return MathUtils.median(definedLengths, Percentile.EstimationType.R_3);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MEAN) {
            return (int) Math.ceil(MathUtils.sum(definedLengths) / (double) definedLengths.length);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MIN) {
            return MathUtils.arrayMin(definedLengths);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MAX) {
            return MathUtils.arrayMax(definedLengths);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.UNDEFINED) {
            return null;
        } else {
            throw new UnsupportedOperationException("Unimplemented insertion summary strategy: " + insertionLengthSummaryStrategy.name());
        }
    }

    protected String collapseIds(final Collection<SVCallRecord> records) {
        Utils.nonNull(records);
        Utils.nonEmpty(records);
        return records.stream().map(SVCallRecord::getId).sorted().collect(Collectors.toList()).get(0);
    }

    /**
     * Gets records likely to have the most accurate breakpoint position. These usually are supported by PE/SR/AS support,
     * whereas depth-only calls tend to be approximate.
     */
    protected Collection<SVCallRecord> getRecordsWithMostPreciseBreakpoints(final Collection<SVCallRecord> items) {
        if (items.stream().allMatch(call -> call.isDepthOnly())) {
            return items;
        } else {
            return items.stream().filter(call -> !call.isDepthOnly()).collect(Collectors.toList());
        }
    }

    /**
     * @param items
     * @return (key, value) entry of (start, end)
     */
    protected Pair<Integer,Integer> collapseInterval(final Collection<SVCallRecord> items) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        final SVCallRecord exampleCall = items.iterator().next();
        if (items.size() > 1) {
            final List<String> contigsA = items.stream().map(SVCallRecord::getContigA).distinct().collect(Collectors.toList());
            Utils.validate(contigsA.size() == 1, "Cannot collapse intervals with multiple position A contigs");
            final List<String> contigsB = items.stream().map(SVCallRecord::getContigB).distinct().collect(Collectors.toList());
            Utils.validate(contigsB.size() == 1, "Cannot collapse intervals with multiple position B contigs");
        }
        final int[] startPositions = items.stream().mapToInt(SVCallRecord::getPositionA).sorted().toArray();
        final int[] endPositions = items.stream().mapToInt(SVCallRecord::getPositionB).sorted().toArray();
        //use the mid value of the sorted list so the start and end represent real breakpoint observations
        final int medianStart = MathUtils.median(startPositions, Percentile.EstimationType.R_1);
        final int medianEnd = MathUtils.median(endPositions, Percentile.EstimationType.R_1);
        final int newStart;
        final int newEnd;
        switch (breakpointSummaryStrategy) {
            case MEDIAN_START_MEDIAN_END:
                newStart = medianStart;
                newEnd = medianEnd;
                break;
            case MIN_START_MAX_END:
                newStart = startPositions[0];
                newEnd = endPositions[endPositions.length - 1];
                break;
            case MAX_START_MIN_END:
                newStart = startPositions[startPositions.length - 1];
                newEnd = endPositions[0];
                break;
            case MEAN_START_MEAN_END:
                newStart = (int)Math.round(MathUtils.sum(startPositions) / (double) startPositions.length);
                newEnd = (int)Math.round(MathUtils.sum(endPositions) / (double) startPositions.length);
                break;
            default:
                throw new UnsupportedOperationException("Unknown breakpoint summary strategy: " + breakpointSummaryStrategy.name());
        }
        if (exampleCall.getType().equals(StructuralVariantType.INS)) {
            // Insertions are a single locus
            return Pair.of(newStart, newStart);
        } else if (exampleCall.getContigA().equals(exampleCall.getContigB())) {
            // Do not let end precede start
            return Pair.of(newStart, Math.max(newStart, newEnd));
        } else {
            // Different contigs, so no constraint on position order
            return Pair.of(newStart, newEnd);
        }
    }

    protected StructuralVariantType collapseTypes(final Collection<SVCallRecord> records) {
        final Set<StructuralVariantType> types = records.stream().map(SVCallRecord::getType).collect(Collectors.toSet());
        if (types.size() == 1) {
            return types.iterator().next();
        }
        if (types.stream().allMatch(GATKSVVariantContextUtils::isCnvType)) {
            return StructuralVariantType.CNV;
        }
        final List<String> typeStrings = types.stream().map(StructuralVariantType::name).collect(Collectors.toList());
        throw new IllegalArgumentException("Incompatible SV types found in cluster: " + String.join(", ", typeStrings));
    }

    protected List<String> collapseAlgorithms(final Collection<SVCallRecord> records) {
        return records.stream()
                .map(SVCallRecord::getAlgorithms)
                .flatMap(Collection::stream)
                .distinct()
                .sorted()
                .collect(Collectors.toList());
    }

    public static final class AlleleCollectionCollapserComparator implements Comparator<Collection<Allele>> {
        @Override
        public int compare(final Collection<Allele> o1, final Collection<Allele> o2) {
            if (o1 == o2) {
                return 0;
            }
            Utils.nonNull(o1);
            Utils.nonNull(o2);
            // Sort on number of alts
            final long numAlt1 = (int) o1.stream().filter(SVCallRecordUtils::isAltAllele).count();
            final long numAlt2 = (int) o2.stream().filter(SVCallRecordUtils::isAltAllele).count();
            if (numAlt1 < numAlt2) {
                return -1;
            } else if (numAlt1 > numAlt2) {
                return 1;
            }
            // Sort on number of calls
            final long numCalled1 = (int) o1.stream().filter(a -> a != null && !a.isNoCall()).count();
            final long numCalled2 = (int) o2.stream().filter(a -> a != null && !a.isNoCall()).count();
            if (numCalled1 > numCalled2) {
                return -1;
            } else if (numCalled1 < numCalled2) {
                return 1;
            }
            final List<Allele> sorted1 = SVCallRecordUtils.sortAlleles(o1);
            final List<Allele> sorted2 = SVCallRecordUtils.sortAlleles(o2);
            // Resort to display strings
            for (int i = 0; i < sorted1.size() && i < sorted2.size(); i++) {
                final int comp = sorted1.get(i).getDisplayString().compareTo(sorted2.get(i).getDisplayString());
                if (comp != 0) {
                    return comp;
                }
            }
            if (o1.size() < o2.size()) {
                return -1;
            } else if (o1.size() > o2.size()) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}
