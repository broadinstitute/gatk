package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
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
public class CanonicalSVCollapser {

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
        MEAN_START_MEAN_END,

        /**
         * Picks a single representative call closest to the overall cluster
         */
        REPRESENTATIVE,

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

    private static final Set<GATKSVVCFConstants.StructuralVariantAnnotationType> SUPPORTED_SV_TYPES = Sets.newHashSet(
            GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
            GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
            GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
            GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
            GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CTX
    );

    /**
     * Comparators used for picking the representative genotype for a given sample
     */
    final Comparator<Genotype> genotypeIsNonRefComparator = (o1, o2) -> {
        final long count1 = Math.min(1, o1.getAlleles().stream().filter(Allele::isNonReference).filter(Allele::isCalled).count());
        final long count2 = Math.min(1, o2.getAlleles().stream().filter(Allele::isNonReference).filter(Allele::isCalled).count());
        return Long.compare(count1, count2);
    };

    final Comparator<Genotype> genotypeNonRefCountComparator = (o1, o2) -> {
        final long count1 = o1.getAlleles().stream().filter(Allele::isNonReference).filter(Allele::isCalled).count();
        final long count2 = o2.getAlleles().stream().filter(Allele::isNonReference).filter(Allele::isCalled).count();
        return Long.compare(count1, count2);
    };

    final Comparator<Genotype> genotypeCalledComparator = (o1, o2) -> {
        final long count1 = o1.getAlleles().stream().filter(Allele::isCalled).count();
        final long count2 = o2.getAlleles().stream().filter(Allele::isCalled).count();
        return Long.compare(count1, count2);
    };

    final Comparator<Genotype> genotypeQualityComparator = (o1, o2) -> {
        final int quality1 = VariantContextGetters.getAttributeAsInt(o1, VCFConstants.GENOTYPE_QUALITY_KEY, 0);
        final int quality2 = VariantContextGetters.getAttributeAsInt(o2, VCFConstants.GENOTYPE_QUALITY_KEY, 0);
        return Integer.compare(quality1, quality2);
    };

    final Comparator<Genotype> genotypeCopyNumberQualityComparator = new Comparator<Genotype>() {
        @Override
        public int compare(Genotype o1, Genotype o2) {
            final int quality1 = VariantContextGetters.getAttributeAsInt(o1, GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT, 0);
            final int quality2 = VariantContextGetters.getAttributeAsInt(o2, GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT, 0);
            return Integer.compare(quality1, quality2);
        }
    };

    final Comparator<Genotype> genotypeCopyNumberComparator = new Comparator<Genotype>() {
        @Override
        public int compare(Genotype o1, Genotype o2) {
            final int expectedQualityNumber1 = VariantContextGetters.getAttributeAsInt(o1, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0);
            final int copyNumber1 = VariantContextGetters.getAttributeAsInt(o1, GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);
            final int expectedQualityNumber2 = VariantContextGetters.getAttributeAsInt(o2, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0);
            final int copyNumber2 = VariantContextGetters.getAttributeAsInt(o2, GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);
            return Double.compare(Math.abs(expectedQualityNumber1 - copyNumber1), Math.abs(expectedQualityNumber2 - copyNumber2));
        }
    };

    private final AltAlleleSummaryStrategy altAlleleSummaryStrategy;
    private final BreakpointSummaryStrategy breakpointSummaryStrategy;
    private final ReferenceSequenceFile reference;
    private final SAMSequenceDictionary dictionary;

    public CanonicalSVCollapser(final ReferenceSequenceFile reference,
                                final AltAlleleSummaryStrategy altAlleleSummaryStrategy,
                                final BreakpointSummaryStrategy breakpointSummaryStrategy) {
        this.reference = Utils.nonNull(reference);
        this.dictionary = reference.getSequenceDictionary();
        this.altAlleleSummaryStrategy = altAlleleSummaryStrategy;
        this.breakpointSummaryStrategy = breakpointSummaryStrategy;
    }

    private static final int distance(final SVCallRecord item, final int newStart, final int newEnd) {
        return Math.abs(item.getPositionA() - newStart) + Math.abs(item.getPositionB() - newEnd);
    }

    protected SVCallRecord getRepresentativeRecord(final Collection<SVCallRecord> items, final int newStart, final int newEnd) {
        return items.stream().sorted(Comparator.comparing(SVCallRecord::getId))
                .min(Comparator.comparing(r -> distance(r, newStart, newEnd))).get();
    }

    public SVCallRecord collapse(final SVClusterEngine.OutputCluster cluster) {
        final List<SVCallRecord> items = cluster.getItems();
        validateRecords(items);

        // Prefer using variants generated with PESR callers, which tend to generate more precise breakpoints
        final Collection<SVCallRecord> mostPreciseCalls = getRecordsWithMostPreciseBreakpoints(items);
        final Pair<Integer, Integer> coordinates = collapseInterval(mostPreciseCalls);
        final int start = coordinates.getKey();
        final int end = coordinates.getValue();
        final SVCallRecord representative = getRepresentativeRecord(mostPreciseCalls, start, end);
        final GATKSVVCFConstants.StructuralVariantAnnotationType type = collapseTypes(items);
        final Integer length = collapseLength(representative, type, start, end);
        final List<String> algorithms = collapseAlgorithms(items);
        final Map<String, Object> attributes = collapseAttributes(representative, items);

        final Boolean strandA = type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV ? null : representative.getStrandA();
        final Boolean strandB = type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV ? null : representative.getStrandB();

        final Allele refAllele = collapseRefAlleles(representative.getContigA(), start);
        final List<Allele> altAlleles;
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.BND) {
            altAlleles = Collections.singletonList(constructBndAllele(strandA, strandB, representative.getContigB(), end, refAllele));
        } else {
            altAlleles = collapseAltAlleles(items);
        }
        final List<Allele> alleles = collapseAlleles(altAlleles, refAllele);
        final List<Genotype> genotypes = harmonizeAltAlleles(altAlleles, collapseAllGenotypes(items, refAllele));

        final Set<String> filters = collapseFilters(items);
        final Double quality = collapseQuality(items);

        return new SVCallRecord(representative.getId(), representative.getContigA(), start, strandA, representative.getContigB(),
                end, strandB, type, representative.getComplexSubtype(), representative.getComplexEventIntervals(),
                length, algorithms, alleles, genotypes, attributes, filters, quality, dictionary);
    }

    protected List<Allele> collapseAlleles(final List<Allele> altAlleles, final Allele refAllele) {
        final List<Allele> alleles = new ArrayList<>(altAlleles.size() + 1);
        alleles.add(refAllele);
        alleles.addAll(altAlleles);
        return alleles;
    }

    protected Set<String> collapseFilters(final List<SVCallRecord> items) {
        return items.stream()
                .map(SVCallRecord::getFilters)
                .flatMap(Collection::stream)
                .collect(Collectors.toSet());
    }

    protected Double collapseQuality(final List<SVCallRecord> items) {
        if (items.size() == 1) {
            return items.get(0).getLog10PError();
        } else {
            return null;
        }
    }

    /**
     * Asserts that the given records are valid for collapsing.
     */
    protected void validateRecords(final Collection<SVCallRecord> records) {
        for (final SVCallRecord r : records) {
            if (!SUPPORTED_SV_TYPES.contains(r.getType())) {
                throw new IllegalArgumentException("Unsupported SV type: " + r.getType());
            }
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
     * Returns the ref allele at the given locus.
     */
    protected Allele collapseRefAlleles(final String contig, final int pos) {
        final byte[] bases = ReferenceUtils.getRefBaseAtPosition(reference, contig, pos);
        Utils.validate(bases != null && bases.length == 1, "Invalid reference locus " + contig + ":" + pos);
        return Allele.create(bases[0], true);
    }

    protected static Allele constructBndAllele(final Boolean strandA, final Boolean strandB, final String contigB,
                                        final int posB, final Allele refAllele) {
        Utils.validateArg(strandA != null, "First breakend strand cannot be null");
        Utils.validateArg(strandB != null, "Second breakend strand cannot be null");
        final String bracket = strandB ? "]" : "[";
        final String str;
        if (strandA) {
            str = refAllele.getBaseString() + bracket + contigB + ":" + posB + bracket;
        } else {
            str = bracket + contigB + ":" + posB + bracket + refAllele.getBaseString();
        }
        return Allele.create(str, false);
    }

    /**
     * Collapses symbolic alleles into a list of representative alleles. Note this supports sub-typed alleles such as
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
            for (final Allele a : altAlleles) {
                if (!validateAltAllele(a)) {
                    throw new IllegalArgumentException("Cannot collapse non-symbolic allele: " + a.getDisplayString());
                }
            }
            // Multiple non-ref alleles need collapsing
            // TODO does not search for subtypes e.g. <DUP:TANDEM>
            int numCnvAlleles = 0;
            int numMultiallelicAlleles = 0;
            // CNVs
            for (final Allele a : altAlleles) {
                if (a.equals(Allele.SV_SIMPLE_CNV)) {
                    numCnvAlleles++;
                    numMultiallelicAlleles++;
                } else if (a.equals(Allele.SV_SIMPLE_DUP) || a.equals(Allele.SV_SIMPLE_DEL)) {
                    numCnvAlleles++;
                }
            }
            if (altAlleles.size() == numCnvAlleles) {
                // Try to match the input alleles if using CNV or DEL/DUP
                if (numMultiallelicAlleles > 0) {
                    return Collections.singletonList(Allele.SV_SIMPLE_CNV);
                } else {
                    return Arrays.asList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP);
                }
            }

            final String[] collapsedAlleleTokens;
            if (altAlleleSummaryStrategy == AltAlleleSummaryStrategy.COMMON_SUBTYPE) {
                collapsedAlleleTokens = collapseAltAllelesCommon(altAlleles);
            } else if (altAlleleSummaryStrategy == AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE) {
                collapsedAlleleTokens = collapseAltAllelesMostSpecific(altAlleles);
            } else {
                throw new UnsupportedOperationException("Unimplemented alt allele summary strategy: " + altAlleleSummaryStrategy.name());
            }
            return Collections.singletonList(Allele.create("<" + String.join(":", collapsedAlleleTokens) + ">", false));
        }
    }

    protected boolean validateAltAllele(final Allele allele) {
        return !allele.isReference() && allele.isSymbolic() && !allele.isBreakpoint() && !allele.isSingleBreakend();
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
                                                final Allele refAllele) {
        return items.stream()
                .map(SVCallRecord::getGenotypes)
                .flatMap(GenotypesContext::stream)
                .collect(Collectors.groupingBy(Genotype::getSampleName))
                .values()
                .stream()
                .map(g -> collapseSampleGenotypes(g, refAllele))
                .collect(Collectors.toList());
    }

    /***
     * Collapses collection of genotypes belonging to a single sample.
     */
    protected Genotype collapseSampleGenotypes(final Collection<Genotype> genotypes,
                                               final Allele refAllele) {

        // Reset attributes and collapse extended attributes
        final Genotype representative = getRepresentativeGenotype(genotypes);
        final GenotypeBuilder builder = new GenotypeBuilder(representative);
        // Replace ref allele with the possibly-new one
        final List<Allele> collapsedAllelesNewRef = representative.getAlleles().stream()
                .map(a -> (a != null && a.isReference()) ? refAllele : a)
                .collect(Collectors.toList());
        builder.alleles(collapsedAllelesNewRef);
        return builder.make();
    }

    /**
     *
     * @param genotypes list of candidate genotypes
     * @return representative genotype
     */
    protected Genotype getRepresentativeGenotype(final Collection<Genotype> genotypes) {

        return genotypes.stream()
                .max(genotypeIsNonRefComparator
                        .thenComparing(genotypeCalledComparator)
                        .thenComparing(genotypeQualityComparator)
                        .thenComparing(genotypeNonRefCountComparator)
                        .thenComparing(genotypeCopyNumberQualityComparator)
                        .thenComparing(genotypeCopyNumberComparator)).get();
    }


    /**
     * Generates genotype alleles, i.e. for the GT field, for CNVs (DEL and/or DUP). Multi-allelics result in
     * {@link Allele#NO_CALL} alleles.
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
        Utils.validateArg(siteAltAlleles.size() <= 2, "No support for variants with over 2 alt alleles");
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
                return Collections.nCopies(expectedCopyNumber, Allele.NO_CALL);
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
        Utils.validateArg(expectedCopyNumber >= 0, "Negative expected copy number");
        if (expectedCopyNumber == 0) {
            return Collections.emptyList();
        } else if (expectedCopyNumber <= copyNumber) {
            // Most common in practice - use faster method
            return Collections.nCopies(expectedCopyNumber, refAllele);
        } else {
            final int numAlt = expectedCopyNumber - copyNumber;
            return makeBiallelicList(Allele.SV_SIMPLE_DEL, refAllele, numAlt, expectedCopyNumber);
        }
    }

    /**
     * Generates genotype alleles for duplication genotypes from the given copy number. Genotypes that cannot be
     * determined unambiguously (e.g. diploid sites) result in {@link Allele#NO_CALL} alleles. Assuming
     * no multi-copy alleles, we return no-call alleles if the apparent copy number exceeds twice the ploidy.
     * @param refAllele  reference allele for the site
     * @param expectedCopyNumber  expected copy number for the genotype
     * @param copyNumber  copy number for the genotype
     * @return  genotype alleles
     */
    public static List<Allele> getDuplicationAllelesFromCopyNumber(final Allele refAllele, final int expectedCopyNumber,
                                                                   final int copyNumber) {
        Utils.nonNull(refAllele);
        Utils.validateArg(expectedCopyNumber >= 0, "Negative expected copy number");
        if (expectedCopyNumber >= copyNumber) {
            // Most common in practice - use faster method
            return Collections.nCopies(expectedCopyNumber, refAllele);
        } else if (expectedCopyNumber == 0) {
            return Collections.emptyList();
        }
        final int numAlt = copyNumber - expectedCopyNumber;
        if (expectedCopyNumber == 1) {
            // Common case on chrY - use faster method
            return Collections.singletonList(Allele.SV_SIMPLE_DUP);
        }
        // Case where we assume no multi-copy alleles or can resolve alleles
        if (numAlt > expectedCopyNumber) {
            throw new IllegalArgumentException("Encountered simple DUP with copy number " + copyNumber + " but the " +
                    "ploidy is only " + expectedCopyNumber);
        }
        return makeBiallelicList(Allele.SV_SIMPLE_DUP, refAllele, Math.min(numAlt, expectedCopyNumber), expectedCopyNumber);
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

    protected Map<String, Object> collapseAttributes(final SVCallRecord representative,
                                                     final Collection<SVCallRecord> items) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        final Map<String, Object> attributes = new HashMap<>();
        for (final Map.Entry<String, Object> entry : representative.getAttributes().entrySet()) {
            attributes.put(entry.getKey(), entry.getValue());
        }
        attributes.put(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, items.stream().map(SVCallRecord::getId).sorted().collect(Collectors.toList()));
        return attributes;
    }

    /***
     * Calculates new SVLEN value.
     * @param representative  representative record
     * @param newType   collapsed sv type
     * @param start    start position of new record
     * @param end    end position of new record
     * @return
     */
    protected final Integer collapseLength(final SVCallRecord representative,
                                           final GATKSVVCFConstants.StructuralVariantAnnotationType newType,
                                           final int start, final int end) {
        Utils.nonNull(representative);
        if (newType == GATKSVVCFConstants.StructuralVariantAnnotationType.INS || newType == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
            return representative.getLength();
        } else if (newType.equals(GATKSVVCFConstants.StructuralVariantAnnotationType.BND) || newType.equals(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX)) {
            return null;
        } else {
            return CoordMath.getLength(start, end);
        }
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
        final int newStart;
        final int newEnd;
        switch (breakpointSummaryStrategy) {
            case MEDIAN_START_MEDIAN_END:
                //use the mid value of the sorted list so the start and end represent real breakpoint observations
                final int medianStart = MathUtils.median(startPositions, Percentile.EstimationType.R_1);
                final int medianEnd = MathUtils.median(endPositions, Percentile.EstimationType.R_1);
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
            case REPRESENTATIVE:
                final SVCallRecord representative = getRepresentativeIntervalItem(items, startPositions, endPositions);
                newStart = representative.getPositionA();
                newEnd = representative.getPositionB();
                break;
            default:
                throw new UnsupportedOperationException("Unknown breakpoint summary strategy: " + breakpointSummaryStrategy.name());
        }
        if (exampleCall.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            // Insertions are represented as a point
            return Pair.of(newStart, newStart);
        } else if (exampleCall.isIntrachromosomal()) {
            // Do not let end precede start
            return Pair.of(newStart, Math.max(newStart, newEnd));
        } else {
            // Different contigs, so no constraint on position order
            return Pair.of(newStart, newEnd);
        }
    }

    private SVCallRecord getRepresentativeIntervalItem(final Collection<SVCallRecord> records,
                                                       final int[] starts,
                                                       final int[] ends) {
        if (records.size() == 1) {
            return records.iterator().next();
        }
        // Favor more common variant with most similar distance
        final Comparator<SVCallRecord> carrierCountComparator = Comparator.comparing(r -> -r.getCarrierGenotypeList().size());
        final Comparator<SVCallRecord> distanceComparator = Comparator.comparing(r -> getDistance(r.getPositionA(), r.getPositionB(), starts, ends));
        final Comparator<SVCallRecord> idComparator = Comparator.comparing(r -> getDistance(r.getPositionA(), r.getPositionB(), starts, ends)); // stabilizes order
        return records.stream().min(
                carrierCountComparator
                        .thenComparing(distanceComparator)
                        .thenComparing(idComparator)).get();
    }

    protected static long getDistance(final int posA,
                                      final int posB,
                                      final int[] starts,
                                      final int[] ends) {
        long d = 0;
        for (int j = 0; j < starts.length; j++) {
            d += FastMath.abs(starts[j] - posA);
        }
        for (int j = 0; j < ends.length; j++) {
            d += FastMath.abs(ends[j] - posB);
        }
        return d;
    }

    protected GATKSVVCFConstants.StructuralVariantAnnotationType collapseTypes(final Collection<SVCallRecord> records) {
        final Set<GATKSVVCFConstants.StructuralVariantAnnotationType> types = records.stream().map(SVCallRecord::getType).collect(Collectors.toSet());
        if (types.size() == 1) {
            return types.iterator().next();
        }
        if (types.stream().allMatch(GATKSVVariantContextUtils::isCnvType)) {
            return GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
        }
        final List<String> typeStrings = types.stream().map(GATKSVVCFConstants.StructuralVariantAnnotationType::name).collect(Collectors.toList());
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
}
