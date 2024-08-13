package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A tool to allow for comparison of two single-sample GVCFs or multi-sample VCFs when some differences are expected
 * input variants should be tagged "actual" or "expected"
 *
 * Example usage:
 * java -jar gatk.jar VCFComparator -R ref.fasta -V:actual actual.vcf -V:expected expected.vcf
 *
 * Caveats:
 * Reference blocks are not checked
 * VCFComparator will ignore version-specific NaN versus empty AS_RAW annotation discrepancies
 * VCFs are both expected to NOT have split multi-alleleics
 */
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Compare two VCFs, as for pipeline test updates",
        oneLineSummary = "Compare two VCFs",
        programGroup = VariantEvaluationProgramGroup.class
)
public class VCFComparator extends MultiVariantWalkerGroupedByOverlap {
    public static final String ALLOW_NEW_STARS_LONG_NAME = "allow-new-stars";
    public static final String IGNORE_QUALS_LONG_NAME = "ignore-quals";
    public static final String DP_CHANGE_ALLOWED_LONG_NAME = "dp-change-allowed";
    public static final String RANK_SUM_CHANGE_ALLOWED_LONG_NAME = "ranksum-change-allowed";
    public static final String POSITIONS_ONLY_LONG_NAME = "positions-only";

    private VariantAnnotatorEngine annotatorEngine;
    private boolean isSingleSample = true;
    private boolean alleleNumberIsDifferent = false;
    private boolean inbreedingCoeffIsDifferent = false;
    private final List<String> annotationsThatVaryWithNoCalls = Arrays.asList(VCFConstants.ALLELE_NUMBER_KEY, GATKVCFConstants.INBREEDING_COEFFICIENT_KEY,
            GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY, GATKVCFConstants.EXCESS_HET_KEY,
            GATKVCFConstants.MLE_ALLELE_COUNT_KEY, GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, VCFConstants.ALLELE_FREQUENCY_KEY,
            VCFConstants.DEPTH_KEY, GATKVCFConstants.QUAL_BY_DEPTH_KEY, GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY,
            GATKVCFConstants.VQS_LOD_KEY, GATKVCFConstants.AS_VQS_LOD_KEY, GATKVCFConstants.AS_FILTER_STATUS_KEY,
            GATKVCFConstants.CULPRIT_KEY, GATKVCFConstants.AS_CULPRIT_KEY); //AC should be the same
    private boolean failOnCompletion = false;

    @Argument(fullName = "warn-on-errors",
            shortName = "warn-on-errors",
            doc = "just emit warnings on errors instead of terminating the run at the first instance. Won't fail unless --finish-before-failing is set",
            optional = true)
    private boolean warnOnError = false;

    @Argument(fullName = "finish-before-failing",
    doc = "emit warnings as necesarry, but perform the whole comparison before failing",
    optional= true)
    private boolean finishBeforeFailing = false;

    @Argument(fullName = "default-ploidy", optional = true, doc = "Overall expected ploidy. Default value is 2.")
    private int defaultPloidy = 2;

    @Argument(fullName = IGNORE_QUALS_LONG_NAME, optional = true, mutex={POSITIONS_ONLY_LONG_NAME})
    private boolean ignoreQuals = false;

    @Argument(fullName = "qual-change-allowed", optional = true)
    private double qualTolerance = 0.001;

    @Argument(fullName = "inbreeding-coeff-change-allowed", optional = true)
    private double inbreedingCoeffTolerance = 0.001;

    @Argument(fullName = "good-qual-threshold", optional = true, doc = "Variants with QUAL at or above this value are considered high quality and should appear in both files")
    private double goodQualThreshold = 100.0;

    @Argument(fullName = DP_CHANGE_ALLOWED_LONG_NAME, optional = true, doc = "Note that this is a signed change (actual - expected)", mutex={POSITIONS_ONLY_LONG_NAME})
    private int dpChange = 0;

    @Argument(fullName = RANK_SUM_CHANGE_ALLOWED_LONG_NAME, optional = true, doc = "Amount of numerical 'noise' allowed in rank sum annotations")
    private double ranksumTolerance = 0.0;

    @Argument(fullName = "likelihood-change-allowed", optional = true)
    private int likelihoodTolerance = 0;

    @Argument(fullName = "ignore-non-ref-data", optional = true)
    private boolean ignoreNonRefData = false;

    @Argument(fullName = "ignore-annotations", optional = true, doc = "Only match on position and alleles, ignoring all INFO annotations", mutex={POSITIONS_ONLY_LONG_NAME})
    private boolean ignoreAnnotations = false;

    @Argument(fullName = "ignore-genotype-annotations", optional = true, doc = "Only match on genotype call, ignoring all FORMAT annotations", mutex={POSITIONS_ONLY_LONG_NAME})
    private boolean ignoreGenotypeAnnotations = false;

    @Argument(fullName = "ignore-genotype-phasing", optional = true, doc = "Only check called genotype alleles, not order or phasing status")
    private boolean ignoreGenotypePhasing = false;

    @Argument(fullName = "ignore-filters", optional = true, doc = "Ignore filter status when comparing variants", mutex={POSITIONS_ONLY_LONG_NAME})
    private boolean ignoreFilters = false;

    @Argument(fullName = "ignore-attribute", optional = true, doc = "Ignore INFO attributes with this key", mutex={POSITIONS_ONLY_LONG_NAME})
    private List<String> ignoreAttributes = new ArrayList<>(5);

    @Argument(fullName = POSITIONS_ONLY_LONG_NAME, optional = true, doc = "Only match on position, ignoring alleles and annotations", mutex={})
    private boolean positionsOnly = false;

    @Argument(fullName = ALLOW_NEW_STARS_LONG_NAME, optional = true, doc = "Allow additional * alleles in actual if there is a corresponding deletion", mutex={POSITIONS_ONLY_LONG_NAME})
    private boolean allowNewStars = false;

    @Argument(fullName = "allow-extra-alleles", optional = true, doc = "Allow extra alleles in actual provided actual is a superset of expected", mutex={ALLOW_NEW_STARS_LONG_NAME, POSITIONS_ONLY_LONG_NAME})
    private boolean allowExtraAlleles = false;

    @Argument(fullName = "allow-missing-stars", optional = true, doc = "Allow missing * in actual if actual has no corresponding deletions")
    private boolean allowMissingStars = false;

    @Argument(fullName = "ignore-star-attributes", optional = true, doc = "Ignore annotation values for spanning deletion (*) alleles")
    private boolean ignoreStarAttributes = false;

    @Argument(fullName = "allow-nan-mismatch", optional = true, doc = "Allow VCF values of missing (.) and NaN to match")
    private boolean allowNanMismatch = false;

    @Argument(fullName = "mute-acceptable-diffs", optional = true, doc = "Suppress warnings or exceptions for differences that are consequences of low quality genotypes or AN discrepancies")
    private boolean muteDiffs = false;

    @Argument(fullName = "ignore-hom-ref-attributes", optional = true, doc = "Skip attribute comparison for homozygous reference genotypes")
    private boolean ignoreHomRefAttributes = false;

    @Argument(fullName = "ignore-dbsnp-ids", optional = true)
    private boolean ignoreDbsnp = false;

    @Advanced
    @Argument(fullName=ReblockGVCF.ANNOTATIONS_TO_KEEP_LONG_NAME, doc="Annotations that are not recognized by GATK that should be checked, otherwise they are removed if ignore-non-ref-data is set.", optional = true)
    private List<String> annotationsToKeep = new ArrayList<>();

    @Hidden
    @Argument(fullName = "output-warnings", optional = true, doc = "Output warnings to a file, useful for testing")
    private GATKPath warningsOutput = null;

    private PrintStream outputStream = null;

    @Override
    public boolean doDictionaryCrossValidation() {
        return false;
    }  //speed things up a bit

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class, AS_StandardAnnotation.class);
    }

    @Override
    public void onTraversalStart() {
        if (getDrivingVariantsFeatureInputs().size() != 2) {
            throw new UserException.BadInput("VCFComparator expects exactly two inputs -- one actual and one expected.");
        }
        final List<FeatureInput<VariantContext>> expected = getDrivingVariantsFeatureInputs().stream().filter(file -> file.getName().equals("expected")).collect(Collectors.toList());
        if (expected.size() != 1) {
            throw new UserException.BadInput("Tool requires exactly one expected input file");
        }

        annotatorEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);

        isSingleSample = getSamplesForVariants().size() == 1 ? true : false;  //if either expected or actual has more than one sample then we're not single-sample

        if (warningsOutput != null) {
            outputStream = new PrintStream(warningsOutput.getOutputStream());
        }
    }

    @Override
    public void apply(final List<VariantContext> variantContexts, final ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
        if (isSingleSample && !variantContexts.get(0).getAlleles().contains(Allele.NON_REF_ALLELE)) {
            throw new UserException.BadInput("Single-sample mode expects two GVCFs with <NON_REF> data for comparison");
        }

        if (variantContexts.size() == 1) {
            final VariantContext vc = variantContexts.get(0);
            //TODO: if there's a variant that overlaps the requested interval the start may not be in the requested intervals
            if (!muteDiffs) {
                throwOrWarn(new UserException("Unmatched variant in " + vc.getSource() + " at position " + vc.getContig() + ":" + vc.getStart()));
            //low coverage sites will have QUAL for singleton hom-vars boosted because of GQ0 hom-refs
            } else if (hasGoodEvidence(vc)) {
                throwOrWarn(new UserException("Unmatched variant in " + vc.getSource() + " at position " + vc.getContig() + ":" + vc.getStart()));
            }
        }

        //here we have all variants overlapping a given position
        //if there's no actual for the expected (at the same position), see if the expected has an overlapping deletion
        for (final VariantContext vc : variantContexts) {
            if (vc.getSource().equals("actual")) {  //there may be more expected than actual variants, but only compare once
                continue;  //TODO: this means we don't report cases where actual has multiple reference blocks spanning expected as unmatched
            }
            //vc is guaranteed to be from expected

            //variantContexts can have a bunch of overlapping variants, but use the start position to compare. Could still have two alleles with same start.
            final List<VariantContext> matches = variantContexts.stream().filter(v -> v.getStart() == vc.getStart()).collect(Collectors.toList());
            //matches includes vc (both actual and expected)
            //TODO: matches.size() can be up to 2 times the number of haplotypes (based on ploidy). So checking only one match won't work for non-diploid organisms
            if (matches.size() == 1 ) {
                // if there's only one match, then it's the expected matching itself, so actual is missing the same vc start
                // if it's a high quality site, throw an error/warning otherwise it's low quality so skip it
                if (isHighQuality(vc) && hasGoodEvidence(vc)
                        && !vc.getGenotype(0).getAlleles().contains(Allele.SPAN_DEL)){
                    throwOrWarn(new UserException("Apparent unmatched high quality variant in " + vc.getSource() + " at " + vc.getContig() + ":" + vc.getStart()));
                } else {
                    return;
                }
            } else {
                if ((!isHighQuality(vc) || (isSingleSample && !GATKVariantContextUtils.genotypeHasConcreteAlt(vc.getGenotype(0).getAlleles())))
                && ignoreHomRefAttributes) {
                    return;
                }
                final VariantContext match;
                // VCFs are expected to NOT have split multi-alleleics
                if (!matches.get(0).getSource().equals(vc.getSource())) {
                    match = matches.get(0);
                } else {
                    match = matches.get(1);
                }

                final List<VariantContext> overlappingDels = variantContexts.stream().filter(v -> v.getStart() < vc.getStart() && v.overlaps(vc))
                        .collect(Collectors.toList());

                final VariantContext expectedTrimmed = trimAlleles(vc, overlappingDels);
                final VariantContext actualTrimmed = trimAlleles(match, overlappingDels);

                //do single-sample GVCF checks, including deletion trimming and dropping
                if (isSingleSample) {
                    try {
                        validateSingleSampleDeletions(vc, match, expectedTrimmed, actualTrimmed, overlappingDels);
                    } catch (UserException e) {
                        throwOrWarn(e);
                    }
                }

                if (positionsOnly) {
                    return;
                }

                //more rigorous checks on annotations and genotypes
                try {
                    checkVariantContextsAreMatching(actualTrimmed, expectedTrimmed, overlappingDels);
                } catch (UserException e) {
                    final boolean hasLowQualityGenotype = expectedTrimmed.getGenotypes().stream().anyMatch(g -> g.getGQ() < 20);
                     if (!muteDiffs || !(alleleNumberIsDifferent || inbreedingCoeffIsDifferent || hasLowQualityGenotype)) {
                        throwOrWarn(e);
                    }
                }

                if (alleleNumberIsDifferent && !muteDiffs) {
                    logger.warn("Observed allele number differed at position " + vc.getContig() + ":" + vc.getStart());
                }
                if (inbreedingCoeffIsDifferent && expectedTrimmed.getGenotypes().stream().anyMatch(g -> g.getGQ() < 20) && !muteDiffs) {
                    logger.warn("Low quality genotype may have caused inbreeding coeff differences at position " + vc.getContig() + ":" + vc.getStart());
                }
            }
        }
    }

    private void validateSingleSampleDeletions(final VariantContext vc, final VariantContext match,
                                               final VariantContext expectedTrimmed, final VariantContext actualTrimmed,
                                               final List<VariantContext> overlappingDels) {
        final Genotype expectedGenotype = expectedTrimmed.getGenotype(0);
        final List<Allele> expectedGenotypeAlleles = expectedGenotype.getAlleles();
        final Genotype actualGenotype = actualTrimmed.getGenotype(0);
        final List<Allele> actualGenotypeAlleles = actualGenotype.getAlleles();

        if (ignoreGenotypePhasing && !actualHasNewAlleles(expectedTrimmed, actualTrimmed)) {
            return;
        }

        // make sure both genotypes are not haploid and do use order here so we can check phasing
        if (!expectedGenotypeAlleles.get(0).equals(actualGenotypeAlleles.get(0)) || ( expectedGenotypeAlleles.size() > 1 && actualGenotypeAlleles.size() > 1 && !expectedGenotypeAlleles.get(1).equals(actualGenotypeAlleles.get(1)))) {
            if (expectedGenotype.isPhased() && expectedGenotypeAlleles.get(1).equals(actualGenotypeAlleles.get(0)) && expectedGenotypeAlleles.get(0).equals(actualGenotypeAlleles.get(1))) {
                throw wrapWithPosition(vc.getContig(), vc.getStart(), new UserException("phasing is swapped. Actual in phaseset "
                        + actualGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, "") + " has "
                        + actualGenotype.getGenotypeString(false)
                        + " expected in phaseset " + expectedGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, "")
                        +" has " + expectedGenotype.getGenotypeString(false)));
            }
            if (overlappingDels.size() == 0) {
                throw wrapWithPosition(vc.getContig(), vc.getStart(), makeGenotypeExceptionFromDifference("called genotype alleles",
                        actualGenotype.getGenotypeString(false),
                        expectedGenotype.getGenotypeString(false)));
            }

            //this is okay if a star got corrected by dropping an upstream hom ref del
            //could be dropped hom ref
            if (overlappingDels.stream().anyMatch(v -> v.getGenotype(0).isHomRef())) {
                try {
                    checkAttributes(vc.getAttributes(), match.getAttributes(), vc.getAlleles(), match.getAlleles(), vc.getPhredScaledQual());
                } catch (UserException e) {
                    throw wrapWithPosition(vc.getContig(), vc.getStart(), new UserException("INFO attributes do not match at " + vc.getContig() + ":" + vc.getStart()));
                }
                try {
                    checkGenotypes(expectedGenotype, actualGenotype);
                } catch (final UserException e) {
                    throw wrapWithPosition(vc.getContig(), vc.getStart(), e);
                }
            } else if (overlappingDels.size() == 0 && expectedGenotypeAlleles.contains(Allele.SPAN_DEL) && !actualGenotypeAlleles.contains(Allele.SPAN_DEL)) {
                return;
            } else {
                if (overlappingDels.stream().anyMatch(v -> v.getSource().equals("actual"))) {  //there should be one overlapping that got trimmed and the other shouldn't overlap anymore
                    if (!allowMissingStars) {
                        throw wrapWithPosition(vc.getContig(), vc.getStart(), new UserException("genotype alleles do not match at spanning deletion site. "
                                +  vc.getSource() + " has " + expectedGenotype.toString() + " and " + match.getSource() + " has "
                                + actualGenotype.toString()));
                    }
                    return;
                } else {
                    final VariantContext overlapper = overlappingDels.get(0);
                    final VariantContext trimmedOverlapper = trimAlleles(overlapper, overlappingDels);
                    if (overlapper.overlaps(vc) && !trimmedOverlapper.overlaps(vc)) {
                        return;
                    }
                }

            }
        }
    }

    private boolean hasGoodEvidence(final VariantContext vc) {
        return (vc.getPhredScaledQual() > goodQualThreshold)
            && (vc.getAttributeAsInt(VCFConstants.DEPTH_KEY,0)/(double)vc.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY,0)) > 5
                && vc.getGenotypes().stream().filter(g -> !g.isHomRef()).mapToInt(Genotype::getGQ).sum() > goodQualThreshold;
    }

    private void checkVariantContextsAreMatching(final VariantContext actual, final VariantContext expected,
                                                 final List<VariantContext> overlappingDels) {
        if (!actual.getContig().equals(expected.getContig())) {
            throw wrapWithPosition(expected.getContig(), expected.getStart(), new UserException("contigs differ for VCs"));
        }

        if (actual.getStart() != expected.getStart()) {
            throw wrapWithPosition(expected.getContig(), expected.getStart(), new UserException("start positions differ for VCs"));
        }

        //don't check end in case we're being lenient about alleles

        //check alleles
        if (actualHasNewAlleles(expected, actual)) {
            try {
                checkAlleles(actual, expected, overlappingDels);
            } catch (final UserException e) {
                throw wrapWithPosition(expected.getContig(), expected.getStart(), e);
            }
        }

        if (!ignoreDbsnp && !actual.getID().equals(expected.getID())) {  //more alleles might mean more dbSNP matches
            throw wrapWithPosition(expected.getContig(), expected.getStart(), new UserException("dbsnp IDs differ for VCs"));
        }

        if (ignoreAnnotations) {
            return;
        }

        if (!ignoreQuals) {
            final double diff = Math.abs(actual.getPhredScaledQual() - expected.getPhredScaledQual());
            if (diff > qualTolerance) {
                throw wrapWithPosition(expected.getContig(), expected.getStart(), new UserException("qual scores differ by " + diff + ", which is more than " + qualTolerance));
            }
        }

        try {
            checkAttributes(actual.getAttributes(), expected.getAttributes(), actual.getAlternateAlleles(), expected.getAlternateAlleles(), expected.getPhredScaledQual());
        } catch (final UserException e) {
            throw wrapWithPosition(expected.getContig() ,expected.getStart() , e);
        }

        if (!alleleNumberIsDifferent) {
            if (!ignoreFilters) {
                if (actual.filtersWereApplied() != expected.filtersWereApplied()) {
                    throw wrapWithPosition(expected.getContig(), expected.getStart(), new UserException(" filters were not applied to both variants"));
                }
                if (!actual.getFilters().equals(expected.getFilters())) {
                    throw wrapWithPosition(expected.getContig(), expected.getStart(), new UserException("variants have different filters: expected has "
                    + (expected.getFilters().isEmpty() ? "PASS" : expected.getFilters()) + " and actual has " + actual.getFilters()));
                }
            }

            for (int i = 0; i < actual.getGenotypes().size(); i++) {
                try {
                    checkGenotypes(actual.getGenotype(0), expected.getGenotype(0));
                } catch (final UserException e) {
                    throw wrapWithPosition(expected.getContig(), expected.getStart(), e);
                }
            }
        }
    }

    private void checkAlleles(final VariantContext actual, final VariantContext expected, final List<VariantContext> overlappingDels) {
        //TODO: this could be more rigorous
        if (ignoreGenotypePhasing && !actualHasNewAlleles(actual, expected)) {
            return;
        }
        if (!allowExtraAlleles && actualHasNewAlleles(actual, expected)) {
            throw new UserException("Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                    + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles());

        } else if (!allowNewStars && hasNewStar(actual, expected)) {
            if (overlappingDels.size() == 0 || !overlappingDels.stream().anyMatch(vc -> vc.getSource().equals("actual"))) {
                throw new UserException("Actual has new unmatched * allele. Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                        + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles());
            }
            //this is a GenomicsDB/CombineGVCFs bug -- there is an overlapping deletion and * should be output, but those tools don't account for multiple haplotypes and upstream variant "ends" the deletion
        } else if (!allowMissingStars && hasMissingStar(actual, expected)) {
            final Set<Allele> remainder = new LinkedHashSet<>(expected.getAlleles());
            remainder.removeAll(actual.getAlleles());
            if (remainder.size() > 1 || !remainder.contains(Allele.SPAN_DEL)) {
                throw new UserException("Actual missing * allele. Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                        + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles());
            }
        } else {
            throw new UserException("Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                    + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles());
        }
    }

    /**
     * Compare VC alternate alleles
     * @param actual
     * @param expected
     * @return
     */
    private boolean actualHasNewAlleles(final VariantContext actual, final VariantContext expected) {
        return !actual.getAlternateAlleles().stream().allMatch(a -> GATKVariantContextUtils.isAlleleInList(
                actual.getReference(), a, expected.getReference(), expected.getAlternateAlleles()));
    }

    private boolean hasMissingStar(final VariantContext actual, final VariantContext expected) {
        return expected.getAlleles().contains(Allele.SPAN_DEL) && !actual.getAlleles().contains(Allele.SPAN_DEL);
    }

    private boolean hasNewStar(final VariantContext actual, final VariantContext expected) {
        return actual.getAlleles().contains(Allele.SPAN_DEL) && !expected.getAlleles().contains(Allele.SPAN_DEL);
    }

    private UserException wrapWithPosition(final String contig, final int start, final UserException e) {
        return new UserException("At position " + contig + ":" +start + " " + e.getMessage());
    }

    @SuppressWarnings("unchecked")
    private void checkAttributes(final Map<String, Object> actual, final Map<String, Object> expected,
                                 final List<Allele> actualAlts, final List<Allele> expectedAlts,
                                 final double expectedQual) {
        final Set<String> expectedKeys = new LinkedHashSet<>(expected.keySet());

        //do a precheck on AN because then we can't expect the rest of the annotations to match
        alleleNumberIsDifferent = actual.containsKey(VCFConstants.ALLELE_NUMBER_KEY) && !isAttributeValueEqual(VCFConstants.ALLELE_NUMBER_KEY, actual.get(VCFConstants.ALLELE_NUMBER_KEY),
                    expected.get(VCFConstants.ALLELE_NUMBER_KEY));

        if (actual.containsKey(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)) {
            try {
                inbreedingCoeffIsDifferent = !isAttributeEqualDoubleSmart(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY,
                        Double.parseDouble(actual.get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY).toString()),
                        Double.parseDouble(expected.get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY).toString()), 0.001);
            } catch (UserException e) {
                inbreedingCoeffIsDifferent = true;
                throw e;
            }
        }

        if (actual.containsKey(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY) && !(actual.get(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY) instanceof List)) {
            try {
                inbreedingCoeffIsDifferent = !isAttributeEqualDoubleSmart(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY,
                        Double.parseDouble(actual.get(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY).toString()),
                        Double.parseDouble(expected.get(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY).toString()), 0.001);
            } catch (final UserException e) {
                inbreedingCoeffIsDifferent = true;
                throw e;
            }
        }

        for (final Map.Entry<String, Object> exp : expected.entrySet()) {
            final Object expectedValue = exp.getValue();
            final String key = exp.getKey();
            if (actual.containsKey(key) && actual.get(key) != null) {
                //we already checked these
                if (key.equals(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)
                        || key.equals(VCFConstants.ALLELE_NUMBER_KEY)
                        || (key.equals(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY) && !(expected.get(key) instanceof List))
                        || ignoreAttributes.contains(key)) {
                    continue;
                }
                final Object actualValue = actual.get(key);
                if (expectedValue instanceof List && actualValue instanceof List || key.equals(GATKVCFConstants.AS_SB_TABLE_KEY)) {
                    // both values are lists, compare element by element
                    final List<? extends Object> expectedList, actualList;
                    if (key.contains("AS_") || key.equals(GATKVCFConstants.AS_SB_TABLE_KEY)) {
                        expectedList = AnnotationUtils.decodeAnyASListWithRawDelim(expectedValue.toString());
                        actualList = AnnotationUtils.decodeAnyASListWithRawDelim(actualValue.toString());
                    } else {
                        expectedList = (List<Object>) expectedValue;
                        actualList = (List<Object>) actualValue;
                    }
                    if (actualList.size() != expectedList.size() && !key.contains("AS_")) {  //TODO: account for subset alleles
                        throw makeVariantExceptionFromDifference("attributes", actual.keySet().toString(),
                                expected.keySet().toString());
                    }
                    final List<? extends Object> expectedListCommon, actualListCommon;
                    if (key.equals(VCFConstants.ALLELE_COUNT_KEY) || key.equals(GATKVCFConstants.MLE_ALLELE_COUNT_KEY)) {
                        if (!isACEqualEnough(key, actualList, expectedList, actualAlts, expectedAlts)) {
                            throw makeVariantExceptionFromDifference(key, actualList.toString(), expectedList.toString());
                        }
                    }
                    //we've already gotten rid of stars... I think -- do we have to sort the alts?  Maybe...  I should be using my own AlleleSpecificAnnotationData!!!
                    final int iterationEnd = ignoreNonRefData && actualAlts.contains(Allele.NON_REF_ALLELE) ? expectedAlts.size()-1 : expectedAlts.size();  //exclusive
                    for (int i = 0; i < iterationEnd; i++) {
                        if (i >= actualList.size() || i >= expectedList.size()) {
                            //TODO: the toStrings here don't do a great job
                            throw makeVariantExceptionFromDifference(key, actualValue.toString(), expectedValue.toString());
                        }
                        if (key.equals(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY)) {
                            final double actualPerAlleleValue = Double.parseDouble(actualList.get(i).toString());
                            final double expectedPerAlleleValue = Double.parseDouble(expectedList.get(i).toString());
                            if (!isAttributeEqualDoubleSmart(key,
                                    actualPerAlleleValue, expectedPerAlleleValue,
                                    inbreedingCoeffTolerance)) {
                                throw makeVariantExceptionFromDifference(key, Double.toString(actualPerAlleleValue), Double.toString(expectedPerAlleleValue));
                            }
                        }
                        if (key.contains("AS_") && key.contains("RankSum")) {
                            final double actualPerAlleleValue, expectedPerAlleleValue;
                            if ((!((String) actualList.get(i)).isEmpty() && ((String) actualList.get(i)).equals("NaN"))) {
                                actualPerAlleleValue = Double.parseDouble(actualList.get(i).toString());
                                expectedPerAlleleValue = Double.parseDouble(expectedList.get(i).toString());
                                if (!isAttributeEqualDoubleSmart(key,
                                        actualPerAlleleValue, expectedPerAlleleValue,
                                        ranksumTolerance)) {
                                    throw makeVariantExceptionFromDifference(key, Double.toString(actualPerAlleleValue), Double.toString(expectedPerAlleleValue));
                                }
                                continue;
                            } else {
                                if (!muteDiffs && !allowNanMismatch) {
                                    logger.warn("GATK version-specific NaN versus empty AS_RAW annotation discrepancy");
                                }
                            }
                        }
                        //this needs some cleanup -- two exceptions??
                       final boolean check;
                        try
                       {check = isAttributeValueEqual(key, actualList.get(i), expectedList.get(i)); }
                       catch (final UserException e){
                           if (!ignoreStarAttributes && !actualAlts.get(i).equals(Allele.SPAN_DEL) && !actualAlts.get(i).equals(Allele.SPAN_DEL)) {
                               throw makeVariantExceptionFromDifference(key, actualList.get(i).toString(), expectedList.get(i).toString());
                           }
                        }
                    }
                } else {  //if not list
                    try {
                        if (key.contains("RankSum") && !key.contains("AS_")) {
                            if (!isAttributeEqualDoubleSmart(key, Double.parseDouble(actualValue.toString()), Double.parseDouble(expectedValue.toString()),
                                    ranksumTolerance)) {
                                throw makeVariantExceptionFromDifference(key, actualValue.toString(), actualValue.toString());
                            }
                        } else {
                            isAttributeValueEqual(key, actualValue, expectedValue);
                        }
                    } catch (final UserException e) {
                        if (annotationsThatVaryWithNoCalls.contains(key) && alleleNumberIsDifferent) {
                            continue;
                        } else {
                            if (key.equals(GATKVCFConstants.QUAL_BY_DEPTH_KEY) || key.equals(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY)) {
                                final double actualDouble = Double.parseDouble(actualValue.toString());
                                final double expectedDouble = Double.parseDouble(expectedValue.toString());
                                final int expectedDP = Integer.parseInt(expected.get(VCFConstants.DEPTH_KEY).toString());
                                if (qualByDepthDifferenceIsAcceptable(actualDouble, expectedDouble, expectedQual, expectedDP)) {  //QD won't match if DP doesn't
                                    if (!muteDiffs) {
                                        logger.warn(key + " difference is within expected tolerances");
                                    }
                                    continue;
                                } else {
                                    final boolean isDepthEqual;
                                    try {
                                        isDepthEqual = isAttributeValueEqual(VCFConstants.DEPTH_KEY, actual.get(VCFConstants.DEPTH_KEY), expected.get(VCFConstants.DEPTH_KEY));
                                    } catch (final UserException e2) {
                                        //check expected QD against expected AS_QD because in a few cases actual is wrong
                                        if (!qualByDepthDifferenceIsAcceptable(Double.parseDouble(actual.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY).toString()),
                                        Double.parseDouble(actual.get(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY).toString()), expectedQual, expectedDP)) {
                                            logger.warn(key + " difference (actual = " + actual.get(key) + " versus expected:"
                                                    + expected.get(key) + " is larger than expected, but so is DP (actual="
                                                    + actual.get(VCFConstants.DEPTH_KEY) + " versus expected:" + expected.get(VCFConstants.DEPTH_KEY) + ")");
                                        }
                                        continue;
                                    }
                                    if (isDepthEqual) {
                                        throw e;
                                    }
                                }
                            }
                        }
                    throw e;
                    }
                }
            }
            expectedKeys.remove(key);
        }
    }

    /**
     * Generate an exception to throw later
     * @param thingToCompare    the key to appear in the exception
     * @param actual
     * @param expected
     * @return
     */
    private UserException makeVariantExceptionFromDifference(final String thingToCompare, final String actual, final String expected) {
        return new UserException("Variant contexts have different " + thingToCompare + ": actual has " + actual
                + " expected has " + expected);
    }

    /**
     * Generate an exception to throw later
     * @param thingToCompare    the key to appear in the exception
     * @param actual
     * @param expected
     * @return
     */
    private UserException makeGenotypeExceptionFromDifference(final String thingToCompare, final String actual, final String expected) {
        return new UserException("Genotypes have different " + thingToCompare + ": actual has " + actual
                + " expected has " + expected);
    }


    /**
     *
     * @param key
     * @param actual    input VariantContext containing key
     * @param expected  input VariantContext containing key
     * @return
     */
    private boolean isAttributeValueEqual(final String key, final Object actual, final Object expected) {
        if (allowNanMismatch && (actual.toString()).equals(".") && expected.toString().equals("NaN")) {
            return true;
        }
        if (!actual.toString().equals(expected.toString())) {
            throw new UserException("Variant contexts have different attribute values for " + key + ": actual has " + actual.toString()
                    + " and expected has " + expected.toString());
        }
        return actual.toString().equals(expected.toString());

    }

    private boolean isAttributeEqualDoubleSmart(final String key, final double actual, final double expected, final double tolerance) {
        final double diff = Math.abs(actual - expected);
        if (diff > tolerance) {
            throw new UserException("Attribute " + key + " has difference " + String.format("%.3f", diff) + ", which is larger difference than allowed delta " + tolerance);
        }
        return true;
    }

    private void checkGenotypes(final Genotype actual, final Genotype expected) {
        if (!actual.getSampleName().equals(expected.getSampleName())) {
            throw new UserException("Sample names do not match");
        }
        if (!CollectionUtils.isEqualCollection(actual.getAlleles(), expected.getAlleles())) {
            throw makeGenotypeExceptionFromDifference("alleles", actual.getAlleles().toString(), expected.getAlleles().toString());
        }

        if (!ignoreGenotypePhasing) {
            if (!actual.getGenotypeString(false).equals(expected.getGenotypeString(false))) {
                throw makeGenotypeExceptionFromDifference("genotype string", actual.getGenotypeString(false), expected.getGenotypeString(false));
            }
            if (actual.isPhased() != expected.isPhased()) {
                throw makeGenotypeExceptionFromDifference("phasing status", Boolean.toString(actual.isPhased()), Boolean.toString(expected.isPhased()));
            }
        }
        if (ignoreHomRefAttributes && actual.isHomRef()) {
            return;
        }
        if (actual.hasDP() != expected.hasDP()) {
            throw makeGenotypeExceptionFromDifference("DP presence", Boolean.toString(actual.hasDP()), Boolean.toString(expected.hasDP()));
        }
        if (actual.hasAD() != expected.hasAD()) {
            throw makeGenotypeExceptionFromDifference("AD presence", Boolean.toString(actual.hasAD()), Boolean.toString(expected.hasAD()));
        }
        if (actual.hasGQ() != expected.hasGQ()) {
            throw makeGenotypeExceptionFromDifference("GQ presence", Boolean.toString(actual.hasGQ()), Boolean.toString(expected.hasGQ()));
        }
        if (!actual.isHomRef()) {
            if (actual.hasGQ() && (Math.abs(actual.getGQ() - expected.getGQ())) > likelihoodTolerance) {
                //allele subsetting can return GQs over 99, which apparently are capped on VCF write
                if (actual.getGQ() < 99 || expected.getGQ() < 99) {
                    throw makeGenotypeExceptionFromDifference("GQ value", Integer.toString(actual.getGQ()), Integer.toString(expected.getGQ()));
                }
            }
        }
        if (ignoreGenotypeAnnotations) {
            return;
        }
        //allow DP difference due to non-ref AD correction
        if (actual.hasDP() && expected.hasDP()) {
            if (actual.getDP() != expected.getDP()) {
                if (dpChange == 0 && expected.getDP() - actual.getDP() != expected.getAD()[expected.getAD().length - 1]) {
                    throw makeGenotypeExceptionFromDifference("DP value", Integer.toString(actual.getDP()), Integer.toString(expected.getDP()));
                }
                if (actual.getDP() - expected.getDP() > dpChange) {
                    throw new UserException("DP difference exceeds allowable tolerance of " + dpChange
                            + ", actual has " + Integer.toString(actual.getDP()) + " expected has " + Integer.toString(expected.getDP()));
                }
                if (actual.hasAD() && (!Arrays.equals(actual.getAD(), expected.getAD()))) {
                    if (dpChange == 0 && actual.getAD()[actual.getAD().length - 1] == expected.getAD()[expected.getAD().length - 1]) {  //non-ref differences are okay, but if there's an AD different that's not non-ref we're in trouble
                        throw makeGenotypeExceptionFromDifference("AD values", Arrays.toString(actual.getAD()), Arrays.toString(expected.getAD()));
                    }
                    for (int i = 0; i < actual.getAD().length - 1; i++) {
                        if (actual.getAD()[i] - expected.getAD()[i] > dpChange) {
                            throw new UserException("AD difference exceeds allowable tolerance of " + dpChange
                                    + ", actual has " + actual.getAD().toString() + " expected has " + expected.getAD().toString());
                        }
                    }
                }
            }

            if (actual.hasPL() != expected.hasPL()) {
                throw makeGenotypeExceptionFromDifference("PL presence", Boolean.toString(actual.hasPL()), Boolean.toString(expected.hasPL()));
            }
            if (actual.hasPL() && expected.hasPL()) {
                final int[] actualPls = actual.getPL();
                final int[] expectedPls = expected.getPL();
                if (actualPls.length != expectedPls.length) {
                    throw new UserException("PL lengths for genotype vary: actual is size " + actualPls.length + " and expected is " + expectedPls.length);
                }
                for (int i = 0; i < actualPls.length; i++) {
                    if (Math.abs(actualPls[i] - expectedPls[i]) > likelihoodTolerance) {
                        throw makeGenotypeExceptionFromDifference("PL value", Arrays.toString(actual.getPL()), Arrays.toString(expected.getPL()));
                    }
                }
            }
        }
    }

    private void throwOrWarn(final UserException e) {
        if (muteDiffs && (alleleNumberIsDifferent || inbreedingCoeffIsDifferent)) {
            return;
        }
        if (warnOnError || finishBeforeFailing) {
            if (warningsOutput != null) {
                outputStream.println(e.getMessage());
            }
            logger.warn("***** " + e.getMessage() + " *****");
            if (finishBeforeFailing) {
                failOnCompletion = true;
            }
        } else {
            throw e;
        }
    }

    //Non-reblocked GVCFs can have alleles that aren't called, so trim them
    private VariantContext trimAlleles(final VariantContext variant, final List<VariantContext> overlappingDels) {
        final Allele ref = variant.getReference();
        final Set<Allele> relevantAlleles = new LinkedHashSet<>();
        for (final Genotype g : variant.getGenotypes()) {
            final List<Allele> gtAlleles = g.getAlleles();
            relevantAlleles.add(ref);
            if (!gtAlleles.contains(Allele.NO_CALL)) {  //in theory there could be half no-calls, but I've never seen one
                if (overlappingDels.size() > 0) {
                    relevantAlleles.addAll(gtAlleles);
                } else relevantAlleles.addAll(gtAlleles.stream().filter(GATKVariantContextUtils::isConcreteAlt).collect(Collectors.toList()));
            }
        }
        final VariantContextBuilder vcBuilder = new VariantContextBuilder(variant);
        final List<Allele> orderedRelevantAlleles = new ArrayList<>(relevantAlleles);
        if (vcBuilder.getGenotypes().size() == 1 && !ignoreNonRefData) {
            orderedRelevantAlleles.add(Allele.NON_REF_ALLELE);
        }
        vcBuilder.alleles(orderedRelevantAlleles);

        //NOTE that we use BEST_MATCH_TO_ORIGINAL for post-reblocked VCFs with no hom ref PLs
        final GenotypesContext gc = AlleleSubsettingUtils.subsetAlleles(variant.getGenotypes(), defaultPloidy,
                variant.getAlleles(), orderedRelevantAlleles, null, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
        vcBuilder.genotypes(gc);

        final Map<String,Object> subsetAnnotations = ReblockGVCF.subsetAnnotationsIfNecessary(annotatorEngine, false, VCFConstants.GENOTYPE_POSTERIORS_KEY, variant, vcBuilder.make(), annotationsToKeep);

        return variant.getGenotype(0).isHomRef() ? vcBuilder.make() : GATKVariantContextUtils.reverseTrimAlleles(vcBuilder.attributes(subsetAnnotations).make());
    }

    /*
     * Some code paths expect significant differences, so focusing on high quality variants is useful.
     */
    private boolean isHighQuality(final VariantContext vc) {
        if (vc.getGenotypes().size() == 1) {
            //genotyping engine returns 0.01 for hom-ref SNPs, probably because of prior
            final Genotype genotype = vc.getGenotype(0);
            if (!genotype.hasPL()) {
                return false;
            }
            if (genotype.isHomRef()) {
                return false;
            }
            if (genotype.getPL()[0] == 0) {
                return false;
            }
            if (genotype.getAlleles().contains(Allele.NON_REF_ALLELE)) {
                return false;
            }
            return vc.getPhredScaledQual() > 0.01 && (MathUtils.arrayMax(genotype.getPL()) > 0);
        } else {
            for (final Genotype genotype : vc.getGenotypes()) {
                if (genotype.isHomRef()) {
                    continue;
                }
                if (passesGnomadAdjCriteria(genotype)) {
                    return true;
                }
            }
            return false;
        }
    }

    /*
     * This is for known significant changes so we can compare only high quality genotypes as defined by these gnomAD adj criteria.
     * https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/
     */
    private boolean passesGnomadAdjCriteria(final Genotype g) {
        if (!g.hasGQ() || g.getGQ() < 20) {
            return false;
        }
        if (!g.hasDP() || g.getDP() < 10) {
            return false;
        }
        if (g.isHet()) {
            if (!g.hasAD()) {
                return false;
            }
            final int[] ad = g.getAD();
            final double alleleBalance;
            if (g.isHetNonRef()) {
                alleleBalance = ad[2] / ((double) ad[1] + ad[2]);

            } else {
                alleleBalance = ad[1] / ((double) ad[0] + ad[1]);
            }
            if (alleleBalance >= 0.2 && alleleBalance <= 0.8) {
                return true;
            }
        }
        return false;
    }

    private boolean isACEqualEnough(final String key, final List<? extends Object> actualList, final List<? extends Object> expectedList,
                                    final List<Allele> actualAlts, final List<Allele> expectedAlts) {
        if (actualList.size() != expectedList.size()) {
            return false;
        }
        //trimming may drop AC=0 alts, but GATK convention is to put them last, so that shouldn't be an issue
        for (int i = 0; i < actualAlts.size(); i++) {
            if (actualAlts.get(i).equals(Allele.SPAN_DEL)) {
                continue;
            }
            if (actualAlts.get(i).equals(expectedAlts.get(i))) {
                final int expectedValue = Integer.parseInt(expectedList.get(i).toString());
                final int actualValue = Integer.parseInt(actualList.get(i).toString());
                //TODO: and....?
            } else {
                throwOrWarn(new UserException("Alleles are not ordered the same"));
                return false;
            }
        }
        return true;
    }

    private boolean qualByDepthWillHaveJitter(final double expectedQual, final int expectedDepth) {
        final double qdEstimate = expectedQual / (double) expectedDepth;
        if (QualByDepth.fixTooHighQD(qdEstimate) != qdEstimate || qdEstimate > 34.9) {  //I don't know why we get different values at 34.9 -- there must be rounding somewhere?
            return true;
        }
        return false;
    }

    private boolean qualByDepthDifferenceIsAcceptable(final double actual, final double expected, final double expectedQual, final int expectedDP) {
        final double diff = Math.abs(expected - actual);
        final double relativeDiff = diff / (expected);
        return expected > 25.0 || relativeDiff < 0.01 || diff < 0.5 ////25 is in the "jitter" zone
                || qualByDepthWillHaveJitter(expectedQual, expectedDP);
    }

    public Object onTraversalSuccess() {
        if (failOnCompletion) {
            throw new UserException("Some comparisons failed.  See stderr log for details.");
        }
        return null;
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
