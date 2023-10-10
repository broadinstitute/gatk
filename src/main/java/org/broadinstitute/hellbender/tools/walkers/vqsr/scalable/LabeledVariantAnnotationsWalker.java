package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiplePassVariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * Base walker for both {@link ExtractVariantAnnotations} and {@link ScoreVariantAnnotations},
 * which enforces identical variant-extraction behavior in both tools via {@link #extractVariantMetadata}.
 *
 * This base implementation covers functionality for {@link ExtractVariantAnnotations}. Namely, it is a single-pass
 * walker, performing the operations:
 *
 *   - nthPassApply(n = 0)
 *      - if variant/alleles pass filters and variant-type/resource-match checks, then:
 *          - add variant/alleles to a {@link LabeledVariantAnnotationsData} collection
 *          - write variant/alleles with labels appended to a sites-only VCF file
 *   - afterNthPass(n = 0)
 *      - write the resulting {@link LabeledVariantAnnotationsData} collection to an HDF5 file
 *
 * This results in the following output:
 *
 *   - an HDF5 file, with the directory structure documented in {@link LabeledVariantAnnotationsData#writeHDF5};
 *     note that the matrix of annotations contains a single row per datum (i.e., per allele, in allele-specific mode,
 *     and per variant otherwise)
 *   - a sites-only VCF file, containing a single line per extracted variant, with labels appended
 *
 * In contrast, the {@link ScoreVariantAnnotations} implementation overrides methods to yield a two-pass walker,
 * performing the operations:
 *
 *   - nthPassApply(n = 0)
 *      - if variant/alleles pass filters and variant-type checks, then:
 *          - add variant/alleles to a {@link LabeledVariantAnnotationsData} collection
 *   - afterNthPass(n = 0)
 *      - write the resulting {@link LabeledVariantAnnotationsData} collection to an HDF5 file
 *      - pass this annotations HDF5 file to a {@link VariantAnnotationsScorer}, which generates and writes scores to an HDF5 file
 *      - read the scores back in and load them into an iterator
 *   - nthPassApply(n = 1)
 *      - if variant/alleles pass filters and variant-type checks (which are identical to the first pass), then:
 *          - draw the corresponding score (or scores, in allele-specific mode) from the iterator
 *          - write the variant (with all alleles, not just those extracted) with the score
 *            (or best score, in allele-specific mode) appended to a VCF file
 *      - else:
 *          - write an unprocessed copy of the variant to a VCF file
 *
 * This results in the following output:
 *
 *   - an HDF5 file, as above
 *   - a VCF file, containing the input variants, with labels, scores, and filters appended/applied for those passing variant-type checks
 */
public abstract class LabeledVariantAnnotationsWalker extends MultiplePassVariantWalker {

    public static final String MODE_LONG_NAME = "mode";
    public static final String IGNORE_FILTER_LONG_NAME = "ignore-filter";
    public static final String IGNORE_ALL_FILTERS_LONG_NAME = "ignore-all-filters";
    public static final String DO_NOT_TRUST_ALL_POLYMORPHIC_LONG_NAME = "do-not-trust-all-polymorphic";
    public static final String RESOURCE_MATCHING_STRATEGY_LONG_NAME = "resource-matching-strategy";
    public static final String OMIT_ALLELES_IN_HDF5_LONG_NAME = "omit-alleles-in-hdf5";
    public static final String DO_NOT_GZIP_VCF_OUTPUT_LONG_NAME = "do-not-gzip-vcf-output";

    public static final String ANNOTATIONS_HDF5_SUFFIX = ".annot.hdf5";

    public static final String RESOURCE_LABEL_INFO_HEADER_LINE_FORMAT_STRING = "This site was labeled as %s according to resources";

    enum ResourceMatchingStrategy {
        START_POSITION, START_POSITION_AND_GIVEN_REPRESENTATION, START_POSITION_AND_MINIMAL_REPRESENTATION
    }

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Prefix for output filenames.")
    String outputPrefix;

    @Argument(
            fullName = StandardArgumentDefinitions.RESOURCE_LONG_NAME,
            doc = "Resource VCFs used to label extracted variants.",
            optional = true)
    private List<FeatureInput<VariantContext>> resources = new ArrayList<>(10);

    @Argument(
            fullName = StandardArgumentDefinitions.ANNOTATION_LONG_NAME,
            shortName = StandardArgumentDefinitions.ANNOTATION_SHORT_NAME,
            doc = "Names of the annotations to extract. Note that a requested annotation may in fact not be present " +
                    "at any extraction site; NaN missing values will be generated for such annotations.",
            minElements = 1)
    List<String> annotationNames = new ArrayList<>();

    @Argument(
            fullName = MODE_LONG_NAME,
            doc = "Variant types to extract.",
            minElements = 1)
    private List<VariantType> variantTypesToExtractList = new ArrayList<>(Arrays.asList(VariantType.SNP, VariantType.INDEL));

    @Argument(
            fullName = IGNORE_FILTER_LONG_NAME,
            doc = "Ignore the specified filter(s) in the input VCF.",
            optional = true)
    private List<String> ignoreInputFilters = new ArrayList<>();

    @Argument(
            fullName = IGNORE_ALL_FILTERS_LONG_NAME,
            doc = "If true, ignore all filters in the input VCF.",
            optional = true)
    private boolean ignoreAllFilters = false;

    // TODO this is a perhaps vestigial argument inherited from VQSR; its impact and necessity could be reevaluated
    @Argument(
            fullName = DO_NOT_TRUST_ALL_POLYMORPHIC_LONG_NAME,
            doc = "If true, do not trust that unfiltered records in the resources contain only polymorphic sites. " +
                    "This may increase runtime if the resources are not sites-only VCFs.",
            optional = true)
    private boolean doNotTrustAllPolymorphic = false;


    @Argument(
            fullName = RESOURCE_MATCHING_STRATEGY_LONG_NAME,
            doc = "The strategy to use for determining whether an input variant is present in a resource " +
                    "in non-allele-specific mode. " +
                    "START_POSITION: Start positions of input and resource variants must match. " +
                    "START_POSITION_AND_GIVEN_REPRESENTATION: The intersection of the sets of input and resource alleles " +
                    "(in their given representations) must also be non-empty. " +
                    "START_POSITION_AND_MINIMAL_REPRESENTATION: The intersection of the sets of input and resource alleles " +
                    "(after converting alleles to their minimal representations) must also be non-empty. " +
                    "This argument has no effect in allele-specific mode, " +
                    "in which the minimal representations of the input and resource alleles must match.",
            optional = true)
    private ResourceMatchingStrategy resourceMatchingStrategy = ResourceMatchingStrategy.START_POSITION;
    @Argument(
            fullName = OMIT_ALLELES_IN_HDF5_LONG_NAME,
            doc = "If true, omit alleles in output HDF5 files in order to decrease file sizes.",
            optional = true
    )
    boolean omitAllelesInHDF5 = false;

    @Argument(
            fullName = DO_NOT_GZIP_VCF_OUTPUT_LONG_NAME,
            doc = "If true, VCF output will not be compressed.",
            optional = true
    )
    boolean doNotGZIPVCFOutput = false;

    private final Set<String> ignoreInputFilterSet = new TreeSet<>();
    Set<VariantType> variantTypesToExtract;
    TreeSet<String> resourceLabels = new TreeSet<>();
    boolean useASAnnotations;

    File outputAnnotationsFile;
    VariantContextWriter vcfWriter;

    LabeledVariantAnnotationsData data;

    @Override
    public void onTraversalStart() {

        ignoreInputFilterSet.addAll(ignoreInputFilters);

        variantTypesToExtract = EnumSet.copyOf(variantTypesToExtractList);

        outputAnnotationsFile = new File(outputPrefix + ANNOTATIONS_HDF5_SUFFIX);
        final String vcfSuffix = doNotGZIPVCFOutput ? ".vcf" : ".vcf.gz";
        final File outputVCFFile = new File(outputPrefix + vcfSuffix);

        // TODO this validation method should perhaps be moved outside of the CNV code
        CopyNumberArgumentValidationUtils.validateOutputFiles(outputAnnotationsFile, outputVCFFile);

        for (final FeatureInput<VariantContext> resource : resources) {
            final TreeSet<String> trackResourceLabels = resource.getTagAttributes().entrySet().stream()
                    .filter(e -> e.getValue().equals("true"))
                    .map(Map.Entry::getKey)
                    .sorted()
                    .collect(Collectors.toCollection(TreeSet::new));
            resourceLabels.addAll(trackResourceLabels);
            logger.info( String.format("Found %s track: labels = %s", resource.getName(), trackResourceLabels));
        }
        resourceLabels.forEach(String::intern); // TODO evaluate if this affects memory usage and remove if not needed

        if (resourceLabels.contains(LabeledVariantAnnotationsData.SNP_LABEL)) {
            throw new UserException.BadInput(String.format("The resource label \"%s\" is reserved for labeling variant types.",
                    LabeledVariantAnnotationsData.SNP_LABEL));
        }

        useASAnnotations = isAlleleSpecificAnnotationRequested();

        if (useASAnnotations && resourceMatchingStrategy != ResourceMatchingStrategy.START_POSITION_AND_MINIMAL_REPRESENTATION) {
            logger.warn(String.format("The %s argument is ignored when allele-specific annotations are requested. The START_POSITION_AND_MINIMAL_REPRESENTATION strategy will be used.",
                    RESOURCE_MATCHING_STRATEGY_LONG_NAME));
            resourceMatchingStrategy = ResourceMatchingStrategy.START_POSITION_AND_MINIMAL_REPRESENTATION;
        }

        data = new LabeledVariantAnnotationsData(annotationNames, resourceLabels, useASAnnotations);
        logger.info(String.format("Using %d annotations %s...", data.getSortedAnnotationNames().size(), data.getSortedAnnotationNames()));

        vcfWriter = createVCFWriter(outputVCFFile);
        vcfWriter.writeHeader(constructVCFHeader(data.getSortedLabels()));

        afterOnTraversalStart();   // perform additional validation, set modes in child tools, etc.
    }

    public void afterOnTraversalStart() {
        // override
    }

    @Override
    protected int numberOfPasses() {
        return 1;
    }

    @Override
    public Object onTraversalSuccess() {
        return null;
    }

    private boolean isAlleleSpecificAnnotationRequested() {
        final Set<String> distinctAnnotationNames = new LinkedHashSet<>(annotationNames);
        final VCFHeader inputHeader = getHeaderForVariants();
        return distinctAnnotationNames.stream().anyMatch(a -> inputHeader.getInfoHeaderLine(a).getCountType() == VCFHeaderLineCount.A);
    }

    static void addExtractedVariantToData(final LabeledVariantAnnotationsData data,
                                          final VariantContext variant,
                                          final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata) {
        data.add(variant,
                metadata.stream().map(Triple::getLeft).collect(Collectors.toList()),
                metadata.stream().map(Triple::getMiddle).collect(Collectors.toList()),
                metadata.stream().map(Triple::getRight).collect(Collectors.toList()));
    }

    void writeExtractedVariantToVCF(final VariantContext variant,
                                    final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata) {
        writeExtractedVariantToVCF(variant,
                metadata.stream().map(Triple::getLeft).flatMap(List::stream).collect(Collectors.toList()),
                metadata.stream().map(Triple::getRight).flatMap(Set::stream).collect(Collectors.toSet()));
    }

    void writeAnnotationsToHDF5() {
        if (data.size() == 0) {
            logger.warn("Found no input variants for extraction. This may be because the specified " +
                    "genomic region contains no input variants of the requested type(s) or, if extracting " +
                    "training labels, because none of the input variants were contained in the resource VCFs " +
                    "or no resource VCFs were provided. The annotations HDF5 file will not be generated.");
            return;
        }
        for (final VariantType variantType : variantTypesToExtract) {
            logger.info(String.format("Extracted annotations for %d variants of type %s.",
                    data.getVariantTypeFlat().stream().mapToInt(t -> t == variantType ? 1 : 0).sum(), variantType));
        }
        for (final String label : data.getSortedLabels()) {
            logger.info(String.format("Extracted annotations for %d variants labeled as %s.",
                    data.isLabelFlat(label).stream().mapToInt(b -> b ? 1 : 0).sum(), label));
        }
        logger.info(String.format("Extracted annotations for %s total records.", data.size()));
        logger.info(String.format("Extracted annotations for %s total variants.", data.flatSize()));

        logger.info("Writing annotations...");
        data.writeHDF5(outputAnnotationsFile, omitAllelesInHDF5);
        logger.info(String.format("Annotations and metadata written to %s.", outputAnnotationsFile.getAbsolutePath()));
    }

    /**
     * Writes a sites-only VCF containing the extracted variants and corresponding labels.
     */
    void writeExtractedVariantToVCF(final VariantContext vc,
                                    final List<Allele> altAlleles,
                                    final Set<String> labels) {
        final List<Allele> alleles = ListUtils.union(Collections.singletonList(vc.getReference()), altAlleles);
        final VariantContextBuilder builder = new VariantContextBuilder(
                vc.getSource(), vc.getContig(), vc.getStart(), vc.getEnd(), alleles);
        labels.forEach(l -> builder.attribute(l, true)); // labels should already be sorted as a TreeSet
        vcfWriter.add(builder.make());
    }

    // modified from VQSR code
    // TODO we're just writing a standard sites-only VCF here, maybe there's a nicer way to do this?
    VCFHeader constructVCFHeader(final List<String> sortedLabels) {
        Set<VCFHeaderLine> hInfo = sortedLabels.stream()
                .map(l -> new VCFInfoHeaderLine(l, 1, VCFHeaderLineType.Flag, String.format(RESOURCE_LABEL_INFO_HEADER_LINE_FORMAT_STRING, l)))
                .collect(Collectors.toCollection(TreeSet::new));
        hInfo.add(GATKVCFHeaderLines.getFilterLine(VCFConstants.PASSES_FILTERS_v4));
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (sequenceDictionary != null) {
            hInfo = VcfUtils.updateHeaderContigLines(hInfo, referenceArguments.getReferencePath(), sequenceDictionary, true);
        }
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        return new VCFHeader(hInfo);
    }

    /**
     * Performs variant-filter and variant-type checks to determine variants/alleles suitable for extraction, and returns
     * a corresponding list of metadata. This method should not be overridden, as it is intended to enforce identical
     * variant-extraction behavior in all child tools. Logic here and below for filtering and determining variant type
     * was retained from VQSR, but has been heavily refactored.
     */
    final List<Triple<List<Allele>, VariantType, TreeSet<String>>> extractVariantMetadata(final VariantContext vc,
                                                                                          final FeatureContext featureContext,
                                                                                          final boolean isExtractUnlabeled) {
        // if variant is filtered, do not consume here
        if (vc == null || !(ignoreAllFilters || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()))) {
            return Collections.emptyList();
        }
        if (!useASAnnotations) {
            // in non-allele-specific mode, get a singleton list of the triple
            // (list of alt alleles passing variant-type and resource-match checks, variant type, set of labels)
            final VariantType variantType = VariantType.getVariantType(vc);
            if (variantTypesToExtract.contains(variantType)) {
                final TreeSet<String> matchingResourceLabels = findMatchingResourceLabels(vc, null, featureContext);
                if (isExtractUnlabeled || !matchingResourceLabels.isEmpty()) {
                    return Collections.singletonList(Triple.of(vc.getAlternateAlleles(), variantType, matchingResourceLabels));
                }
            }
        } else {
            // in allele-specific mode, get a list containing the triples
            // (singleton list of alt allele, variant type, set of labels)
            // corresponding to alt alleles that pass variant-type and resource-match checks
            return vc.getAlternateAlleles().stream()
                    .filter(a -> !GATKVCFConstants.isSpanningDeletion(a))
                    .filter(a -> variantTypesToExtract.contains(VariantType.getAlleleSpecificVariantType(vc, a)))
                    .map(a -> Triple.of(Collections.singletonList(a), VariantType.getAlleleSpecificVariantType(vc, a),
                            findMatchingResourceLabels(vc, a, featureContext)))
                    .filter(t -> isExtractUnlabeled || !t.getRight().isEmpty())
                    .collect(Collectors.toList());
        }
        // if variant-type and resource-match checks failed, return an empty list
        return Collections.emptyList();
    }

    /**
     * @param altAllele     {@code null} if non-allele-specific mode ({@code useASAnnotations} is false)
     */
    private TreeSet<String> findMatchingResourceLabels(final VariantContext vc,
                                                       final Allele altAllele,
                                                       final FeatureContext featureContext) {
        final TreeSet<String> matchingResourceLabels = new TreeSet<>();
        for (final FeatureInput<VariantContext> resource : resources) {
            final List<VariantContext> resourceVCs = featureContext.getValues(resource, featureContext.getInterval().getStart());
            for (final VariantContext resourceVC : resourceVCs) {
                // we should have set resourceMatchingStrategy = ResourceMatchingStrategy.START_POSITION_AND_MINIMAL_REPRESENTATION if useASAnnotations is true
                if (isMatchingVariant(vc, resourceVC, altAllele, !doNotTrustAllPolymorphic, resourceMatchingStrategy)) {
                    resource.getTagAttributes().entrySet().stream()
                            .filter(e -> e.getValue().equals("true"))
                            .map(Map.Entry::getKey)
                            .forEach(matchingResourceLabels::add);
                }
            }
        }
        return matchingResourceLabels;
    }

    /**
     * @param altAllele     {@code null} if non-allele-specific mode ({@code useASAnnotations} is false)
     */
    private static boolean isMatchingVariant(final VariantContext vc,
                                             final VariantContext resourceVC,
                                             final Allele altAllele,
                                             final boolean trustAllPolymorphic,
                                             final ResourceMatchingStrategy resourceMatchingStrategy) {
        if (resourceVC != null && resourceVC.isNotFiltered() && resourceVC.isVariant() && VariantType.checkVariantType(vc, resourceVC) &&
                (trustAllPolymorphic || !resourceVC.hasGenotypes() || resourceVC.isPolymorphicInSamples())) { // this is the check originally performed by VQSR
            switch (resourceMatchingStrategy) {
                case START_POSITION:
                    return true;
                case START_POSITION_AND_GIVEN_REPRESENTATION:
                    // we further require that at least one alt allele is present in the resource alt alleles, but don't reconcile representations
                    return !Sets.intersection(Sets.newHashSet(vc.getAlternateAlleles()), Sets.newHashSet(resourceVC.getAlternateAlleles())).isEmpty();
                case START_POSITION_AND_MINIMAL_REPRESENTATION:
                    // we further require that at least one alt allele is present in the resource alt alleles, and do reconcile representations
                    try {
                        if (altAllele == null) {
                            // non-allele-specific mode
                            return vc.getAlternateAlleles().stream()
                                    .anyMatch(alt -> GATKVariantContextUtils.isAlleleInList(vc.getReference(), alt, resourceVC.getReference(), resourceVC.getAlternateAlleles()));
                        }
                        // allele-specific mode
                        return GATKVariantContextUtils.isAlleleInList(vc.getReference(), altAllele, resourceVC.getReference(), resourceVC.getAlternateAlleles());
                    } catch (final IllegalStateException e) {
                        throw new IllegalStateException("Reference allele mismatch at position " + resourceVC.getContig() + ':' + resourceVC.getStart() + ": ", e);
                    }
                default:
                    throw new GATKException.ShouldNeverReachHereException("Unknown ResourceMatchingStrategy.");
            }
        }
        return false;
    }
}