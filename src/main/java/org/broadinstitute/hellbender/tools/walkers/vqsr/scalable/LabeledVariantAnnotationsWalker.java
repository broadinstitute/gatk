package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiplePassVariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;
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
 *      - if variant/alleles pass filters and variant-type/overlapping-resource checks, then:
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
 *   - a VCF file, containing the input variants, with labels and scores appended for those passing variant-type checks TODO + truth-sensitivity scores + filters applied?
 */
@CommandLineProgramProperties(
        // TODO
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public abstract class LabeledVariantAnnotationsWalker extends MultiplePassVariantWalker {

    static final String ANNOTATIONS_HDF5_SUFFIX = ".annot.hdf5";
    private static final String VCF_SUFFIX = ".vcf.gz";

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output prefix.")
    String outputPrefix;

    /**
     * TODO
     */
    @Argument(
            fullName = StandardArgumentDefinitions.RESOURCE_LONG_NAME,
            doc = "") // TODO
    private List<FeatureInput<VariantContext>> resources = new ArrayList<>(10);

    @Argument(
            fullName = "mode",
            doc = "Variant types to extract")
    private List<VariantType> variantTypesToExtractList = new ArrayList<>(Arrays.asList(VariantType.SNP, VariantType.INDEL));

    /**
     * TODO fix up
     * Extract per-allele annotations.
     * Annotations should be specified using their full names with AS_ prefix.
     * Non-allele-specific annotations will be applied to all alleles.
     */
    @Argument(
            fullName = "use-allele-specific-annotations",
            doc = "If specified, attempt to use the allele-specific versions of the specified annotations.",
            optional = true)
    boolean useASAnnotations = false;

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(
            fullName = StandardArgumentDefinitions.ANNOTATION_LONG_NAME,
            shortName = "A",
            doc = "The names of the annotations to extract.",
            minElements = 1)
    List<String> annotationNames = new ArrayList<>();

    @Argument(
            fullName = "ignore-filter",
            doc = "If specified, use variants marked as filtered by the specified filter name in the input VCF file.",
            optional = true)
    private List<String> ignoreInputFilters = new ArrayList<>();

    @Argument(
            fullName = "ignore-all-filters",
            doc = "If specified, ignore all input filters.",
            optional = true)
    private boolean ignoreAllFilters = false;

    @Advanced
    @Argument(
            fullName = "trust-all-polymorphic",
            doc = "Trust that unfiltered records in the resources contain only polymorphic sites to decrease runtime.",
            optional = true)
    private boolean trustAllPolymorphic = false;

    @Advanced
    @Argument(
            fullName = "omit-alleles-in-hdf5"
    )
    boolean omitAllelesInHDF5 = false;

    private final Set<String> ignoreInputFilterSet = new TreeSet<>();
    Set<VariantType> variantTypesToExtract;
    TreeSet<String> resourceLabels = new TreeSet<>();

    File outputAnnotationsFile;
    VariantContextWriter vcfWriter;

    LabeledVariantAnnotationsData data;

    @Override
    public void onTraversalStart() {

        // TODO generally clean up validation
        // TODO validate annotation names and AS mode

        ignoreInputFilterSet.addAll(ignoreInputFilters);

        variantTypesToExtract = EnumSet.copyOf(variantTypesToExtractList);

        outputAnnotationsFile = new File(outputPrefix + ANNOTATIONS_HDF5_SUFFIX);
        final File outputVCFFile = new File(outputPrefix + VCF_SUFFIX);

        for (final File outputFile : Arrays.asList(outputAnnotationsFile, outputVCFFile)) {
            if ((outputFile.exists() && !outputFile.canWrite()) ||
                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
            }
        }

        // TODO extract some of this here and in findOverlappingResourceLabels
        for (final FeatureInput<VariantContext> resource : resources) {
            final TreeSet<String> trackResourceLabels = resource.getTagAttributes().entrySet().stream()
                    .filter(e -> e.getValue().equals("true"))
                    .map(Map.Entry::getKey)
                    .sorted()
                    .collect(Collectors.toCollection(TreeSet::new));
            resourceLabels.addAll(trackResourceLabels);
            logger.info( String.format("Found %s track: labels = %s", resource.getName(), trackResourceLabels));
        }
        resourceLabels.forEach(String::intern); // TODO not sure if this is necessary?

        // TODO let child tools specify required labels
        if (!resourceLabels.contains(LabeledVariantAnnotationsData.TRAINING_LABEL)) {
            throw new CommandLineException(
                    "No training set found! Please provide sets of known polymorphic loci marked with the training=true feature input tag. " +
                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
        }

        if (!resourceLabels.contains(LabeledVariantAnnotationsData.TRUTH_LABEL)) {
            throw new CommandLineException(
                    "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true feature input tag. " +
                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
        }

        data = new LabeledVariantAnnotationsData(annotationNames, resourceLabels, useASAnnotations);

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

    // TODO maybe clean up all this Triple and metadata business with a class?
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

    void writeAnnotationsToHDF5AndClearData() {
        if (data.size() == 0) {
            throw new GATKException("None of the specified input variants were present in the resource VCFs.");
        }
        for (final VariantType variantType : variantTypesToExtract) {
            logger.info(String.format("Extracted annotations for %d variants of type %s.",
                    data.getVariantTypeFlat().stream().mapToInt(t -> t == variantType ? 1 : 0).sum(), variantType));
        }
        for (final String label : data.getSortedLabels()) {
            logger.info(String.format("Extracted annotations for %d variants labeled as %s.",
                    data.isLabelFlat(label).stream().mapToInt(b -> b ? 1 : 0).sum(), label));
        }
        logger.info(String.format("Extracted annotations for %s total variants.", data.size()));

        logger.info("Writing annotations...");
        data.writeHDF5(outputAnnotationsFile, omitAllelesInHDF5);
        logger.info(String.format("Annotations and metadata written to %s.", outputAnnotationsFile.getAbsolutePath()));

        data.clear();
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
        Set<VCFHeaderLine> hInfo = getDefaultToolVCFHeaderLines();
        // TODO extract
        hInfo.addAll(sortedLabels.stream()
                .map(l -> new VCFInfoHeaderLine(l, 1, VCFHeaderLineType.Flag, String.format("This site was labeled as %s according to resources", l)))
                .collect(Collectors.toList()));
        hInfo.add(GATKVCFHeaderLines.getFilterLine(VCFConstants.PASSES_FILTERS_v4));
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        hInfo = VcfUtils.updateHeaderContigLines(hInfo, null, sequenceDictionary, true);
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
            // (list of alt alleles passing variant-type and overlapping-resource checks, variant type, set of labels)
            final VariantType variantType = VariantType.getVariantType(vc);
            if (variantTypesToExtract.contains(variantType)) {
                final TreeSet<String> overlappingResourceLabels = findOverlappingResourceLabels(vc, null, null, featureContext);
                if (isExtractUnlabeled || !overlappingResourceLabels.isEmpty()) {
                    return Collections.singletonList(Triple.of(vc.getAlternateAlleles(), variantType, overlappingResourceLabels));
                }
            }
        } else {
            // in allele-specific mode, get a list containing the triples
            // (singleton list of alt allele, variant type, set of labels)
            // corresponding to alt alleles that pass variant-type and overlapping-resource checks
            return vc.getAlternateAlleles().stream()
                    .filter(a -> !GATKVCFConstants.isSpanningDeletion(a))
                    .filter(a -> variantTypesToExtract.contains(VariantType.getVariantType(vc, a)))
                    .map(a -> Triple.of(Collections.singletonList(a), VariantType.getVariantType(vc, a),
                            findOverlappingResourceLabels(vc, vc.getReference(), a, featureContext)))
                    .filter(t -> isExtractUnlabeled || !t.getRight().isEmpty())
                    .collect(Collectors.toList());
        }
        // if variant-type and overlapping-resource checks failed, return an empty list
        return Collections.emptyList();
    }

    private TreeSet<String> findOverlappingResourceLabels(final VariantContext vc,
                                                          final Allele refAllele,
                                                          final Allele altAllele,
                                                          final FeatureContext featureContext) {
        final TreeSet<String> overlappingResourceLabels = new TreeSet<>();
        for (final FeatureInput<VariantContext> resource : resources) {
            final List<VariantContext> resourceVCs = featureContext.getValues(resource, featureContext.getInterval().getStart());
            for (final VariantContext resourceVC : resourceVCs) {
                if (useASAnnotations && !doAllelesMatch(refAllele, altAllele, resourceVC)) {
                    continue;
                }
                if (isValidVariant(vc, resourceVC, trustAllPolymorphic)) {
                    resource.getTagAttributes().entrySet().stream()
                            .filter(e -> e.getValue().equals("true"))
                            .map(Map.Entry::getKey)
                            .forEach(overlappingResourceLabels::add);
                }
            }
        }
        return overlappingResourceLabels;
    }

    private static boolean isValidVariant(final VariantContext vc,
                                          final VariantContext resourceVC,
                                          final boolean trustAllPolymorphic) {
        return resourceVC != null && resourceVC.isNotFiltered() && resourceVC.isVariant() && VariantType.checkVariantType(vc, resourceVC) &&
                (trustAllPolymorphic || !resourceVC.hasGenotypes() || resourceVC.isPolymorphicInSamples());
    }

    private static boolean doAllelesMatch(final Allele refAllele,
                                          final Allele altAllele,
                                          final VariantContext resourceVC) {
        if (altAllele == null) {
            return true;
        }
        try {
            return GATKVariantContextUtils.isAlleleInList(refAllele, altAllele, resourceVC.getReference(), resourceVC.getAlternateAlleles());
        } catch (final IllegalStateException e) {
            throw new IllegalStateException("Reference allele mismatch at position " + resourceVC.getContig() + ':' + resourceVC.getStart() + " : ", e);
        }
    }
}