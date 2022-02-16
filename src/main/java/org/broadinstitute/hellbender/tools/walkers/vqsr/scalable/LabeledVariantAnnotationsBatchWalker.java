package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
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
 * TODO
 */
@CommandLineProgramProperties(
        // TODO
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public class LabeledVariantAnnotationsBatchWalker extends MultiVariantWalker {

    private static final String ANNOTATIONS_HDF5_SUFFIX = ".annot.hdf5";

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output prefix.")
    String outputPrefix;

    /**
     * Any set of VCF files to use as lists of training or truth sites.
     * Training - The program builds the model using input variants that overlap with these training sites.
     * Truth - The program uses these truth sites to determine where to set the cutoff in sensitivity.
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
     * Extract per-allele annotations.
     * Annotations should be specified using their full names with AS_ prefix.
     * Non-allele-specific annotations will be applied to all alleles.
     */
    @Argument(
            fullName = "use-allele-specific-annotations",
            doc = "If specified, attempt to use the allele-specific versions of the specified annotations.",
            optional = true)
    private boolean useASAnnotations = false;

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(
            fullName = StandardArgumentDefinitions.ANNOTATION_LONG_NAME,
            shortName = "A",
            doc = "The names of the annotations to extract.",
            minElements = 1)
    private List<String> annotationNames = new ArrayList<>();

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

    // TODO document and validate batchSize * number of annotations <= HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX
    @Argument(
            fullName = "batch-size",
            minValue = 1,
            optional = true
    )
    private int batchSize = 100000;

    // TODO document and validate batchSize * number of annotations <= HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX
    @Advanced
    @Argument(
            fullName = "omit-alleles-in-hdf5"
    )
    private boolean omitAllelesInHDF5 = false;

    Set<VariantType> variantTypesToExtract;
    File outputAnnotationsFile;
    private final Set<String> ignoreInputFilterSet = new TreeSet<>();

    LabeledVariantAnnotationsData dataBatch;

    int batchIndex = 0;

    // TODO document, make enum (extract labeled vs. extract all)
    public boolean isExtractVariantsNotOverlappingResources() {
        return true;
    }

    public void beforeOnTraversalStart() {
        // override
    }

    @Override
    public void onTraversalStart() {

        beforeOnTraversalStart();

        // TODO validate annotation names and AS mode

        variantTypesToExtract = EnumSet.copyOf(variantTypesToExtractList);

        outputAnnotationsFile = new File(outputPrefix + ANNOTATIONS_HDF5_SUFFIX);

        for (final File outputFile : Collections.singletonList(outputAnnotationsFile)) {
            if ((outputFile.exists() && !outputFile.canWrite()) ||
                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
            }
        }

        if (ignoreInputFilters != null) {
            ignoreInputFilterSet.addAll(ignoreInputFilters);
        }

        final Set<String> resourceLabels = new TreeSet<>();
        for (final FeatureInput<VariantContext> resource : resources) {
            final TreeSet<String> trackResourceLabels = resource.getTagAttributes().entrySet().stream()
                    .filter(e -> e.getValue().equals("true"))
                    .map(Map.Entry::getKey)
                    .sorted()
                    .collect(Collectors.toCollection(TreeSet::new));
            resourceLabels.addAll(trackResourceLabels);
            logger.info( String.format("Found %s track: labels = %s", resource.getName(), trackResourceLabels));
        }
        resourceLabels.forEach(String::intern);

//        if (!resourceLabels.contains(LabeledVariantAnnotationsData.TRAINING_LABEL)) {
//            throw new CommandLineException(
//                    "No training set found! Please provide sets of known polymorphic loci marked with the training=true feature input tag. " +
//                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
//        }
//
//        if (!resourceLabels.contains(LabeledVariantAnnotationsData.TRUTH_LABEL)) {
//            throw new CommandLineException(
//                    "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true feature input tag. " +
//                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
//        }

        dataBatch = new LabeledVariantAnnotationsData(annotationNames, resourceLabels, batchSize, useASAnnotations);
    }

    @Override
    public void apply(final VariantContext vc,
                      final ReadsContext readsContext,
                      final ReferenceContext ref,
                      final FeatureContext featureContext) {
        addVariantToBatchIfItPassesExtractionChecks(vc, featureContext);

        if (dataBatch.size() == batchSize) {
            consumeBatch();
        }
    }

    void consumeBatch() {
        dataBatch.writeBatchToHDF5(outputAnnotationsFile, batchIndex, omitAllelesInHDF5);
        doBatchWork();
        dataBatch.clear();
        batchIndex++;
    }

    void doBatchWork() {
        // override
    }

    @Override
    public Object onTraversalSuccess() {

        // final batch
        if (dataBatch.size() > 0) {
            consumeBatch();
        }

        afterOnTraversalSuccess();

        return null;
    }

    public void afterOnTraversalSuccess() {
        // override
    }

    // code here and below for filtering and determining variant type was essentially retained from VQSR
    private void addVariantToBatchIfItPassesExtractionChecks(final VariantContext vc,
                                                             final FeatureContext featureContext) {
        // TODO dump unneeded VariantContext info if in Extract
        if (vc == null) {
            return;
        }
        if (!(ignoreAllFilters || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()))) {
            return;
        }
        if (!useASAnnotations) {
            final VariantType variantType = VariantType.getVariantType(vc);
            if (variantTypesToExtract.contains(variantType)) {
                final Set<String> overlappingResourceLabels = findOverlappingResourceLabels(vc, null, null, featureContext);
                if (isExtractVariantsNotOverlappingResources() || !overlappingResourceLabels.isEmpty()) {
                    dataBatch.add(vc,
                            vc.getAlternateAlleles(),
                            Collections.singletonList(variantType),
                            Collections.singletonList(overlappingResourceLabels));
                }
            }
        } else {
            final int numAltAlleles = vc.getAlternateAlleles().size();
            final List<Allele> extractedAltAllelePerDatum = new ArrayList<>(numAltAlleles);
            final List<VariantType> extractedVariantTypePerDatum = new ArrayList<>(numAltAlleles);
            final List<Set<String>> extractedLabelsPerDatum = new ArrayList<>(numAltAlleles);
            for (final Allele altAllele : vc.getAlternateAlleles()) {
                if (GATKVCFConstants.isSpanningDeletion(altAllele)) {
                    continue;
                }
                final VariantType variantType = VariantType.getVariantType(vc, altAllele);
                if (variantTypesToExtract.contains(variantType)) {
                    final Set<String> overlappingResourceLabels = findOverlappingResourceLabels(vc, vc.getReference(), altAllele, featureContext);
                    if (isExtractVariantsNotOverlappingResources() || !overlappingResourceLabels.isEmpty()) {
                        extractedAltAllelePerDatum.add(altAllele);
                        extractedVariantTypePerDatum.add(variantType);
                        extractedLabelsPerDatum.add(overlappingResourceLabels);
                    }
                }
            }
            dataBatch.add(vc,
                    extractedAltAllelePerDatum,
                    extractedVariantTypePerDatum,
                    extractedLabelsPerDatum);
        }
    }

    private Set<String> findOverlappingResourceLabels(final VariantContext vc,
                                                      final Allele refAllele,
                                                      final Allele altAllele,
                                                      final FeatureContext featureContext) {
        final Set<String> overlappingResourceLabels = new TreeSet<>();
        for (final FeatureInput<VariantContext> resource : resources) {
            final List<VariantContext> resourceVCs = featureContext.getValues(resource, featureContext.getInterval().getStart());
            for (final VariantContext resourceVC : resourceVCs) {
                if (useASAnnotations && !doAllelesMatch(refAllele, altAllele, resourceVC)) {
                    continue;
                }
                if (isValidVariant(vc, resourceVC, trustAllPolymorphic)) {
                    overlappingResourceLabels.addAll(resource.getTagAttributes().entrySet().stream()
                            .filter(e -> e.getValue().equals("true"))
                            .map(Map.Entry::getKey)
                            .collect(Collectors.toSet()));
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