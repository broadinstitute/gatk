package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

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
public class VariantAnnotationWalker extends MultiVariantWalker {

    /**
     * Any set of VCF files to use as lists of training or truth sites.
     * Training - The program builds the model using input variants that overlap with these training sites.
     * Truth - The program uses these truth sites to determine where to set the cutoff in sensitivity.
     */
    @Argument(
            fullName = StandardArgumentDefinitions.RESOURCE_LONG_NAME,
            doc = "") // TODO
    private List<FeatureInput<VariantContext>> resource = new ArrayList<>();

    @Argument(
            fullName = "mode",
            shortName = "mode",
            doc = "Variant type to extract")
    public VariantTypeMode mode = VariantTypeMode.SNP;

    /**
     * Extract per-allele annotations.
     * Annotations should be specified using their full names with AS_ prefix.
     * Non-allele-specific (scalar) annotations will be applied to all alleles.
     */
    @Argument(
            fullName = "use-allele-specific-annotations",
            shortName = "AS",
            doc = "If specified, attempt to use the allele-specific versions of the specified annotations.",
            optional = true)
    public boolean useASannotations = false;

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(
            fullName = "use-annotation",
            shortName = "an",
            doc = "The names of the annotations to extract.",
            minElements = 1)
    private List<String> useAnnotations = new ArrayList<>();

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

    VariantDataManager dataManager;
    boolean isExtractTrainingAndTruthOnly;
    private final Set<String> ignoreInputFilterSet = new TreeSet<>();
    private final List<ImmutablePair<VariantContext, FeatureContext>> variantsAtLocus = new ArrayList<>(10);

    public void beforeOnTraversalStart() {
        // override
    }

    @Override
    public void onTraversalStart() {

        beforeOnTraversalStart();

        dataManager = new VariantDataManager(new ArrayList<>(useAnnotations), useASannotations, trustAllPolymorphic);

        if (ignoreInputFilters != null) {
            ignoreInputFilterSet.addAll(ignoreInputFilters);
        }

        for (final FeatureInput<VariantContext> variantFeature : resource) {
            dataManager.addVariantSet(new VariantSet(variantFeature));
        }

        if (!dataManager.checkHasTrainingSet()) {
            throw new CommandLineException(
                    "No training set found! Please provide sets of known polymorphic loci marked with the training=true feature input tag. " +
                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
        }

        if (!dataManager.checkHasTruthSet()) {
            throw new CommandLineException(
                    "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true feature input tag. " +
                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
        }
    }

    @Override
    public void apply(final VariantContext vc,
                      final ReadsContext readsContext,
                      final ReferenceContext ref,
                      final FeatureContext featureContext) {
        // Queue up all variant/featureContext pairs that share a start locus, and defer
        // processing until the start-position or contig changes. In practice this will
        // rarely queue up more than a single variant at a time.
        if (variantLocusChanged(vc)) {
            consumeQueuedVariants();
        }
        variantsAtLocus.add(new ImmutablePair<>(vc, featureContext));
    }

    // Check to see if the start locus for this variant is different from the
    // ones in the queue.
    private boolean variantLocusChanged(final VariantContext nextVC) {
        if (variantsAtLocus.isEmpty()) {
            return false;
        } else {
            final ImmutablePair<VariantContext, FeatureContext> previous = variantsAtLocus.get(0);
            return (nextVC.getStart() != previous.left.getStart() || !nextVC.getContig().equals(previous.left.getContig()));
        }
    }

    private void consumeQueuedVariants() {
        variantsAtLocus.forEach(v -> addVariantDatum(v.left, v.right));
        variantsAtLocus.clear();
    }

    // TODO check presence in training/truth sets outside of addDatum method
    private void addVariantDatum(final VariantContext vc,
                                 final FeatureContext context) {
        if (vc != null && (ignoreAllFilters || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()))) {
            if (VariantDataManager.checkVariationClass(vc, mode) && !useASannotations) {
                dataManager.addDatum(context, vc, null, null, isExtractTrainingAndTruthOnly);
            } else if (useASannotations) {
                for (final Allele allele : vc.getAlternateAlleles()) {
                    if (!GATKVCFConstants.isSpanningDeletion(allele) && VariantDataManager.checkVariationClass(vc, allele, mode)) {
                        //note that this may not be the minimal representation for the ref and alt allele
                        dataManager.addDatum(context, vc, vc.getReference(), allele, isExtractTrainingAndTruthOnly);
                    }
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {

        consumeQueuedVariants(); // finish processing any queued variants

        afterTraversalSuccess();

        return true;
    }

    public void afterTraversalSuccess() {
        // override
    }

    void writeAnnotationsHDF5(final File file) {
        try (final HDF5File hdf5File = new HDF5File(file, HDF5File.OpenMode.CREATE)) { // TODO allow appending
            IOUtils.canReadFile(hdf5File.getFile());

            hdf5File.makeStringArray("/data/annotation_names", dataManager.getAnnotationKeys().toArray(new String[0]));
            hdf5File.makeDoubleMatrix("/data/annotations", dataManager.getData().stream().map(vd -> vd.annotations).toArray(double[][]::new));
            hdf5File.makeDoubleArray("/data/is_transition", dataManager.getData().stream().mapToDouble(vd -> vd.isTransition ? 1 : 0).toArray());
            hdf5File.makeDoubleArray("/data/is_training", dataManager.getData().stream().mapToDouble(vd -> vd.atTrainingSite ? 1 : 0).toArray());
            hdf5File.makeDoubleArray("/data/is_truth", dataManager.getData().stream().mapToDouble(vd -> vd.atTruthSite ? 1 : 0).toArray());
        } catch (final RuntimeException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                    exception, file.getAbsolutePath()));
        }
    }
}