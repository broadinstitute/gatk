package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * TODO
 */
@CommandLineProgramProperties(
        summary = "Build a recalibration model to score variant quality for filtering purposes",
        oneLineSummary = "Build a recalibration model to score variant quality for filtering purposes",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public class ScikitLearnVariantTrain extends MultiVariantWalker {

    @ArgumentCollection
    final private ScikitLearnVariantTrainArgumentCollection VTAC = new ScikitLearnVariantTrainArgumentCollection();

    /////////////////////////////
    // Inputs
    /////////////////////////////

    /**
     * Any set of VCF files to use as lists of training, truth, or known sites.
     * Training - The program builds the Gaussian mixture model using input variants that overlap with these training sites.
     * Truth - The program uses these truth sites to determine where to set the cutoff in VQSLOD sensitivity.
     * Known - The program only uses known sites for reporting purposes (to indicate whether variants are already known or
     * novel). They are not used in any calculations by the algorithm itself.
     * Bad - A database of known bad variants can be used to supplement the set of worst ranked variants (compared to the
     * Gaussian mixture model) that the program selects from the data to model "bad" variants.
     */
    @Argument(fullName=StandardArgumentDefinitions.RESOURCE_LONG_NAME,
            doc="A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth sets are required to run)",
            optional=false)
    private List<FeatureInput<VariantContext>> resource = new ArrayList<>();

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output recal file used by ApplyVQSR", optional=false)
    private GATKPath output;

    @Argument(fullName="tranches-file", doc="The output tranches file used by ApplyVQSR", optional=false)
    private File TRANCHES_FILE; // not GATKPath since this name must be accessible to R code

    /////////////////////////////
    // Additional Command Line Arguments
    /////////////////////////////

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(fullName="use-annotation",
            shortName="an",
            doc="The names of the annotations which should used for calculations",
            optional=false)
    private List<String> USE_ANNOTATIONS = new ArrayList<>();

    /**
     * Add truth sensitivity slices through the call set at the given values. The default values are 100.0, 99.9, 99.0, and 90.0
     * which will result in 4 estimated tranches in the final call set: the full set of calls (100% sensitivity at the accessible
     * sites in the truth set), a 99.9% truth sensitivity tranche, along with progressively smaller tranches at 99% and 90%.
     * Note: You must pass in each tranche as a separate value (e.g. -tranche 100.0 -tranche 99.9).
     */
    @Argument(fullName="truth-sensitivity-tranche",
            shortName="tranche",
            doc="The levels of truth sensitivity at which to slice the data. (in percent, that is 1.0 for 1 percent)",
            optional=true)
    private List<Double> TS_TRANCHES = new ArrayList<Double>(Arrays.asList(100.0, 99.9, 99.0, 90.0));

    /**
     * For this to work properly, the --ignore-filter argument should also be applied to the ApplyVQSR command.
     */
    @Argument(fullName="ignore-filter",
            doc="If specified, the variant recalibrator will also use variants marked as filtered by the specified filter name in the input VCF file",
            optional=true)
    private List<String> IGNORE_INPUT_FILTERS = new ArrayList<>();

    @Argument(fullName="ignore-all-filters",
            doc="If specified, the variant recalibrator will ignore all input filters. Useful to rerun the VQSR from a filtered output file.",
            optional=true)
    private boolean IGNORE_ALL_FILTERS = false;

    /**
     *  TODO
     */
    @Argument(fullName="output-model",
            doc="If specified, the variant recalibrator will output the VQSR model to this file path.",
            optional=true)
    private GATKPath outputModelPath = null;

    @Argument(fullName="output-scores")
    private File outputScoresFile = null;

    @Argument(fullName="python-script")
    private File pythonScriptFile = null;

    /**
     * This argument is intended to be used in a more complicated VQSR scheme meant for very large WGS callsets that
     * require a prohibitive amount of memory for classic VQSR. Given that training data is downsampled once it exceeds
     * --max-num-training-data, reading in additional data to build the model only serves to consume resources. However,
     * with this argument the output recal file will also be downsampled. The recommended VQSR procedure when using this
     * argument is to run VariantRecalibrator once with sampling and designate an --output-model file. Then
     * VariantRecalibrator can be run a second time scattered using the -scatterTranches argument and that file as an
     * --input-model.  The scattered recal files can be gathered with the GatherVcfs tool and the scattered tranches can
     * be gathered with the GatherTranches tool.
     *
     */
    @Argument(fullName="sample-every-Nth-variant",
            shortName = "sample-every",
            doc="If specified, the variant recalibrator will use (and output) only a subset of variants consisting of every Nth variant where N is specified by this argument; for use with --output-model -- see argument details",
            optional=true)
    @Hidden
    private int sampleMod = 1;

    /**
     * Add VQSLOD slices through the call set at the given values, to be used with the -scatterTranches argument. The
     * default values span from -10 to +10 at varying resolution. The resolution output here will affect the accuracy
     * with which the gathered tranches are able to match the requested truth sensitivity levels.
     *
     * Alternative lists of tranche values should only be used for testing.
     */
    @Hidden
    @Argument(fullName="vqslod-tranche",
            doc="The levels of VQSLOD at which to slice the data.",
            optional=true)
    private List<Double> VQSLOD_TRANCHES = new ArrayList<>(1000);
    {
        for (double i=10.0; i>5; i-=0.1) {
            VQSLOD_TRANCHES.add(i);
        }
        for (double i=5.0; i>-5; i-=0.01) {
            VQSLOD_TRANCHES.add(i);
        }
        for (double i=-5.0; i>-10; i-=0.1) {
            VQSLOD_TRANCHES.add(i);
        }
    };

    /**
     * The statistical model being built by this tool may fail due to simple statistical sampling
     * issues. Rather than dying immediately when the initial model fails, this argument allows the
     * tool to restart with a different random seed and try to build the model again. The first
     * successfully built model will be kept.
     *
     * Note that the most common underlying cause of model building failure is that there is insufficient data to
     * build a really robust model. This argument provides a workaround for that issue but it is
     * preferable to provide this tool with more data (typically by including more samples or more territory)
     * in order to generate a more robust model.
     */
    @Advanced
    @Argument(fullName="max-attempts",
            doc="Number of attempts to build a model before failing",
            optional=true)
    @VisibleForTesting
    protected int max_attempts = 1;

    /////////////////////////////
    // Debug Arguments
    /////////////////////////////
    @Advanced
    @Argument(fullName = "trust-all-polymorphic",
            doc = "Trust that all the input training sets' unfiltered records contain only polymorphic sites to drastically speed up the computation.",
            optional=true)
    private boolean TRUST_ALL_POLYMORPHIC = false;

    @VisibleForTesting
    protected List<Integer> annotationOrder = null;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private ScikitLearnVariantTrainDataManager dataManager;
    private VariantContextWriter recalWriter;
    private PrintStream tranchesStream;
    private final Set<String> ignoreInputFilterSet = new TreeSet<>();
    final private ArrayList<VariantDatum> reduceSum = new ArrayList<>(2000);
    final private List<ImmutablePair<VariantContext, FeatureContext>> variantsAtLocus = new ArrayList<>();
    private long counter = 0;

    //---------------------------------------------------------------------------------------------------------------
    //
    // onTraversalStart
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void onTraversalStart() {

        PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");

        IOUtils.canReadFile(pythonScriptFile);
        dataManager = new ScikitLearnVariantTrainDataManager( new ArrayList<>(USE_ANNOTATIONS), VTAC);

        if ( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( IGNORE_INPUT_FILTERS );
        }

        try {
            tranchesStream = new PrintStream(TRANCHES_FILE);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(TRANCHES_FILE, e);
        }

        for (final FeatureInput<VariantContext> variantFeature : resource ) {
            dataManager.addTrainingSet( new TrainingSet( variantFeature ) );
        }

        if ( !dataManager.checkHasTrainingSet() ) {
            throw new CommandLineException(
                    "No training set found! Please provide sets of known polymorphic loci marked with the training=true feature input tag. For example, -resource hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }
        if ( !dataManager.checkHasTruthSet() ) {
            throw new CommandLineException(
                    "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true feature input tag. For example, -resource hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }

        //TODO: this should be refactored/consolidated as part of
        // https://github.com/broadinstitute/gatk/issues/2112
        // https://github.com/broadinstitute/gatk/issues/121,
        // https://github.com/broadinstitute/gatk/issues/1116 and
        // Initialize VCF header lines
        Set<VCFHeaderLine> hInfo = getDefaultToolVCFHeaderLines();
        VariantRecalibrationUtils.addVQSRStandardHeaderLines(hInfo);
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (hasReference()) {
            hInfo = VcfUtils.updateHeaderContigLines(
                    hInfo, referenceArguments.getReferencePath(), sequenceDictionary, true);
        }
        else if (null != sequenceDictionary) {
            hInfo = VcfUtils.updateHeaderContigLines(hInfo, null, sequenceDictionary, true);
        }

        recalWriter = createVCFWriter(output);
        recalWriter.writeHeader( new VCFHeader(hInfo) );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // apply
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        // Queue up all variant/featureContext pairs that share a start locus, and defer
        // processing until the start-position or contig changes. In practice this will
        // rarely queue up more than a single variant at a time.
        if (variantLocusChanged(vc)) {
            consumeQueuedVariants();
        }
        if (counter % sampleMod == 0)
            variantsAtLocus.add(new ImmutablePair<>(vc, featureContext));
        counter++;

    }

    // Check to see if the start locus for this variant is different from the
    // ones in the queue.
    private boolean variantLocusChanged(final VariantContext nextVC) {
        if (variantsAtLocus.isEmpty()) {
            return false;
        }
        else {
            final ImmutablePair<VariantContext, FeatureContext> previous = variantsAtLocus.get(0);
            return (nextVC.getStart() != previous.left.getStart() || !nextVC.getContig().equals(previous.left.getContig()));
        }
    }

    private void consumeQueuedVariants() {
        variantsAtLocus.forEach(v -> addVariantDatum(v.left, v.right));
        variantsAtLocus.clear();
    }

    private void addVariantDatum(final VariantContext vc, final FeatureContext context) {
        if( vc != null && ( IGNORE_ALL_FILTERS || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()) ) ) {
            if( GMMVariantTrainDataManager.checkVariationClass( vc, VTAC.MODE ) && !VTAC.useASannotations) {
                addDatum(reduceSum, true, context, vc, null, null);
            }
            else if( VTAC.useASannotations ) {
                for (final Allele allele : vc.getAlternateAlleles()) {
                    if (!GATKVCFConstants.isSpanningDeletion(allele) && GMMVariantTrainDataManager.checkVariationClass(vc, allele, VTAC.MODE)) {
                        //note that this may not be the minimal representation for the ref and alt allele
                        addDatum(reduceSum, true, context, vc, vc.getReference(), allele);
                    }
                }
            }
        }
    }

    /**
     * add a datum representing a variant site (or allele) to the data in {@code variants}, which represents the callset to be recalibrated
     * @param variants is modified by having a new VariantDatum added to it
     */
    private void addDatum(
            final ArrayList<VariantDatum> variants,
            final boolean isInput,
            final FeatureContext featureContext,
            final VariantContext vc,
            final Allele refAllele,
            final Allele altAllele) {
        final VariantDatum datum = new VariantDatum();

        // Populate the datum with lots of fields from the VariantContext, unfortunately the VC is too big so we just
        // pull in only the things we absolutely need.
        datum.referenceAllele = refAllele;
        datum.alternateAllele = altAllele;
        dataManager.decodeAnnotations(datum, vc, true);

        // non-deterministic because order of calls depends on load of machine
        datum.loc = (isInput ? new SimpleInterval(vc) : null);

        datum.originalQual = vc.getPhredScaledQual();
        datum.isSNP = vc.isSNP() && vc.isBiallelic();
        datum.isTransition = datum.isSNP && GATKVariantContextUtils.isTransition(vc);
        datum.isAggregate = !isInput;

        // Loop through the training data sets and if they overlap this locus (and allele, if applicable) then update
        // the prior and training status appropriately. The locus used to find training set variants is retrieved
        // by parseTrainingSets from the FeatureContext argument.
        dataManager.parseTrainingSets(featureContext, vc, datum, TRUST_ALL_POLYMORPHIC);
        final double priorFactor = QualityUtils.qualToProb(datum.prior);
        datum.prior = Math.log10(priorFactor) - Math.log10(1.0 - priorFactor);

        variants.add(datum);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // on traversal success
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Object onTraversalSuccess() {

        consumeQueuedVariants(); // finish processing any queued variants

        dataManager.setData(reduceSum);

        final String rawAnnotationsOutput = output.toString().endsWith(".recal") ? output.toString().split(".recal")[0] : output.toString();
        final File rawAnnotationsFile = new File(rawAnnotationsOutput + ".annot.raw.hdf5");
        writeAnnotationsHDF5(rawAnnotationsFile);

        dataManager.normalizeData(true, annotationOrder); // Each data point is now (x - mean) / standard deviation

        final String annotationsOutput = output.toString().endsWith(".recal") ? output.toString().split(".recal")[0] : output.toString();
        final File annotationsFile = new File(annotationsOutput + ".annot.hdf5");
        writeAnnotationsHDF5(annotationsFile);

        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                new Resource(pythonScriptFile.getAbsolutePath(), ScikitLearnVariantTrain.class),
                null,
                composePythonArguments(rawAnnotationsFile, annotationsFile, VTAC.hyperparametersJSONFile, outputScoresFile));

        if (pythonProcessOutput.getExitValue() != 0) {
            throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
        }

        logger.info(String.format("Scores written to %s.", outputScoresFile.getAbsolutePath()));

        try (final HDF5File outputScoresFileHDF5File = new HDF5File(outputScoresFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(outputScoresFileHDF5File.getFile());
            final double[] scores = outputScoresFileHDF5File.readDoubleArray("/scores");
            dataManager.setScores(dataManager.getData(), scores);
        } catch (final RuntimeException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s: %s",
                    outputScoresFile.getAbsolutePath(), exception));
        }

        // Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
        final int nCallsAtTruth = TrancheManager.countCallsAtTruth(dataManager.getData(), Double.NEGATIVE_INFINITY);
        final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric(nCallsAtTruth);
        final List<? extends Tranche> tranches = TrancheManager.findTranches(dataManager.getData(), TS_TRANCHES, metric, VTAC.MODE);
        tranchesStream.print(TruthSensitivityTranche.printHeader());
        tranchesStream.print(Tranche.tranchesString(tranches));

        logger.info("Writing out recalibration table...");
        dataManager.writeOutRecalibrationTable(recalWriter, getBestAvailableSequenceDictionary());

        return true;
    }

    @Override
    public void closeTool(){
        if (recalWriter != null) {
            recalWriter.close();
        }
        if (tranchesStream != null) {
            tranchesStream.close();
        }
    }

    public void writeAnnotationsHDF5(final File file) {
        try (final HDF5File hdf5File = new HDF5File(file, HDF5File.OpenMode.CREATE)) { // TODO allow appending
            IOUtils.canReadFile(hdf5File.getFile());

            hdf5File.makeStringArray("/data/annotation_names", dataManager.getAnnotationKeys().stream().toArray(String[]::new));
            hdf5File.makeDoubleMatrix("/data/annotations", dataManager.getData().stream().map(vd -> vd.annotations).toArray(double[][]::new));
            hdf5File.makeDoubleArray("/data/is_training", dataManager.getData().stream().mapToDouble(vd -> vd.atTrainingSite ? 1 : 0).toArray());
            hdf5File.makeDoubleArray("/data/is_truth", dataManager.getData().stream().mapToDouble(vd -> vd.atTruthSite ? 1 : 0).toArray());
            hdf5File.makeDoubleArray("/data/is_anti_training", dataManager.getData().stream().mapToDouble(vd -> vd.atAntiTrainingSite ? 1 : 0).toArray());

            logger.info(String.format("Annotations written to %s.", file.getAbsolutePath()));
        } catch (final RuntimeException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                    exception, file.getAbsolutePath()));
        }
    }

    private static List<String> composePythonArguments(final File rawAnnotationsFile,
                                                       final File annotationsFile,
                                                       final File hyperparametersJSONFile,
                                                       final File outputScoresFile) {
        try {
            return new ArrayList<>(Arrays.asList(
                    "--raw_annotations_file=" + rawAnnotationsFile.getCanonicalPath(),
                    "--annotations_file=" + annotationsFile.getCanonicalPath(),
                    "--hyperparameters_json_file=" + hyperparametersJSONFile.getCanonicalPath(),
                    "--output_scores_file=" + outputScoresFile.getCanonicalPath()));
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Encountered exception resolving canonical file paths: %s", e));
        }
    }
}