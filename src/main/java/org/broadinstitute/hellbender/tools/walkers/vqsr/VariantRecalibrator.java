package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.common.annotations.VisibleForTesting;
import Jama.Matrix;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.ExpandingArrayList;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import org.apache.commons.lang3.tuple.ImmutablePair;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Build a recalibration model to score variant quality for filtering purposes
 *
 * <p>
 * Note: this tool only accepts a single input variant file (unlike GATK3, which accepted multiple
 * input variant files).
 * </p>
 *
 * <p>
 * The purpose of variant recalibration is to assign a well-calibrated probability to each variant call in a call set.
 * You can then create highly accurate call sets by filtering based on this single estimate for the accuracy of each call.
 * The approach taken by variant quality score recalibration is to develop a continuous, covarying estimate of the relationship
 * between SNP call annotations (such as QD, MQ, and ReadPosRankSum, for example) and the probability that a SNP is a true genetic
 * variant versus a sequencing or data processing artifact. This model is determined adaptively based on "true sites" provided
 * as input, typically HapMap 3 sites and those sites found to be polymorphic on the Omni 2.5M SNP chip array (in humans). This
 * adaptive error model can then be applied to both known and novel variation discovered in the call set of interest to evaluate
 * the probability that each call is real. The score that gets added to the INFO field of each variant is called the VQSLOD. It is
 * the log odds of being a true variant versus being false under the trained Gaussian mixture model.
 * </p>
 *
 * <p>
 * This tool performs the first pass in a two-stage process called VQSR; the second pass is performed by the
 * <a href='https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php'>ApplyRecalibration</a> tool.
 * In brief, the first pass consists of creating a Gaussian mixture model by looking at the distribution of annotation
 * values over a high quality subset of the input call set, and then scoring all input variants according to the model.
 * The second pass consists of filtering variants based on score cutoffs identified in the first pass.
 *</p>
 *
 * <p>VQSR is probably the hardest part of the Best Practices to get right, so be sure to read the
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=39'>method documentation</a>,
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=1259'>parameter recommendations</a> and
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=2805'>tutorial</a> to really understand what these
 * tools and how to use them for best results on your own data.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 * <li>The input raw variants to be recalibrated. These variant calls must be annotated with the annotations that will be
 * used for modeling. If the calls come from multiple samples, they must have been obtained by joint calling the samples,
 * either directly (running HaplotypeCaller on all samples together) or via the GVCF workflow (HaplotypeCaller with -ERC
 * GVCF per-sample then GenotypeGVCFs on the resulting gVCFs) which is more scalable.</li>
 * <li>Known, truth, and training sets to be used by the algorithm. See the method documentation for more details.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 * <li>A recalibration table file that will be used by the ApplyRecalibration tool.</li>
 * <li>A tranches file which shows various metrics of the recalibration callset for slices of the data.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <p>Recalibrating SNPs in exome data:</p>
 * <pre>
 * ./gatk-launch JAVA_OPTS=-Xmx4g \
 *   VariantRecalibrator \
 *   --variant raw_variants.vcf \
 *   --resource hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
 *   --resource omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
 *   --resource 1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf
 *   --resource dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_135.b37.vcf \
 *   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
 *   -mode SNP \
 *   --recalFile output.recal \
 *   -tranchesFile output.tranches \
 *   --rscriptFile output.plots.R
 * </pre>
 *
 * <h3>Allele-specfic usage</h3>
 * <pre>
 * ./gatk-launch JAVA_OPTS=-Xmx4g \
 *   VariantRecalibrator \
 *   --variant raw_variants.withASannotations.vcf \
 *   -AS \
 *   --resource hapmap,known=false,training=true,truth=true,prior=15.0:hapmap_3.3.b37.sites.vcf \
 *   --resource omni,known=false,training=true,truth=false,prior=12.0:1000G_omni2.5.b37.sites.vcf \
 *   --resource 1000G,known=false,training=true,truth=false,prior=10.0:1000G_phase1.snps.high_confidence.vcf
 *   --resource dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_135.b37.vcf \
 *   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
 *   -mode SNP \
 *   --recalFile output.AS.recal \
 *   -tranchesFile output.AS.tranches \
 *   --rscriptFile output.plots.AS.R
 * </pre>
 * The input VCF must have been produced using allele-specific annotations in HaplotypeCaller.
 * Note that each allele will have a separate line in the output .recal file with its own VQSLOD and culprit that will be
 * transferred to the final VCF in ApplyRecalibration.
 *
 * <h3>Caveats</h3>
 *
 * <ul>
 * <li>SNPs and indels must be recalibrated in separate runs (but it is not necessary to separate them into different
 * files). Mixed records are treated as indels.</li>
 * <li>The values used in the example above are only meant to show how the command lines are composed.
 * They are not meant to be taken as specific recommendations of values to use in your own work, and they may be
 * different from the values cited elsewhere in our documentation. For the latest and greatest recommendations on
 * how to set parameter values for you own analyses, please read the Best Practices section of the documentation,
 * especially the <a href='https://www.broadinstitute.org/gatk/guide/article?id=1259'>FAQ document</a> on VQSR parameters.</li>
 * <li>Whole genomes and exomes take slightly different parameters, so make sure you adapt your commands accordingly! See
 * the documents linked above for details.</li>
 * <li>If you work with small datasets (e.g. targeted capture experiments or small number of exomes), you will run into
 * problems. Read the docs linked above for advice on how to deal with those issues.</li>
 * <li>In order to create the model reporting plots Rscript needs to be in your environment PATH (this is the scripting
 * version of R, not the interactive version). See <a target="r-project" href="http://www.r-project.org">http://www.r-project.org</a>
 * for more info on how to download and install R.</li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = "Build a recalibration model to score variant quality for filtering purposes",
        oneLineSummary = "Build a recalibration model to score variant quality for filtering purposes",
        programGroup = VariantProgramGroup.class
)
public class VariantRecalibrator extends MultiVariantWalker {

    private static final String PLOT_TRANCHES_RSCRIPT = "plot_Tranches.R";

    @ArgumentCollection
    private VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();

    /////////////////////////////
    // Inputs
    /////////////////////////////

    /**
     * These additional calls should be unfiltered and annotated with the error covariates that are intended to be used for modeling.
     */
    @Argument(fullName="aggregate",
            shortName = "aggregate", doc="Additional raw input variants to be used in building the model",
            optional=true)
    private List<FeatureInput<VariantContext>> aggregate = new ArrayList<>();

    /**
     * Any set of VCF files to use as lists of training, truth, or known sites.
     * Training - The program builds the Gaussian mixture model using input variants that overlap with these training sites.
     * Truth - The program uses these truth sites to determine where to set the cutoff in VQSLOD sensitivity.
     * Known - The program only uses known sites for reporting purposes (to indicate whether variants are already known or
     * novel). They are not used in any calculations by the algorithm itself.
     * Bad - A database of known bad variants can be used to supplement the set of worst ranked variants (compared to the
     * Gaussian mixture model) that the program selects from the data to model "bad" variants.
     */
    @Argument(fullName="resource",
            shortName = "resource",
            doc="A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth sets are required to run)",
            optional=false)
    private List<FeatureInput<VariantContext>> resource = new ArrayList<>();

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output recal file used by ApplyRecalibration", optional=false)
    private String output;

    @Argument(fullName="tranches_file", shortName="tranchesFile", doc="The output tranches file used by ApplyRecalibration", optional=false)
    private String TRANCHES_FILE;

    /////////////////////////////
    // Additional Command Line Arguments
    /////////////////////////////
    /**
     * The expected transition / transversion ratio of true novel variants in your targeted region (whole genome, exome, specific
     * genes), which varies greatly by the CpG and GC content of the region. See expected Ti/Tv ratios section of the GATK best
     * practices documentation (http://www.broadinstitute.org/gatk/guide/best-practices) for more information.
     * Normal values are 2.15 for human whole genome values and 3.2 for human whole exomes. Note
     * that this parameter is used for display purposes only and isn't used anywhere in the algorithm!
     */
    @Argument(fullName="target_titv",
            shortName="titv",
            doc="The expected novel Ti/Tv ratio to use when calculating FDR tranches and for display on the optimization curve output figures. (approx 2.15 for whole genome experiments). ONLY USED FOR PLOTTING PURPOSES!",
            optional=true)
    private double TARGET_TITV = 2.15;

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(fullName="use_annotation",
            shortName="an",
            doc="The names of the annotations which should used for calculations",
            optional=false)
    private List<String> USE_ANNOTATIONS = new ArrayList<>();

    /**
     * Add truth sensitivity slices through the call set at the given values. The default values are 100.0, 99.9, 99.0, and 90.0
     * which will result in 4 estimated tranches in the final call set: the full set of calls (100% sensitivity at the accessible
     * sites in the truth set), a 99.9% truth sensitivity tranche, along with progressively smaller tranches at 99% and 90%.
     */
    @Argument(fullName="TStranche",
            shortName="tranche",
            doc="The levels of truth sensitivity at which to slice the data. (in percent, that is 1.0 for 1 percent)",
            optional=true)
    private List<Double> TS_TRANCHES = new ArrayList<Double>(Arrays.asList(100.0, 99.9, 99.0, 90.0));

    /**
     * For this to work properly, the -ignoreFilter argument should also be applied to the ApplyRecalibration command.
     */
    @Argument(fullName="ignore_filter",
            shortName="ignoreFilter",
            doc="If specified, the variant recalibrator will also use variants marked as filtered by the specified filter name in the input VCF file",
            optional=true)
    private List<String> IGNORE_INPUT_FILTERS = new ArrayList<>();

    @Argument(fullName="ignore_all_filters",
            shortName="ignoreAllFilters",
            doc="If specified, the variant recalibrator will ignore all input filters. Useful to rerun the VQSR from a filtered output file.",
            optional=true)
    private boolean IGNORE_ALL_FILTERS = false;

    @Argument(fullName="rscript_file", shortName="rscriptFile", doc="The output rscript file generated by the VQSR to aid in visualization of the input data and learned model", optional=true)
    private String RSCRIPT_FILE = null;

    /**
     *  This GATKReport gives information to describe the VQSR model fit. Normalized means for the positive model are
     *  concatenated as one table and negative model normalized means as another table. Covariances are also concatenated
     *  for positive and negative models, respectively. Tables of annotation means and standard deviations are provided
     *  to help describe the normalization. The model fit report can be read in with our R gsalib package. Individual
     *  model Gaussians can be subset by the value in the "Gaussian" column if desired.
     */
    @Argument(fullName="output_model",
            shortName = "outputModel",
            doc="If specified, the variant recalibrator will output the VQSR model fit to the file specified by -modelFile or to stdout",
            optional=true)
    private boolean outputModel = false;

    @Argument(fullName="model_file",
            shortName = "modelFile",
            doc="A GATKReport containing the positive and negative model fits",
            optional=true)
    private String modelReport = null;

    @Hidden
    @Argument(fullName="replicate",
            shortName="replicate",
            doc="Used to debug the random number generation inside the VQSR. Do not use.",
            optional=true)
    private int REPLICATE = 200;

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
    @Argument(fullName="max_attempts",
            shortName = "max_attempts",
            doc="Number of attempts to build a model before failing",
            optional=true)
    @VisibleForTesting
    protected int max_attempts = 1;

    /////////////////////////////
    // Debug Arguments
    /////////////////////////////
    @Advanced
    @Argument(fullName = "trustAllPolymorphic",
            shortName = "allPoly",
            doc = "Trust that all the input training sets' unfiltered records contain only polymorphic sites to drastically speed up the computation.",
            optional=true)
    private boolean TRUST_ALL_POLYMORPHIC = false;

    // Temporary argument for validation of GATK4 implementation against GATKs results:
    @Advanced
    @Argument(fullName="gatk3Compatibility", shortName="gatk3",
            doc="Set initial random number generator state for an exact match against GATK3 results.",
            optional=true)
    private boolean gatk3Compatibility = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantDataManager dataManager;
    private VariantContextWriter recalWriter;
    private PrintStream tranchesStream;
    private ArrayList<Double> replicate = new ArrayList<>(REPLICATE * 2);
    private final Set<String> ignoreInputFilterSet = new TreeSet<>();
    private final VariantRecalibratorEngine engine = new VariantRecalibratorEngine( VRAC );
    private ExpandingArrayList<VariantDatum> reduceSum = new ExpandingArrayList<>(2000);
    private List<ImmutablePair<VariantContext, FeatureContext>> variantsAtLocus = new ArrayList<>();

    //---------------------------------------------------------------------------------------------------------------
    //
    // onTraversalStart
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void onTraversalStart() {

        if (gatk3Compatibility) {
            // Temporary argument for validation of GATK4 implementation against GATK3 results:
            // Reset the RNG and draw a single int to align the RNG initial state with that used
            // by GATK3 to allow comparison of results with GATK3
            Utils.resetRandomGenerator();
            Utils.getRandomGenerator().nextInt();
        }
        dataManager = new VariantDataManager( new ArrayList<>(USE_ANNOTATIONS), VRAC );

        if (RSCRIPT_FILE != null && !RScriptExecutor.RSCRIPT_EXISTS)
            Utils.warnUser(logger, String.format(
                    "Rscript not found in environment path. %s will be generated but PDF plots will not.",
                    RSCRIPT_FILE));

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( IGNORE_INPUT_FILTERS );
        }

        try {
            tranchesStream = new PrintStream(TRANCHES_FILE);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(TRANCHES_FILE, e);
        }

        for( FeatureInput<VariantContext> variantFeature : resource ) {
            dataManager.addTrainingSet( new TrainingSet( variantFeature ) );
        }

        if( !dataManager.checkHasTrainingSet() ) {
            throw new UserException.CommandLineException(
                    "No training set found! Please provide sets of known polymorphic loci marked with the training=true feature input tag. For example, -resource hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }
        if( !dataManager.checkHasTruthSet() ) {
            throw new UserException.CommandLineException(
                    "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true feature input tag. For example, -resource hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf" );
        }

        //TODO: this should be refactored/consolidated as part of
        // https://github.com/broadinstitute/gatk/issues/2112
        // https://github.com/broadinstitute/gatk/issues/121,
        // https://github.com/broadinstitute/gatk/issues/1116 and
        // Initialize VCF header lines
        Set<VCFHeaderLine> hInfo = new HashSet<>();
        hInfo.add(new VCFHeaderLine("source", this.getClass().getSimpleName()));
        VariantRecalibrationUtils.addVQSRStandardHeaderLines(hInfo);
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (hasReference()) {
            hInfo = VcfUtils.updateHeaderContigLines(
                    hInfo, referenceArguments.getReferenceFile(), sequenceDictionary, true);
        }
        else if (null != sequenceDictionary) {
            hInfo = VcfUtils.updateHeaderContigLines(hInfo, null, sequenceDictionary, true);
        }

        recalWriter = createVCFWriter(new File(output));
        recalWriter.writeHeader( new VCFHeader(hInfo) );

        for( int iii = 0; iii < REPLICATE * 2; iii++ ) {
            replicate.add(Utils.getRandomGenerator().nextDouble());
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // apply
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {
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
        }
        else {
            ImmutablePair<VariantContext, FeatureContext> previous = variantsAtLocus.get(0);
            return (nextVC.getStart() != previous.left.getStart() || !nextVC.getContig().equals(previous.left.getContig()));
        }
    }

    private void consumeQueuedVariants() {
        variantsAtLocus.forEach(v -> addVariantDatum(v.left, true, v.right));
        if (!aggregate.isEmpty()) {
            // use the first featureContext in the queue for the aggregate resources
            addOverlappingAggregateVariants(aggregate, false, variantsAtLocus.get(0).getRight());
        }
        variantsAtLocus.clear();
    }

    /**
     * Find overlapping variants and pull out the necessary information to create the VariantDatum
     * @param aggregateInputs the input sources to search within
     * @param isInput   is this the driving variant input (true) or an aggregate input ?
     * @param context   the FeatureContext from the apply call
     * @return  a list of VariantDatums, can be empty
     */
    private void addOverlappingAggregateVariants(
            final List<FeatureInput<VariantContext>> aggregateInputs,
            final boolean isInput,
            final FeatureContext context ) {
        if( aggregateInputs == null ) { throw new IllegalArgumentException("aggregateInputs cannot be null."); }
        if( context == null ) { throw new IllegalArgumentException("context cannot be null."); }

        for( final VariantContext vc : context.getValues(aggregateInputs, context.getInterval().getStart()) ) {
            addVariantDatum(vc, isInput, context);
        }
    }

    private void addVariantDatum( VariantContext vc, final boolean isInput, final FeatureContext context ) {
        if( vc != null && ( IGNORE_ALL_FILTERS || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()) ) ) {
            if( VariantDataManager.checkVariationClass( vc, VRAC.MODE ) && !VRAC.useASannotations) {
                addDatum(reduceSum, isInput, context, vc, null, null);
            }
            else if( VRAC.useASannotations ) {
                for (final Allele allele : vc.getAlternateAlleles()) {
                    if ( allele == Allele.SPAN_DEL ) {
                        continue;
                    }
                    if ( VariantDataManager.checkVariationClass(vc, allele, VRAC.MODE )) {
                        addDatum(reduceSum, isInput, context, vc, vc.getReference(), allele);
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
            final ExpandingArrayList<VariantDatum> variants,
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

        for (int i = 1; i <= max_attempts; i++) {
            try {
                dataManager.setData(reduceSum);
                dataManager.normalizeData(); // Each data point is now (x - mean) / standard deviation

                // Generate the positive model using the training data and evaluate each variant
                final List<VariantDatum> positiveTrainingData = dataManager.getTrainingData();

                final GaussianMixtureModel goodModel = engine.generateModel(positiveTrainingData, VRAC.MAX_GAUSSIANS);
                engine.evaluateData(dataManager.getData(), goodModel, false);

                // Generate the negative model using the worst performing data and evaluate each variant contrastively
                final List<VariantDatum> negativeTrainingData = dataManager.selectWorstVariants();
                final GaussianMixtureModel badModel = engine.generateModel(negativeTrainingData,
                        Math.min(VRAC.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL, VRAC.MAX_GAUSSIANS));
                dataManager.dropAggregateData(); // Don't need the aggregate data anymore so let's free up the memory
                engine.evaluateData(dataManager.getData(), badModel, true);

                if (badModel.failedToConverge || goodModel.failedToConverge) {
                    throw new UserException(
                            "NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider " + (badModel.failedToConverge ? "raising the number of variants used to train the negative model (via --minNumBadVariants 5000, for example)." : "lowering the maximum number of Gaussians allowed for use in the model (via --maxGaussians 4, for example)."));
                }

                if (outputModel) {
                    GATKReport report = writeModelReport(goodModel, badModel, USE_ANNOTATIONS);
                    try(final PrintStream modelReportStream = new PrintStream(modelReport)) {
                        report.print(modelReportStream);
                    } catch (FileNotFoundException e) {
                        throw new UserException.CouldNotCreateOutputFile("File: (" + modelReport + ")", e);
                    }
                }

                engine.calculateWorstPerformingAnnotation(dataManager.getData(), goodModel, badModel);

                // Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
                final int nCallsAtTruth = TrancheManager.countCallsAtTruth(dataManager.getData(), Double.NEGATIVE_INFINITY);
                final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric(nCallsAtTruth);
                final List<Tranche> tranches = TrancheManager.findTranches(dataManager.getData(), TS_TRANCHES, metric, VRAC.MODE);
                tranchesStream.print(Tranche.tranchesString(tranches));

                logger.info("Writing out recalibration table...");
                dataManager.writeOutRecalibrationTable(recalWriter, getBestAvailableSequenceDictionary());
                if (RSCRIPT_FILE != null) {
                    logger.info("Writing out visualization Rscript file...");
                    createVisualizationScript(dataManager.getRandomDataForPlotting(
                            1000,
                            positiveTrainingData,
                            negativeTrainingData,
                            dataManager.getEvaluationData()),
                            goodModel,
                            badModel,
                            0.0,
                            dataManager.getAnnotationKeys().toArray(new String[USE_ANNOTATIONS.size()]));
                }

                if (VRAC.MODE == VariantRecalibratorArgumentCollection.Mode.INDEL) {
                    // Print out an info message to make it clear why the tranches plot is not generated
                    logger.info("Tranches plot will not be generated since we are running in INDEL mode");
                } else {
                    // Execute the RScript command to plot the table of truth values
                    RScriptExecutor executor = new RScriptExecutor();
                    executor.addScript(new Resource(PLOT_TRANCHES_RSCRIPT, VariantRecalibrator.class));
                    executor.addArgs(new File(TRANCHES_FILE).getAbsoluteFile(), TARGET_TITV);
                    // Print out the command line to make it clear to the user what is being executed and how one might modify it
                    logger.info("Executing: " + executor.getApproximateCommandLine());
                    executor.exec();
                }
                return true;
            }
            catch (Exception e) {
                //TODO: see https://github.com/broadgsa/gatk-protected/pull/19
                // If we retain this loop, it needs to be restructured. Currently, any existing state (including data
                // already written to output streams) is retained for the next attempt. The output streams should
                // really be re-created and all state should be reset. Also, because we're catching Exception, this
                // will retry even if we have bad input, such as poorly annotated training data that results in an
                // IllegalArgumentException.
                if (i == max_attempts) {
                    throw e;
                } else {
                    logger.info(String.format(
                            "Exception occurred on attempt %d of %d. Trying again. Message was: '%s'",
                            i,
                            max_attempts,
                            e.getMessage()));
                }
            }
        }

        return false;
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

    protected GATKReport writeModelReport(
            final GaussianMixtureModel goodModel,
            final GaussianMixtureModel badModel,
            List<String> annotationList) {
        final String formatString = "%.3f";
        final GATKReport report = new GATKReport();

        if (dataManager != null) {  //for unit test
            final double[] meanVector = dataManager.getMeanVector();
            GATKReportTable annotationMeans = makeVectorTable(
                    "AnnotationMeans",
                    "Mean for each annotation, used to normalize data",
                    dataManager.annotationKeys,
                    meanVector,
                    "Mean",
                    formatString);
            report.addTable(annotationMeans);

            final double[] varianceVector = dataManager.getVarianceVector();  //"varianceVector" is actually stdev
            GATKReportTable annotationVariances = makeVectorTable(
                    "AnnotationStdevs",
                    "Standard deviation for each annotation, used to normalize data",
                    dataManager.annotationKeys,
                    varianceVector,
                    "Standard deviation",
                    formatString);
            report.addTable(annotationVariances);
        }

        //The model and Gaussians don't know what the annotations are, so get them from this class
        //VariantDataManager keeps the annotation in the same order as the argument list
        GATKReportTable positiveMeans = makeMeansTable(
                "PositiveModelMeans", "Vector of annotation values to describe the (normalized) mean for each Gaussian in the positive model",
                annotationList,
                goodModel,
                formatString);
        report.addTable(positiveMeans);

        GATKReportTable positiveCovariance = makeCovariancesTable(
                "PositiveModelCovariances", "Matrix to describe the (normalized) covariance for each Gaussian in the positive model; covariance matrices are joined by row",
                annotationList,
                goodModel,
                formatString);
        report.addTable(positiveCovariance);

        //do the same for the negative model means
        GATKReportTable negativeMeans = makeMeansTable(
                "NegativeModelMeans", "Vector of annotation values to describe the (normalized) mean for each Gaussian in the negative model",
                annotationList, badModel, formatString);
        report.addTable(negativeMeans);

        GATKReportTable negativeCovariance = makeCovariancesTable(
                "NegativeModelCovariances",
                "Matrix to describe the (normalized) covariance for each Gaussian in the negative model; covariance matrices are joined by row",
                annotationList,
                badModel,
                formatString);
        report.addTable(negativeCovariance);

        return report;
    }

    protected GATKReportTable makeVectorTable(
            final String tableName,
            final String tableDescription,
            final List<String> annotationList,
            final double[] perAnnotationValues,
            final String columnName,
            final String formatString) {
        GATKReportTable vectorTable = new GATKReportTable(
                tableName, tableDescription, annotationList.size(), GATKReportTable.Sorting.DO_NOT_SORT);
        vectorTable.addColumn("Annotation", "%s");
        vectorTable.addColumn(columnName, formatString);
        for (int i = 0; i < perAnnotationValues.length; i++) {
            vectorTable.addRowIDMapping(annotationList.get(i), i, true);
            vectorTable.set(i, 1, perAnnotationValues[i]);
        }
        return vectorTable;
    }

    private GATKReportTable makeMeansTable(
            final String tableName,
            final String tableDescription,
            final List<String> annotationList,
            final GaussianMixtureModel model,
            final String formatString) {
        GATKReportTable meansTable = new GATKReportTable(
                tableName, tableDescription, annotationList.size(), GATKReportTable.Sorting.DO_NOT_SORT);
        meansTable.addColumn("Gaussian", "");
        for (final String annotationName : annotationList) {
            meansTable.addColumn(annotationName, formatString);
        }
        final List<MultivariateGaussian> modelGaussians = model.getModelGaussians();
        for (int i = 0; i < modelGaussians.size(); i++) {
            final MultivariateGaussian gaussian = modelGaussians.get(i);
            final double[] meanVec = gaussian.mu;
            if (meanVec.length != annotationList.size())
                throw new IllegalStateException("Gaussian mean vector does not have the same size as the list of annotations");
            meansTable.addRowIDMapping(i, i, true);
            for (int j = 0; j < annotationList.size(); j++)
                meansTable.set(i, annotationList.get(j), meanVec[j]);
        }
        return meansTable;
    }

    private GATKReportTable makeCovariancesTable(
            final String tableName,
            final String tableDescription,
            final List<String> annotationList,
            final GaussianMixtureModel model,
            final String formatString) {
        GATKReportTable modelCovariances = new GATKReportTable(
                tableName, tableDescription, annotationList.size()+2, GATKReportTable.Sorting.DO_NOT_SORT); //+2 is for Gaussian and Annotation columns
        modelCovariances.addColumn("Gaussian", "");
        modelCovariances.addColumn("Annotation", "");
        for (final String annotationName : annotationList) {
            modelCovariances.addColumn(annotationName, formatString);
        }
        final List<MultivariateGaussian> modelGaussians = model.getModelGaussians();
        for (int i = 0; i < modelGaussians.size(); i++) {
            final MultivariateGaussian gaussian = modelGaussians.get(i);
            final Matrix covMat = gaussian.sigma;
            if (covMat.getRowDimension() != annotationList.size() || covMat.getColumnDimension() != annotationList.size())
                throw new IllegalStateException("Gaussian covariance matrix does not have the same size as the list of annotations");
            for (int j = 0; j < annotationList.size(); j++) {
                modelCovariances.set(j + i * annotationList.size(), "Gaussian", i);
                modelCovariances.set(j + i * annotationList.size(), "Annotation", annotationList.get(j));
                for (int k = 0; k < annotationList.size(); k++) {
                    modelCovariances.set(j + i * annotationList.size(), annotationList.get(k), covMat.get(j, k));

                }
            }
        }
        return modelCovariances;
    }

    //TODO: does this R code have to be embedded here?
    private void createVisualizationScript(
            final List<VariantDatum> randomData,
            final GaussianMixtureModel goodModel,
            final GaussianMixtureModel badModel,
            final double lodCutoff,
            final String[] annotationKeys ) {
        PrintStream stream;
        try {
            stream = new PrintStream(RSCRIPT_FILE);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(RSCRIPT_FILE, e);
        }

        // We make extensive use of the ggplot2 R library: http://had.co.nz/ggplot2/
        stream.println("library(ggplot2)");
        // For compactPDF in R 2.13+
        stream.println("library(tools)");
        // For graphical functions R 2.14.2+
        stream.println("library(grid)");

        createArrangeFunction( stream );

        stream.println("outputPDF <- \"" + RSCRIPT_FILE + ".pdf\"");
        stream.println("pdf(outputPDF)"); // Unfortunately this is a huge pdf file, BUGBUG: need to work on reducing the file size

        for(int iii = 0; iii < annotationKeys.length; iii++) {
            for( int jjj = iii + 1; jjj < annotationKeys.length; jjj++) {
                logger.info( "Building " + annotationKeys[iii] + " x " + annotationKeys[jjj] + " plot...");

                final List<VariantDatum> fakeData = new ExpandingArrayList<>();
                double minAnn1 = 100.0, maxAnn1 = -100.0, minAnn2 = 100.0, maxAnn2 = -100.0;
                for( final VariantDatum datum : randomData ) {
                    minAnn1 = Math.min(minAnn1, datum.annotations[iii]);
                    maxAnn1 = Math.max(maxAnn1, datum.annotations[iii]);
                    minAnn2 = Math.min(minAnn2, datum.annotations[jjj]);
                    maxAnn2 = Math.max(maxAnn2, datum.annotations[jjj]);
                }
                // Create a fake set of data which spans the full extent of these two annotation dimensions in order
                // to calculate the model PDF projected to 2D
                final double NUM_STEPS = 60.0;
                for(double ann1 = minAnn1; ann1 <= maxAnn1; ann1+= (maxAnn1 - minAnn1) / NUM_STEPS) {
                    for(double ann2 = minAnn2; ann2 <= maxAnn2; ann2+= (maxAnn2 - minAnn2) / NUM_STEPS) {
                        final VariantDatum datum = new VariantDatum();
                        datum.prior = 0.0;
                        datum.annotations = new double[randomData.get(0).annotations.length];
                        datum.isNull = new boolean[randomData.get(0).annotations.length];
                        for(int ann=0; ann< datum.annotations.length; ann++) {
                            datum.annotations[ann] = 0.0;
                            datum.isNull[ann] = true;
                        }
                        datum.annotations[iii] = ann1;
                        datum.annotations[jjj] = ann2;
                        datum.isNull[iii] = false;
                        datum.isNull[jjj] = false;
                        fakeData.add(datum);
                    }
                }

                engine.evaluateData( fakeData, goodModel, false );
                engine.evaluateData( fakeData, badModel, true );

                stream.print("surface <- c(");
                for( final VariantDatum datum : fakeData ) {
                    stream.print(String.format("%.4f, %.4f, %.4f, ",
                            dataManager.denormalizeDatum(datum.annotations[iii], iii),
                            dataManager.denormalizeDatum(datum.annotations[jjj], jjj),
                            Math.min(4.0, Math.max(-4.0, datum.lod))));
                }
                stream.println("NA,NA,NA)");
                stream.println("s <- matrix(surface,ncol=3,byrow=T)");

                stream.print("data <- c(");
                for( final VariantDatum datum : randomData ) {
                    stream.print(String.format("%.4f, %.4f, %.4f, %d, %d,",
                            dataManager.denormalizeDatum(datum.annotations[iii], iii),
                            dataManager.denormalizeDatum(datum.annotations[jjj], jjj),
                            (datum.lod < lodCutoff ? -1.0 : 1.0),
                            (datum.atAntiTrainingSite ? -1 : (datum.atTrainingSite ? 1 : 0)), (datum.isKnown ? 1 : -1)));
                }
                stream.println("NA,NA,NA,NA,1)");
                stream.println("d <- matrix(data,ncol=5,byrow=T)");

                final String surfaceFrame = "sf." + annotationKeys[iii] + "." + annotationKeys[jjj];
                final String dataFrame = "df." + annotationKeys[iii] + "." + annotationKeys[jjj];

                stream.println(surfaceFrame + " <- data.frame(x=s[,1], y=s[,2], lod=s[,3])");
                stream.println(dataFrame + " <- data.frame(x=d[,1], y=d[,2], retained=d[,3], training=d[,4], novelty=d[,5])");
                stream.println("dummyData <- " + dataFrame + "[1,]");
                stream.println("dummyData$x <- NaN");
                stream.println("dummyData$y <- NaN");
                stream.println("p <- ggplot(data=" + surfaceFrame + ", aes(x=x, y=y)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("p1 = p +ggtitle(\"model PDF\") + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + geom_tile(aes(fill = lod)) + scale_fill_gradient(high=\"green\", low=\"red\", space=\"rgb\")");
                stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=retained, alpha=I(1/7),legend=FALSE) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=retained),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p2 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(name=\"outcome\", high=\"black\", low=\"red\",breaks=c(-1,1),guide=\"legend\",labels=c(\"filtered\",\"retained\"))");
                stream.println("p <- qplot(x,y,data="+ dataFrame + "["+dataFrame+"$training != 0,], color=training, alpha=I(1/7)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=training),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p3 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(high=\"green\", low=\"purple\",breaks=c(-1,1),guide=\"legend\", labels=c(\"neg\", \"pos\"))");
                stream.println("p <- qplot(x,y,data=" + dataFrame + ", color=novelty, alpha=I(1/7)) +theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))");
                stream.println("q <- geom_point(aes(x=x,y=y,color=novelty),data=dummyData, alpha=1.0, na.rm=TRUE)");
                stream.println("p4 = p + q + labs(x=\""+ annotationKeys[iii] +"\", y=\""+ annotationKeys[jjj] +"\") + scale_colour_gradient(name=\"novelty\", high=\"blue\", low=\"red\",breaks=c(-1,1),guide=\"legend\", labels=c(\"novel\",\"known\"))");
                stream.println("arrange(p1, p2, p3, p4, ncol=2)");
            }
        }
        stream.println("dev.off()");

        stream.println("if (exists(\"compactPDF\")) {");
        stream.println("compactPDF(outputPDF)");
        stream.println("}");

        stream.close();

        // Execute Rscript command to generate the clustering plots
        RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new File(RSCRIPT_FILE));
        logger.info("Executing: " + executor.getApproximateCommandLine());
        executor.exec();
    }

    // The Arrange function is how we place the 4 model plots on one page
    // from http://gettinggeneticsdone.blogspot.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html
    private void createArrangeFunction( final PrintStream stream ) {
        stream.println("vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)");
        stream.println("arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {");
        stream.println("dots <- list(...)");
        stream.println("n <- length(dots)");
        stream.println("if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}");
        stream.println("if(is.null(nrow)) { nrow = ceiling(n/ncol)}");
        stream.println("if(is.null(ncol)) { ncol = ceiling(n/nrow)}");
        stream.println("grid.newpage()");
        stream.println("pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )");
        stream.println("ii.p <- 1");
        stream.println("for(ii.row in seq(1, nrow)){");
        stream.println("ii.table.row <- ii.row ");
        stream.println("if(as.table) {ii.table.row <- nrow - ii.table.row + 1}");
        stream.println("for(ii.col in seq(1, ncol)){");
        stream.println("ii.table <- ii.p");
        stream.println("if(ii.p > n) break");
        stream.println("print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))");
        stream.println("ii.p <- ii.p + 1");
        stream.println("}");
        stream.println("}");
        stream.println("}");
    }
}
