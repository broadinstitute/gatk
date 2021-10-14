package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.math3.linear.MatrixUtils;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixture;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import org.apache.commons.lang3.tuple.ImmutablePair;

import java.io.*;
import java.util.*;

/**
 * Build a recalibration model to score variant quality for filtering purposes
 *
 * <p>This tool performs the first pass in a two-stage process called Variant Quality Score Recalibration (VQSR).
 * Specifically, it builds the model that will be used in the second step to actually filter variants. This model
 * attempts to describe the relationship between variant annotations (such as QD, MQ and ReadPosRankSum, for example)
 * and the probability that a variant is a true genetic variant versus a sequencing or data processing artifact. It is
 * developed adaptively based on "true sites" provided as input, typically HapMap sites and those sites found to be
 * polymorphic on the Omni 2.5M SNP chip array (in humans). This adaptive error model can then be applied to both known
 * and novel variation discovered in the call set of interest to evaluate the probability that each call is real. The
 * result is a score called the VQSLOD that gets added to the INFO field of each variant. This score is the log odds of
 * being a true variant versus being false under the trained Gaussian mixture model. </p>
 *
 * <h4>Summary of the VQSR procedure</h4>
 * <p>The purpose of variant recalibration is to assign a well-calibrated probability to each variant call in a call set.
 * These probabilities can then be used to filter the variants with a greater level of accuracy and flexibility than
 * can typically be achieved by traditional hard-filter (filtering on individual annotation value thresholds). The first
 * pass consists of building a model that describes how variant annotation values co-vary with the truthfulness of
 * variant calls in a training set, and then scoring all input variants according to the model. The second pass simply
 * consists of specifying a target sensitivity value (which corresponds to an empirical VQSLOD cutoff) and applying
 * filters to each variant call according to their ranking. The result is a VCF file in which variants have been
 * assigned a score and filter status.</p>
 *
 * <p>VQSR is probably the hardest part of the Best Practices to get right, so be sure to read the
 * <a href='https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-'>method documentation</a> and
 * <a href='https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering'>tutorial</a> to really understand what these
 * tools do and how to use them for best results on your own data.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *      <li>The input variants to be recalibrated. These variant calls must be annotated with the annotations that will be
 * used for modeling. If the calls come from multiple samples, they must have been obtained by joint calling the samples,
 * either directly (running HaplotypeCaller on all samples together) or via the GVCF workflow (HaplotypeCaller with -ERC
 * GVCF per-sample then GenotypeGVCFs on the resulting gVCFs) which is more scalable.</li>
 *      <li>Known, truth, and training sets to be used by the algorithm. See the method documentation linked above for
 * more details.</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 * <li>A recalibration table file that will be used by the ApplyVQSR tool.</li>
 * <li>A tranches file that shows various metrics of the recalibration callset for slices of the data.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Recalibrating SNPs in exome data</h4>
 * <pre>
 * gatk VariantRecalibrator \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V input.vcf.gz \
 *   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.sites.vcf.gz \
 *   --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.sites.vcf.gz \
 *   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
 *   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
 *   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
 *   -mode SNP \
 *   -O output.recal \
 *   --tranches-file output.tranches \
 *   --rscript-file output.plots.R
 * </pre>
 *
 * <h4>Allele-specific version of the SNP recalibration (beta)</h4>
 * <pre>
 * gatk VariantRecalibrator \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V input.vcf.gz \
 *   -AS \
 *   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.sites.vcf.gz \
 *   --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.sites.vcf.gz \
 *   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
 *   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
 *   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
 *   -mode SNP \
 *   -O output.AS.recal \
 *   --tranches-file output.AS.tranches \
 *   --rscript-file output.plots.AS.R
 * </pre>
 * <p>Note that to use the allele-specific (AS) mode, the input VCF must have been produced using allele-specific
 * annotations in HaplotypeCaller. Note also that each allele will have a separate line in the output recalibration
 * file with its own VQSLOD and `culprit`, which will be transferred to the final VCF by the ApplyVQSR tool.
 *
 * <h3>Caveats</h3>
 *
 * <ul>
 * <li>The values used in the example above are only meant to show how the command lines are composed.
 * They are not meant to be taken as specific recommendations of values to use in your own work, and they may be
 * different from the values cited elsewhere in our documentation. For the latest and greatest recommendations on
 * how to set parameter values for your own analyses, please read the Best Practices section of the documentation,
 * especially the <a href='https://software.broadinstitute.org/gatk/guide/article?id=1259'>FAQ document</a> on VQSR parameters.</li>
 * <li>Whole genomes and exomes take slightly different parameters, so make sure you adapt your commands accordingly! See
 * the documents linked above for details.</li>
 * <li>If you work with small datasets (e.g. targeted capture experiments or small number of exomes), you will run into
 * problems. Read the docs linked above for advice on how to deal with those issues.</li>
 * <li>In order to create the model reporting plots, the Rscript executable needs to be in your environment PATH
 * (this is the scripting version of R, not the interactive version).
 * See <a target="r-project" href="http://www.r-project.org">http://www.r-project.org</a> for more information on how
 * to download and install R.</li>
 * </ul>
 *
 * <h3>Additional notes</h3>
 * <ul>
 *     <li>This tool only accepts a single input variant file unlike earlier version of GATK, which accepted multiple
 *     input variant files.</li>
 *     <li>SNPs and indels must be recalibrated in separate runs, but it is not necessary to separate them into different
 * files. See the tutorial linked above for an example workflow. Note that mixed records are treated as indels.</li>
 *     <li></li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = "Build a recalibration model to score variant quality for filtering purposes",
        oneLineSummary = "Build a recalibration model to score variant quality for filtering purposes",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public class VariantRecalibrator extends MultiVariantWalker {

    @ArgumentCollection
    final private VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();

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
     * The expected transition / transversion ratio of true novel variants in your targeted region (whole genome, exome, specific
     * genes), which varies greatly by the CpG and GC content of the region. See expected Ti/Tv ratios section of the GATK best
     * practices documentation (https://software.broadinstitute.org/gatk/guide/best-practices) for more information.
     * Normal values are 2.15 for human whole genome values and 3.2 for human whole exomes. Note
     * that this parameter is used for display purposes only and isn't used anywhere in the algorithm!
     */
    @Argument(fullName="target-titv",
            shortName="titv",
            doc="The expected novel Ti/Tv ratio to use when calculating FDR tranches and for display on the optimization curve output figures. (approx 2.15 for whole genome experiments). ONLY USED FOR PLOTTING PURPOSES!",
            optional=true)
    private double TARGET_TITV = 2.15;

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
     *  This GATKReport gives information to describe the VQSR model fit. Normalized means for the positive model are
     *  concatenated as one table and negative model normalized means as another table. Covariances are also concatenated
     *  for positive and negative models, respectively. Tables of annotation means and standard deviations are provided
     *  to help describe the normalization. The model fit report can be read in with our R gsalib package. Individual
     *  model Gaussians can be subset by the value in the "Gaussian" column if desired.
     */
    @Argument(fullName="output-model",
            doc="If specified, the variant recalibrator will output the VQSR model to this file path.",
            optional=true)
    private GATKPath outputModel = null;

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
    private VariantDataManager dataManager;
    private VariantContextWriter recalWriter;
    private PrintStream tranchesStream;
    final private ArrayList<Double> replicate = new ArrayList<>(REPLICATE * 2);
    private final Set<String> ignoreInputFilterSet = new TreeSet<>();
    private final VariantRecalibratorEngine engine = new VariantRecalibratorEngine( VRAC );
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

        dataManager = new VariantDataManager( new ArrayList<>(USE_ANNOTATIONS), VRAC );

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

        for ( int iii = 0; iii < REPLICATE * 2; iii++ ) {
            replicate.add(Utils.getRandomGenerator().nextDouble());
        }
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

    private void addVariantDatum(final VariantContext vc, final boolean isInput, final FeatureContext context ) {
        if( vc != null && ( IGNORE_ALL_FILTERS || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()) ) ) {
            if( VariantDataManager.checkVariationClass( vc, VRAC.MODE ) && !VRAC.useASannotations) {
                addDatum(reduceSum, isInput, context, vc, null, null);
            }
            else if( VRAC.useASannotations ) {
                for (final Allele allele : vc.getAlternateAlleles()) {
                    if (!GATKVCFConstants.isSpanningDeletion(allele) && VariantDataManager.checkVariationClass(vc, allele, VRAC.MODE)) {
                        //note that this may not be the minimal representation for the ref and alt allele
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
        dataManager.normalizeData(true, annotationOrder); // Each data point is now (x - mean) / standard deviation

        final List<VariantDatum> positiveTrainingData = dataManager.getTrainingData();
        final double[][] positiveTrainingDataArray = positiveTrainingData.stream().map(vd -> vd.annotations).toArray(double[][]::new);

        final int nFeatures = USE_ANNOTATIONS.size();
        final double[] meanPrior = new double[nFeatures];
        Arrays.fill(meanPrior, 0.);
        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();

        final BayesianGaussianMixture bgmm = new BayesianGaussianMixture.Builder()
                .nComponents(VRAC.MAX_GAUSSIANS)
                .maxIter(VRAC.MAX_ITERATIONS)
                .nInit(max_attempts)
                .initMethod(BayesianGaussianMixture.InitMethod.K_MEANS_PLUS_PLUS)
                .weightConcentrationPrior(VRAC.DIRICHLET_PARAMETER)
                .meanPrior(meanPrior)
                .degreesOfFreedomPrior(nFeatures)
                .covariancePrior(covariancePrior)
                .seed(1)
                .warmStart(true)
                .verboseInterval(1)
                .build();
        bgmm.fit(positiveTrainingDataArray);

        System.out.println("weights: " + bgmm.getWeights());
        System.out.println("meanPrecision: " + bgmm.getMeanPrecision());
        System.out.println("means: " + bgmm.getMeans());
        System.out.println("precisionsCholesky: " + bgmm.getPrecisionsCholesky());
        System.out.println("covariances: " + bgmm.getCovariances());
        System.out.println("degreesOfFreedom: " + bgmm.getDegreesOfFreedom());

        final double[][] data = dataManager.getData().stream().map(vd -> vd.annotations).toArray(double[][]::new);
        final double[] scores = bgmm.scoreSamples(data);
//        dumpScores(scores, output + ".scores.tsv");

        dataManager.dropAggregateData(); // Don't need the aggregate data anymore so let's free up the memory
        dataManager.setScores(dataManager.getData(), scores);

        // Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
        final int nCallsAtTruth = TrancheManager.countCallsAtTruth(dataManager.getData(), Double.NEGATIVE_INFINITY);
        final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric(nCallsAtTruth);
        final List<? extends Tranche> tranches = TrancheManager.findTranches(dataManager.getData(), TS_TRANCHES, metric, VRAC.MODE);
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

    private static void dumpScores(double[] scores, String output) {
        try (FileWriter fos = new FileWriter(output);
             PrintWriter dos = new PrintWriter(fos)) {

            for (int i = 0; i < scores.length; i++)
            {
                dos.print(scores[i] + "\n");
            }
        } catch (IOException e) {
            System.out.println("Error printing scores.");
        }
    }
}