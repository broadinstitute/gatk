package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.common.annotations.VisibleForTesting;
import Jama.Matrix;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
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
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportColumn;
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
 * <a href='https://software.broadinstitute.org/gatk/guide/article?id=39'>method documentation</a>,
 * <a href='https://software.broadinstitute.org/gatk/guide/article?id=1259'>parameter recommendations</a> and
 * <a href='https://software.broadinstitute.org/gatk/guide/article?id=2805'>tutorial</a> to really understand what these
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

    private static final String PLOT_TRANCHES_RSCRIPT = "plot_Tranches.R";

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

    @Argument(fullName="rscript-file",doc="The output rscript file generated by the VQSR to aid in visualization of the input data and learned model", optional=true)
    private File RSCRIPT_FILE = null;

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
     *  The filename for a VQSR model fit to use to recalibrate the input variants. This model should be generated using
     *  a previous VariantRecalibration run with the --output-model argument.
     */
    @Argument(fullName="input-model",
            doc="If specified, the variant recalibrator will read the VQSR model from this file path.",
            optional=true)
    private GATKPath inputModel = null;

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
     *  This argument must be specified if VariantRecalibrator is run scattered because the tranch output format changes.
     *  The scattered recal files can be gathered with the GatherVcfs tool and the scattered tranches can be gathered
     *  with the GatherTranches tool. See the description of the -sampleEvery argument for more information on running
     *  scattered VariantRecalibrator.
     */

    @Argument(fullName="output-tranches-for-scatter",
            doc="Output tranches in a format appropriate to running VariantRecalibrator in scatter-gather",
            optional = true)
    @Hidden
    private boolean scatterTranches = false;

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
    private GATKReportTable nmcTable;
    private GATKReportTable nmmTable;
    private GATKReportTable nPMixTable;
    private GATKReportTable pmcTable;
    private GATKReportTable pmmTable;
    private GATKReportTable pPMixTable;
    private int numAnnotations;
    private RScriptExecutor rScriptExecutor;

    //---------------------------------------------------------------------------------------------------------------
    //
    // onTraversalStart
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void onTraversalStart() {

        dataManager = new VariantDataManager( new ArrayList<>(USE_ANNOTATIONS), VRAC );

        if (RSCRIPT_FILE != null) {
            rScriptExecutor = new RScriptExecutor();
            if(!rScriptExecutor.externalExecutableExists()) {
                Utils.warnUser(logger, String.format(
                        "Rscript not found in environment path. %s will be generated but PDF plots will not.",
                        RSCRIPT_FILE));
            }
        }

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

        if (inputModel != null) { // Load GMM from a file
            logger.info("Loading model from:" + inputModel);
            try (final InputStream is = inputModel.getInputStream()) {
                final GATKReport reportIn = new GATKReport(is);

                // Read all the tables
                nmcTable = reportIn.getTable("NegativeModelCovariances");
                nmmTable = reportIn.getTable("NegativeModelMeans");
                nPMixTable = reportIn.getTable("BadGaussianPMix");
                pmcTable = reportIn.getTable("PositiveModelCovariances");
                pmmTable = reportIn.getTable("PositiveModelMeans");
                pPMixTable = reportIn.getTable("GoodGaussianPMix");
                final GATKReportTable anMeansTable = reportIn.getTable("AnnotationMeans");
                final GATKReportTable anStDevsTable = reportIn.getTable("AnnotationStdevs");

                orderAndValidateAnnotations(anMeansTable, dataManager.annotationKeys);
                numAnnotations = annotationOrder.size();

                final Map<String, Double> anMeans = getMapFromVectorTable(anMeansTable);
                final Map<String, Double> anStdDevs = getMapFromVectorTable(anStDevsTable);
                dataManager.setNormalization(anMeans, anStdDevs);
            } catch (IOException e) {
                throw new UserException.CouldNotReadInputFile("File: (" + inputModel + ")", e);
            }
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

    /**
     * Order and validate annotations according to the annotations in the serialized model
     * Annotations on the command line must be the same as those in the model report or this will throw an exception.
     * Sets the {@code annotationOrder} list to map from command line order to the model report's order.
     * n^2 because we typically use 7 or less annotations.
     * @param annotationTable GATKReportTable of annotations read from the serialized model file
     */
    protected void orderAndValidateAnnotations(final GATKReportTable annotationTable, final List<String> annotationKeys){
        annotationOrder = new ArrayList<Integer>(annotationKeys.size());

        for (int i = 0; i < annotationTable.getNumRows(); i++){
            String serialAnno = (String)annotationTable.get(i, "Annotation");
            for (int j = 0; j < annotationKeys.size(); j++) {
                if (serialAnno.equals( annotationKeys.get(j))){
                    annotationOrder.add(j);
                }
            }
        }

        if(annotationOrder.size() != annotationTable.getNumRows() || annotationOrder.size() != annotationKeys.size()) {
            throw new CommandLineException( "Annotations specified on the command line do not match annotations in the model report." );
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

        for (int i = 1; i <= max_attempts; i++) {
            try {
                dataManager.setData(reduceSum);
                dataManager.normalizeData(inputModel == null, annotationOrder); // Each data point is now (x - mean) / standard deviation

                final GaussianMixtureModel goodModel;
                final GaussianMixtureModel badModel;

                final List<VariantDatum> positiveTrainingData = dataManager.getTrainingData();
                final List<VariantDatum> negativeTrainingData;

                if (inputModel != null) {  // GMMs were loaded from a file
                    logger.info("Using serialized GMMs from file...");
                    goodModel = GMMFromTables(pmmTable, pmcTable, pPMixTable, numAnnotations, positiveTrainingData.size());
                    engine.evaluateData(dataManager.getData(), goodModel, false);
                    negativeTrainingData = dataManager.selectWorstVariants();
                    badModel = GMMFromTables(nmmTable, nmcTable, nPMixTable, numAnnotations, negativeTrainingData.size());
                } else { // Generate the GMMs from scratch
                    // Generate the positive model using the training data and evaluate each variant
                    goodModel = engine.generateModel(positiveTrainingData, VRAC.MAX_GAUSSIANS);
                    engine.evaluateData(dataManager.getData(), goodModel, false);
                    if (goodModel.failedToConverge) {
                        if (outputModel != null) {
                            final GATKReport report = writeModelReport(goodModel, null, USE_ANNOTATIONS);
                            saveModelReport(report, outputModel);
                        }
                        throw new UserException.VQSRPositiveModelFailure("Positive training model failed to converge.  One or more annotations " +
                                "(usually MQ) may have insufficient variance.  Please consider lowering the maximum number" +
                                " of Gaussians allowed for use in the model (via --max-gaussians 4, for example).");
                    }
                    // Generate the negative model using the worst performing data and evaluate each variant contrastively
                    negativeTrainingData = dataManager.selectWorstVariants();
                    badModel = engine.generateModel(negativeTrainingData,
                            Math.min(VRAC.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL, VRAC.MAX_GAUSSIANS));
                    if (badModel.failedToConverge) {
                        throw new UserException.VQSRNegativeModelFailure(
                                "NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe." +
                                        " Please consider raising the number of variants used to train the negative model " +
                                        "(via --minimum-bad-variants 5000, for example).");
                    }
                }

                dataManager.dropAggregateData(); // Don't need the aggregate data anymore so let's free up the memory
                engine.evaluateData(dataManager.getData(), badModel, true);

                if (outputModel != null) {
                    final GATKReport report = writeModelReport(goodModel, badModel, USE_ANNOTATIONS);
                    saveModelReport(report, outputModel);
                }

                engine.calculateWorstPerformingAnnotation(dataManager.getData(), goodModel, badModel);


                // Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
                final int nCallsAtTruth = TrancheManager.countCallsAtTruth(dataManager.getData(), Double.NEGATIVE_INFINITY);
                final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric(nCallsAtTruth);
                if ( !scatterTranches ) {
                    final List<? extends Tranche> tranches = TrancheManager.findTranches(dataManager.getData(), TS_TRANCHES, metric, VRAC.MODE);
                    tranchesStream.print(TruthSensitivityTranche.printHeader());
                    tranchesStream.print(Tranche.tranchesString(tranches));
                }
                else {
                    final List<? extends Tranche> tranches = TrancheManager.findVQSLODTranches(dataManager.getData(), VQSLOD_TRANCHES, metric, VRAC.MODE);
                    tranchesStream.print(VQSLODTranche.printHeader());
                    tranchesStream.print(Tranche.tranchesString(tranches));
                }

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
                } else if (scatterTranches) {
                    //skip R plots for scattered tranches because the format is different and the R code parses them
                    logger.info("Tranches plot will not be generated since we are running in scattered mode");
                } else if (RSCRIPT_FILE != null) { //we don't use the RSCRIPT_FILE for tranches, but here it's an indicator if we're setup to run R
                    // Execute the RScript command to plot the table of truth values
                    rScriptExecutor.addScript(new Resource(PLOT_TRANCHES_RSCRIPT, VariantRecalibrator.class));
                    rScriptExecutor.addArgs(TRANCHES_FILE.getAbsoluteFile(), TARGET_TITV);
                    // Print out the command line to make it clear to the user what is being executed and how one might modify it
                    logger.info("Executing: " + rScriptExecutor.getApproximateCommandLine());
                    rScriptExecutor.exec();
                }
                return true;
            }
            catch (final Exception e) {
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

    /**
     * Rebuild a Gaussian Mixture Model from gaussian means and co-variates stored in a GATKReportTables
     * @param muTable           Table of Gaussian means
     * @param sigmaTable        Table of Gaussian co-variates
     * @param pmixTable         Table of PMixLog10 values
     * @param numAnnotations    Number of annotations, i.e. Dimension of the annotation space in which the Gaussians live
     * @return  a GaussianMixtureModel whose state reflects the state recorded in the tables.
     */
    protected GaussianMixtureModel GMMFromTables(final GATKReportTable muTable,
         final GATKReportTable sigmaTable,
         final GATKReportTable pmixTable,
         final int numAnnotations,
         final int numVariants){
         final List<MultivariateGaussian> gaussianList = new ArrayList<>();

        int curAnnotation = 0;
        for (final GATKReportColumn reportColumn : muTable.getColumnInfo() ) {
            if (!reportColumn.getColumnName().equals("Gaussian")) {
                for (int row = 0; row < muTable.getNumRows(); row++) {
                    if (gaussianList.size() <= row) {
                        final MultivariateGaussian mg = new MultivariateGaussian(numVariants, numAnnotations);
                        gaussianList.add(mg);
                    }
                    gaussianList.get(row).mu[curAnnotation] = (Double) muTable.get(row, reportColumn.getColumnName());
                }
                curAnnotation++;
            }
        }

        for (final GATKReportColumn reportColumn : pmixTable.getColumnInfo() ) {
            if (reportColumn.getColumnName().equals("pMixLog10")) {
                for (int row = 0; row < pmixTable.getNumRows(); row++) {
                    gaussianList.get(row).pMixtureLog10 = (Double) pmixTable.get(row, reportColumn.getColumnName());
                }
            }
        }

        int curJ = 0;
        for (final GATKReportColumn reportColumn : sigmaTable.getColumnInfo() ) {
            if (reportColumn.getColumnName().equals("Gaussian")) continue;
            if (reportColumn.getColumnName().equals("Annotation")) continue;

            for (int row = 0; row < sigmaTable.getNumRows(); row++) {
                final int curGaussian = row / numAnnotations;
                final int curI = row % numAnnotations;
                final double curVal = (Double) sigmaTable.get(row, reportColumn.getColumnName());
                gaussianList.get(curGaussian).sigma.set(curI, curJ, curVal);

            }
            curJ++;

        }

        return new GaussianMixtureModel(gaussianList, VRAC.SHRINKAGE, VRAC.DIRICHLET_PARAMETER, VRAC.PRIOR_COUNTS);

    }

    private Map<String, Double> getMapFromVectorTable(final GATKReportTable vectorTable) {
        final Map<String, Double> dataMap = new HashMap<>();

        //do a row-major traversal
        for (int i = 0; i < vectorTable.getNumRows(); i++) {
            dataMap.put((String) vectorTable.get(i, 0), (Double) vectorTable.get(i, 1));
        }
        return dataMap;
    }

    protected GATKReport writeModelReport(
            final GaussianMixtureModel goodModel,
            final GaussianMixtureModel badModel,
            final List<String> annotationList) {
        final String formatString = "%.16E";
        final GATKReport report = new GATKReport();

        if (dataManager != null) {  //for unit test
            final double[] meanVector = dataManager.getMeanVector();
            final GATKReportTable annotationMeans = makeVectorTable(
                    "AnnotationMeans",
                    "Mean for each annotation, used to normalize data",
                    dataManager.annotationKeys,
                    meanVector,
                    "Mean",
                    formatString);
            report.addTable(annotationMeans);

            final double[] varianceVector = dataManager.getVarianceVector();  //"varianceVector" is actually stdev
            final GATKReportTable annotationVariances = makeVectorTable(
                    "AnnotationStdevs",
                    "Standard deviation for each annotation, used to normalize data",
                    dataManager.annotationKeys,
                    varianceVector,
                    "StandardDeviation",  //column header must be one word
                    formatString);
            report.addTable(annotationVariances);
        }

        final List<String> gaussianStrings = new ArrayList<>();
        final double[] pMixtureLog10s = new double[goodModel.getModelGaussians().size()];
        int idx = 0;

        for( final MultivariateGaussian gaussian : goodModel.getModelGaussians() ) {
            pMixtureLog10s[idx] = gaussian.pMixtureLog10;
            gaussianStrings.add(Integer.toString(idx++) );
        }

        final GATKReportTable goodPMix = makeVectorTable("GoodGaussianPMix", "Pmixture log 10 used to evaluate model", gaussianStrings, pMixtureLog10s, "pMixLog10", formatString, "Gaussian");
        report.addTable(goodPMix);

        if (badModel != null) {
            gaussianStrings.clear();
            final double[] pMixtureLog10sBad = new double[badModel.getModelGaussians().size()];
            idx = 0;

            for (final MultivariateGaussian gaussian : badModel.getModelGaussians()) {
                pMixtureLog10sBad[idx] = gaussian.pMixtureLog10;
                gaussianStrings.add(Integer.toString(idx++));
            }
            final GATKReportTable badPMix = makeVectorTable("BadGaussianPMix", "Pmixture log 10 used to evaluate model", gaussianStrings, pMixtureLog10sBad, "pMixLog10", formatString, "Gaussian");
            report.addTable(badPMix);
        }

        //The model and Gaussians don't know what the annotations are, so get them from this class
        //VariantDataManager keeps the annotation in the same order as the argument list
        final GATKReportTable positiveMeans = makeMeansTable(
                "PositiveModelMeans", "Vector of annotation values to describe the (normalized) mean for each Gaussian in the positive model",
                annotationList,
                goodModel,
                formatString);
        report.addTable(positiveMeans);

        final GATKReportTable positiveCovariance = makeCovariancesTable(
                "PositiveModelCovariances", "Matrix to describe the (normalized) covariance for each Gaussian in the positive model; covariance matrices are joined by row",
                annotationList,
                goodModel,
                formatString);
        report.addTable(positiveCovariance);

        if (badModel != null) {
            //do the same for the negative model means
            final GATKReportTable negativeMeans = makeMeansTable(
                    "NegativeModelMeans", "Vector of annotation values to describe the (normalized) mean for each Gaussian in the negative model",
                    annotationList, badModel, formatString);
            report.addTable(negativeMeans);

            final GATKReportTable negativeCovariance = makeCovariancesTable(
                    "NegativeModelCovariances",
                    "Matrix to describe the (normalized) covariance for each Gaussian in the negative model; covariance matrices are joined by row",
                    annotationList,
                    badModel,
                    formatString);
            report.addTable(negativeCovariance);
        }

        return report;
    }

    private static void saveModelReport(final GATKReport report, final GATKPath modelPath) {
        try (final PrintStream modelReportStream = new PrintStream(modelPath.getOutputStream())) {
            if (modelReportStream.checkError()) {
                throw new IOException("checkError failure condition in PrintStream constructor");
            }
            report.print(modelReportStream);
            if (modelReportStream.checkError()) {
                throw new IOException("checkError failure condition in output writing");
            }
        } catch (final Exception e) {
            throw new UserException.CouldNotCreateOutputFile(modelPath, "Exception writing to report output. ", e);
        }
    }

    protected GATKReportTable makeVectorTable(final String tableName,
              final String tableDescription,
              final List<String> annotationList,
              final double[] perAnnotationValues,
              final String columnName,
              final String formatString) {
        return makeVectorTable(tableName,
                tableDescription,
                annotationList,
                perAnnotationValues,
                columnName,
                formatString,
                "Annotation");
    }

private GATKReportTable makeVectorTable(final String tableName,
            final String tableDescription,
            final List<String> annotationList,
            final double[] perAnnotationValues,
            final String columnName,
            final String formatString,
            final String firstColumn) {
        final GATKReportTable vectorTable = new GATKReportTable(tableName,
                tableDescription,
                annotationList.size(),
                GATKReportTable.Sorting.DO_NOT_SORT);
        vectorTable.addColumn(firstColumn, "");
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
        final GATKReportTable meansTable = new GATKReportTable(
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
        final GATKReportTable modelCovariances = new GATKReportTable(
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
        final PrintStream stream;
        try {
            stream = new PrintStream(RSCRIPT_FILE);
        } catch (final FileNotFoundException e ) {
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

                final List<VariantDatum> fakeData = new ArrayList<>();
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
        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(RSCRIPT_FILE);
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
