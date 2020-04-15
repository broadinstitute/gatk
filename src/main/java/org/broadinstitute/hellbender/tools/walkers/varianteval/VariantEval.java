package org.broadinstitute.hellbender.tools.walkers.varianteval;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.*;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.AlleleFrequency.*;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.IntervalStratification;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 *
 *
 * <p>
 * Given a variant callset, it is common to calculate various quality control metrics. These metrics include the number of
 * raw or filtered SNP counts; ratio of transition mutations to transversions; concordance of a particular sample's calls
 * to a genotyping chip; number of s per sample; etc. Furthermore, it is often useful to stratify these metrics
 * by various criteria like functional class (missense, nonsense, silent), whether the site is CpG site, the amino acid
 * degeneracy of the site, etc. VariantEval facilitates these calculations in two ways: by providing several built-in
 * evaluation and stratification modules, and by providing a framework that permits the easy development of new evaluation
 * and stratification modules.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more variant sets to evaluate plus any number of comparison sets.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Evaluation tables detailing the results of the eval modules which were applied.
 * For example:
 * </p>
 * <pre>
 * output.eval.grp:
 * ##:GATKReport.v0.1 CountVariants : Counts different classes of variants in the sample
 * CountVariants  CompFeatureInput  CpG      EvalFeatureInput  JexlExpression  Novelty  nProcessedLoci  nCalledLoci  nRefLoci  nVariantLoci  variantRate ...
 * CountVariants  dbsnp             CpG      eval              none            all      65900028        135770       0         135770        0.00206024  ...
 * CountVariants  dbsnp             CpG      eval              none            known    65900028        47068        0         47068         0.00071423  ...
 * CountVariants  dbsnp             CpG      eval              none            novel    65900028        88702        0         88702         0.00134601  ...
 * CountVariants  dbsnp             all      eval              none            all      65900028        330818       0         330818        0.00502000  ...
 * CountVariants  dbsnp             all      eval              none            known    65900028        120685       0         120685        0.00183133  ...
 * CountVariants  dbsnp             all      eval              none            novel    65900028        210133       0         210133        0.00318866  ...
 * CountVariants  dbsnp             non_CpG  eval              none            all      65900028        195048       0         195048        0.00295976  ...
 * CountVariants  dbsnp             non_CpG  eval              none            known    65900028        73617        0         73617         0.00111710  ...
 * CountVariants  dbsnp             non_CpG  eval              none            novel    65900028        121431       0         121431        0.00184265  ...
 * ...
 * </pre>
 * </p>
 *
 * <h3>Usage examples</h3>
 * <pre>
 * gatk VariantEval \
 *   -R reference.fasta \
 *   -O output.eval.grp \
 *   --eval set1:set1.vcf \
 *   --eval set2:set2.vcf \
 *   [--comp comp.vcf]
 * </pre>
 *
 * Count Mendelian violations for each family in a callset with multiple families (and provided pedigree)
 * <pre>
 * gatk VariantEval \
 *   -R reference.fasta \
 *   -O output.MVs.byFamily.table \
 *   --eval multiFamilyCallset.vcf \
 *   -no-ev -noST \
 *   -ST Family \
 *   -EV MendelianViolationEvaluator
 * </pre>
 *
 * <h3>Caveat</h3>
 *
 * <p>Some stratifications and evaluators are incompatible with each other due to their respective memory requirements,
 * such as AlleleCount and VariantSummary, or Sample and VariantSummary. If you specify such a combination, the program
 * will output an error message and ask you to disable one of these options. We do not currently provide an exhaustive
 * list of incompatible combinations, so we recommend trying out combinations that you are interested in on a dummy
 * command line, to rapidly ascertain whether it will work or not.</p>
 *
 */
@CommandLineProgramProperties(
        summary = "Given a variant callset, it is common to calculate various quality control metrics. These metrics include the number of " +
                "raw or filtered SNP counts; ratio of transition mutations to transversions; concordance of a particular sample's calls " +
                "to a genotyping chip; number of singletons per sample; etc. Furthermore, it is often useful to stratify these metrics " +
                "by various criteria like functional class (missense, nonsense, silent), whether the site is CpG site, the amino acid " +
                "degeneracy of the site, etc. VariantEval facilitates these calculations in two ways: by providing several built-in " +
                "evaluation and stratification modules, and by providing a framework that permits the easy development of new evaluation " +
                "and stratification modules.",
        oneLineSummary = "General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class VariantEval extends MultiVariantWalker {
    public static final String IS_SINGLETON_KEY = "ISSINGLETON";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    protected File outFile;

    /**
     * The variant file(s) to evaluate.
     */
    @Argument(fullName="eval", shortName = "eval", doc="Input evaluation file(s)", optional=false)
    public List<FeatureInput<VariantContext>> evals;

    /**
     * The variant file(s) to compare against.
     */
    @Argument(fullName = StandardArgumentDefinitions.COMPARISON_LONG_NAME, shortName = StandardArgumentDefinitions.COMPARISON_SHORT_NAME, doc="Input comparison file(s)", optional=true)
    public List<FeatureInput<VariantContext>> compsProvided = new ArrayList<>();
    private List<FeatureInput<VariantContext>> comps = new ArrayList<>();

    /**
     * dbSNP comparison VCF.  By default, the dbSNP file is used to specify the set of "known" variants.
     * Other sets can be specified with the -known-name (--known_names) argument.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * Some analyses want to count overlap not with dbSNP (which is in general very open) but
     * actually want to itemize their overlap specifically with a set of gold standard sites
     * such as HapMap, OMNI, or the gold standard indels.  This argument provides a mechanism
     * for communicating which file to use
     */
    @Argument(fullName="gold-standard", shortName = "gold", doc="Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison", optional=true)
    public FeatureInput<VariantContext> goldStandard = null;

    /**
     * Note that the --list argument requires a fully resolved and correct command-line to work.
     */
    @Argument(fullName="list", shortName="ls", doc="List the available eval modules and exit", optional=true)
    protected Boolean LIST = false;

    // Partitioning the data arguments
    @Argument(shortName="select", doc="One or more stratifications to use when evaluating the data", optional=true)
    protected ArrayList<String> SELECT_EXPS = new ArrayList<String>();

    @Argument(shortName="select-name", doc="Names to use for the list of stratifications (must be a 1-to-1 mapping)", optional=true)
    protected ArrayList<String> SELECT_NAMES = new ArrayList<String>();

    @Argument(fullName="sample", shortName="sn", doc="Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context", optional=true)
    protected Set<String> SAMPLE_EXPRESSIONS = new TreeSet<>();

    @Argument(fullName = StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName = StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for determining the population \"founders\"", optional=true)
    private GATKPath pedigreeFile;

    /**
     * List of feature tracks to be used for specifying "known" variants other than dbSNP.
     */
    @Argument(shortName="known-name", doc="Name of feature bindings containing variant sites that should be treated as known when splitting eval features into known and novel subsets", optional=true)
    protected Set<String> KNOWN_NAMES = new HashSet<String>();
    List<FeatureInput<VariantContext>> knowns = new ArrayList<>();

    // Stratification arguments
    @Argument(fullName="stratification-module", shortName="ST", doc="One or more specific stratification modules to apply to the eval track(s) (in addition to the standard stratifications, unless -noS is specified)", optional=true)
    protected List<String> STRATIFICATIONS_TO_USE = new ArrayList<>();

    @Argument(fullName="do-not-use-all-standard-stratifications", shortName="no-st", doc="Do not use the standard stratification modules by default (instead, only those that are specified with the -S option)", optional=true)
    protected Boolean NO_STANDARD_STRATIFICATIONS = false;

    /**
     * See the -list argument to view available modules.
     */
    @Argument(fullName="eval-module", shortName="EV", doc="One or more specific eval modules to apply to the eval track(s) (in addition to the standard modules, unless -no-ev is specified)", optional=true)
    protected List<String> MODULES_TO_USE = new ArrayList<>();

    @Argument(fullName="do-not-use-all-standard-modules", shortName="no-ev", doc="Do not use the standard modules by default (instead, only those that are specified with the -EV option)", optional=true)
    protected Boolean NO_STANDARD_MODULES = false;

    @Argument(fullName="min-phase-quality", shortName="mpq", doc="Minimum phasing quality", optional=true)
    protected double MIN_PHASE_QUALITY = 10.0;

    @Argument(shortName="mvq", fullName="mendelian-violation-qual-threshold", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation. Default is 50.", optional=true)
    protected double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 50;

    @Argument(shortName="ploidy", fullName="sample-ploidy", doc="Per-sample ploidy (number of chromosomes per sample)", optional=true)
    protected int ploidy = GATKVariantContextUtils.DEFAULT_PLOIDY;

    @Argument(fullName="ancestral-alignments", shortName="aa", doc="Fasta file with ancestral alleles", optional=true)
    private File ancestralAlignmentsFile = null;

    @Argument(fullName="require-strict-allele-match", shortName="strict", doc="If provided only comp and eval tracks with exactly matching reference and alternate alleles will be counted as overlapping", optional=true)
    private boolean requireStrictAlleleMatch = false;

    @Argument(fullName="keep-ac0", shortName="keep-ac0", doc="If provided, modules that track polymorphic sites will not require that a site have AC > 0 when the input eval has genotypes", optional=true)
    protected boolean keepSitesWithAC0 = false;

    @Hidden
    @Argument(fullName="num-samples", doc="If provided, modules that track polymorphic sites will not require that a site have AC > 0 when the input eval has genotypes", optional=true)
    private int numSamplesFromArgument = 0;

    /**
     * If true, VariantEval will treat -eval 1 -eval 2 as separate tracks from the same underlying
     * variant set, and evaluate the union of the results.  Useful when you want to do -eval chr1.vcf -eval chr2.vcf etc.
     */
    @Argument(fullName="merge-evals", shortName="merge-evals", doc="If provided, all -eval tracks will be merged into a single eval track", optional=true)
    public boolean mergeEvals = false;

    /**
     * File containing tribble-readable features for the IntervalStratificiation
     */
    @Argument(fullName="strat-intervals", shortName="strat-intervals", doc="File containing tribble-readable features for the IntervalStratificiation", optional=true)
    public FeatureInput<Feature> intervalsFile = null;

    /**
     * File containing tribble-readable features containing known CNVs.  For use with VariantSummary table.
     */
    @Argument(fullName="known-cnvs", shortName="known-cnvs", doc="File containing tribble-readable features describing a known list of copy number variants", optional=true)
    public FeatureInput<Feature> knownCNVsFile = null;

    protected StratifyingScale AFScale = StratifyingScale.LINEAR;
    protected boolean useCompAFStratifier = false;

    @Override
    protected MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MultiVariantInputArgumentCollection() {
            private static final long serialVersionUID = 1L;

            @Override
            public List<GATKPath> getDrivingVariantPaths() {
                //driving variants will be determined by initializeDrivingVariants()
                return Collections.emptyList();
            }
        };
    }

    @Override
    protected void initializeDrivingVariants() {
        getDrivingVariantsFeatureInputs().addAll(evals);
        if (dbsnp.dbsnp != null) {
            getDrivingVariantsFeatureInputs().add(dbsnp.dbsnp);
        }

        getDrivingVariantsFeatureInputs().addAll(compsProvided);

        super.initializeDrivingVariants();
    }

    // Variables
    private Set<SortableJexlVCMatchExp> jexlExpressions = new TreeSet<SortableJexlVCMatchExp>();

    private boolean isSubsettingSamples;
    private Set<String> sampleNamesForEvaluation = new LinkedHashSet<String>();
    private Set<String> familyNamesForEvaluation = new LinkedHashSet<String>();
    private Set<String> sampleNamesForStratification = new LinkedHashSet<String>();
    private Set<String> familyNamesForStratification = new LinkedHashSet<String>();

    // important stratifications
    private boolean byFilterIsEnabled = false;
    private boolean perSampleIsEnabled = false;
    private boolean perFamilyIsEnabled = false;

    // Public constants
    final private static String ALL_SAMPLE_NAME = "all";
    final private static String ALL_FAMILY_NAME = "all";

    // Utility class
    private final VariantEvalUtils variantEvalUtils = new VariantEvalUtils(this);

    // Ancestral alignments
    private ReferenceSequenceFile ancestralAlignments = null;

    // The set of all possible evaluation contexts
    StratificationManager<VariantStratifier, EvaluationContext> stratManager;

    private SampleDB sampleDB = null;

    // maintain the mapping of FeatureInput to name used in output file
    Map<FeatureInput<VariantContext>, String> inputToNameMap = new HashMap<>();

    /**
     * Initialize the stratifications, evaluations, evaluation contexts, and reporting object
     */
    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);

        // Just list the modules, and exit quickly.
        if (LIST) { variantEvalUtils.listModulesAndExit(); }

        sampleDB = SampleDB.createSampleDBFromPedigreeAndDataSources(pedigreeFile, getSamplesForVariants(), PedigreeValidationType.STRICT);

        comps.addAll(compsProvided);
        compsProvided.forEach(x -> inputToNameMap.put(x, x.hasUserSuppliedName() ? x.getName() : StandardArgumentDefinitions.COMPARISON_SHORT_NAME));
        if ( dbsnp.dbsnp != null ) {
            comps.add(dbsnp.dbsnp);
            inputToNameMap.put(dbsnp.dbsnp, "dbsnp");
            knowns.add(dbsnp.dbsnp);
        }

        evals.forEach(x -> inputToNameMap.put(x, x.hasUserSuppliedName() ? x.getName() : "eval"));

        // Set up set of additional knowns
        for ( FeatureInput<VariantContext> compInput : comps ) {
            if ( KNOWN_NAMES.contains(getNameForInput(compInput)))
                knowns.add(compInput);
        }

        // Now that we have all the inputs categorized, determine the sample list from the eval inputs.
        Map<String, VCFHeader> vcfInputs = new HashMap<>();
        evals.forEach(x -> vcfInputs.put(x.getName(), (VCFHeader)getHeaderForFeatures(x)));

        Set<String> vcfSamples = new HashSet<>();
        vcfInputs.forEach((k,v) -> vcfSamples.addAll(v.getSampleNamesInOrder()));

        // Load the sample list, using an intermediate tree set to sort the samples
        final Set<String> allSampleNames = new HashSet<>(vcfSamples);
        sampleNamesForEvaluation.addAll(new TreeSet<>(SAMPLE_EXPRESSIONS.isEmpty() ? vcfSamples : Utils.filterCollectionByExpressions(vcfSamples, SAMPLE_EXPRESSIONS, false)));

        isSubsettingSamples = ! sampleNamesForEvaluation.containsAll(allSampleNames);
        familyNamesForEvaluation.addAll(sampleDB.getFamilyIDs());

        //If stratifying by sample name, assign a stratification for each sample we're evaluating (based on commandline args)...
        if (STRATIFICATIONS_TO_USE.contains("Sample") ) {
            sampleNamesForStratification.addAll(sampleNamesForEvaluation);
        }
        //...and also a stratification for the sum over all samples
        sampleNamesForStratification.add(ALL_SAMPLE_NAME);

        //If stratifying by sample name, assign a stratification for each family...
        if ( STRATIFICATIONS_TO_USE.contains("Family") ) {
            familyNamesForStratification.addAll(familyNamesForEvaluation);
        }
        //...and also a stratification for the sum over all families
        familyNamesForStratification.add(ALL_FAMILY_NAME);

        // Initialize select expressions
        for (VariantContextUtils.JexlVCMatchExp jexl : VariantContextUtils.initializeMatchExps(SELECT_NAMES, SELECT_EXPS)) {
            SortableJexlVCMatchExp sjexl = new SortableJexlVCMatchExp(jexl.name, jexl.exp);
            jexlExpressions.add(sjexl);
        }

        // Initialize the set of stratifications and evaluations to use
        // The list of stratifiers and evaluators to use
        final List<VariantStratifier> stratificationObjects = variantEvalUtils.initializeStratificationObjects(NO_STANDARD_STRATIFICATIONS, STRATIFICATIONS_TO_USE);
        final Set<Class<? extends VariantEvaluator>> evaluationClasses = variantEvalUtils.initializeEvaluationObjects(NO_STANDARD_MODULES, MODULES_TO_USE);

        checkForIncompatibleEvaluatorsAndStratifiers(stratificationObjects, evaluationClasses);

        for ( VariantStratifier vs : stratificationObjects ) {
            if ( vs.getName().equals("Filter") )
                byFilterIsEnabled = true;
            else if ( vs.getName().equals("Sample") )
                perSampleIsEnabled = true;
            else if ( vs.getName().equals("Family"))
                perFamilyIsEnabled = true;
        }

        if (perSampleIsEnabled && perFamilyIsEnabled)
            throw new CommandLineException.BadArgumentValue("ST", "Variants cannot be stratified by sample and family at the same time");

        if (perFamilyIsEnabled && sampleDB.getTrios().isEmpty())
            throw new CommandLineException.BadArgumentValue("ST", "Cannot stratify by family without *.ped file");


        if ( intervalsFile != null ) {
            boolean fail = true;
            for ( final VariantStratifier vs : stratificationObjects ) {
                if ( vs.getClass().equals(IntervalStratification.class) )
                    fail = false;
            }
            if ( fail )
                throw new CommandLineException.BadArgumentValue("ST", "stratIntervals argument provided but -ST IntervalStratification not provided");
        }

        // Initialize the evaluation contexts
        createStratificationStates(stratificationObjects, evaluationClasses);

        // Load ancestral alignments
        if (ancestralAlignmentsFile != null) {
            try {
                ancestralAlignments = new IndexedFastaSequenceFile(ancestralAlignmentsFile.toPath());
            } catch (FileNotFoundException e) {
                throw new GATKException(String.format("The ancestral alignments file, '%s', could not be found", ancestralAlignmentsFile.getAbsolutePath()));
            }
        }

        assertThatTerritoryIsSpecifiedIfNecessary();
    }

    private void assertThatTerritoryIsSpecifiedIfNecessary() {
        final Set<String> evaluatorsWhichRequireTerritory = stratManager.values()
                .stream()
                .flatMap(ctx -> ctx.getVariantEvaluators().stream())
                .filter(Objects::nonNull)
                .filter(VariantEvaluator::requiresTerritoryToBeSpecified)
                .map(VariantEvaluator::getSimpleName)
                .collect(Collectors.toSet());
        if(!evaluatorsWhichRequireTerritory.isEmpty() && getTraversalIntervals() == null){
            throw new UserException("You specified evaluators which require a covered territory to be specified.  " +
                    "\nPlease specify intervals or a reference file or disable all of the following evaluators:" +
                    evaluatorsWhichRequireTerritory.stream()
                            .collect(Collectors.joining(", ")));
        }
    }

    private void checkForIncompatibleEvaluatorsAndStratifiers( final List<VariantStratifier> stratificationObjects,
                                                             Set<Class<? extends VariantEvaluator>> evaluationClasses) {
        for ( final VariantStratifier vs : stratificationObjects ) {
            for ( Class<? extends VariantEvaluator> ec : evaluationClasses )
                if ( vs.getIncompatibleEvaluators().contains(ec) )
                    throw new CommandLineException.BadArgumentValue("ST and ET",
                            "The selected stratification " + vs.getName() +
                                    " and evaluator " + ec.getSimpleName() +
                                    " are incompatible due to combinatorial memory requirements." +
                                    " Please disable one");
        }
    }

    final void createStratificationStates(final List<VariantStratifier> stratificationObjects, final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        final List<VariantStratifier> strats = new ArrayList<VariantStratifier>(stratificationObjects);
        stratManager = new StratificationManager<>(strats);

        logger.info("Creating " + stratManager.size() + " combinatorial stratification states");
        for ( int i = 0; i < stratManager.size(); i++ ) {
            EvaluationContext ec = createEvaluationContext(evaluationObjects);
            stratManager.set(i, ec);
        }
    }

    /**
     * Create the EvaluationContext (new instance) for the provided set of VariantEvaluators.
     *
     * @param evaluationObjects The list of VariantEvaluator classes
     * @return The EvaluationContext for this set of VariantEvaluator classes
     */
    protected EvaluationContext createEvaluationContext(final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        return new EvaluationContext(this, evaluationObjects);
    }

    private class PositionAggregator {
        private SimpleInterval i = null;

        private ReadsContext readsContext;
        private ReferenceContext referenceContext;
        private FeatureContext featureContext;

        private void addVariant(VariantContext vc, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
            if (i == null || !vc.getContig().equals(i.getContig()) || vc.getStart() != i.getStart()) {
                callDoApply();

                i = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd());
                this.readsContext = readsContext;
                this.referenceContext = referenceContext;
                this.featureContext = featureContext;
            }
            else if (vc.getEnd() > i.getEnd()) {
                //expand region
                i = new SimpleInterval(i.getContig(), i.getStart(), vc.getEnd());

                this.readsContext = new ReadsContext(this.readsContext, i);
                this.referenceContext = new ReferenceContext(this.referenceContext, i);
                this.featureContext = new FeatureContext(this.featureContext, i);
            }
        }

        public void callDoApply(){
            if (i != null) {
                doApply(this.readsContext, this.referenceContext, this.featureContext);
                i = null;
            }
        }

        public void onComplete() {
            callDoApply();
        }
    }

    final PositionAggregator aggr = new PositionAggregator();

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        aggr.addVariant(variant, readsContext, referenceContext, featureContext);
    }

    public String getNameForInput(FeatureInput<VariantContext> input) {
        return inputToNameMap.get(input);
    }

    /**
     * This will get called once per site where a variant is present in any input
     */
    public void doApply(ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        HashMap<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> evalVCs = variantEvalUtils.bindVariantContexts(referenceContext, featureContext, evals, byFilterIsEnabled, true, perSampleIsEnabled, perFamilyIsEnabled, mergeEvals);
        HashMap<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> compVCs = variantEvalUtils.bindVariantContexts(referenceContext, featureContext, comps, byFilterIsEnabled, false, false, false, false);

        // for each eval track
        for ( final FeatureInput<VariantContext> evalInput : evals ) {
            final Map<String, Collection<VariantContext>> emptyEvalMap = Collections.emptyMap();
            final Map<String, Collection<VariantContext>> evalSet = evalVCs.containsKey(evalInput) ? evalVCs.get(evalInput) : emptyEvalMap;

            Set<String> statificationLevels;

            // for each sample stratifier
            if (perFamilyIsEnabled)
                statificationLevels = familyNamesForStratification;
            else
                statificationLevels = sampleNamesForStratification;
            for ( final String stratLevelName : statificationLevels ) {
                Collection<VariantContext> evalSetBySample = evalSet.get(stratLevelName);

                if ( evalSetBySample == null ) {
                    evalSetBySample = new HashSet<>(1);
                    evalSetBySample.add(null);
                }

                // for each eval in the track
                for ( VariantContext eval : evalSetBySample ) {
                    String aastr = (ancestralAlignments == null) ? null : new String(ancestralAlignments.getSubsequenceAt(eval.getContig(), eval.getStart(), eval.getEnd()).getBases());

                    // deal with ancestral alleles if requested
                    if ( eval != null && aastr != null ) {
                        eval = new VariantContextBuilder(eval).attribute("ANCESTRALALLELE", aastr).make();
                    }

                    String evalName = getNameForInput(evalInput);

                    // for each comp track
                    for ( final FeatureInput<VariantContext> compInput : comps ) {
                        processComp(referenceContext, readsContext, featureContext, eval, evalName, compInput, stratLevelName, compVCs, evalSetBySample);
                    }

                    if (comps.isEmpty()) {
                        processComp(referenceContext, readsContext, featureContext, eval, evalName, null, stratLevelName, compVCs, evalSetBySample);
                    }
                }
            }

            if ( mergeEvals ) break; // stop processing the eval tracks
        }
    }

    private void processComp(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext eval, String evalName, FeatureInput<VariantContext> compInput, String stratLevelName, HashMap<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> compVCs, Collection<VariantContext> evalSetBySample) {
        String compName = getNameForInput(compInput);

        // no sample stratification for comps
        final HashMap<String, Collection<VariantContext>> compSetHash = compInput == null ? null : compVCs.get(compInput);
        final Collection<VariantContext> compSet = (compSetHash == null || compSetHash.isEmpty()) ? Collections.<VariantContext>emptyList() : compVCs.get(compInput).values().iterator().next();

        // find the comp
        final VariantContext comp = findMatchingComp(eval, compSet);

        Collection<EvaluationContext> contextsForStratification;
        if (perFamilyIsEnabled)
            contextsForStratification = getEvaluationContexts(referenceContext, readsContext, featureContext, eval, evalName, comp, compName, null, stratLevelName);
        else {
            String familyID;
            if (stratLevelName.equals("all"))
                familyID = "all";
            else
                familyID = sampleDB.getSample(stratLevelName).getFamilyID();
            contextsForStratification = getEvaluationContexts(referenceContext, readsContext, featureContext, eval, evalName, comp, compName, stratLevelName, familyID);
        }
        for ( EvaluationContext nec : contextsForStratification ) {

            // eval against the comp
            synchronized (nec) {
                nec.apply(referenceContext, readsContext, featureContext, comp, eval);
            }

            // eval=null against all comps of different type that aren't bound to another eval
            for ( VariantContext otherComp : compSet ) {
                if ( otherComp != comp && ! compHasMatchingEval(otherComp, evalSetBySample) ) {
                    synchronized (nec) {
                        nec.apply(referenceContext, readsContext, featureContext, otherComp, null);
                    }
                }
            }
        }
    }

    /**
     * Given specific eval and comp VCs and the sample name, return an iterable
     * over all of the applicable state keys.
     *
     * this code isn't structured yet for efficiency.  Here we currently are
     * doing the following inefficient algorithm:
     *
     * for each strat:
     *   get list of relevant states that eval and comp according to strat
     *   add this list of states to a list of list states
     *
     * then
     *
     * ask the strat manager to look up all of the keys associated with the combinations
     * of these states.  For example, suppose we have a single variant S.  We have active
     * strats EvalFeatureInput, CompFeatureInput, and Novelty.  We produce a list that looks like:
     *
     *   L = [[Eval], [Comp], [All, Novel]]
     *
     * We then go through the strat manager tree to produce the keys associated with these states:
     *
     *   K = [0, 1] where EVAL x COMP x ALL = 0 and EVAL x COMP x NOVEL = 1
     *
     * It's clear that a better
     *
     *
     * @param referenceContext
     * @param readsContext
     * @param featureContext
     * @param eval
     * @param evalName
     * @param comp
     * @param compName
     * @param sampleName
     * @return
     */
    protected Collection<EvaluationContext> getEvaluationContexts(final ReferenceContext referenceContext,
                                                                  final ReadsContext readsContext,
                                                                  final FeatureContext featureContext,
                                                                  final VariantContext eval,
                                                                  final String evalName,
                                                                  final VariantContext comp,
                                                                  final String compName,
                                                                  final String sampleName,
                                                                  final String familyName) {
        final List<List<Object>> states = new LinkedList<List<Object>>();
        for ( final VariantStratifier vs : stratManager.getStratifiers() ) {
            states.add(vs.getRelevantStates(referenceContext, readsContext, featureContext, comp, compName, eval, evalName, sampleName, familyName));
        }
        return stratManager.values(states);
    }


    private boolean compHasMatchingEval(final VariantContext comp, final Collection<VariantContext> evals) {
        // find all of the matching comps
        for ( final VariantContext eval : evals ) {
            if ( eval != null && doEvalAndCompMatch(comp, eval, requireStrictAlleleMatch) != EvalCompMatchType.NO_MATCH )
                return true;
        }

        // nothing matched
        return false;
    }

    private enum EvalCompMatchType { NO_MATCH, STRICT, LENIENT }

    private EvalCompMatchType doEvalAndCompMatch(final VariantContext eval, final VariantContext comp, boolean requireStrictAlleleMatch) {
        if ( comp.getType() == VariantContext.Type.NO_VARIATION || eval.getType() == VariantContext.Type.NO_VARIATION )
            // if either of these are NO_VARIATION they are LENIENT matches
            return EvalCompMatchType.LENIENT;

        if ( comp.getType() != eval.getType() )
            return EvalCompMatchType.NO_MATCH;

        // find the comp which matches both the reference allele and alternate allele from eval
        final Allele altEval = eval.getAlternateAlleles().size() == 0 ? null : eval.getAlternateAllele(0);
        final Allele altComp = comp.getAlternateAlleles().size() == 0 ? null : comp.getAlternateAllele(0);
        if ((altEval == null && altComp == null) || (altEval != null && altEval.equals(altComp) && eval.getReference().equals(comp.getReference())))
            return EvalCompMatchType.STRICT;
        else
            return requireStrictAlleleMatch ? EvalCompMatchType.NO_MATCH : EvalCompMatchType.LENIENT;
    }

    private VariantContext findMatchingComp(final VariantContext eval, final Collection<VariantContext> comps) {
        // if no comps, return null
        if ( comps == null || comps.isEmpty() )
            return null;

        // if no eval, return any comp
        if ( eval == null )
            return comps.iterator().next();

        // find all of the matching comps
        VariantContext lenientMatch = null;
        for ( final VariantContext comp : comps ) {
            switch ( doEvalAndCompMatch(comp, eval, requireStrictAlleleMatch) ) {
                case STRICT:
                    return comp;
                case LENIENT:
                    if ( lenientMatch == null ) lenientMatch = comp;
                    break;
                case NO_MATCH:
                    // do nothing
            }
        }

        // nothing matched, just return lenientMatch, which might be null
        return lenientMatch;
    }

    @Override
    public Object onTraversalSuccess() {
        aggr.onComplete();

        logger.info("Finalizing variant report");
        
        // go through the evaluations and finalize them
        for ( final EvaluationContext nec : stratManager.values() )
            for ( final VariantEvaluator ve : nec.getVariantEvaluators() )
                ve.finalizeEvaluation();

        //send data to MetricsCollection
        CompOverlap compOverlap = null;
        IndelSummary indelSummary = null;
        CountVariants countVariants = null;
        MultiallelicSummary multiallelicSummary = null;
        TiTvVariantEvaluator tiTvVariantEvaluator = null;
        MetricsCollection metricsCollection = null;
        for(final EvaluationContext nec: stratManager.values()) {
            for(final VariantEvaluator ve : nec.getVariantEvaluators()) {
                if (ve instanceof CompOverlap)
                    compOverlap = (CompOverlap) ve;
                else if (ve instanceof IndelSummary)
                    indelSummary = (IndelSummary) ve;
                else if (ve instanceof CountVariants)
                    countVariants = (CountVariants) ve;
                else if (ve instanceof MultiallelicSummary)
                    multiallelicSummary = (MultiallelicSummary) ve;
                else if (ve instanceof TiTvVariantEvaluator)
                    tiTvVariantEvaluator = (TiTvVariantEvaluator) ve;
                else if (ve instanceof MetricsCollection)
                    metricsCollection = (MetricsCollection) ve;
            }

            if(metricsCollection != null)
                metricsCollection.setData(compOverlap.concordantRate, indelSummary.n_SNPs, countVariants.nSNPs, indelSummary.n_indels, multiallelicSummary.nIndels, indelSummary.insertion_to_deletion_ratio, countVariants.insertionDeletionRatio, tiTvVariantEvaluator.tiTvRatio);
        }

        try (PrintStream out = IOUtils.makePrintStreamMaybeGzipped(new GATKPath(outFile.getAbsolutePath()))) {
            VariantEvalReportWriter.writeReport(out, stratManager, stratManager.getStratifiers(), stratManager.get(0).getVariantEvaluators());
        }
        catch(IOException e) {
            throw new UserException.CouldNotCreateOutputFile(e.getMessage(), e);
        }

        return null;
    }

    // Accessors
    public Logger getLogger() { return logger; }

    public double getMinPhaseQuality() { return MIN_PHASE_QUALITY; }

    public int getSamplePloidy() { return ploidy; }
    public double getMendelianViolationQualThreshold() { return MENDELIAN_VIOLATION_QUAL_THRESHOLD; }

    public static String getAllSampleName() { return ALL_SAMPLE_NAME; }
    public static String getAllFamilyName() { return ALL_FAMILY_NAME; }

    public List<FeatureInput<VariantContext>> getKnowns() { return knowns; }

    public List<FeatureInput<VariantContext>> getEvals() { return evals; }

    public boolean isSubsettingToSpecificSamples() { return isSubsettingSamples; }
    public Set<String> getSampleNamesForEvaluation() { return sampleNamesForEvaluation; }

    public Set<String> getFamilyNamesForEvaluation() { return familyNamesForEvaluation; }

    public int getNumberOfSamplesForEvaluation() {
        if (sampleNamesForEvaluation!= null &&  !sampleNamesForEvaluation.isEmpty())
            return sampleNamesForEvaluation.size();
        else {
            return numSamplesFromArgument;
        }

    }
    public Set<String> getSampleNamesForStratification() { return sampleNamesForStratification; }

    public Set<String> getFamilyNamesForStratification() { return familyNamesForStratification; }

    public List<FeatureInput<VariantContext>> getComps() { return comps; }

    public Set<SortableJexlVCMatchExp> getJexlExpressions() { return jexlExpressions; }


    public StratifyingScale getAFScale() { return AFScale; }
    public boolean getCompAFStratifier() { return useCompAFStratifier; }


    public Set<String> getContigNames() {
        final TreeSet<String> contigs = new TreeSet<>();
        for( final SAMSequenceRecord r :  getSequenceDictionaryForDrivingVariants().getSequences()) {
            contigs.add(r.getSequenceName());
        }
        return contigs;
    }

    public boolean ignoreAC0Sites() {
        return ! keepSitesWithAC0;
    }

    public SampleDB getSampleDB() {
        return sampleDB;
    }

    /**
     * If an evaluator calls this method it must override {@link VariantEvaluator#requiresTerritoryToBeSpecified()} to return true.
     * @return either the size of the interval list given to the tool or the size of the reference given to the tool
     */
    public long getnProcessedLoci() {
        if(getTraversalIntervals() == null){
            throw new GATKException("BUG: One of the evaluators used should have overriden requiresTerritoryToBeSpecified, please report this to the developers." +
                    "\nEvaluators: " + stratManager.values()
                    .stream()
                    .flatMap(evaluator -> evaluator.getVariantEvaluators().stream())
                    .map(VariantEvaluator::getSimpleName)
                    .sorted()
                    .distinct()
                    .collect(Collectors.joining(", ")));
        }
        return getTraversalIntervals().stream().mapToLong(SimpleInterval::size).sum();
    }

    public FeatureInput<Feature> getKnownCNVsFile() {
        return knownCNVsFile;
    }

    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        // Reason this override is needed:
        // Because getDefaultReadFilters() is nearly always included (see GATKTool), this pulls in SampleReadFilter,
        // which has an argument named "sample", which conflicts with the local argument.  we get an exception like:
        // org.broadinstitute.barclay.argparser.CommandLineException$CommandLineParserInternalException: [sample, sample] has already been used.
        // at org.broadinstitute.barclay.argparser.CommandLineArgumentParser.handleArgumentAnnotation(CommandLineArgumentParser.java:1002)
        // therefore override the method and return an empty list.
        return Collections.emptyList();
    }
}
