package org.broadinstitute.hellbender.tools.walkers.varianteval;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.*;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.*;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.reflections.Reflections;

import javax.annotation.Nullable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This class allows other classes to replicate the behavior of VariantEval
 *
 * Usage:
 * -Pass the genotype args into the constructor, which will the initialize the engine completely
 */
public class VariantEvalEngine {
    public static final String IS_SINGLETON_KEY = "ISSINGLETON";

    private final VariantEvalArgumentCollection variantEvalArgs;

    private final Logger logger = LogManager.getLogger(VariantEvalEngine.class);

    private final SAMSequenceDictionary samSequenceDictionaryForDrivingVariants;
    private final List<SimpleInterval> traversalIntervals;
    private final FeatureManager features;

    private final static Map<String, Class<? extends VariantStratifier>> stratifierClasses;
    private final static Set<String> standardStratificationNames;
    private final static Set<String> requiredStratificationNames;

    private final static Map<String, Class<? extends VariantEvaluator>> evaluatorClasses;
    private final static Set<String> standardEvaluatorNames;

    static {
        stratifierClasses = new HashMap<>();
        standardStratificationNames = new HashSet<>();
        requiredStratificationNames = new HashSet<>();

        Reflections reflectionsStrat = new Reflections(VariantStratifier.class.getPackage().getName());
        Set<Class<? extends VariantStratifier>> allClasses = reflectionsStrat.getSubTypesOf(VariantStratifier.class);
        for (Class<? extends VariantStratifier> clazz : allClasses) {
            stratifierClasses.put(clazz.getSimpleName(), clazz);

            if (StandardStratification.class.isAssignableFrom(clazz)) {
                standardStratificationNames.add(clazz.getSimpleName());
            }

            if (RequiredStratification.class.isAssignableFrom(clazz)) {
                requiredStratificationNames.add(clazz.getSimpleName());
            }
        }

        evaluatorClasses = new HashMap<>();
        standardEvaluatorNames= new HashSet<>();

        Reflections reflectionsEval = new Reflections(VariantEvaluator.class.getPackage().getName());
        Set<Class<? extends VariantEvaluator>> allEvalClasses = reflectionsEval.getSubTypesOf(VariantEvaluator.class);
        for (Class<? extends VariantEvaluator> clazz : allEvalClasses) {
            evaluatorClasses.put(clazz.getSimpleName(), clazz);

            if (StandardEval.class.isAssignableFrom(clazz)) {
                standardEvaluatorNames.add(clazz.getSimpleName());
            }
        }
    }

    // Ancestral alignments
    private ReferenceSequenceFile ancestralAlignments = null;

    // The set of all possible evaluation contexts
    private StratificationManager<VariantStratifier, EvaluationContext> stratManager;

    private SampleDB sampleDB = null;

    private List<FeatureInput<VariantContext>> knowns = new ArrayList<>();

    // maintain the mapping of FeatureInput to name used in output file
    private Map<FeatureInput<VariantContext>, String> inputToNameMap = new HashMap<>();

    // Variables
    private Set<SortableJexlVCMatchExp> jexlExpressions = new TreeSet<>();

    private boolean isSubsettingSamples;
    private Set<String> sampleNamesForEvaluation = new LinkedHashSet<>();
    private Set<String> familyNamesForEvaluation = new LinkedHashSet<>();
    private Set<String> sampleNamesForStratification = new LinkedHashSet<>();
    private Set<String> familyNamesForStratification = new LinkedHashSet<>();

    // important stratifications
    private boolean byFilterIsEnabled = false;
    private boolean perSampleIsEnabled = false;
    private boolean perFamilyIsEnabled = false;

    private AlleleFrequency.StratifyingScale AFScale = AlleleFrequency.StratifyingScale.LINEAR;
    private boolean useCompAFStratifier = false;

    // maintain the mapping of source name (from VC) to FeatureInput name
    private Map<String, FeatureInput<VariantContext>> drivingVariantSourceMap;

    // No args constructor for unit testing only
    @VisibleForTesting
    protected VariantEvalEngine() {
        this.variantEvalArgs = new VariantEvalArgumentCollection();
        this.samSequenceDictionaryForDrivingVariants = null;
        this.traversalIntervals = null;
        this.features = null;
    }

    public VariantEvalEngine(VariantEvalArgumentCollection variantEvalArgs, FeatureManager features, List<SimpleInterval> traversalIntervals, SAMSequenceDictionary samSequenceDictionaryForDrivingVariants, @Nullable Collection<String> samples) {
        this.variantEvalArgs = variantEvalArgs;
        this.samSequenceDictionaryForDrivingVariants = samSequenceDictionaryForDrivingVariants;
        this.traversalIntervals = traversalIntervals;
        this.features = features;

        // Cache map of source name -> FeatureInput
        drivingVariantSourceMap = new HashMap<>();
        variantEvalArgs.getFeatureInputsForDrivingVariants().forEach(x -> drivingVariantSourceMap.put(x.getName(), x));

        validateAndInitialize(samples);
    }

    /**
     * Initialize the stratifications, evaluations, evaluation contexts, and reporting object
     */
    protected void validateAndInitialize(@Nullable Collection<String> samples) {
        sampleDB = SampleDB.createSampleDBFromPedigreeAndDataSources(variantEvalArgs.pedigreeFile, samples, PedigreeValidationType.STRICT);

        variantEvalArgs.comps.addAll(variantEvalArgs.compsProvided);
        variantEvalArgs.compsProvided.forEach(comp -> inputToNameMap.put(comp, comp.hasUserSuppliedName() ? comp.getName() : StandardArgumentDefinitions.COMPARISON_SHORT_NAME));
        if ( variantEvalArgs.dbsnp.dbsnp != null ) {
            variantEvalArgs.comps.add(variantEvalArgs.dbsnp.dbsnp);
            inputToNameMap.put(variantEvalArgs.dbsnp.dbsnp, "dbsnp");
            knowns.add(variantEvalArgs.dbsnp.dbsnp);
        }

        variantEvalArgs.evals.forEach(eval -> inputToNameMap.put(eval, eval.hasUserSuppliedName() ? eval.getName() : "eval"));

        // Set up set of additional knowns. dbSNP was addressed above, so use compsProvided, not comps
        for ( FeatureInput<VariantContext> compInput : variantEvalArgs.compsProvided ) {
            if (variantEvalArgs.knownNames.contains(getNameForInput(compInput)))
                knowns.add(compInput);
        }

        // Now that we have all the inputs categorized, determine the sample list from the eval inputs.
        Map<String, VCFHeader> vcfInputs = new HashMap<>();
        variantEvalArgs.evals.forEach(eval -> vcfInputs.put(eval.getName(), (VCFHeader)features.getHeader(eval)));

        Set<String> vcfSamples = new HashSet<>();
        vcfInputs.forEach((k,v) -> vcfSamples.addAll(v.getSampleNamesInOrder()));

        // Load the sample list, using an intermediate tree set to sort the samples
        final Set<String> allSampleNames = new HashSet<>(vcfSamples);
        sampleNamesForEvaluation.addAll(new TreeSet<>(variantEvalArgs.sampleExpressions.isEmpty() ? vcfSamples : Utils.filterCollectionByExpressions(vcfSamples, variantEvalArgs.sampleExpressions, false)));

        isSubsettingSamples = ! sampleNamesForEvaluation.containsAll(allSampleNames);
        familyNamesForEvaluation.addAll(sampleDB.getFamilyIDs());

        //If stratifying by sample name, assign a stratification for each sample we're evaluating (based on commandline args)...
        if (variantEvalArgs.stratificationsToUse.contains("Sample") ) {
            sampleNamesForStratification.addAll(sampleNamesForEvaluation);
        }
        //...and also a stratification for the sum over all samples
        sampleNamesForStratification.add(VariantEvalArgumentCollection.ALL_SAMPLE_NAME);


        //If stratifying by sample name, assign a stratification for each family...
        if ( variantEvalArgs.stratificationsToUse.contains("Family") ) {
            familyNamesForStratification.addAll(familyNamesForEvaluation);
        }
        //...and also a stratification for the sum over all families
        familyNamesForStratification.add(VariantEvalArgumentCollection.ALL_FAMILY_NAME);

        // Initialize select expressions
        for (VariantContextUtils.JexlVCMatchExp jexl : VariantContextUtils.initializeMatchExps(variantEvalArgs.selectNames, variantEvalArgs.selectExps)) {
            SortableJexlVCMatchExp sjexl = new SortableJexlVCMatchExp(jexl.name, jexl.exp);
            jexlExpressions.add(sjexl);
        }

        // Initialize the set of stratifications and evaluations to use
        // The list of stratifiers and evaluators to use
        final List<VariantStratifier> stratificationObjects = initializeStratificationObjects(variantEvalArgs.noStandardStratifications, variantEvalArgs.stratificationsToUse);
        final Set<Class<? extends VariantEvaluator>> evaluationClasses = initializeEvaluationObjects(variantEvalArgs.noStandardModules, variantEvalArgs.modulesToUse);

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


        if ( variantEvalArgs.intervalsFile != null ) {
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
        if (variantEvalArgs.ancestralAlignmentsFile != null) {
            try {
                ancestralAlignments = new IndexedFastaSequenceFile(variantEvalArgs.ancestralAlignmentsFile.toPath());
            } catch (FileNotFoundException e) {
                throw new GATKException(String.format("The ancestral alignments file, '%s', could not be found", variantEvalArgs.ancestralAlignmentsFile.getAbsolutePath()));
            }
        }

        assertThatTerritoryIsSpecifiedIfNecessary();
    }

    public String getNameForInput(FeatureInput<VariantContext> input) {
        return inputToNameMap.get(input);
    }

    private void assertThatTerritoryIsSpecifiedIfNecessary() {
        final Set<String> evaluatorsWhichRequireTerritory = stratManager.values()
                .stream()
                .flatMap(ctx -> ctx.getVariantEvaluators().stream())
                .filter(Objects::nonNull)
                .filter(VariantEvaluator::requiresTerritoryToBeSpecified)
                .map(VariantEvaluator::getSimpleName)
                .collect(Collectors.toSet());
        if (!evaluatorsWhichRequireTerritory.isEmpty() && traversalIntervals == null){
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

    protected void createStratificationStates(final List<VariantStratifier> stratificationObjects, final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        final List<VariantStratifier> strats = new ArrayList<VariantStratifier>(stratificationObjects);
        stratManager = new StratificationManager<>(strats);

        logger.info("Creating " + stratManager.size() + " combinatorial stratification states");
        for ( int i = 0; i < stratManager.size(); i++ ) {
            EvaluationContext ec = createEvaluationContext(evaluationObjects);
            stratManager.set(i, ec);
        }
    }

    public void finalizeReport(File outFile) {
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
        for (final EvaluationContext nec: stratManager.values()) {
            for (final VariantEvaluator ve : nec.getVariantEvaluators()) {
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

            if (metricsCollection != null)
                metricsCollection.setData(compOverlap.concordantRate, indelSummary.n_SNPs, countVariants.nSNPs, indelSummary.n_indels, multiallelicSummary.nIndels, indelSummary.insertion_to_deletion_ratio, countVariants.insertionDeletionRatio, tiTvVariantEvaluator.tiTvRatio);
        }

        try (PrintStream out = IOUtils.makePrintStreamMaybeGzipped(new GATKPath(outFile.getAbsolutePath()))) {
            VariantEvalReportWriter.writeReport(out, stratManager, stratManager.getStratifiers(), stratManager.get(0).getVariantEvaluators());
        }
        catch(IOException e) {
            throw new UserException.CouldNotCreateOutputFile(e.getMessage(), e);
        }
    }

    private Map<FeatureInput<VariantContext>, List<VariantContext>> groupVariantsByFeatureInput(final List<VariantContext> variants) {
        final Map<FeatureInput<VariantContext>, List<VariantContext>> byFeatureInput = new HashMap<>();
        variants.forEach(vc -> byFeatureInput.compute(drivingVariantSourceMap.get(vc.getSource()),
                (k, v) -> {
                    final List<VariantContext> variantList = v == null ? new ArrayList<>() : v;
                    variantList.add(vc);
                    return variantList;
                }
        ));
        return byFeatureInput;
    }

    public void apply(final List<VariantContext> variantContexts, final ReferenceContext referenceContext) {
        final Map<FeatureInput<VariantContext>, List<VariantContext>> variantMap = groupVariantsByFeatureInput(variantContexts);

        final List<VariantContext> allEvals = new ArrayList<>();
        for (FeatureInput<VariantContext> eval : variantEvalArgs.evals) {
            if (variantMap.containsKey(eval)) {
                allEvals.addAll(variantMap.get(eval));
            }
        }

        final List<VariantContext> allComps = new ArrayList<>();
        if (variantEvalArgs.comps != null) {
            for (FeatureInput<VariantContext> comp : variantEvalArgs.comps) {
                if (variantMap.containsKey(comp)) {
                    allComps.addAll(variantMap.get(comp));
                }
            }
        }

        final SimpleInterval interval = allEvals.isEmpty() ? new SimpleInterval(variantContexts.get(0).getContig(), variantContexts.get(0).getStart(), variantContexts.get(0).getStart()) : generateContextInterval(allEvals);
        final FeatureContext featureContext = new FeatureContext(features, interval);

        final Map<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> evalVCs = allEvals.isEmpty() ? Collections.emptyMap() : bindVariantContexts(variantMap, variantEvalArgs.evals, byFilterIsEnabled, true, perSampleIsEnabled, perFamilyIsEnabled, variantEvalArgs.mergeEvals);
        final Map<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> compVCs = allComps.isEmpty() ? Collections.emptyMap() : bindVariantContexts(variantMap, variantEvalArgs.comps, byFilterIsEnabled, false, false, false, false);

        final VariantEvalContext variantEvalContext = new VariantEvalContext(referenceContext, featureContext, variantMap, this);

        // for each eval track
        for ( final FeatureInput<VariantContext> evalInput : variantEvalArgs.evals ) {
            final Map<String, Collection<VariantContext>> evalSet = evalVCs.containsKey(evalInput) ? evalVCs.get(evalInput) : Collections.emptyMap();

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
                    for ( final FeatureInput<VariantContext> compInput : variantEvalArgs.comps ) {
                        processComp(variantEvalContext, eval, evalName, compInput, stratLevelName, compVCs, evalSetBySample);
                    }

                    if (variantEvalArgs.comps.isEmpty()) {
                        processComp(variantEvalContext, eval, evalName, null, stratLevelName, compVCs, evalSetBySample);
                    }
                }
            }

            if ( variantEvalArgs.mergeEvals ) break; // stop processing the eval tracks
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

    public boolean isSubsettingToSpecificSamples() { return isSubsettingSamples; }

    public Set<String> getSampleNamesForEvaluation() { return sampleNamesForEvaluation; }

    public Set<String> getFamilyNamesForEvaluation() { return familyNamesForEvaluation; }

    public int getNumberOfSamplesForEvaluation() {
        if (sampleNamesForEvaluation!= null &&  !sampleNamesForEvaluation.isEmpty())
            return sampleNamesForEvaluation.size();
        else {
            return variantEvalArgs.numSamplesFromArgument;
        }
    }
    public Set<String> getSampleNamesForStratification() { return sampleNamesForStratification; }

    public Set<String> getFamilyNamesForStratification() { return familyNamesForStratification; }

    public Set<SortableJexlVCMatchExp> getJexlExpressions() { return jexlExpressions; }


    public AlleleFrequency.StratifyingScale getAFScale() { return AFScale; }
    public boolean getCompAFStratifier() { return useCompAFStratifier; }

    public SampleDB getSampleDB() {
        return sampleDB;
    }

    public List<FeatureInput<VariantContext>> getKnowns() {
        return knowns;
    }

    /**
     * If an evaluator calls this method it must override {@link VariantEvaluator#requiresTerritoryToBeSpecified()} to return true.
     * @return either the size of the interval list given to the tool or the size of the reference given to the tool
     */
    public long getnProcessedLoci() {
        if (traversalIntervals == null){
            throw new GATKException("BUG: One of the evaluators used should have overriden requiresTerritoryToBeSpecified, please report this to the developers." +
                    "\nEvaluators: " + stratManager.values()
                    .stream()
                    .flatMap(evaluator -> evaluator.getVariantEvaluators().stream())
                    .map(VariantEvaluator::getSimpleName)
                    .sorted()
                    .distinct()
                    .collect(Collectors.joining(", ")));
        }
        return traversalIntervals.stream().mapToLong(SimpleInterval::size).sum();
    }

    private boolean compHasMatchingEval(final VariantContext comp, final Collection<VariantContext> evals) {
        // find all of the matching comps
        for ( final VariantContext eval : evals ) {
            if ( eval != null && doEvalAndCompMatch(comp, eval, variantEvalArgs.requireStrictAlleleMatch) != EvalCompMatchType.NO_MATCH )
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
            switch ( doEvalAndCompMatch(comp, eval, variantEvalArgs.requireStrictAlleleMatch) ) {
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

    private SimpleInterval generateContextInterval(List<VariantContext> variantContexts) {
        final int maxEnd = variantContexts.stream().map(VariantContext::getEnd).max(Integer::compareTo).get();

        return new SimpleInterval(variantContexts.get(0).getContig(), variantContexts.get(0).getStart(), maxEnd);
    }

    private void processComp(VariantEvalContext vec, VariantContext eval, String evalName, FeatureInput<VariantContext> compInput, String stratLevelName, Map<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> compVCs, Collection<VariantContext> evalSetBySample) {
        String compName = getNameForInput(compInput);

        // no sample stratification for comps
        final HashMap<String, Collection<VariantContext>> compSetHash = compInput == null ? null : compVCs.get(compInput);
        final Collection<VariantContext> compSet = (compSetHash == null || compSetHash.isEmpty()) ? Collections.<VariantContext>emptyList() : compVCs.get(compInput).values().iterator().next();

        // find the comp
        final VariantContext comp = findMatchingComp(eval, compSet);

        Collection<EvaluationContext> contextsForStratification;
        if (perFamilyIsEnabled)
            contextsForStratification = getEvaluationContexts(vec, eval, evalName, comp, compName, null, stratLevelName);
        else {
            String familyID;
            if (stratLevelName.equals("all"))
                familyID = "all";
            else
                familyID = sampleDB.getSample(stratLevelName).getFamilyID();
            contextsForStratification = getEvaluationContexts(vec, eval, evalName, comp, compName, stratLevelName, familyID);
        }
        for ( EvaluationContext nec : contextsForStratification ) {

            // eval against the comp
            synchronized (nec) {
                nec.apply(vec, comp, eval);
            }

            // eval=null against all comps of different type that aren't bound to another eval
            for ( VariantContext otherComp : compSet ) {
                if ( otherComp != comp && ! compHasMatchingEval(otherComp, evalSetBySample) ) {
                    synchronized (nec) {
                        nec.apply(vec, otherComp, null);
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
     * @param vec
     * @param eval
     * @param evalName
     * @param comp
     * @param compName
     * @param sampleName
     * @return
     */
    protected Collection<EvaluationContext> getEvaluationContexts(final VariantEvalContext vec,
                                                                  final VariantContext eval,
                                                                  final String evalName,
                                                                  final VariantContext comp,
                                                                  final String compName,
                                                                  final String sampleName,
                                                                  final String familyName) {
        final List<List<Object>> states = new LinkedList<>();
        for ( final VariantStratifier vs : stratManager.getStratifiers() ) {
            states.add(vs.getRelevantStates(vec, comp, compName, eval, evalName, sampleName, familyName));
        }
        return stratManager.values(states);
    }

    protected List<String> getModulesToUse() {
        return Collections.unmodifiableList(variantEvalArgs.modulesToUse);
    }

    /**
     * For a list of track names, bind the variant contexts to a trackName->sampleName->VariantContext mapping.
     * Additional variant contexts per sample are automatically generated and added to the map unless the sample name
     * matches the ALL_SAMPLE_NAME constant.
     *
     * @return the mapping of track to VC list that should be populated
     */
    public HashMap<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>>
    bindVariantContexts(Map<FeatureInput<VariantContext>, List<VariantContext>> variantMap,
                        List<FeatureInput<VariantContext>> tracks,
                        boolean byFilter,
                        boolean subsetBySample,
                        boolean trackPerSample,
                        boolean trackPerFamily,
                        boolean mergeTracks) {
        HashMap<FeatureInput<VariantContext>, HashMap<String, Collection<VariantContext>>> bindings = new HashMap<>();

        FeatureInput<VariantContext> firstTrack = tracks.isEmpty() ? null : tracks.get(0);
        for (FeatureInput<VariantContext> track : tracks) {
            HashMap<String, Collection<VariantContext>> mapping = new HashMap<>();

            if (variantMap.containsKey(track)) {
                //Note: these are limiting to only those w/ the same start, as was the GATK3 behavior.
                for (VariantContext vc : variantMap.get(track)) {

                    // First, filter the VariantContext to represent only the samples for evaluation
                    VariantContext vcsub = vc;

                    if ((subsetBySample) && vc.hasGenotypes())
                        vcsub = getSubsetOfVariantContext(vc, getSampleNamesForEvaluation());

                    //always add a mapping for all samples together
                    if ((byFilter || !vcsub.isFiltered())) {
                        addMapping(mapping, VariantEvalArgumentCollection.ALL_SAMPLE_NAME, vcsub);
                    }

                    // Now, if stratifying, split the subsetted vc per sample and add each as a new context
                    if (vc.hasGenotypes() && trackPerSample) {
                        for (String sampleName : getSampleNamesForEvaluation()) {
                            final VariantContext samplevc = getSubsetOfVariantContext(vc, sampleName);

                            if (byFilter || !samplevc.isFiltered()) {
                                addMapping(mapping, sampleName, samplevc);
                            }
                        }
                    } else if (vc.hasGenotypes() && trackPerFamily) {
                        for (final String familyName : getFamilyNamesForEvaluation()) {
                            Set<String> familyMemberNames;
                            //if the current stratification family name is "all", then add all the families to the VC for evaluation here
                            if (familyName.equals(VariantEvalArgumentCollection.ALL_FAMILY_NAME)) {
                                familyMemberNames = getSampleNamesForEvaluation();
                            } else {
                                familyMemberNames = getSampleDB().getFamily(familyName).stream().map(Sample::getID).collect(Collectors.toSet());
                            }
                            final VariantContext samplevc = getSubsetOfVariantContext(vc, familyMemberNames);

                            if (byFilter || !samplevc.isFiltered()) {
                                addMapping(mapping, familyName, samplevc);
                            }
                        }
                    }
                }
            }

            if (mergeTracks && bindings.containsKey(firstTrack)) {
                // go through each binding of sample -> value and add all of the bindings from this entry
                HashMap<String, Collection<VariantContext>> firstMapping = bindings.get(firstTrack);
                for (Map.Entry<String, Collection<VariantContext>> elt : mapping.entrySet()) {
                    Collection<VariantContext> firstMappingSet = firstMapping.get(elt.getKey());
                    if (firstMappingSet != null) {
                        firstMappingSet.addAll(elt.getValue());
                    } else {
                        firstMapping.put(elt.getKey(), elt.getValue());
                    }
                }
            } else {
                bindings.put(track, mapping);
            }
        }

        return bindings;
    }

    private void addMapping(HashMap<String, Collection<VariantContext>> mappings, String sample, VariantContext vc) {
        if (!mappings.containsKey(sample))
            mappings.put(sample, new ArrayList<>(1));
        mappings.get(sample).add(vc);
    }

    /**
     * Subset a VariantContext to a single sample
     *
     * @param vc         the VariantContext object containing multiple samples
     * @param sampleName the sample to pull out of the VariantContext
     * @return a new VariantContext with just the requested sample
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, String sampleName) {
        return getSubsetOfVariantContext(vc, Collections.singleton(sampleName));
    }

    /**
     * Subset a VariantContext to a set of samples
     *
     * @param vc          the VariantContext object containing multiple samples
     * @param sampleNames the samples to pull out of the VariantContext
     * @return a new VariantContext with just the requested samples
     */
    public VariantContext getSubsetOfVariantContext(VariantContext vc, Set<String> sampleNames) {
        // if we want to preserve AC0 sites as polymorphic we need to not rederive alleles
        final boolean deriveAlleles = variantEvalArgs.ignoreAC0Sites();
        return ensureAnnotations(vc, vc.subContextFromSamples(sampleNames, deriveAlleles));
    }

    public VariantContext ensureAnnotations(final VariantContext vc, final VariantContext vcsub) {
        final int originalAlleleCount = vc.getHetCount() + 2 * vc.getHomVarCount();
        final int newAlleleCount = vcsub.getHetCount() + 2 * vcsub.getHomVarCount();
        final boolean isSingleton = originalAlleleCount == newAlleleCount && newAlleleCount == 1;
        final boolean hasChrCountAnnotations = vcsub.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) &&
                vcsub.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) &&
                vcsub.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY);

        if ( ! isSingleton && hasChrCountAnnotations ) {
            // nothing to update
            return vcsub;
        } else {
            // have to do the work
            VariantContextBuilder builder = new VariantContextBuilder(vcsub);

            if ( isSingleton )
                builder.attribute(IS_SINGLETON_KEY, true);

            if ( ! hasChrCountAnnotations )
                VariantContextUtils.calculateChromosomeCounts(builder, true);

            return builder.make();
        }
    }

    /**
     * Initialize required, standard and user-specified stratification objects
     *
     * @param noStandardStrats  don't use the standard stratifications
     * @param modulesToUse      the list of stratification modules to use
     * @return set of stratifications to use
     */
    public List<VariantStratifier> initializeStratificationObjects(boolean noStandardStrats, List<String> modulesToUse) {
        TreeSet<VariantStratifier> strats = new TreeSet<>();
        Set<String> stratsToUse = new HashSet<>(requiredStratificationNames);

        // By default, use standard stratification modules.
        if (!noStandardStrats) {
            stratsToUse.addAll(standardStratificationNames);
        }

        // Now add the user-selected modules
        stratsToUse.addAll(modulesToUse);

        // Instantiate the stratifications
        for (String module : stratsToUse) {
            if (!stratifierClasses.containsKey(module)) {
                throw new CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            if (stratifierClasses.containsKey(module)) {
                Class<? extends VariantStratifier> c = stratifierClasses.get(module);

                VariantStratifier vs = createVariantStratifier(c);
                strats.add(vs);
            }
        }

        return new ArrayList<>(strats);
    }

    public VariantStratifier createVariantStratifier(Class<? extends VariantStratifier> clazz) {
        return createClass(clazz);
    }

    public VariantEvaluator createVariantEvaluator(Class<? extends VariantEvaluator> clazz) {
        return createClass(clazz);
    }

    private <T> T createClass(Class<? extends T> clazz) {
        try {
            return clazz.getDeclaredConstructor(VariantEvalEngine.class).newInstance(this);
        }
        catch (final InstantiationException | IllegalAccessException | InvocationTargetException | NoSuchMethodException e) {
            if (e.getCause() instanceof CommandLineException) {
                throw (CommandLineException)e.getCause();
            }

            throw new GATKException("Problem making an instance of " + clazz + " Do check that the class has a constructor that accepts VariantEvalEngine", e);
        }
    }

    /**
     * Initialize required, standard and user-specified evaluation objects
     *
     * @param noStandardEvals don't use the standard evaluations
     * @param modulesToUse    the list of evaluation modules to use
     * @return set of evaluations to use
     */
    public Set<Class<? extends VariantEvaluator>> initializeEvaluationObjects(boolean noStandardEvals, List<String> modulesToUse) {
        Set<String> evalsToUse = new TreeSet<>(modulesToUse);

        // By default, use standard eval modules.
        if (!noStandardEvals) {
            evalsToUse.addAll(standardEvaluatorNames);
        }

        // Get the specific classes provided.
        Set<Class<? extends VariantEvaluator>> evals = new HashSet<>();
        for (String module : evalsToUse) {
            if (!evaluatorClasses.containsKey(module)) {
                throw new CommandLineException("Module " + module + " could not be found; please check that you have specified the class name correctly");
            }

            evals.add(evaluatorClasses.get(module));
        }

        //add MetricsCollection if required modules are included
        if (evals.contains(evaluatorClasses.get("CompOverlap")) && evals.contains(evaluatorClasses.get("IndelSummary")) && evals.contains(evaluatorClasses.get("TiTvVariantEvaluator")) && evals.contains(evaluatorClasses.get("CountVariants")) && evals.contains(evaluatorClasses.get("MultiallelicSummary")) )
            evals.add(evaluatorClasses.get("MetricsCollection"));

        return evals;
    }

    public static Map<String, Class<? extends VariantStratifier>> getStratifierClasses() {
        return Collections.unmodifiableMap(stratifierClasses);
    }

    public static Set<String> getStandardStratificationNames() {
        return Collections.unmodifiableSet(standardStratificationNames);
    }

    public static Set<String> getRequiredStratificationNames() {
        return Collections.unmodifiableSet(requiredStratificationNames);
    }

    public static Map<String, Class<? extends VariantEvaluator>> getEvaluatorClasses() {
        return Collections.unmodifiableMap(evaluatorClasses);
    }

    public static Set<String> getStandardEvaluatorNames() {
        return Collections.unmodifiableSet(standardEvaluatorNames);
    }

    public VariantEvalArgumentCollection getVariantEvalArgs() {
        return variantEvalArgs;
    }

    public void setAFScale(AlleleFrequency.StratifyingScale AFScale) {
        this.AFScale = AFScale;
    }

    public void setUseCompAFStratifier(boolean useCompAFStratifier) {
        this.useCompAFStratifier = useCompAFStratifier;
    }

    public final <T extends Feature> Object getHeaderForFeatures( final FeatureInput<T> featureDescriptor ) {
        return features.getHeader(featureDescriptor);
    }

    public SAMSequenceDictionary getSequenceDictionaryForDrivingVariants() {
        return samSequenceDictionaryForDrivingVariants;
    }

    @VisibleForTesting
    protected StratificationManager<VariantStratifier, EvaluationContext> getStratManager() {
        return stratManager;
    }
}
