package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import net.minidev.json.JSONValue;
import net.minidev.json.parser.ParseException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.samples.Trio;

import java.io.*;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.*;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.jetbrains.annotations.NotNull;

import static org.apache.commons.math3.util.FastMath.*;


/**
 * -Extract matrix of properties for each variant.
 * -Also extract, numVariants x numTrios x 3 tensors of allele count and genotype quality.
 * -These data will be used to train a variant filter based on min GQ (and stratified by other variant properties) that
 * maximizes the admission of variants with Mendelian inheritance pattern while omitting non-Mendelian variants.
 * -Derived class must implement abstract method trainFilter()
 */
public abstract class MinGqVariantFilterBase extends VariantWalker {
    @Argument(fullName=StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME,
              shortName=StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME,
              doc="Pedigree file, necessary for \"Train\" mode, ignored in \"Filter\" mode.", optional=true)
    public GATKPath pedigreeFile = null;

    @Argument(fullName=StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              doc="Output VCF produced in \"Filter\" mode, ignored in \"Train\" mode.", optional=true)
    public GATKPath outputFile = null;

    @Argument(fullName="model-file", shortName="m",
              doc="Path to saved pre-existing filter model. In \"Filter\" mode, this is a mandatory argument."
                 +" In \"Train\" mode, this is an optional argument for additional training based on an existing model.")
    public GATKPath modelFile = null;

    @Argument(fullName="properties-table-file", optional=true,
            doc="Path to save properties table as JSON for analysis outside of this program. Only has an effect in TRAIN mode.")
    public GATKPath propertiesTableFile = null;

    @Argument(fullName="truth-file", shortName="t", optional=true,
              doc="Path to JSON file with truth data. Keys are sample IDs and values are objects with key \"good\""
                 +" corresponding to a list of known true variant IDs, and key \"bad\" corresponding to a list of known"
                 +" bad variant IDs")
    public GATKPath truthFile = null;

    enum RunMode {Train, Filter}
    @Argument(fullName="mode", doc="Mode of operation: either \"Train\" or \"Filter\"")
    public RunMode runMode;

    @Argument(fullName="keep-homvar", shortName="kh", doc="Keep homvar variants even if their GQ is less than min-GQ", optional=true)
    public boolean keepHomvar = false;

    @Argument(fullName="train-homref", doc="Do not train xgboost classifier with HOMREF genotypes", optional=true)
    public boolean trainHomref = true;

    @Argument(fullName="keep-multiallelic", shortName="km", doc="Keep multiallelic variants even if their GQ is less than min-GQ", optional=true)
    public boolean keepMultiallelic = true;

    @Argument(fullName="use-copy-number-calls", doc="If true, attempt to use copy number info if genotype is NO CALL. "+
                       "This also entails keeping track of copy number call quality (RD_GQ) and reconciling that with GQ",
                       optional=true)
    public boolean useCopyNumberCalls = false;

    @Argument(fullName="keep-non-var", shortName="kn", doc="In \"Filter\" mode: if true, keep variants with no HET or HOMVAR GTs (all REF or no-call), otherwise omit those variants", optional=true)
    public boolean keepNonVar = false;

    @Argument(fullName="validation-proportion", shortName="vp", doc="Proportion of variants to set aside for cross-validation",
              optional=true, minValue=0.0, maxValue=1.0)
    public double validationProportion = 0.2;

    @Argument(fullName="report-min-gq-filter-threshold", shortName="rf", optional=true, minValue=0.0, maxValue=1.0,
              doc="Add \"" + EXCESSIVE_MIN_QUALITY_FILTER_KEY + "\" to FILTER if the proportion of samples with calls "
                  + "filtered due to low quality is greater than this threshold.")
    public double reportMinGqFilterThreshold = 0.005;

    @Argument(fullName="truth-weight", shortName="tw", optional=true, minValue=0.0, maxValue=1.0,
            doc="Weight for truth data in loss function.")
    public static double truthWeight = 0.5;

    @Argument(fullName="progress-verbosity", shortName="p",
              doc="Level of verbosity for progress log", optional = true)
    public int progressVerbosity = 2;

    static final String minSamplesToEstimateAlleleFrequencyKey = "min-samples-to-estimate-allele-frequency";
    @Argument(fullName=minSamplesToEstimateAlleleFrequencyKey, shortName="ms", optional=true,
              doc="If the VCF does not have allele frequency, estimate it from the sample population if there are at least this many samples. Otherwise throw an exception.")
    public int minSamplesToEstimateAlleleFrequency = 100;

    @Argument(fullName="genome-track", shortName="gt", optional=true)
    final List<String> genomeTrackFiles = new ArrayList<>();

    @Argument(fullName="min-scaled-logits", shortName="sl", optional=true)
    short minScaledLogits = probToScaledLogits(0.5);
    short FAILING_GQ;

    @Argument(fullName="target-precision", doc="Desired minimum precision to achieve before maximizing recall", optional=true)
    public static double targetPrecision = 0.90;

    @Argument(fullName="min-truth-variants-to-estimate-recall",
              doc="Minimum number of variants with truth data to estimate recall", optional=true)
    public static int minVariantsToEstimateRecall = 50;

    @Argument(fullName="min-variants-for-scaling-weight",
            doc="The weight for each variant category is inversely-proportional to number of variants, to prevent " +
                "drowning out rare variant types. Place a lower limit on the inverse scaling, so a few ultra-rare " +
                " variants don't cause the objective to misbehave", optional=true, minValue=1)
    public static int minVariantsForScalingWeight = 100;

    @Argument(fullName="max-inheritance-af", optional=true,
        doc="Max allele frequency where inheritance will be considered truthful. AF in range (mIAf, 1.0 - mIAf) have no score")
    public static double maxInheritanceAf = 0.05;

    @Argument(fullName="large-af-weight-penalty", optional=true,
        doc="Multiplicative penalty for learning weights based on inheritance for variants with large allele frequency.")
    public static double largeAfWeightPenalty = 1e-6;

    @Argument(fullName="strict-mendelian", optional=true,
            doc="If true, apply normal rules for mendelian inheritance; if false, allow confusion between HET and HOMVAR"
    )
    public static boolean strictMendelian = true;

    @Argument(fullName="min-bin-weight", optional=true,
            doc="Minimum value of weight for a given SV category bin"
    )
    public static double minBinWeight = 0.01;

    @Argument(fullName="scaled-logit-property", optional=true,
              doc="Name of genotype property where logit-scaled variant quality will be stored"
    )
    public static String scaledLogitProperty = "SL";

    @Argument(fullName="error-on-no-trios", optional=true, doc="Throw exception if there are no trios in training mode")
    public static Boolean errorOnNoTrios = true;

    List<TrackOverlapDetector> trackOverlapDetectors = null;

    static final Map<String, double[]> propertyBinsMap = new LinkedHashMap<>();
    // numVariants array of property-category bins
    protected int[] propertyBins = null;
    int numPropertyBins;
    protected String[] propertyBinLabels = null;
    protected List<String> propertyBinLabelColumns = null;
    protected double[] propertyBinInheritanceWeights = null;
    protected double[] propertyBinTruthWeights = null;
    protected float[] propertyBinMinGqWeights = null;
    protected float[] propertyBinGoodInheritanceWeights = null;
    protected float[] propertyBinBadInheritanceWeights = null;
    protected float[] propertyBinGoodTruthWeights = null;
    protected float[] propertyBinBadTruthWeights = null;
    protected boolean[] propertyBinIsLargeAlleleFraction = null;

    protected List<String> sampleIds = null; // numSamples list of IDs for samples that will be used for training
    protected int[][] trioSampleIndices = null; // numTrios x 3 matrix of sample indices (paternal, maternal, child) for each trio
    protected Set<Integer> allTrioSampleIndices; // set of all sample indices that are in a trio
    protected final Map<String, Integer> sampleIndices = new HashMap<>(); // map from sampleId to numerical index, for samples used in training
    protected final PropertiesTable propertiesTable = new PropertiesTable();
    // direct links to important properties, so they don't need to be looked up in the table each time they are used:
    protected PropertiesTable.ShortMatProperty sampleVariantCallQualities = null;
    protected PropertiesTable.ByteMatProperty sampleVariantAlleleCounts = null;
    protected PropertiesTable.ByteMatProperty sampleVariantNoCallCounts = null;
    protected List<String> variantIds = new ArrayList<>();

    protected Random randomGenerator = Utils.getRandomGenerator();

    private VariantContextWriter vcfWriter = null;

    private int numVariants;
    private int numSamples;
    private int numTrios;

    private Set<String> vcfHeaderIds = null;
    private static final float PROB_EPS = 1.0e-3F;
    private static final String NO_CALL_COUNTS_KEY = "NO_CALL_COUNTS";
    private static final String SVLEN_KEY = "SVLEN";
    private static final String EV_KEY = "EV";
    private static final String CONC_ST_KEY = "CONC_ST";
    private static final String FILTER_KEY = "FILTER";
    private static final String ALGORITHMS_KEY = "ALGORITHMS";
    private static final String EVIDENCE_KEY = "EVIDENCE";
    private static final String STATUS_KEY = "STATUS";
    private static final String NON_REF_GENOTYPE_CONCORDANCE_KEY = "NON_REF_GENOTYPE_CONCORDANCE";
    private static final String VAR_PPV_KEY = "VAR_PPV";
    private static final String VAR_SENSITIVITY_KEY = "VAR_SENSITIVITY";
    private static final String TRUTH_AF_KEY = "TRUTH_AF";
    private static final String MIN_QUALITY_KEY = "MIN" + scaledLogitProperty;
    private static final String EXCESSIVE_MIN_QUALITY_FILTER_KEY = "LOW_QUALITY";
    private static final String MULTIALLELIC_FILTER = "MULTIALLELIC";
    private static final String GOOD_VARIANT_TRUTH_KEY = "good_variant_ids";
    private static final String BAD_VARIANT_TRUTH_KEY = "bad_variant_ids";
    private static final String CALL_QUALITY_KEY = "CALL_QUALITY";
    protected static final int FATHER_IND = 0;
    protected static final int MOTHER_IND = 1;
    protected static final int CHILD_IND = 2;

    private static final String IS_COPY_NUMBER_CALL_KEY = "IS_CN_CALL";

    private static final String PE_GQ_KEY = "PE_GQ";
    private static final String RD_CN_KEY = "RD_CN";
    private static final String RD_GQ_KEY = "RD_GQ";
    private static final String SR_GQ_KEY = "SR_GQ";
    private static final short MISSING_GQ_VAL = -1;


    // properties used to gather main matrix / tensors during apply function
    // map from variant ID to array of known good/bad sample indices for that variant
    private Map<String, Set<Integer>> goodVariantSampleIndices = null;
    private Map<String, Set<Integer>> badVariantSampleIndices = null;


    // train/validation split indices
    private int[] trainingIndices;
    private int[] validationIndices;

    // properties for calculating f1 score or estimating its pseudo-derivatives
    private MinGq[] perVariantOptimalMinGq = null;
    protected double optimalProportionOfSampleVariantsPassing;

    protected final int getNumVariants() { return propertiesTable.getNumRows(); }
    protected final int getNumSamples() { return numSamples; }
    protected final int getNumTrios() { return numTrios; }
    protected final int getNumProperties() { return propertiesTable.getNumProperties(); }
    protected final int[] getTrainingIndices() { return trainingIndices; }
    protected final int[] getValidationIndices() { return validationIndices; }

    // stats on tool actions, input/output VCFs
    private int numFilterableGenotypes;
    private int numFilteredGenotypes;
    private int numInputVar;
    private int numInputRef;
    private int numInputNoCall;
    private int numOutputVar;
    private int numOutputRef;
    private int numOutputNoCall;
    private int numDroppedVariants;

    /**
     * Entry-point function to initialize the samples database from input data
     */
    private void getPedTrios() {
        final Set<Trio> pedTrios;
        if(pedigreeFile == null) {
            pedTrios = new HashSet<>();
            if(errorOnNoTrios) {
                throw new UserException.BadInput(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME
                        + " must be specified in \"TRAIN\" mode");
            } else {
                System.err.println(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME + " was not specified in \"Train\" mode");
            }
        } else {
            final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
            final Set<String> vcfSamples = new HashSet<>(getHeaderForVariants().getGenotypeSamples());
            // get the trios from the pedigree file, keeping only those that are fully present in the input VCF
            pedTrios = sampleDBBuilder.getFinalSampleDB()
                    .getTrios()
                    .stream()
                    .filter(trio -> vcfSamples.contains(trio.getPaternalID()) && vcfSamples.contains(trio.getMaternalID())
                            && vcfSamples.contains(trio.getChildID()))
                    .collect(Collectors.toSet());
            if (pedTrios.isEmpty()) {
                if(errorOnNoTrios) {
                    throw new UserException.BadInput(
                            "The pedigree file (" + pedigreeFile + ") does not contain any trios that are present in the input VCF "
                                    + "(" + drivingVariantFile + ")"
                    );
                } else {
                    System.err.println(
                    "The pedigree file (" + pedigreeFile + ") does not contain any trios that are present in the input VCF "
                            + "(" + drivingVariantFile + ")"
                    );
                }
            }
        }
        numTrios = pedTrios.size();
        if(numTrios == 0) {
            truthWeight = 1.0;
        }
        // collect ped trios into training samples
        pedTrios.stream()
                .flatMap(trio -> Stream.of(trio.getPaternalID(), trio.getMaternalID(), trio.getChildID()))
                .forEach(this::addSampleIndex);
        // keep track of the samples that are in trios (those will always have "truth" through inheritance)
        trioSampleIndices = pedTrios.stream()
                .map(
                    trio -> Stream.of(trio.getPaternalID(), trio.getMaternalID(), trio.getChildID())
                            .mapToInt(sampleIndices::get)
                            .toArray()
                )
                .toArray(int[][]::new);

        allTrioSampleIndices = pedTrios.stream()
                .flatMap(trio -> Stream.of(trio.getPaternalID(), trio.getMaternalID(), trio.getChildID()))
                .map(sampleIndices::get)
                .collect(Collectors.toSet());

    }

    private void loadVariantTruthData() {
        if(truthFile == null) {
            System.out.println("No truth file specified. Not using truth data.");
            return;
        } else {
            System.out.println("Loading truth data from " + truthFile);
        }
        final JSONObject jsonObject;
        try (final InputStream fileInputStream = truthFile.getInputStream()){
            jsonObject = (JSONObject) JSONValue.parseWithException(fileInputStream);
        } catch (IOException | ParseException ioException) {
            throw new GATKException("Unable to parse JSON from inputStream", ioException);
        }
        final Set<String> vcfSamples = new HashSet<>(getHeaderForVariants().getGenotypeSamples());
        final Map<String, Set<String>> goodVariantSampleIds = new HashMap<>();
        final Map<String, Set<String>> badVariantSampleIds = new HashMap<>();
        for(final Map.Entry<String, Object> sampleTruthEntry : jsonObject.entrySet()) {
            final String sampleId = sampleTruthEntry.getKey();
            if(!vcfSamples.contains(sampleId)) {
                continue; // don't bother storing truth data that doesn't overlap the input VCF
            }
            final JSONObject sampleTruth = (JSONObject)sampleTruthEntry.getValue();
            for(final Object variantIdObj : (JSONArray)sampleTruth.get(GOOD_VARIANT_TRUTH_KEY)) {
                final String variantId = (String)variantIdObj;
                if(goodVariantSampleIds.containsKey(variantId)) {
                    goodVariantSampleIds.get(variantId).add(sampleId);
                } else {
                    goodVariantSampleIds.put(variantId, new HashSet<>(Collections.singleton(sampleId)) );
                }
            }
            for(final Object variantIdObj : (JSONArray)sampleTruth.get(BAD_VARIANT_TRUTH_KEY)) {
                final String variantId = (String)variantIdObj;
                if(badVariantSampleIds.containsKey(variantId)) {
                    badVariantSampleIds.get(variantId).add(sampleId);
                } else {
                    badVariantSampleIds.put(variantId, new HashSet<>(Collections.singleton(sampleId)));
                }
            }
        }
        if(goodVariantSampleIds.isEmpty() && badVariantSampleIds.isEmpty()) {
            System.out.println("Truth file specified (" + truthFile + "), but no samples/variants overlap with input VCF ("
                               + drivingVariantFile + "). Not using truth data.");
            goodVariantSampleIndices = null;
            badVariantSampleIndices = null;
            return;
        }
        // Add any new samples that have truth but are not in trios file to the training samples
        goodVariantSampleIds.values().stream().flatMap(Collection::stream).forEach(this::addSampleIndex);
        badVariantSampleIds.values().stream().flatMap(Collection::stream).forEach(this::addSampleIndex);

        // Finally replace the sets of sample IDs with sets of sample indices. We'll use this data structure because
        // it's smaller.
        goodVariantSampleIndices = new HashMap<>();
        goodVariantSampleIds.forEach(
                (variantId, sampleIds) -> goodVariantSampleIndices.put(
                        variantId, sampleIds.stream().mapToInt(sampleIndices::get).boxed().collect(Collectors.toSet())
                    )
        );
        badVariantSampleIndices = new HashMap<>();
        badVariantSampleIds.forEach(
                (variantId, sampleIds) -> badVariantSampleIndices.put(
                        variantId, sampleIds.stream().mapToInt(sampleIndices::get).boxed().collect(Collectors.toSet())
                )
        );
    }

    private void addSampleIndex(final String sampleId) {
        if(!sampleIndices.containsKey(sampleId)) {
            sampleIndices.put(sampleId, sampleIndices.size());
        }
    }

    private void setSampleIds() {
        sampleIds = sampleIndices.entrySet()
            .stream()
            .sorted(Comparator.comparingInt(Map.Entry::getValue))
            .map(Map.Entry::getKey)
            .collect(Collectors.toList());
        numSamples = sampleIndices.size();
        if(numSamples == 0) {
            throw new IllegalArgumentException(
                "No " + (runMode == RunMode.Train ? "trainable" : "filterable") + " samples were found"
            );
        }
    }

    private void initializeVcfWriter() {
        List<Options> options = new ArrayList<>();
        options.add(Options.INDEX_ON_THE_FLY);
        vcfWriter = GATKVariantContextUtils.createVCFWriter(
                outputFile.toPath(),
                getHeaderForVariants().getSequenceDictionary(),
                createOutputVariantMD5,
                options.toArray(new Options[0]));
        final Set<VCFHeaderLine> hInfo = new LinkedHashSet<>(getHeaderForVariants().getMetaDataInInputOrder());
        final String filterableVariant = (keepMultiallelic ? "biallelic " : "") +
                                         (keepHomvar ? "non-homvar" : "") +
                                         "variant";
        hInfo.add(new VCFFormatHeaderLine(scaledLogitProperty, 1, VCFHeaderLineType.Integer,
                               LOGIT_SCALE + " times the logits that a genotype is correct"));
        hInfo.add(new VCFInfoHeaderLine(MIN_QUALITY_KEY, 1, VCFHeaderLineType.Integer,
                            "Minimum passing GQ for each " + filterableVariant));
        hInfo.add(new VCFFilterHeaderLine(EXCESSIVE_MIN_QUALITY_FILTER_KEY,
                               "More than " + (100 * reportMinGqFilterThreshold) + "% of sample GTs were masked as no-call GTs due to low GQ"));
        final VCFHeader vcfHeader = new VCFHeader(hInfo, getHeaderForVariants().getGenotypeSamples());
        vcfWriter.setHeader(vcfHeader);
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void onTraversalStart() {
        vcfHeaderIds = Stream.concat(
            getHeaderForVariants().getInfoHeaderLines().stream().map(VCFCompoundHeaderLine::getID),
            getHeaderForVariants().getFormatHeaderLines().stream().map(VCFCompoundHeaderLine::getID)
        ).collect(Collectors.toSet());
        // get failing GQ by converting min passing scaled logits to GQ, then decreasing by 1
        FAILING_GQ = (short)(probToPhred(scaledLogitsToP(minScaledLogits)) - 1);
        // set propertyBinsMap bins
        propertyBinsMap.put(SVLEN_KEY, new double[] {50.0, 500.0, 5000.0});
        propertyBinsMap.put(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {maxInheritanceAf, 1.0 - maxInheritanceAf});

        loadTrainedModel();  // load model and saved properties stats
        trackOverlapDetectors = genomeTrackFiles.stream()  // load genome tracks
            .map(TrackOverlapDetector::new)
            .collect(Collectors.toList());

        numInputVar = 0;
        numInputRef = 0;
        numInputNoCall = 0;
        if(runMode == RunMode.Train) {
            getPedTrios();  // get trios from pedigree file
            loadVariantTruthData(); // load variant truth data from JSON file
            // at this point, included samples will only be those that are in a trio or have some variant truth data
            setSampleIds(); // set data organizing training samples
        } else {
            getHeaderForVariants().getGenotypeSamples().forEach(this::addSampleIndex); // use all samples in VCF
            setSampleIds();
            initializeVcfWriter();  // initialize vcfWriter and write header
            numFilterableGenotypes = 0;
            numFilteredGenotypes = 0;
            numOutputVar = 0;
            numOutputRef = 0;
            numOutputNoCall = 0;
            numDroppedVariants = 0;
            propertiesTable.setNumAllocatedRows(1);
            System.out.println("minScaledLogits=" + minScaledLogits);
        }
    }

    private boolean getIsMultiallelic(final VariantContext variantContext) {
        return variantContext.getNAlleles() > 2 || variantContext.getFilters().contains(MULTIALLELIC_FILTER);
    }

    private boolean getVariantIsTrainable(final VariantContext variantContext,
                                          final byte[] sampleAlleleCounts,
                                          final byte[] sampleNoCallCounts) {
        if(keepMultiallelic && getIsMultiallelic(variantContext)) {
            return false;  // can't filter it, so can't train
        }

        // check if any of the samples have truth data and can be filtered. If so, the variant is trainable.
        // first check if any member of trios are filterable
        final boolean inheritanceTrainable = Arrays.stream(trioSampleIndices).anyMatch(
                trioIndices -> {  // trainable if any trio is filterable: all samples present, no no-calls, at least one
                    //                                      non-ref sample, and at least one filterable allele count
                    if(Arrays.stream(trioIndices).anyMatch(sampleIndex -> sampleNoCallCounts[sampleIndex] > 0)) {
                        return false; // trio has no-calls
                    }
                    final int fatherAc = sampleAlleleCounts[trioIndices[FATHER_IND]];
                    final int motherAc = sampleAlleleCounts[trioIndices[MOTHER_IND]];
                    final int childAc = sampleAlleleCounts[trioIndices[CHILD_IND]];
                    return (fatherAc > 0 || motherAc > 0 || childAc > 0) &&
                            (alleleCountIsFilterable(fatherAc) || alleleCountIsFilterable(motherAc) ||
                                    alleleCountIsFilterable(childAc));
                }
        );

        if(inheritanceTrainable) {
            return true;
        }
        // next check if there are any samples in the truth set that are filterable
        if(goodVariantSampleIndices != null &&
           !getFilterableTruthSampleIndices(variantContext, goodVariantSampleIndices, sampleAlleleCounts,
                                            sampleNoCallCounts).isEmpty()) {
           //return true;
           return badVariantSampleIndices != null &&
                   !getFilterableTruthSampleIndices(variantContext, badVariantSampleIndices, sampleAlleleCounts,
                           sampleNoCallCounts).isEmpty();
        }
//        return badVariantSampleIndices != null &&
//              !getFilterableTruthSampleIndices(variantContext, badVariantSampleIndices, sampleAlleleCounts,
//                                               sampleNoCallCounts).isEmpty();
        return false;
    }

    /**
     * form numTrios x 3 matrix of allele counts for specified variantIndex (paternal, maternal, child)
     */
    protected int[][] getTrioAlleleCountsMatrix(final int variantIndex) {
        return Arrays.stream(trioSampleIndices)
            .map(
                trioIndices -> Arrays.stream(trioIndices)
                                .map(sampleIndex -> sampleVariantAlleleCounts.getAsInt(variantIndex, sampleIndex))
                                .toArray()
            )
            .toArray(int[][]::new);
    }

    /**
     * form numTrios x 3 matrix of no-call counts for specified variantIndex (paternal, maternal, child)
     */
    protected int[][] getTrioNoCallCountsMatrix(final int variantIndex) {
        return Arrays.stream(trioSampleIndices)
                .map(
                        trioIndices -> Arrays.stream(trioIndices)
                                .map(sampleIndex -> sampleVariantNoCallCounts.getAsInt(variantIndex, sampleIndex))
                                .toArray()
                )
                .toArray(int[][]::new);
    }

    /**
     * form numTrios x 3 matrix of genotype qualities for specified variantIndex (paternal, maternal, child)
     */
    protected short[][] getTrioCallQualitiesMatrix(final int variantIndex) {
        final short[] sampleCallQualities =  sampleVariantCallQualities.values[variantIndex];
        return Arrays.stream(trioSampleIndices)
            .map(trioIndices -> new short[] {sampleCallQualities[trioIndices[0]],
                                             sampleCallQualities[trioIndices[1]],
                                             sampleCallQualities[trioIndices[2]]}
            )
            .toArray(short[][]::new);
    }


    protected boolean getSampleVariantIsFilterable(final int variantIndex, final int sampleIndex) {
        return genotypeIsFilterable(
                sampleVariantAlleleCounts.getAsInt(variantIndex, sampleIndex),
                sampleVariantNoCallCounts.getAsInt(variantIndex, sampleIndex)
        );
    }

    /**
     * A genotype is filterable if its allele count is filterable and it isn't fully no-call
     */
    protected boolean genotypeIsFilterable(final int alleleCount, final int noCallCount) {
        return noCallCount < 2 && alleleCountIsFilterable(alleleCount);
    }

    /**
     * An allele count is filterable if it's HOMREF or HET, or if it's HOMVAR and keepHomvar is not true,
     */
    protected boolean alleleCountIsFilterable(final int alleleCount) {
        return !keepHomvar || alleleCount < 2;
    }


    protected Set<Integer> getInheritanceTrainableSampleIndices(final int variantIndex) {
        final byte[] sampleAlleleCounts = sampleVariantAlleleCounts.values[variantIndex];
        final byte[] sampleNoCallCounts = sampleVariantNoCallCounts.values[variantIndex];
        return Arrays.stream(trioSampleIndices)
                .filter(
                        trio -> trioPasses(
                                sampleAlleleCounts[trio[FATHER_IND]], sampleAlleleCounts[trio[MOTHER_IND]],
                                sampleAlleleCounts[trio[CHILD_IND]],
                                sampleNoCallCounts[trio[FATHER_IND]], sampleNoCallCounts[trio[MOTHER_IND]],
                                sampleNoCallCounts[trio[CHILD_IND]]
                        )
                )
                .flatMapToInt(Arrays::stream)
                .boxed()
                .collect(Collectors.toSet());
    }

    protected Set<Integer> getGoodSampleIndices(final int variantIndex) {
        return goodVariantSampleIndices.getOrDefault(variantIds.get(variantIndex), Collections.emptySet());
    }

    protected Set<Integer> getBadSampleIndices(final int variantIndex) {
        return badVariantSampleIndices.getOrDefault(variantIds.get(variantIndex), Collections.emptySet());
    }


    private List<Set<Integer>> filterableButNotTrainable = null;
    /**
     * A sample Variant is trainable if a) it is filterable, and b) there is truth data (inheritance or known good/bad)
     */
    protected boolean getSampleVariantIsTrainable(final int variantIndex, final int sampleIndex) {
        if(filterableButNotTrainable == null) {
            filterableButNotTrainable = IntStream.range(0, numVariants)
                .mapToObj(
                    vIndex -> {
                        final Set<Integer> inheritanceTrainable = getInheritanceTrainableSampleIndices(vIndex);
                        final Set<Integer> goodSampleIndices = getGoodSampleIndices(vIndex);
                        final Set<Integer> badSampleIndices = getBadSampleIndices(vIndex);
                        return IntStream.range(0, numSamples)
                            .filter(sIndex -> getSampleVariantIsFilterable(vIndex, sIndex))
                            .filter(sIndex -> !inheritanceTrainable.contains(sIndex))
                            .filter(sIndex -> !goodSampleIndices.contains(sIndex))
                            .filter(sIndex -> !badSampleIndices.contains(sIndex))
                            .boxed()
                            .collect(Collectors.toSet());
                    }
                )
                .collect(Collectors.toList());
        }
        if(getSampleVariantIsFilterable(variantIndex, sampleIndex)) {
            // don't train on data with any no-call counts, or with no training data
            // return !filterableButNotTrainable.get(variantIndex).contains(sampleIndex);
            return trainHomref || sampleVariantAlleleCounts.getAsInt(variantIndex, sampleIndex) > 0;
        } else {
            return false;
        }
    }

    private void getTrackProperties(final TrackOverlapDetector trackOverlapDetector,
                                    final VariantContext variantContext) {
        final List<SimpleInterval> genomeTrackOverlapLocations =
            SimpleSvInterval.streamFrom(variantContext)
                    .flatMap(SimpleSvInterval::streamGenomeTrackOverlapLocations)
                    .distinct()
                    .collect(Collectors.toList());

        final double[] overlaps;
        try {
            if (trackOverlapDetector.hasOther()) { // This track has paired intervals
                overlaps = genomeTrackOverlapLocations.stream()
                    .mapToDouble(location -> trackOverlapDetector.getPrimaryOverlapFraction(location) +
                                             trackOverlapDetector.getOtherOverlapfraction(location))
                    .toArray();
                // check for spanning by streaming over all possible pairs
                final boolean spans = IntStream.range(0, genomeTrackOverlapLocations.size() - 1)
                    .anyMatch(
                        i -> IntStream.range(i + 1, genomeTrackOverlapLocations.size()).anyMatch(
                            j -> trackOverlapDetector.spansPrimaryAndOther(
                                genomeTrackOverlapLocations.get(i),
                                genomeTrackOverlapLocations.get(j)
                            )
                        )
                    );
                propertiesTable.append(
                        trackOverlapDetector.getName() + "_spans",
                        spans
                );
            } else {  // this track has simple intervals
                overlaps = genomeTrackOverlapLocations.stream()
                        .mapToDouble(trackOverlapDetector::getPrimaryOverlapFraction)
                        .toArray();
            }
        } catch(IllegalArgumentException illegalArgumentException) {
            throw new IllegalArgumentException(
                "Error getting overlap for " + variantContext.getID() + " and " + trackOverlapDetector.getName(),
                illegalArgumentException
            );
        }
        propertiesTable.append(
            trackOverlapDetector.getName() + "_min",
            Arrays.stream(overlaps).min().orElse(0.0)
        );
        propertiesTable.append(
            trackOverlapDetector.getName() + "_max",
            Arrays.stream(overlaps).max().orElse(0.0)
        );
    }

    private static short getGenotypeAttributeAsShort(final Genotype genotype, final String key, Short defaultValue) {
        if(key.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
            return (short)genotype.getGQ();
        }
        Object x = genotype.getExtendedAttribute(key);
        if (x == null ||
            x == VCFConstants.MISSING_VALUE_v4 ||
            (x instanceof java.lang.Character && (char)x == VCFConstants.NO_CALL_ALLELE)) {
            if(defaultValue == null) {
                throw new IllegalArgumentException("Genotype is missing value of " + key);
            } else {
                return defaultValue;
            }
        }
        if ( x instanceof Short ) return (Short)x;
        try {
            //noinspection ConstantConditions
            return Short.parseShort((String) x); // throws an exception if this isn't a string
        } catch(ClassCastException classCastException) {
            throw new GATKException("Unable to extract value " + x + " for key " + key);
        }
    }

    @SuppressWarnings("SameParameterValue")
    private short[] getGenotypeAttributeAsShort(final Iterable<Genotype> sampleGenotypes, final String attributeKey,
                                                final Short missingAttributeValue) {
        short[] values = new short[sampleIds.size()];
        int index = 0;
        for(final Genotype genotype : sampleGenotypes) {
            values[index] = getGenotypeAttributeAsShort(genotype, attributeKey, missingAttributeValue);
            ++index;
        }
        return values;
    }

    private static Set<String> getVcfPropertyAsStringSet(Object rawProperty) {
        if ( rawProperty == null ||
                rawProperty == VCFConstants.EMPTY_ALLELE ||
                (rawProperty instanceof java.lang.Character && (char)rawProperty == VCFConstants.NO_CALL_ALLELE) ) {
            return Collections.emptySet();
        } else {
            // Throws an exception if this isn't a String
            if (rawProperty instanceof String) {
                // Remove any weird brackets, and split by commas
                return Arrays.stream(((String)rawProperty)
                        .replaceAll("[\\[\\] ]", "")
                        .split(",")
                ).collect(Collectors.toSet());
            } else {
                throw new IllegalArgumentException("Value was not a String");
            }
        }
    }

    private static boolean isMissing(final List<String> values) {
        return values == null || values.isEmpty() || (
                values.size() == 1 && (values.get(0) == null || values.get(0).equals(VCFConstants.EMPTY_ALLELE))
        );
    }

    private static Set<String> getInfoAttributeAsStringSet(final VariantContext variantContext, final String key) {
        try {
            final List<String> values = variantContext.getAttributeAsStringList(key, null);
            return isMissing(values) ? Collections.emptySet() : new HashSet<>(values);
        } catch(IllegalArgumentException illegalArgumentException) {
            throw new IllegalArgumentException(
                "Getting " + key + " from INFO field: " + illegalArgumentException.getMessage(),
                illegalArgumentException
            );
        }
    }

    private static Set<String> getGenotypeAttributeAsStringSet(final Genotype genotype, final String key) {
        Object rawProperty = genotype.getExtendedAttribute(key);
        try {
            return getVcfPropertyAsStringSet(rawProperty);
        } catch(IllegalArgumentException illegalArgumentException) {
            throw new IllegalArgumentException(
                "Getting " + key + " for sample " + genotype.getSampleName() + ": " + illegalArgumentException.getMessage(),
                illegalArgumentException
           );
        }
    }

    @SuppressWarnings("SameParameterValue")
    private Set<String>[] getGenotypeAttributeAsStringSet(final Iterable<Genotype> sampleGenotypes,
                                                          final String key) {
        @SuppressWarnings("unchecked")
        final Set<String>[] values = (Set<String>[]) new Set<?>[numSamples];
        int index = 0;
        for(final Genotype genotype : sampleGenotypes) {
            values[index] = getGenotypeAttributeAsStringSet(genotype, key);
            ++index;
        }

        return values;
    }

    private Set<Integer> getFilterableTruthSampleIndices(final VariantContext variantContext,
                                                         final Map<String, Set<Integer>> truthSampleIndices,
                                                         final byte[] sampleAlleleCounts,
                                                         final byte[] sampleNoCallCounts) {
        if(truthSampleIndices.containsKey(variantContext.getID())) {
            // Get sample IDs of known bad variants that are filterable
            return truthSampleIndices.get(variantContext.getID())
                    .stream()
                    .filter(sampleIndex -> genotypeIsFilterable(sampleAlleleCounts[sampleIndex],
                                                                sampleNoCallCounts[sampleIndex]))
                    .collect(Collectors.toSet());
        }
        return Collections.emptySet();
    }


    private class AlleleCountsGetter {
        final byte nonRefCounts;
        final byte numCalledCounts;
        final byte noCallCounts;
        final short callQuality;
        final boolean isCopyNumberCall;

        AlleleCountsGetter(final Genotype genotype) {
            byte numAlleles = 0;
            byte noCallCounts = 0;
            byte nonRefCounts = 0;
            short callQuality = MISSING_GQ_VAL;
            for(final Allele allele : genotype.getAlleles()) {
                // Do normal counting of alleles from genotype field
                ++numAlleles;
                if(allele.isNoCall()) {
                    ++noCallCounts;
                } else if(!allele.isReference()) {
                    ++nonRefCounts;
                }
            }

            boolean isCopyNumberCall = false;
            if(useCopyNumberCalls && genotype.isNoCall()) {
                // Check if maybe this is a copy number call
                // Get the value from this key, if it's available
                final short rd_cn = getGenotypeAttributeAsShort(genotype, RD_CN_KEY, (short)-1);
                if(rd_cn != -1) { // got copy number info, use it
                    nonRefCounts = (byte)(rd_cn == 2 ? 0 : 1);
                    noCallCounts = 0;
                    numAlleles = 2;
                    isCopyNumberCall = true;
                    callQuality = getGenotypeAttributeAsShort(genotype, RD_GQ_KEY, MISSING_GQ_VAL);
                }
            }
            if(!isCopyNumberCall) {
                callQuality = getGenotypeAttributeAsShort(genotype, VCFConstants.GENOTYPE_QUALITY_KEY, MISSING_GQ_VAL);
            }
            this.nonRefCounts = nonRefCounts;
            this.numCalledCounts = (byte)(numAlleles - noCallCounts);
            this.noCallCounts = noCallCounts;
            this.callQuality = callQuality;
            this.isCopyNumberCall = isCopyNumberCall;
        }
    }

    /**
     * Accumulate properties for variant matrix, and allele counts, genotype quality for trio tensors
     */
    @Override
    public void apply(VariantContext variantContext, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {
        try {
            wrapped_apply(variantContext);
        } catch(Exception exception) {
            throw new IllegalArgumentException(
                "Error processing variant " + variantContext.getID() + ": " + exception.getMessage(),
                exception
            );
        }
    }

    private void wrapped_apply(VariantContext variantContext) {
        /////////////////////////////////////// Process allele counts and quality //////////////////////////////////////
        int numCalledAlleles = 0;
        int numNonRefAlleles = 0;
        int numVariantInputVar = 0;
        int numVariantInputNoCall = 0;
        int numVariantInputRef = 0;
        final byte[] sampleAlleleCounts = new byte[numSamples];
        final byte[] sampleNoCallCounts = new byte[numSamples];
        final boolean[] isCopyNumberCall = new boolean[numSamples];
        final short[] callQuality = new short[numSamples];
        int sampleIndex = 0;
        for(final Genotype genotype : variantContext.getGenotypesOrderedBy(sampleIds)) {
            final AlleleCountsGetter alleleCountsGetter = new AlleleCountsGetter(genotype);
            sampleAlleleCounts[sampleIndex] = alleleCountsGetter.nonRefCounts;
            sampleNoCallCounts[sampleIndex] = alleleCountsGetter.noCallCounts;
            isCopyNumberCall[sampleIndex] = alleleCountsGetter.isCopyNumberCall;
            callQuality[sampleIndex] = alleleCountsGetter.callQuality;
            ++sampleIndex;

            if(alleleCountsGetter.nonRefCounts > 0) {
                ++numVariantInputVar;
            } else if(alleleCountsGetter.noCallCounts > 0) {
                ++numVariantInputNoCall;
            } else {
                ++numVariantInputRef;
            }
            numNonRefAlleles += alleleCountsGetter.nonRefCounts;
            numCalledAlleles += alleleCountsGetter.numCalledCounts;
        }
        numInputVar += numVariantInputVar;
        numInputNoCall += numVariantInputNoCall;
        numInputRef += numVariantInputRef;

        if(runMode == RunMode.Train && !getVariantIsTrainable(variantContext, sampleAlleleCounts, sampleNoCallCounts)) {
            // no need to train on unfilterable variants
            return;
        }
        ++numVariants;
        /////////////////////////////////////////// Get values from FORMAT /////////////////////////////////////////////
        propertiesTable.append(VCFConstants.ALLELE_COUNT_KEY, sampleAlleleCounts);
        propertiesTable.append(NO_CALL_COUNTS_KEY, sampleNoCallCounts);
        if(useCopyNumberCalls) {
            // if using copy number calls, add the reconciled call quality, and whether this is a copy number call
            propertiesTable.append(IS_COPY_NUMBER_CALL_KEY, isCopyNumberCall);
            propertiesTable.append(CALL_QUALITY_KEY, callQuality);
        }
        final Iterable<Genotype> sampleGenotypes = variantContext.getGenotypesOrderedBy(sampleIds);
        // Append various qualities
        Stream.of(VCFConstants.GENOTYPE_QUALITY_KEY, PE_GQ_KEY, SR_GQ_KEY, RD_GQ_KEY)
                .filter(filterKey -> vcfHeaderIds.contains(filterKey))
                .forEach(
                    genotypeAttributeKey -> propertiesTable.append(
                        genotypeAttributeKey,
                        getGenotypeAttributeAsShort(sampleGenotypes, genotypeAttributeKey, MISSING_GQ_VAL)
                    )
                );
        // also get EV, set of per-genotype evidence, CONC_ST, genotype concordance contingency state
        Stream.of(EV_KEY, CONC_ST_KEY)
                .filter(formatKey -> vcfHeaderIds.contains(formatKey))
                .forEach(
                    formatKey -> propertiesTable.append(formatKey, getGenotypeAttributeAsStringSet(sampleGenotypes, formatKey))
                );
        //////////////////////////////////////// Get or estimate allele frequency //////////////////////////////////////
        float alleleFrequency = (float)variantContext.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.0);
        if(alleleFrequency <= 0) {
            if(variantContext.getNSamples() <= minSamplesToEstimateAlleleFrequency) {
                throw new GATKException("VCF does not have " + VCFConstants.ALLELE_FREQUENCY_KEY + " annotated or enough samples to estimate it ("
                                        + minSamplesToEstimateAlleleFrequencyKey + "=" + minSamplesToEstimateAlleleFrequency + " but there are "
                                        + variantContext.getNSamples() + " samples)");
            }
            // VCF not annotated with allele frequency, guess it from allele counts. If somehow we have a variant with
            // no called alleles, just set alleleFrequency to 0. There's nothing to do with it anyway, so it hardly
            // matters
            alleleFrequency = numCalledAlleles > 0 ? numNonRefAlleles / (float) numCalledAlleles : 0F;
        }
        propertiesTable.append(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFrequency);
        ///////////////////////////////////////////////// Get FILTER ///////////////////////////////////////////////////
        propertiesTable.append(FILTER_KEY, variantContext.getFilters());
        ////////////////////////////////////////////// Get INFO fields /////////////////////////////////////////////////
        // SVTYPE is mandatory
        final String svType = variantContext.getAttributeAsString(VCFConstants.SVTYPE, null);
        if(svType == null) {
            throw new GATKException("Missing " + VCFConstants.SVTYPE + " for variant " + variantContext.getID());
        }
        propertiesTable.append(VCFConstants.SVTYPE, svType);
        // Sometimes BNDs have no SVLEN, in that case set SVLEN to -1
        propertiesTable.append(SVLEN_KEY, SimpleSvInterval.getOrInferSvLen(variantContext));
        // get INFO fields that are sets of strings:
        Stream.of(ALGORITHMS_KEY, EVIDENCE_KEY, STATUS_KEY)
                .filter(infoKey -> vcfHeaderIds.contains(infoKey))
                .forEach(
                    infoKey -> propertiesTable.append(infoKey, getInfoAttributeAsStringSet(variantContext, infoKey))
                );
        // add INFO fields that are floats
        Stream.of(NON_REF_GENOTYPE_CONCORDANCE_KEY, VAR_PPV_KEY, VAR_SENSITIVITY_KEY, TRUTH_AF_KEY)
                .filter(infoKey -> vcfHeaderIds.contains(infoKey))
                .forEach(
                    infoKey -> propertiesTable.append(infoKey,
                                                     (float)variantContext.getAttributeAsDouble(infoKey, Double.NaN))
                );

        /////////////////////////////////////////  Get genome track overlaps ///////////////////////////////////////////
        for(final TrackOverlapDetector trackOverlapDetector : trackOverlapDetectors) {
            getTrackProperties(trackOverlapDetector, variantContext);
        }

        //////////////////////// Wrap up: store data if training, filter variant context if filtering //////////////////
        if(runMode == RunMode.Train) {
            variantIds.add(variantContext.getID());
        } else {
            if(sampleVariantCallQualities == null) {
                variantIds.add(variantContext.getID());
                sampleVariantCallQualities = useCopyNumberCalls ?
                    (PropertiesTable.ShortMatProperty) propertiesTable.get(CALL_QUALITY_KEY) :
                    (PropertiesTable.ShortMatProperty) propertiesTable.get(VCFConstants.GENOTYPE_QUALITY_KEY);
                sampleVariantAlleleCounts = (PropertiesTable.ByteMatProperty) propertiesTable.get(VCFConstants.ALLELE_COUNT_KEY);
                sampleVariantNoCallCounts = (PropertiesTable.ByteMatProperty) propertiesTable.get(NO_CALL_COUNTS_KEY);
                propertiesTable.validateAndFinalize();
            } else {
                variantIds.set(0, variantContext.getID()); // only one variant ID at a time in context when filtering
                propertiesTable.oneHot();
            }

            final VariantContext filteredVariantContext = filterVariantContext(variantContext);
            if(filteredVariantContext != null) {
                vcfWriter.add(filteredVariantContext);
            }

            propertiesTable.clearRows();
        }
    }

    Genotype[] filterGenotypes = null;
    float[] variantPropertiesForFilterVariantContext = null;
    boolean[] sampleVariantFilterableForFilterVariantContext = null;
    private VariantContext filterVariantContext(final VariantContext variantContext) {
        final int numProperties = getNumProperties();
        final int numSamples;
        if(filterGenotypes == null) {
            numSamples = variantContext.getNSamples();
            filterGenotypes = new Genotype[numSamples];
            variantPropertiesForFilterVariantContext = new float[numSamples * numProperties];
            sampleVariantFilterableForFilterVariantContext = new boolean[numSamples];
        } else {
            numSamples = filterGenotypes.length;
            if(numSamples != variantContext.getNSamples()) {
                throw new IllegalArgumentException("At " + variantContext.getID() + ": number of samples changed from "
                                                   + numSamples + " to " + variantContext.getNSamples());
            }
        }

        // filter all samples for this variant context in a batch to avoid DMatrix creation / prediction overhead
        final boolean maybeFilterable = !(keepMultiallelic && getIsMultiallelic(variantContext));
        final short[] adjustedGq, adjustedLogits;
        if(maybeFilterable) {
            int numFilterableSamples = 0;
            int sampleIndex = 0;
            for (final Genotype genotype : variantContext.getGenotypes()) {
                filterGenotypes[sampleIndex] = genotype;
                sampleVariantFilterableForFilterVariantContext[sampleIndex] =
                        getSampleVariantIsFilterable(0, sampleIndex);
                if (sampleVariantFilterableForFilterVariantContext[sampleIndex]) {
                    ++numFilterableSamples;
                }
                ++sampleIndex;
            }
            numFilterableGenotypes += numFilterableSamples;
            final short[][] adjustedQuality = adjustedQualityBatch();
            adjustedGq = adjustedQuality[0];
            adjustedLogits = adjustedQuality[1];
        } else {
            adjustedGq = null;
            adjustedLogits = null;
            int sampleIndex = 0;
            for(final Genotype genotype : variantContext.getGenotypes()) {
                filterGenotypes[sampleIndex] = genotype;
                sampleVariantFilterableForFilterVariantContext[sampleIndex] = false;
                ++sampleIndex;
            }
        }

        int numFiltered = 0;
        int predictSampleIndex = 0;
        int numVariantOutputVar = 0;
        for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
            final Genotype genotype = filterGenotypes[sampleIndex];
            final short callQuality;
            final short scaledLogits;
            if(sampleVariantFilterableForFilterVariantContext[sampleIndex]) {
                //noinspection ConstantConditions
                callQuality = adjustedGq[predictSampleIndex];
                //noinspection ConstantConditions
                scaledLogits = adjustedLogits[predictSampleIndex];
                ++predictSampleIndex;
            } else {
                // Still need to compute scaledLogits. Also: don't allow no-calls to have a passing quality
                callQuality = maybeFilterable && filterGenotypes[sampleIndex].isNoCall() ?
                    (short)Integer.min(sampleVariantCallQualities.values[0][sampleIndex], FAILING_GQ) :
                    sampleVariantCallQualities.values[0][sampleIndex];
                scaledLogits = probToScaledLogits(1.0 - phredToProb(callQuality));
            }
            final boolean needsFilter = sampleVariantFilterableForFilterVariantContext[sampleIndex] &&
                                        scaledLogits < minScaledLogits;
            final Genotype filteredGenotype = getFilteredGenotype(genotype, callQuality, scaledLogits, needsFilter);
            if(needsFilter) {
                ++numFiltered;
            }

            filterGenotypes[sampleIndex] = filteredGenotype;
            final AlleleCountsGetter alleleCountsGetter = new AlleleCountsGetter(filteredGenotype);
            if (alleleCountsGetter.nonRefCounts > 0) {
                ++numVariantOutputVar;
            } else if (alleleCountsGetter.noCallCounts > 0) {
                ++numOutputNoCall;
            } else {
                ++numOutputRef;
            }
        }
        numOutputVar += numVariantOutputVar;
        numFilteredGenotypes += numFiltered;
        if(maybeFilterable && numVariantOutputVar == 0 && !keepNonVar) {
            ++numDroppedVariants;
            return null;  // This variant has no called non-ref alleles. Drop it!
        }

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variantContext)
            .genotypes(filterGenotypes).attribute(MIN_QUALITY_KEY, minScaledLogits);
        if(numFiltered > variantContext.getNSamples() * reportMinGqFilterThreshold) {
            variantContextBuilder.filter(EXCESSIVE_MIN_QUALITY_FILTER_KEY);
        }

        return variantContextBuilder.make();
    }

    private Genotype getFilteredGenotype(final Genotype inputGenotype, final short callQuality,
                                         final short scaledLogits, final boolean needsFilter) {
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(inputGenotype)
                .attribute(scaledLogitProperty, scaledLogits);
        if(needsFilter) {
            // Set GT and copy-number calls to no-call
            genotypeBuilder.alleles(GATKVariantContextUtils.noCallAlleles(inputGenotype.getPloidy()));
            if(inputGenotype.hasExtendedAttribute(RD_CN_KEY)) {
                genotypeBuilder.attribute(RD_CN_KEY, VCFConstants.NO_CALL_ALLELE);
            }
        }
        if(inputGenotype.hasGQ()) {
            genotypeBuilder.GQ(callQuality);
        }
        if(inputGenotype.hasExtendedAttribute(RD_GQ_KEY)) {
            genotypeBuilder.attribute(RD_GQ_KEY, callQuality);
        }
        return genotypeBuilder.make();
    }


    private boolean notCloseEnough(final String strRep, final double val,
                                   @SuppressWarnings("SameParameterValue") final double tol) {
        return FastMath.abs(Double.parseDouble(strRep) - val) > FastMath.abs(tol * val);
    }

    static String padWidth(final String baseString, final int width) {
        return baseString.length() >= width ?
            baseString :
            baseString + String.join("", Collections.nCopies(width - baseString.length(), " "));
    }

    static private String edgeValueToString(final double binEdge, final int precision) {
        return String.format("%." + precision + "f", binEdge);
    }

    static private String[] edgeValuesToStrings(final double[] binEdges, final int precision) {
        return Arrays.stream(binEdges)
                .mapToObj(binEdge -> edgeValueToString(binEdge, precision))
                .toArray(String[]::new);
    }

    static private boolean edgeStringsRepeat(final String[] edgeStrings) {
        return IntStream.range(0, edgeStrings.length - 1)
            .anyMatch(i -> edgeStrings[i].equals(edgeStrings[i + 1]));
    }

    private String[] getPropertyBinNames( final PropertiesTable.Property property, final double[] bins) {
        double minVal = bins[0];
        double maxVal = bins[bins.length - 1];
        for(int variantIndex = 0; variantIndex < property.getNumRows(); ++variantIndex) {
            final float val = property.getAsFloat(variantIndex, 0, false);
            if(val < minVal) {
                minVal = val;
            } else if(val > maxVal) {
                maxVal = val;
            }
        }
        final boolean addLeftEdge = minVal < bins[0];
        final boolean addRightEdge = maxVal > bins[bins.length - 1];
        final double[] fullBinEdges = new double[bins.length + (addLeftEdge ? 1 : 0) + (addRightEdge ? 1: 0)];
        System.arraycopy(bins, 0, fullBinEdges, addLeftEdge ? 1 : 0, bins.length);
        if(addLeftEdge) {
            fullBinEdges[0] = minVal;
        }
        if(addRightEdge) {
            fullBinEdges[fullBinEdges.length - 1] = maxVal;
        }

        // increase precision of string representation of all bin edges until there is no confusion
        int precision = 0;
        String[] binEdgeNames = edgeValuesToStrings(fullBinEdges, precision);
        while(edgeStringsRepeat(binEdgeNames)) {
            ++precision;
            binEdgeNames = edgeValuesToStrings(fullBinEdges, precision);
        }
        // increase precision of string representation of individual bin edges until they are accurate enough
        for(int i = 0; i < binEdgeNames.length; ++i) {
            int precision_i = precision;
            while (notCloseEnough(binEdgeNames[i], fullBinEdges[i], 0.1)) {
                ++precision_i;
                binEdgeNames[i] = edgeValueToString(fullBinEdges[i], precision_i);
            }
        }
        // form the bin names, keeping track of the width
        int maxWidth = property.name.length();
        String[] propertyBinNames = new String[binEdgeNames.length - 1];
        for(int i = 0; i < propertyBinNames.length; ++i) {
            propertyBinNames[i] = binEdgeNames[i] + "-" + binEdgeNames[i + 1];
            maxWidth = max(maxWidth, propertyBinNames[i].length());
        }
        // pad bin names so they are all the same length, for easy table viewing
        for(int i = 0; i < propertyBinNames.length; ++i) {
            propertyBinNames[i] = padWidth(propertyBinNames[i], maxWidth);
        }

        return propertyBinNames;
    }

    private void setPropertyBins() {
        propertyBins = new int[numVariants]; // guaranteed to be all zeros by Java spec

        // Bin variants by SVTYPE
        propertyBinLabelColumns = new ArrayList<>(1 + propertyBinsMap.size());

        propertyBinLabelColumns.add(VCFConstants.SVTYPE);
        final PropertiesTable.StringArrProperty svType =
            (PropertiesTable.StringArrProperty) propertiesTable.get(VCFConstants.SVTYPE);
        final List<String> allSvTypes = svType.getAllLabels();

        final Map<String, Integer> binNames = IntStream.range(0, allSvTypes.size())
            .boxed()
            .collect(Collectors.toMap(allSvTypes::get, i -> i));
        IntStream.range(0, numVariants).forEach(variantIndex ->
            propertyBins[variantIndex] = binNames.get(svType.getAsString(variantIndex))
        );

        // Successively refine bins by each element of propertyBinsMap
        for(Map.Entry<String, double[]> propertyBinsEntry : propertyBinsMap.entrySet()) {
            final String propertyName = propertyBinsEntry.getKey();
            final double[] bins = propertyBinsEntry.getValue();

            // get the old bin names as an ordered list, so that existing propertyBins map correctly to oldBinNames
            final List<String> oldBinNames = binNames.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByValue())
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
            binNames.clear();

            // get the property values, and strings representing the bins we'll use
            final PropertiesTable.Property property = propertiesTable.get(propertyName);
            final String[] propertyBinNames = getPropertyBinNames(property, bins);
            // add property name to the header with extra padding if needed
            propertyBinLabelColumns.add(padWidth(propertyName, propertyBinNames[0].length()));

            IntStream.range(0, numVariants).forEach(variantIndex -> {
                // for each variant, find which of the new bins it should index into
                int propertyBin = Arrays.binarySearch(bins, property.getAsFloat(variantIndex, 0, false));
                if(propertyBin < 0) { propertyBin = ~propertyBin; }
                // form the new bin name by adding a new column with the appropriate binned value
                final String newBinName = oldBinNames.get(propertyBins[variantIndex]) + "\t" + propertyBinNames[propertyBin];
                if(binNames.containsKey(newBinName)) {  // this bin exists, just get its index
                    propertyBins[variantIndex] = binNames.get(newBinName);
                } else {  // this is a new bin, create it and get its index
                    final int newBin = binNames.size();
                    binNames.put(newBinName, newBin);
                    propertyBins[variantIndex] = newBin;
                }
            });
        }

        numPropertyBins = binNames.size();
        propertyBinLabels = new String[numPropertyBins];
        binNames.forEach((key, value) -> propertyBinLabels[value] = key);

        // Finally, want the property bin descriptions to be in sorted order for easier reading out fitting results
        final int[] sortInds = IntStream.range(0, numPropertyBins)
            .boxed()
            .sorted(Comparator.comparing(i -> propertyBinLabels[i]))
            .mapToInt(Integer::new)
            .toArray();
        propertyBinLabels = IntStream.range(0, numPropertyBins)
            .mapToObj(i -> propertyBinLabels[sortInds[i]])
            .toArray(String[]::new);

        final int[] unsortInds = new int[numPropertyBins];
        IntStream.range(0, numPropertyBins).forEach(i -> unsortInds[sortInds[i]] = i);
        propertyBins = Arrays.stream(propertyBins)
            .map(i -> unsortInds[i])
            .toArray();
    }

    private void setPropertyBinIsLargeAlleleFraction() {
        propertyBinIsLargeAlleleFraction = new boolean[numPropertyBins];
        final boolean[] isSet = new boolean[numPropertyBins];
        int numSet = 0;
        for(int variantIndex = 0; variantIndex < numVariants; ++variantIndex) {
            final int propertyBin = propertyBins[variantIndex];
            if(!isSet[propertyBin]) {
                final float alleleFrequency = getAlleleFrequency(variantIndex);
                propertyBinIsLargeAlleleFraction[propertyBin] = alleleFrequency > maxInheritanceAf;
                isSet[propertyBin] = true;
                ++numSet;
                if(numSet == numPropertyBins) {
                    return;
                }
            }
        }
    }

    private IntStream streamFilterableGq(final int variantIndex) {
        return IntStream.range(0, numSamples)
                .filter(
                    sampleIndex -> getSampleVariantIsFilterable(variantIndex, sampleIndex)
                )
                .map(sampleIndex -> sampleVariantCallQualities.getAsInt(variantIndex, sampleIndex));
    }

    private IntStream streamFilterableGq() {
        return IntStream.range(0, numVariants).flatMap(this::streamFilterableGq);
    }

    protected Stream<MinGq> getCandidateMinGqs(final int variantIndex) {
        final Stream.Builder<Short> hetBuilder = Stream.builder();
        final Stream.Builder<Short> homvarBuilder = Stream.builder();
        final byte[] sampleAlleleCounts = sampleVariantAlleleCounts.values[variantIndex];
        final short[] sampleGqs = sampleVariantCallQualities.values[variantIndex];
        for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
            final byte alleleCount = sampleAlleleCounts[sampleIndex];
            if(alleleCount <= 1) {
                hetBuilder.add(sampleGqs[sampleIndex]);
            } else if(!keepHomvar) {
                homvarBuilder.add(sampleGqs[sampleIndex]);
            }
        }
        final List<Short> candidateHetGq = hetBuilder.build().sorted().distinct().collect(Collectors.toList());
        final List<Short> candidateHomVarGq = homvarBuilder.build().sorted().distinct().collect(Collectors.toList());
        if(candidateHetGq.isEmpty()) {
            if(candidateHomVarGq.isEmpty()) {
                return null; // no candidate minGq, this variant can't really be filtered
            } else {  // can filter some HOMVAR GQs, so just make trivial HET filter
                candidateHetGq.add(Short.MIN_VALUE);
            }
        } else {
            // consider filtering out all HET variants by adding 1 to highest filterable GQ value
            candidateHetGq.add((short) (candidateHetGq.get(candidateHetGq.size() - 1) + 1));
        }
        if(!candidateHomVarGq.isEmpty()) {
            // consider filtering out all HOMVAR variants by adding 1 to highest filterable GQ value
            candidateHomVarGq.add((short) (candidateHomVarGq.get(candidateHomVarGq.size() - 1) + 1));
        }
        return candidateHetGq.stream()
            .map(
                hetGq -> {
                    final List<MinGq> hetList = candidateHomVarGq.stream()
                        .filter(homVarGq -> homVarGq >= hetGq)
                        .map(homVarGq -> new MinGq(hetGq, homVarGq))
                        .collect(Collectors.toList());
                    if(hetList.isEmpty()) {
                        hetList.add(new MinGq(hetGq, hetGq));
                    }
                    return hetList;
                }
            )
            .flatMap(Collection::stream);
    }

    private boolean isMendelian(final int[] trioAc, final int[] trioNoCalls) {
        return isMendelian(trioAc[FATHER_IND], trioAc[MOTHER_IND], trioAc[CHILD_IND],
                           trioNoCalls[FATHER_IND], trioNoCalls[MOTHER_IND], trioNoCalls[CHILD_IND]);
    }

    private boolean isMendelian(final int fatherAc, final int motherAc, final int childAc,
                                final int fatherNoCalls, final int motherNoCalls, final int childNoCalls) {
        if(fatherNoCalls > 0 || motherNoCalls > 0 || childNoCalls > 0 ||
           (fatherAc == 0 && motherAc == 0 && childAc == 0)) {
            // don't try to figure out what happened if there are no-calls, don't count trios with all HOMREF
            return false;
        } else if(strictMendelian) {
            // child allele counts should not exhibit de-novo mutations nor be missing inherited homvar
            final int maxAc = (fatherAc > 0 ? 1 : 0) + (motherAc > 0 ? 1 : 0);
            final int minAc = fatherAc / 2 + motherAc / 2;
            return (minAc <= childAc) && (childAc <= maxAc);
        } else {
            // child allele counts should not exhibit de-novo mutations
            final int maxAc = fatherAc > 0 || motherAc > 0 ? 2 : 0;
            return childAc <= maxAc;
        }
    }

    private boolean trioPasses(final int[] trioAc, final int[] trioNoCalls) {
        return trioPasses(trioAc[FATHER_IND], trioAc[MOTHER_IND], trioAc[CHILD_IND],
                          trioNoCalls[FATHER_IND], trioNoCalls[MOTHER_IND], trioNoCalls[CHILD_IND]);
    }

    private boolean trioPasses(final int fatherAc, final int motherAc, final int childAc,
                               final int fatherNoCalls, final int motherNoCalls, final int childNoCalls) {
        return (fatherAc > 0 || motherAc > 0 || childAc > 0) &&
                fatherNoCalls == 0 && motherNoCalls == 0 && childNoCalls == 0;
    }


    final void setSampleIndexToPredictionIndex(final Integer[] sampleIndexToPredictionIndex, final int variantIndex) {
        int predictionIndex = 0;
        for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
            if(getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                sampleIndexToPredictionIndex[sampleIndex] = predictionIndex;
                ++predictionIndex;
            } else {
                sampleIndexToPredictionIndex[sampleIndex] = null;
            }
        }
    }

    /**
     * Convert from dLoss / dProb to dLoss / dLogits (assumes prob = sigma(logits))
     * @param prob: probability predicted by classifier
     * @param d1LossDProb:  1st derivative of loss with respect to prob
     * @param d2LossDProb: 2nd derivative of loss with respect to prob
     * @param weight: weight of this prediction
     * @param d1LossDLogits: array to hold 1st derivative of loss with respect to logits
     * @param d2LossDLogits: array to hold 2nd derivative of loss with respect to logits
     * @param predictionIndex: index into arrays of dLoss/dLogits arrays
     */
    protected void setLossDerivs(final double prob, final double d1LossDProb, final double d2LossDProb,
                                 final double weight, final float[] d1LossDLogits, final float[] d2LossDLogits,
                                 final int predictionIndex) {
        final double d1ProbDLogit = prob * (1.0 - prob) / LOGIT_SCALE;
        final double d2ProbDLogit = d1ProbDLogit * (1.0 - 2.0 * prob) / LOGIT_SCALE;
        final double scaledWeight = weight * 4.0 * LOGIT_SCALE; // help scale derivatives to be close to 1.0
        d1LossDLogits[predictionIndex] += (float)(scaledWeight * d1LossDProb * d1ProbDLogit);
        d2LossDLogits[predictionIndex] += (float)(scaledWeight * (d2LossDProb * d1ProbDLogit * d1ProbDLogit
                                                                  + d1LossDProb * d2ProbDLogit));
    }

    private float getAlleleFrequency(final int variantIndex) {
        return propertiesTable.get(VCFConstants.ALLELE_FREQUENCY_KEY).getAsFloat(variantIndex, 0);
    }

    protected FilterSummary getFilterSummary(final MinGq minGq, final int variantIndex, final String label) {
        final int[][] trioAlleleCountsMatrix = getTrioAlleleCountsMatrix(variantIndex);
        final int[][] trioNoCallCountsMatrix = getTrioNoCallCountsMatrix(variantIndex);
        final short[][] trioCallQualitiesMatrix = getTrioCallQualitiesMatrix(variantIndex);
        final float variantAlleleFrequency = getAlleleFrequency(variantIndex);
        final int propertyBin = propertyBins[variantIndex];
        final double variantTruthWeight = propertyBinTruthWeights[propertyBin];
        final double variantInheritanceWeight = propertyBinInheritanceWeights[propertyBin];
        return getFilterSummary(minGq, trioAlleleCountsMatrix, trioNoCallCountsMatrix,
                                trioCallQualitiesMatrix, variantAlleleFrequency, new GoodBadGqs(variantIndex),
                                variantTruthWeight, variantInheritanceWeight, label);
    }

    protected BinnedFilterSummaries getBinnedFilterSummary(final MinGq minGq, final int variantIndex) {
        final int propertyBin = propertyBins[variantIndex];
        return new BinnedFilterSummaries(getFilterSummary(minGq, variantIndex, propertyBinLabels[propertyBin]),
                                         propertyBin, numPropertyBins);
    }

    private void filterTrioCall(final int[] alleleCounts, final int[] noCallCounts, final int sampleIndex,
                                final boolean callPasses) {
        if (!callPasses && alleleCountIsFilterable(alleleCounts[sampleIndex])) {
            alleleCounts[sampleIndex] = 0;
            noCallCounts[sampleIndex] = 2;
        }
    }

    private void filterTrioCallsByMinGq(final int[] alleleCounts, final int[] noCallCounts, final short[] trioGqs,
                                        final MinGq minGq) {
        for(int indexInTrio = 0; indexInTrio < 3; ++indexInTrio) {
            final boolean callPasses = trioGqs[indexInTrio] >= (alleleCounts[indexInTrio] <= 1 ? minGq.minGqHet :
                                                                minGq.minGqHomVar);
            filterTrioCall(alleleCounts, noCallCounts, indexInTrio, callPasses);
        }
    }

    private int countTrioVariantPassing(final int[] trioAcs, final int[] trioNoCalls) {
        int numPassing = 0;
        for(int i = 0; i < 3; ++i) {
            if(trioAcs[i] > 0 && trioNoCalls[i] == 0) {
                ++numPassing;
            }
        }
        return numPassing;
    }

    protected FilterSummary getFilterSummary(final MinGq minGq,
                                             final int[][] trioAlleleCountsMatrix, final int[][] trioNoCallCountsMatrix,
                                             final short[][] trioCallQualitiesMatrix, final float alleleFrequency,
                                             GoodBadGqs goodBadGqs,
                                             final double variantTruthWeight, final double variantInheritanceWeight,
                                             final String label) {
        long numVariantsPassed = 0;
        long numVariants = 0;
        long numMendelianTrios = 0;
        long numMendelianTriosPassed = 0;
        long numTriosPassed = 0;
        long numTruePositives = 0;
        long numFalsePositives = 0;
        long numFalseNegatives = 0;
        long numTrueNegatives = 0;
        for (final int gq : goodBadGqs.goodHetGqs) {
            if(gq >= minGq.minGqHet) {
                ++numTruePositives;
            } else {
                ++numFalseNegatives;
            }
        }
        if(!keepHomvar) {
            for (final int gq : goodBadGqs.goodHomVarGqs) {
                if(gq >= minGq.minGqHomVar) {
                    ++numTruePositives;
                } else {
                    ++numFalseNegatives;
                }
            }
        }
        for(final int gq : goodBadGqs.badHetGqs) {
            if(gq >= minGq.minGqHet) {
                ++numFalsePositives;
            } else {
                ++numTrueNegatives;
            }
        }
        if(!keepHomvar) {
            for(final int gq : goodBadGqs.badHomVarGqs) {
                if(gq >= minGq.minGqHomVar) {
                    ++numFalsePositives;
                } else {
                    ++numTrueNegatives;
                }
            }
        }

        final int[] trioAc = new int[3];
        final int[] trioNoCalls = new int[3];
        for (int trioIndex = 0; trioIndex < numTrios; ++trioIndex) {
            System.arraycopy(trioAlleleCountsMatrix[trioIndex], 0, trioAc, 0, 3);
            System.arraycopy(trioNoCallCountsMatrix[trioIndex], 0, trioNoCalls, 0, 3);
            final short[] trioGq = trioCallQualitiesMatrix[trioIndex];
            final int numVariantsTrio = countTrioVariantPassing(trioAc, trioNoCalls);
            numVariants += numVariantsTrio;
            final boolean trioIsMendelian = isMendelian(trioAc, trioNoCalls);
            if(trioIsMendelian) {
                ++numMendelianTrios;
            }

            filterTrioCallsByMinGq(trioAc, trioNoCalls, trioGq, minGq);
            numVariantsPassed += countTrioVariantPassing(trioAc, trioNoCalls);
            if(trioPasses(trioAc, trioNoCalls)) {
                ++numTriosPassed;
                if(trioIsMendelian) {
                    ++numMendelianTriosPassed;
                }
            }
        }

        return new FilterSummary(minGq, numMendelianTriosPassed, numMendelianTrios, numTriosPassed, numVariants,
                                 numVariantsPassed, numTruePositives, numFalsePositives, numFalseNegatives,
                                 numTrueNegatives,alleleFrequency > maxInheritanceAf,
                                 variantTruthWeight, variantInheritanceWeight, label);
    }

    protected FilterSummary getFilterSummary(final boolean[] samplePasses, final int variantIndex, final int propertyBin) {
        long numVariantsPassed = 0;
        long numVariants = 0;
        long numMendelianTrios = 0;
        long numMendelianTriosPassed = 0;
        long numTriosPassed = 0;
        long numTruePositives = 0;
        long numFalsePositives = 0;
        long numFalseNegatives = 0;
        long numTrueNegatives = 0;

        final Set<Integer> goodSampleIndices = getGoodSampleIndices(variantIndex);
        final Set<Integer> badSampleIndices = getBadSampleIndices(variantIndex);
        // copy these arrays so that a) I can make them ints without taking a huge amount of memory, and
        //                           b) I don't have to make two versions of filterTrioCall()
        final Integer[] sampleIndexToPredictionIndex = new Integer[numSamples];
        setSampleIndexToPredictionIndex(sampleIndexToPredictionIndex, variantIndex);
        for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
            final Integer predictionIndex = sampleIndexToPredictionIndex[sampleIndex];
            if(predictionIndex == null) {
                continue;
            }
            if(goodSampleIndices.contains(sampleIndex)) {
                if(samplePasses[predictionIndex]) {
                    ++numTruePositives;
                } else {
                    ++numFalseNegatives;
                }
            } else if(badSampleIndices.contains(sampleIndex)) {
                if(samplePasses[predictionIndex]) {
                    ++numFalsePositives;
                } else {
                    ++numTrueNegatives;
                }
            }
        }

        final byte[] sampleAlleleCounts = sampleVariantAlleleCounts.values[variantIndex];
        final byte[] sampleNoCallCounts = sampleVariantNoCallCounts.values[variantIndex];
        final int[] trioAc = new int[3];
        final int[] trioNoCalls = new int[3];
        for(int trioIndex = 0; trioIndex < numTrios; ++trioIndex) {
            final int[] sampleIndices = trioSampleIndices[trioIndex];

            // copy trio allele counts / no-call counts into temporary buffers (so they can be filtered cleanly)
            // copy samplePasses just to keep things in the same form
            for(int indexInTrio = 0; indexInTrio < 3; ++indexInTrio) {
                final int sampleIndex = sampleIndices[indexInTrio];
                trioAc[indexInTrio] = sampleAlleleCounts[sampleIndex];
                trioNoCalls[indexInTrio] = sampleNoCallCounts[sampleIndex];
            }
            numVariants += countTrioVariantPassing(trioAc, trioNoCalls);
            final boolean trioIsMendelian = isMendelian(trioAc, trioNoCalls);
            if(trioIsMendelian) {
                ++numMendelianTrios;
            }

            for(int indexInTrio = 0; indexInTrio < 3; ++indexInTrio) {
                final int sampleIndex = sampleIndices[indexInTrio];
                final Integer predictionIndex = sampleIndexToPredictionIndex[sampleIndex];
                filterTrioCall(trioAc, trioNoCalls, indexInTrio,
                        predictionIndex == null || samplePasses[predictionIndex]);
            }

            numVariantsPassed += countTrioVariantPassing(trioAc, trioNoCalls);
            if(trioPasses(trioAc, trioNoCalls)) {
                ++numTriosPassed;
                if(trioIsMendelian) {
                    ++numMendelianTriosPassed;
                }
            }
        }

        final String label = propertyBinLabels[propertyBin];
        final double variantTruthWeight = propertyBinTruthWeights[propertyBin];
        final double variantInheritanceWeight = propertyBinInheritanceWeights[propertyBin];
        final boolean isLargeAlleleFrequency = propertyBinIsLargeAlleleFraction[propertyBin];
        final MinGq minGq = new MinGq(minScaledLogits, minScaledLogits);
        return new FilterSummary(minGq, numMendelianTriosPassed, numMendelianTrios, numTriosPassed, numVariants,
                numVariantsPassed, numTruePositives, numFalsePositives, numFalseNegatives, numTrueNegatives,
                isLargeAlleleFrequency, variantTruthWeight, variantInheritanceWeight, label);
    }

    protected BinnedFilterSummaries getBinnedFilterSummary(final boolean[] samplePasses, final int variantIndex) {
        final int propertyBin = propertyBins[variantIndex];
        return new BinnedFilterSummaries(getFilterSummary(samplePasses, variantIndex, propertyBin),
                                         propertyBin, numPropertyBins);
    }

    static protected class FilterQuality extends FilterSummary implements Comparable<FilterQuality> {
        final FilterLoss loss;

        FilterQuality(final FilterSummary filterSummary) {
            super(filterSummary);
            this.loss = new FilterLoss(filterSummary);
        }

        @Override
        public int compareTo(final MinGqVariantFilterBase.@NotNull FilterQuality other) {
            final int lossCompare = this.loss.compareTo(other.loss);
            // for two equivalent-loss filters, take the more permissive one
            return lossCompare == 0 ? this.minGq.compareTo(other.minGq) : lossCompare;
        }

        static final FilterQuality EMPTY = new FilterQuality(FilterSummary.EMPTY);
    }

    class GoodBadGqs {
        // not worth fighting Java's inability to stream short here
        final int[] goodHetGqs;
        final int[] badHetGqs;
        final int[] goodHomVarGqs;
        final int[] badHomVarGqs;

        GoodBadGqs(final int variantIndex) {
            final int[] emptyGqs = new int[0];
            final Set<Integer> goodSampleIndices = getGoodSampleIndices(variantIndex);
            final Set<Integer> badSampleIndices = getBadSampleIndices(variantIndex);
            final byte[] sampleAlleleCounts = sampleVariantAlleleCounts.values[variantIndex];
            final short[] sampleCallQualities = sampleVariantCallQualities.values[variantIndex];
            this.goodHetGqs = goodSampleIndices == null ?
                emptyGqs :
                goodSampleIndices.stream()
                    .filter(sampleIndex -> sampleAlleleCounts[sampleIndex] == 1)
                    .mapToInt(sampleIndex -> sampleCallQualities[sampleIndex])
                    .toArray();
            this.badHetGqs = badSampleIndices == null ?
                emptyGqs :
                badSampleIndices.stream()
                    .filter(sampleIndex -> sampleAlleleCounts[sampleIndex] == 1)
                    .mapToInt(sampleIndex -> sampleCallQualities[sampleIndex])
                    .toArray();
            this.goodHomVarGqs = keepHomvar || goodSampleIndices == null ?
                emptyGqs :
                goodSampleIndices.stream()
                    .filter(sampleIndex -> sampleAlleleCounts[sampleIndex] > 1)
                    .mapToInt(sampleIndex -> sampleCallQualities[sampleIndex])
                    .toArray();
            this.badHomVarGqs = keepHomvar || badSampleIndices == null ?
                emptyGqs :
                badSampleIndices.stream()
                    .filter(sampleIndex -> sampleAlleleCounts[sampleIndex] > 1)
                    .mapToInt(sampleIndex -> sampleCallQualities[sampleIndex])
                    .toArray();
        }
    }


    protected FilterQuality getOptimalVariantMinGq(final int variantIndex, final FilterSummary backgroundFilterSummary) {
        final Stream<MinGq> candidateMinGqs = getCandidateMinGqs(variantIndex);
        if(candidateMinGqs == null) {
            // minGq doesn't matter for this row, so return previous optimal filter or trivial filter
            return backgroundFilterSummary == null ? FilterQuality.EMPTY : new FilterQuality(backgroundFilterSummary);
        }

        final int[][] trioAlleleCountsMatrix = getTrioAlleleCountsMatrix(variantIndex);
        final int[][] trioNoCallCountsMatrix = getTrioNoCallCountsMatrix(variantIndex);
        final short[][] trioCallQualitiesMatrix = getTrioCallQualitiesMatrix(variantIndex);
        final float variantAlleleFrequency = getAlleleFrequency(variantIndex);
        final String label = propertyBinLabels[propertyBins[variantIndex]];
        if(backgroundFilterSummary == null) {
            // doing optimization only considering loss of each individual variant
            return candidateMinGqs
                    .parallel()
                    .map(minGq -> new FilterQuality(
                            getFilterSummary(
                                minGq, trioAlleleCountsMatrix, trioNoCallCountsMatrix, trioCallQualitiesMatrix,
                                variantAlleleFrequency, new GoodBadGqs(variantIndex), truthWeight,
                                 1.0 - truthWeight, label
                            )
                        ))
                    .min(FilterQuality::compareTo)
                    .orElseThrow(RuntimeException::new);
        } else {
            // doing optimization considering overall loss
            return candidateMinGqs
                    .parallel()
                    .map(minGq -> new FilterQuality(
                            getFilterSummary(
                                minGq, trioAlleleCountsMatrix, trioNoCallCountsMatrix, trioCallQualitiesMatrix,
                                variantAlleleFrequency, new GoodBadGqs(variantIndex), truthWeight,
                                 1.0 - truthWeight, label
                            )
                            .add(backgroundFilterSummary)
                    ))
                    .min(FilterQuality::compareTo)
                    .orElseThrow(RuntimeException::new);
        }
    }

    private static double assignWeightFromLoss(final double loss, final boolean isLargeAlleleFrequency) {
        final double baseWeight = Double.isFinite(loss) ?
                                    FastMath.max(1.0 - 2 * loss, minBinWeight) :
                                    0.0;
        return isLargeAlleleFrequency ? largeAfWeightPenalty * baseWeight : baseWeight;
    }

    private void scalePropertyBinWeights(
        final int propertyBinIndex, final double[] rawBinWeights, final double rawTotalWeight,
        final boolean isTruthWeight, final int numTrainableSampleVariants,
        final float[] goodWeights, final float[] badWeights
    ) {
        // scale weights so that:
        //    -overall average weight across variants is 1.0
        //    -each bin has weight proportional to evidence strength (inheritance or truth)
        //    -each bin has weight inversely proportional to number of variants in that bin, so that rare variant types
        //     are not swamped in the optimization
        //    -passing variants and failing variants have their weights balanced
        int numGoodBin = 0;
        int numBadBin = 0;
        int numTrainable = 0;
        for (final int variantIndex : trainingIndices) {
            final Set<Integer> trainableSampleIndices;
            if(isTruthWeight) {
                trainableSampleIndices = new HashSet<>(
                        getGoodSampleIndices(variantIndex)
                );
                trainableSampleIndices.addAll(
                        getBadSampleIndices(variantIndex)
                );
            } else {
                trainableSampleIndices = getInheritanceTrainableSampleIndices(variantIndex);
            }
            numTrainable += trainableSampleIndices.size();
            if(propertyBins[variantIndex] != propertyBinIndex) {
                continue;
            }
            for(final int trainableSampleIndex : trainableSampleIndices) {
                if(samplePasses(variantIndex, trainableSampleIndex)) {
                    ++numGoodBin;
                } else {
                    ++numBadBin;
                }
            }
        }

        // ws_b = c * wr_b / T_b
        //totSVW = sum_b T_b * ws_b = sum_b c * wr_b = c rawTotalWeight
        // meanSVW = totSVW / sum_b T_b = 1.0
        // c * rawTotalWeight = sum_b T_b
        // c = NT / rawTotalWeight

        final double rawBinWeight = rawBinWeights[propertyBinIndex];
        final double weightScale = isTruthWeight ? truthWeight : 1.0 - truthWeight;
        final double averageBinWeight = rawTotalWeight / numPropertyBins;
        final int numBinVariants = FastMath.max(numGoodBin + numBadBin, 1);
        final double averageBinVariants = numTrainableSampleVariants / (float)numPropertyBins;

//        final double overallScaledWeight = weightScale * rawBinWeight / averageBinWeight
//                                         * averageBinVariants / FastMath.max(numBinVariants, minVariantsForScalingWeight);
        final double overallScaledWeight = weightScale * rawBinWeight * averageBinVariants
                                          / FastMath.max(numBinVariants, minVariantsForScalingWeight);
        rawBinWeights[propertyBinIndex] = overallScaledWeight;
        final double goodWeight, badWeight;
        if(numGoodBin == 0 || numBadBin == 0) {  // so unbalanced don't attempt scaling, just evenly distribute weight
            goodWeight = overallScaledWeight;
            badWeight = overallScaledWeight;
        } else {
            badWeight = overallScaledWeight * numBinVariants / numBadBin / 2.0;
            goodWeight = overallScaledWeight * numBinVariants / numGoodBin / 2.0;
        }
        goodWeights[propertyBinIndex] = (float)goodWeight;
        badWeights[propertyBinIndex] = (float)badWeight;
    }

    private void setPerVariantOptimalMinGq() {
        // Get initial optimal filter qualities, optimizing each variant separately
        // Collect total summary stats, store min GQ
        final List<List<Integer>> indicesList = IntStream.range(0, numPropertyBins)
            .mapToObj(propertyBin -> new ArrayList<Integer>())
            .collect(Collectors.toList());
        IntStream.range(0, numVariants)
            .forEach(index -> indicesList.get(propertyBins[index]).add(index));

        perVariantOptimalMinGq = new MinGq[numVariants];
        propertyBinInheritanceWeights = new double[numPropertyBins];
        propertyBinTruthWeights = new double[numPropertyBins];

        FilterSummary unbinnedFilterSummary = FilterSummary.EMPTY;
        double totalInheritWeight = 0;
        double totalTruthWeight = 0;
        for(int propertyBin = 0; propertyBin < numPropertyBins; ++propertyBin) {
            final List<Integer> indices = indicesList.get(propertyBin);
            final FilterQuality binFilterQuality = setBinOptimalMinGq(indices, propertyBin);
            unbinnedFilterSummary = unbinnedFilterSummary.add(binFilterQuality);
            // set weight for truth and inheritance for this bin based on the loss / performance
            final double inheritanceLoss = binFilterQuality.loss.inheritanceLoss;
            propertyBinInheritanceWeights[propertyBin] = assignWeightFromLoss(inheritanceLoss,
                                                                              binFilterQuality.isLargeAlleleFrequency);
            totalInheritWeight += propertyBinInheritanceWeights[propertyBin];
            final double truthLoss = binFilterQuality.loss.truthLoss;
            propertyBinTruthWeights[propertyBin] = assignWeightFromLoss(truthLoss, false);
            totalTruthWeight += propertyBinTruthWeights[propertyBin];
        }

        if(progressVerbosity > 0) {
            System.out.format("Done setting optimal minGQ.\n");
            System.out.format("totalInheritWeight = %f, totalTruthWeight = %f\n",
                              totalInheritWeight, totalTruthWeight);
        }

        // scale weights so that:
        //    -overall average weight across variants is 1.0
        //    -each bin has weight proportional to evidence strength (inheritance or truth)
        //    -each bin has weight inversely proportional to number of variants in that bin, so that rare variant types
        //     are not swamped in the optimization
        //    -passing variants and failing variants have their weights balanced
        propertyBinMinGqWeights = new float[numPropertyBins];
        propertyBinGoodTruthWeights = new float[numPropertyBins];
        propertyBinBadTruthWeights = new float[numPropertyBins];
        propertyBinGoodInheritanceWeights = new float[numPropertyBins];
        propertyBinBadInheritanceWeights = new float[numPropertyBins];
        final int numTrainableSampleVariants = getNumTrainableSampleVariants(trainingIndices);

        for(int propertyBin = 0; propertyBin < numPropertyBins; ++propertyBin) {
            long onlyMinGqTrainable = 0;
            for (final int variantIndex : trainingIndices) {
                if(propertyBins[variantIndex] != propertyBin) {
                    continue;
                }
                final Set<Integer> trainableSampleIndices = getInheritanceTrainableSampleIndices(variantIndex);
                        trainableSampleIndices.addAll(
                                getGoodSampleIndices(variantIndex)
                        );
                        trainableSampleIndices.addAll(
                                getBadSampleIndices(variantIndex)
                        );
                onlyMinGqTrainable += getNumTrainableSamples(variantIndex) - trainableSampleIndices.size();
            }
            propertyBinMinGqWeights[propertyBin] = (float)minBinWeight * numTrainableSampleVariants / numPropertyBins
                                                    / FastMath.max(onlyMinGqTrainable, 1);

            scalePropertyBinWeights(propertyBin, propertyBinInheritanceWeights, totalInheritWeight,false,
                                    numTrainableSampleVariants, propertyBinGoodInheritanceWeights,
                                    propertyBinBadInheritanceWeights);
            scalePropertyBinWeights(propertyBin, propertyBinTruthWeights, totalTruthWeight,true,
                                    numTrainableSampleVariants, propertyBinGoodTruthWeights, propertyBinBadTruthWeights);
        }

        if(progressVerbosity > 0) {
            System.out.println("got propertyBin weights");
            final int labelWidth = Arrays.stream(propertyBinLabels).mapToInt(String::length).max()
                .orElseThrow(RuntimeException::new);
            final String labelFormat = "%" + labelWidth + "s";
            System.out.format(String.join("\t", propertyBinLabelColumns) +
                                 "\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\n",
                              "minGqWeight", "inh weight", "+inh weight", "-inh weight",
                                        "true weight", "+true weight", "-true weight");
            for(int propertyBin = 0; propertyBin < numPropertyBins; ++propertyBin) {
                System.out.format(labelFormat + "\t%9.4g\t%9.4g\t%9.4g\t%9.4g\t%9.4g\t%9.4g\t%9.4g\n",
                    propertyBinLabels[propertyBin], propertyBinMinGqWeights[propertyBin],
                    propertyBinInheritanceWeights[propertyBin], propertyBinGoodInheritanceWeights[propertyBin],
                    propertyBinBadInheritanceWeights[propertyBin],
                    propertyBinTruthWeights[propertyBin], propertyBinGoodTruthWeights[propertyBin],
                    propertyBinBadTruthWeights[propertyBin]
                );
            }
        }

        // compute the proportion of filterable variants that passed the optimal minGQ filter
        final long numFilterable = IntStream.range(0, numVariants)
                .mapToLong(
                        variantIndex -> IntStream.range(0, numSamples)
                                .filter(sampleIndex -> getSampleVariantIsFilterable(variantIndex, sampleIndex))
                                .count()
                )
                .sum();
        if(progressVerbosity > 0) {
            System.out.format("%d filterable variant x sample genotypes\n", numFilterable);
        }

        optimalProportionOfSampleVariantsPassing = unbinnedFilterSummary.numVariantsPassed
                                                 / (double)unbinnedFilterSummary.numVariants;
        if(progressVerbosity > 0) {
            System.out.format("Optimal proportion of variants passing: %.3f\n", optimalProportionOfSampleVariantsPassing);
        }
    }

    FilterQuality setBinOptimalMinGq(final List<Integer> indices, final int propertyBin) {
        System.out.println(
            "Optimizing minGqs for\t" + String.join("\t", propertyBinLabelColumns) + "\n" +
            "                     \t" + propertyBinLabels[propertyBin]
        );
        System.out.println();
        FilterSummary overallSummary = FilterSummary.EMPTY;
        final FilterSummary[] filterSummaries = new FilterSummary[indices.size()];
        for(int i = 0; i < indices.size(); ++i) {
            final int variantIndex = indices.get(i);
            final FilterQuality greedyFilter = getOptimalVariantMinGq(variantIndex, null);
            filterSummaries[i] = greedyFilter;
            overallSummary = overallSummary.add(greedyFilter);
            perVariantOptimalMinGq[variantIndex] = greedyFilter.minGq;
        }

        // Iteratively improve filters, optimizing for OVERALL loss
        boolean anyImproved = true;
        int numIterations = 0;
        FilterLoss previousLoss = new FilterLoss(overallSummary);
        System.out.println("           \t" + FilterLoss.getHeader(propertyBinLabelColumns));
        System.out.println("Iteration " + numIterations + "\t" + previousLoss);
        while(anyImproved) {
            ++numIterations;
            anyImproved = false;
            for(int i = 0; i < indices.size(); ++i) {
                final int variantIndex = indices.get(i);
                final FilterSummary previousFilter = filterSummaries[i];
                final FilterSummary backgroundFilter = overallSummary.subtract(previousFilter);
                final FilterQuality greedyFilter = getOptimalVariantMinGq(variantIndex, backgroundFilter);

                if(greedyFilter.loss.lt(previousLoss)) {
                    anyImproved = true;
                    overallSummary = greedyFilter;
                    perVariantOptimalMinGq[variantIndex] = greedyFilter.minGq;
                    filterSummaries[i] = greedyFilter.subtract(backgroundFilter);
                } else if(greedyFilter.loss.gt(previousLoss)) {
                    throw new GATKException(
                            "Loss increased. This is a bug!\n" +
                            "previous:" + previousLoss + "\n" +
                            "current: " + greedyFilter.loss
                    );
                }
                previousLoss = greedyFilter.loss;
            }
            System.out.println("Iteration " + numIterations + "\t" + previousLoss);
        }
        displayHistogram(keepHomvar ? "Optimal minGq histogram" : "Optimal minGqHet histogram",
                indices.stream().mapToInt(index -> perVariantOptimalMinGq[index].minGqHet),
                true);
        if(!keepHomvar) {
            displayHistogram("Optimal minGqHomVar histogram",
                    indices.stream().mapToInt(index -> perVariantOptimalMinGq[index].minGqHomVar),
                    true);
        }
        return new FilterQuality(overallSummary);
    }

    protected int getNumTrainableSamples(final int variantIndex) {
        return (int)IntStream.range(0, getNumSamples())
                .filter(sampleIndex -> getSampleVariantIsTrainable(variantIndex, sampleIndex))
                .count();
    }
    /**
     * Get number of rows, account for the fact that unfilterable (e.g. already HOMREF) samples will not be used
     */
    private final Map<int[], Integer> numTrainableSampleVariantsMap = new HashMap<>();
    protected int getNumTrainableSampleVariants(final int[] variantIndices) {
        if(!numTrainableSampleVariantsMap.containsKey(variantIndices)) {
            final long numTrainableSampleVariants = Arrays.stream(variantIndices)
                    .mapToLong(this::getNumTrainableSamples)
                    .sum();
            if(numTrainableSampleVariants > Integer.MAX_VALUE) {
                throw new GATKException("This data has " + numTrainableSampleVariants + " trainable samples x variants,"
                                        + " which is more than Integer.MAX_VALUE. It is not possible to make large"
                                        + " enough arrays for training");
            }
            numTrainableSampleVariantsMap.put(variantIndices, (int)numTrainableSampleVariants);
            return (int)numTrainableSampleVariants;
        }
        return numTrainableSampleVariantsMap.get(variantIndices);
    }

    protected boolean samplePasses(final MinGq minGq, final short gq, final byte alleleCount) {
        return gq >= (alleleCount > 1 ? minGq.minGqHomVar : minGq.minGqHet);
    }

    protected boolean samplePasses(final int variantIndex, final int sampleIndex) {
        return samplePasses(perVariantOptimalMinGq[variantIndex],
                            sampleVariantCallQualities.values[variantIndex][sampleIndex],
                            sampleVariantAlleleCounts.values[variantIndex][sampleIndex]);
    }

    protected int fillSamplePasses(final int variantIndex, int flatIndex, final boolean[] samplePasses) {
        final short[] sampleGqs = sampleVariantCallQualities.values[variantIndex];
        final byte[] sampleAlleleCounts = sampleVariantAlleleCounts.values[variantIndex];
        final MinGq minGq = perVariantOptimalMinGq[variantIndex];
        for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
            if(getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                samplePasses[flatIndex] = samplePasses(minGq, sampleGqs[sampleIndex], sampleAlleleCounts[sampleIndex]);
                ++flatIndex;
            }
        }
        return flatIndex;
    }

    private final Map<int[], boolean[]> sampleVariantTruthMap = new HashMap<>();
    protected boolean[] getSampleVariantTruth(final int[] variantIndices) {
        if(sampleVariantTruthMap.containsKey(variantIndices)) {
            return sampleVariantTruthMap.get(variantIndices);
        }
        final int numRows = getNumTrainableSampleVariants(variantIndices);
        final boolean[] sampleVariantTruth = new boolean[numRows];

        int flatIndex = 0;
        for(final int variantIndex : variantIndices) {
            flatIndex = fillSamplePasses(variantIndex, flatIndex, sampleVariantTruth);
        }
        sampleVariantTruthMap.put(variantIndices, sampleVariantTruth);
        return sampleVariantTruth;
    }

    @SuppressWarnings("SameParameterValue")
    protected void displayHistogram(final String description, final IntStream intStream, boolean binValues) {
        final Map<Integer, Integer> rawValuesMap = new HashMap<>();
        intStream.forEach(gq -> {
            if (rawValuesMap.containsKey(gq)) {
                rawValuesMap.put(gq, 1 + rawValuesMap.get(gq));
            } else {
                rawValuesMap.put(gq, 1);
            }
        });
        if(rawValuesMap.size() == 0) {
            System.out.println(description + ": no data");
        }
        final int minGqValue = rawValuesMap.keySet().stream().min(Integer::compareTo).orElseThrow(RuntimeException::new);
        final int maxGqValue = rawValuesMap.keySet().stream().max(Integer::compareTo).orElseThrow(RuntimeException::new);

        final Map<Integer, Integer> displayValuesMap;
        if(binValues) {
            displayValuesMap = new HashMap<>();
            rawValuesMap.forEach((gq, numGq) -> {
                final int binGq;
                if (gq == 0) {
                    binGq = gq;
                } else {
                    final int magnitude = (int) Math.pow(10.0, Math.floor(Math.log10(Math.abs(gq))));
                    binGq = magnitude * (gq / magnitude);
                }

                if (displayValuesMap.containsKey(binGq)) {
                    displayValuesMap.put(binGq, numGq + displayValuesMap.get(binGq));
                } else {
                    displayValuesMap.put(binGq, numGq);
                }
            });
        } else {
            displayValuesMap = rawValuesMap;
        }
        System.out.println(description + ":");
        System.out.println("min=" + minGqValue + ", max=" + maxGqValue);
        displayValuesMap.keySet()
            .stream()
            .sorted()
            .forEach(minGq -> System.out.println(minGq + ": " + displayValuesMap.get(minGq)));
    }

    @SuppressWarnings("SameParameterValue")
    protected void displayHistogram(final String description, final DoubleStream doubleStream, final int numBins,
                                    final double minValue, final double maxValue) {
        final double safety = 100;
        final double valueEps = safety * FastMath.max(FastMath.ulp(FastMath.abs(minValue)),
                                                      FastMath.ulp(FastMath.abs(maxValue)));
        final double numBinsEps = safety * FastMath.ulp((double)numBins);
        final double valueOffset = minValue - valueEps;
        final double binScale = (numBins - numBinsEps) / (maxValue - valueOffset);
        final long[] valueBins = new long [numBins];
        double[] valueRange = new double[] {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
        doubleStream.forEach(value -> {
            final int bin = (int)FastMath.floor((value - valueOffset) * binScale);
            try {
                ++valueBins[bin];
            } catch(ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
                final String errStr = String.format("bad bin: min: %f, max: %f, scale: %f, offset: %f, value: %f\n",
                                                    minValue, maxValue, binScale, valueOffset, value);
                throw new GATKException(errStr, arrayIndexOutOfBoundsException);
            }
            if(value < valueRange[0]) {
                valueRange[0] = value;
            }
            if(value > valueRange[1]) {
                valueRange[1] = value;
            }
        });
        final long numValues = Arrays.stream(valueBins).sum();
        if(valueRange[0] > valueRange[1]) {
            System.out.println(description + ": no data");
        } else {
            System.out.format("%s: %d values\n", description, numValues);
            System.out.format("\tactual range: [%f, %f]\n", valueRange[0], valueRange[1]);
            System.out.println("\t low - high    %");
            double high = minValue;
            for(int bin = 0; bin < numBins; ++bin) {
                double low = high;
                high = low + 1.0 / binScale;
                System.out.format("\t%.2f - %.2f   %.1f\n", low, high, valueBins[bin] * 100.0 / numValues);
            }
        }
    }

    protected void displayPercentiles(final String description, final DoubleStream doubleStream) {
        final int[] numNans = new int[1];
        final double[] values = doubleStream
            .filter(v -> {if(Double.isNaN(v)) { ++numNans[0]; return false;} else return true;})
            .sorted()
            .toArray();
        if(values.length == 0) {
            if(numNans[0] == 0) {
                System.out.println(description + ": no data");
            } else {
                System.out.format("%s: %d values 100%% are NaN\n", description, numNans[0]);
            }
        } else if(values[0] == values[values.length - 1]) {
            System.out.format("%s: %d NaNs, %d non-NaN values 100%% equal %.3g\n",
                              description, numNans[0], values.length, values[0]);
        } else {
            System.out.format("%s: %d NaNs, %d non-NaN values\n", description, numNans[0], values.length);
            System.out.format("\tmin: %.3g\n", values[0]);
            for(final int percentile : new int[]{25, 50, 75}) {
                final int percentileIndex = (int)FastMath.round((percentile / 100.0) * (values.length - 1));
                System.out.format("\t%d%%: %.3g\n", percentile, values[percentileIndex]);
            }
            System.out.format("\tmax: %.3g\n", values[values.length - 1]);
        }
    }

    void printPropertiesValuesSummary() {
        final int nameWidth = FastMath.max(
            "propertyName".length(),
            StreamSupport.stream(propertiesTable.spliterator(), false)
                .mapToInt(prop -> prop.name.length())
                .max().orElse(0)
        );
        final String nameFormat = "%"+ nameWidth + "s";
        System.out.format("index\t" + nameFormat + "\tnRows\tnColumns\t%10s\t%10s\n",
                          "propertyName", "baseline", "scale");
        int idx = 0;
        for(final PropertiesTable.Property property : propertiesTable) {
            System.out.format("%d\t" + nameFormat + "\t%d\t%8d\t%10.5g\t%10.5g\n",
                              idx, property.name, property.getNumRows(), property.getNumColumns(),
                                        property.getBaseline(), property.getScale());
            ++idx;
        }
    }

    void printPropertiesDebugInfo() {
        System.out.println("########################################");
        System.out.println("numVariants: " + getNumVariants());
        System.out.println("numSamples: " + getNumSamples());
        System.out.println("numTrios: " + getNumTrios());
        System.out.println("numProperties: " + getNumProperties());
        System.out.format("Input VCF had %.1f variants per sample.\n", numInputVar / (double)numSamples);
        printPropertiesValuesSummary();

        final Map<String, List<String>> labelsEncoding = propertiesTable.getLabelsEncoding();
        for(final Map.Entry<String, List<String>> entry : labelsEncoding.entrySet()) {
            System.out.println(entry.getKey() + ":");
            int idx = 0;
            for(final String label : entry.getValue()) {
                System.out.println(idx + "\t" + label);
                ++idx;
            }
        }

        displayHistogram("Filterable alleles Gq histogram:", streamFilterableGq(),true);

        System.out.println("########################################");
    }


    protected FilterLoss getTrainingLoss(final float[] pSampleVariantIsGood, final float[] d1Loss, final float[] d2Loss,
                                         final int[] variantIndices) {
        final boolean[] sampleVariantIsGood = getSampleVariantTruth(variantIndices);
        // zero out derivative arrays
        for(int predictIndex = 0; predictIndex < pSampleVariantIsGood.length; ++predictIndex) {
            d1Loss[predictIndex] = 0F;
            d2Loss[predictIndex] = 0F;
        }

        int predictIndex = 0;
        FilterLoss loss = FilterLoss.EMPTY;
        for(final int variantIndex : variantIndices) {
            final Set<Integer> inheritanceTrainableSampleIndices = getInheritanceTrainableSampleIndices(variantIndex);
            final Set<Integer> truthTrainableSampleIndices = new HashSet<>(
                    getGoodSampleIndices(variantIndex)
            );
            truthTrainableSampleIndices.addAll(
                    getBadSampleIndices(variantIndex)
            );
            final int propertyBin = propertyBins[variantIndex];
            final float minGqWeight = propertyBinMinGqWeights[propertyBin];
            final float halfMinGqWeight = minGqWeight / 2F;
            final float goodTruthWeight = halfMinGqWeight + propertyBinGoodTruthWeights[propertyBin];
            final float badTruthWeight = halfMinGqWeight + propertyBinBadTruthWeights[propertyBin];
            final float goodInheritanceWeight = halfMinGqWeight + propertyBinGoodInheritanceWeights[propertyBin];
            final float badInheritanceWeight = halfMinGqWeight + propertyBinBadInheritanceWeights[propertyBin];
            final float truthWeight = halfMinGqWeight + (float)propertyBinTruthWeights[propertyBin];
            final float inheritanceWeight = halfMinGqWeight + (float)propertyBinInheritanceWeights[propertyBin];

            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                final double p = pSampleVariantIsGood[predictIndex];
                final boolean sampleIsGood = sampleVariantIsGood[predictIndex];
                final double delta, deltaLoss, d1, d2;
                if(sampleIsGood) {
                    delta = FastMath.max(PROB_EPS, p);
                    deltaLoss = -FastMath.log(delta);
                    d1 = -1.0 / delta;
                    d2 = -d1 / delta;
                } else {
                    delta = FastMath.max(PROB_EPS, 1.0 - p);
                    deltaLoss = -FastMath.log(delta);
                    d1 = 1.0 / delta;
                    d2 = d1 / delta;
                }

                final double sampleInheritanceWeight, sampleTruthWeight;
                if(inheritanceTrainableSampleIndices.contains(sampleIndex)) {
                    sampleInheritanceWeight = sampleIsGood ? goodInheritanceWeight : badInheritanceWeight;
                    //sampleInheritanceWeight = inheritanceWeight;
                } else {
                    sampleInheritanceWeight = halfMinGqWeight;
                }
                if(truthTrainableSampleIndices.contains(sampleIndex)) {
                    sampleTruthWeight = sampleIsGood ? goodTruthWeight : badTruthWeight;
                    //sampleTruthWeight = truthWeight;
                } else {
                    sampleTruthWeight = halfMinGqWeight;
                }
                final double derivWeight = sampleInheritanceWeight + sampleTruthWeight;

                setLossDerivs(p, d1, d2, derivWeight, d1Loss, d2Loss, predictIndex);
                final FilterLoss sampleVariantLoss =
                    new FilterLoss(deltaLoss, deltaLoss, sampleInheritanceWeight, sampleTruthWeight, null);
                loss = FilterLoss.add(loss, sampleVariantLoss);
                ++predictIndex;
            }
        }
        return loss;
    }

    protected FilterLoss getTrainingLoss(final float[] pSampleVariantIsGood, final int[] variantIndices) {
        final boolean[] sampleVariantIsGood = getSampleVariantTruth(variantIndices);

        int predictIndex = 0;
        FilterLoss loss = FilterLoss.EMPTY;
        for(final int variantIndex : variantIndices) {
            final Set<Integer> inheritanceTrainableSampleIndices = getInheritanceTrainableSampleIndices(variantIndex);
            final Set<Integer> truthTrainableSampleIndices = new HashSet<>(
                    getGoodSampleIndices(variantIndex)
            );
            truthTrainableSampleIndices.addAll(
                    getBadSampleIndices(variantIndex)
            );
            final int propertyBin = propertyBins[variantIndex];
            final float minGqWeight = propertyBinMinGqWeights[propertyBin];
            final float halfMinGqWeight = minGqWeight / 2F;
            final float goodTruthWeight = halfMinGqWeight + propertyBinGoodTruthWeights[propertyBin];
            final float badTruthWeight = halfMinGqWeight + propertyBinBadTruthWeights[propertyBin];
            final float goodInheritanceWeight = halfMinGqWeight + propertyBinGoodInheritanceWeights[propertyBin];
            final float badInheritanceWeight = halfMinGqWeight + propertyBinBadInheritanceWeights[propertyBin];
            final float truthWeight = halfMinGqWeight + (float)propertyBinTruthWeights[propertyBin];
            final float inheritanceWeight = halfMinGqWeight + (float)propertyBinInheritanceWeights[propertyBin];

            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                final double p = pSampleVariantIsGood[predictIndex];
                final boolean sampleIsGood = sampleVariantIsGood[predictIndex];
                final double deltaLoss = -FastMath.log(FastMath.max(PROB_EPS, sampleIsGood ? p : 1.0 - p));

                final double sampleInheritanceWeight, sampleTruthWeight;
                if(inheritanceTrainableSampleIndices.contains(sampleIndex)) {
                    sampleInheritanceWeight = sampleIsGood ? goodInheritanceWeight : badInheritanceWeight;
                    //sampleInheritanceWeight = inheritanceWeight;
                } else {
                    sampleInheritanceWeight = halfMinGqWeight;
                }
                if(truthTrainableSampleIndices.contains(sampleIndex)) {
                    sampleTruthWeight = sampleIsGood ? goodTruthWeight : badTruthWeight;
                    //sampleTruthWeight = truthWeight;
                } else {
                    sampleTruthWeight = halfMinGqWeight;
                }

                final FilterLoss sampleVariantLoss =
                        new FilterLoss(deltaLoss, deltaLoss, sampleInheritanceWeight, sampleTruthWeight, null);
                loss = FilterLoss.add(loss, sampleVariantLoss);
                ++predictIndex;
            }
        }
        return loss;
    }

    private boolean[][] getSampleVariantPasses(final float[] pSampleVariantIsGood, final int[] variantIndices) {
        final boolean[][] sampleVariantPasses = new boolean[variantIndices.length][];
        int flatIndex = 0;
        for(int row = 0; row < variantIndices.length; ++row) {
            final int variantIndex = variantIndices[row];
            final boolean[] samplePasses = new boolean[getNumTrainableSamples(variantIndex)];
            sampleVariantPasses[row] = samplePasses;
            for(int col = 0; col < samplePasses.length; ++col) {
                samplePasses[col] = pSampleVariantIsGood[flatIndex] >= 0.5F;
                ++flatIndex;
            }
        }
        return sampleVariantPasses;
    }

    protected FilterLoss getLoss(final float[] pSampleVariantIsGood, final int[] variantIndices) {
        final boolean[][] sampleVariantPasses = getSampleVariantPasses(pSampleVariantIsGood, variantIndices);

        return new FilterLoss(
            IntStream.range(0, variantIndices.length)
                .mapToObj(i -> getBinnedFilterSummary(sampleVariantPasses[i], variantIndices[i]))
                .reduce(BinnedFilterSummaries::add).orElse(BinnedFilterSummaries.EMPTY),
            propertyBinLabelColumns
        );
    }

    protected FilterLoss getLoss(final MinGq[] minGq, final int[] variantIndices) {
        if(minGq.length != variantIndices.length) {
            throw new GATKException(
                    "Length of minGq (" + minGq.length + ") does not match length of variantIndices (" + variantIndices.length + ")"
            );
        }
        return new FilterLoss(
            IntStream.range(0, minGq.length)
                .mapToObj(i -> getBinnedFilterSummary(minGq[i], variantIndices[i]))
                .reduce(BinnedFilterSummaries::add).orElse(BinnedFilterSummaries.EMPTY),
            propertyBinLabelColumns
        );
    }

    private void displayNaiveLoss(final int[] variantIndices, final String name) {
        final MinGq[] minGq = Arrays.stream(variantIndices).mapToObj(i -> new MinGq((short)0, (short)0)).toArray(MinGq[]::new);
        final FilterLoss loss = getLoss(minGq, variantIndices);
        System.out.format("\tNaive loss for %s\n\t\t%s\n", name, loss.toString().replaceAll("\n", "\n\t\t"));
    }

    private void displayNaiveLoss() {
        if(progressVerbosity > 0) {
            displayNaiveLoss(trainingIndices, "training");
            displayNaiveLoss(validationIndices, "validation");
        }
    }

    private void setTrainingAndValidationIndices() {
        final int numValidationIndices = (int)round(validationProportion * numVariants);
        final List<Integer> shuffleIndices = IntStream.range(0, numVariants).boxed().collect(Collectors.toList());
        Collections.shuffle(shuffleIndices, randomGenerator);

        validationIndices = shuffleIndices.subList(0, numValidationIndices).stream()
            .sorted().mapToInt(Integer::intValue).toArray();
        trainingIndices = shuffleIndices.subList(numValidationIndices, numVariants).stream()
            .sorted().mapToInt(Integer::intValue).toArray();
    }

    private void saveTrainedModel() {
        try (final OutputStream outputStream = modelFile.getOutputStream()) {
            final OutputStream unclosableOutputStream = new FilterOutputStream(outputStream) {
                @Override
                public void close() {
                    // don't close the stream in one of the subroutines
                }
            };
            propertiesTable.saveDataEncoding(unclosableOutputStream);
            saveModel(unclosableOutputStream);
        } catch(IOException ioException) {
            throw new GATKException("Error saving modelFile " + modelFile, ioException);
        }
    }

    private void savePropertiesTable() {
        if(propertiesTableFile == null) {
            return;
        }
        try (final OutputStream outputStream = propertiesTableFile.getOutputStream()) {
            final OutputStream unclosableOutputStream = new FilterOutputStream(outputStream) {
                @Override
                public void close() {
                    // don't close the stream in one of the subroutines
                }
            };
            propertiesTable.save(unclosableOutputStream);
        } catch(IOException ioException) {
            throw new GATKException("Error saving propertiesTable " + propertiesTableFile, ioException);
        }
    }

    private void loadTrainedModel() {
        if(modelFile == null || !Files.exists(modelFile.toPath())) {
            if(runMode == RunMode.Filter) {
                throw new UserException("mode=FILTER, but trained model file was not provided.");
            }
            return;
        }
        try (final InputStream inputStream = modelFile.getInputStream()) {
            final InputStream unclosableInputStream = new FilterInputStream(inputStream) {
                @Override
                public void close() {
                    // don't close the stream in one of the subroutines
                }
            };
            propertiesTable.loadDataEncoding(unclosableInputStream);
            loadModel(unclosableInputStream );
        } catch (Exception exception) {
            throw new GATKException("Error loading modelFile " + modelFile + " (malformed file?)", exception);
        }
        System.out.println("loadTrainedModel complete");
    }

    private byte[] modelCheckpoint = null;

    protected void saveModelCheckpoint() {
        final ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
        saveModel(outputStream);
        modelCheckpoint = outputStream.toByteArray();
    }

    protected void loadModelCheckpoint() {
        final ByteArrayInputStream inputStream = new ByteArrayInputStream(modelCheckpoint);
        loadModel(inputStream);
    }

    protected void clearUnneededProperties() {
        // clear properties that are only used for training (e.g. keep allele counts, which are used for scoring too)
        // as those properties have all been copied into DMatrices, and aren't needed anymore
        final Set<String> neededProperties = new HashSet<>(Arrays.asList(
            VCFConstants.ALLELE_COUNT_KEY, NO_CALL_COUNTS_KEY
        ));

        if(progressVerbosity > 0) {
            System.out.println("clearUnneededProperties()");
        }
        for (final String propertyName : propertiesTable.getPropertyNames()) {
            if(!neededProperties.contains(propertyName)) {
                if(progressVerbosity > 0) {
                    System.out.format("\t%s", propertyName);
                }
                propertiesTable.get(propertyName).clearAllocatedRows();
            }
        }
        // We're doing this to try to prevent out-of-memory crashes, so garbage collect and wait for a few seconds to
        // hopefully get the maximum benefit
        System.gc();
        try {
            Thread.sleep(5);
        } catch (InterruptedException e) {
            // Not worth crashing over, just warn and continue
            System.out.println("Unable to sleep after garbage collection");
        }
        if(progressVerbosity > 0) {
            System.out.println("\tdone.");
        }
    }

    protected abstract boolean needsNormalizedProperties();
    protected abstract int predictBatch(final float[] outputProbabilities);
    protected abstract void trainFilter();
    protected abstract void saveModel(final OutputStream outputStream);
    protected abstract void loadModel(final InputStream inputStream);

    protected short probToPhred(final double p) {
        return p >= 1 ?
                Short.MAX_VALUE :
                (short) FastMath.floor(1. - 10. * FastMath.log10(2. * (1. - p)));
    }

    static final double PHRED_COEF = FastMath.log(10.0) / 10.0;

    protected double phredToProb(final double phred) {
        return 1.0 - 0.5 * FastMath.exp(-PHRED_COEF * (1.0 - phred));
    }

    private float[] outputProbabilitiesForAdjustedGq = null;
    private short[][] outputQualitiesForAdjustedGq = null;
    private short[] outputGqForAdjustedGq = null;
    private short[] outputLogitsForAdjustedGq = null;
    private short[][] adjustedQualityBatch() {
        if(outputProbabilitiesForAdjustedGq == null) {
            outputProbabilitiesForAdjustedGq = new float[numSamples];
            outputGqForAdjustedGq = new short[numSamples];
            outputLogitsForAdjustedGq = new short[numSamples];
            outputQualitiesForAdjustedGq = new short[][] {outputGqForAdjustedGq, outputLogitsForAdjustedGq};
        }
        final int numRows = predictBatch(outputProbabilitiesForAdjustedGq);
        for(int row = 0; row < numRows; ++row) {
            outputGqForAdjustedGq[row] = probToPhred(outputProbabilitiesForAdjustedGq[row]);
            outputLogitsForAdjustedGq[row] = probToScaledLogits(outputProbabilitiesForAdjustedGq[row]);
        }
        return outputQualitiesForAdjustedGq;
    }

    // near p=0.5, each scaled logit is ~ 0.1% change in likelihood
    static final double LOGIT_SCALE = 1.0 / FastMath.log(0.501 / 0.499);

    protected short probToScaledLogits(final double p) {
        return p == 0 ? Short.MIN_VALUE :
                        p == 1 ? Short.MAX_VALUE :
                                 (short)FastMath.floor(LOGIT_SCALE * FastMath.log(p / (1.0 - p)));
    }

    protected float scaledLogitsToP(final float scaled_logits) {
        return (float)(1.0 / (1.0 + FastMath.exp(-scaled_logits / LOGIT_SCALE)));
    }

    @Override
    public Object onTraversalSuccess() {
        if(runMode == RunMode.Train) {
            if(numVariants == 0) {
                throw new GATKException("No variants contained in vcf: " + drivingVariantFile);
            }

            setPropertyBins();
            setPropertyBinIsLargeAlleleFraction();
            propertiesTable.validateAndFinalize();
            sampleVariantCallQualities = useCopyNumberCalls ?
                    (PropertiesTable.ShortMatProperty) propertiesTable.get(CALL_QUALITY_KEY) :
                    (PropertiesTable.ShortMatProperty) propertiesTable.get(VCFConstants.GENOTYPE_QUALITY_KEY);
            sampleVariantAlleleCounts = (PropertiesTable.ByteMatProperty) propertiesTable.get(VCFConstants.ALLELE_COUNT_KEY);
            sampleVariantNoCallCounts = (PropertiesTable.ByteMatProperty) propertiesTable.get(NO_CALL_COUNTS_KEY);

            printPropertiesDebugInfo();

            savePropertiesTable();
            setTrainingAndValidationIndices();
            setPerVariantOptimalMinGq();
            displayNaiveLoss();

            trainFilter();
            saveTrainedModel();
        } else {
            System.out.println("Filter summary:");
            System.out.println("\tFiltered " + numFilteredGenotypes + " of " + numFilterableGenotypes +
                               " filterable genotypes in " + numVariants + " variants x " + numSamples + " samples");
            System.out.format("\t\t = %.1f%% of genotypes filtered.\n",
                              100.0 * numFilteredGenotypes / (double)(numFilterableGenotypes));
            System.out.format("\tOmitted %d trivial variants (%.1f%% of total)",
                              numDroppedVariants, numDroppedVariants * 100.0 / numVariants);
            System.out.format("\tInput VCF had %.1f variants per sample\n",
                              numInputVar / (double)numSamples);
            System.out.format("\t\t(%d ref, %d non-ref, %d no-call)\n", numInputRef, numInputVar, numInputNoCall);
            System.out.format("\tOutput VCF had %.1f variants per sample\n",
                             numOutputVar / (double)numSamples);
            System.out.format("\t\t(%d ref, %d non-ref, %d no-call)\n", numOutputRef, numOutputVar, numOutputNoCall);
        }
        return null;
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
            System.out.println("Closed " + outputFile);
        }
    }
}
