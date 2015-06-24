package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.*;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.tools.walkers.bqsr.ReadRecalibrationInfo;
import org.broadinstitute.hellbender.tools.walkers.bqsr.RecalibrationEngine;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.*;

/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various covariates
 * (such as read group, reported quality score, machine cycle, and nucleotide context).
 *
 * Based off of BaseRecalibrator, but changed to be able to work in Dataflow.
 *
 * Expected usage:
 * -ship the command line to Dataflow workers
 * -create a BaseCalibratorUprooted there,
 * -feed it the right inputs (use the CalibrationTablesBuilder wrapper
 *  to help)
 * -extract the outputs, merge them, and save.
 * -use a BaseRecalibratorUprooted to format the outputs into a report.
 */

@CommandLineProgramProperties(
        usage = "First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        usageShort = "Generates recalibration table (do not run this one from the cmd line)",
        programGroup = ReadProgramGroup.class
)
public final class BaseRecalibratorWorker {
    final protected static Logger logger = LogManager.getLogger(BaseRecalibratorWorker.class);

    @ArgumentCollection(doc="all the command line arguments for BQSR and its covariates")
    private BaseRecalibrationArgumentCollection BRAC;


    // --------------------------------------------------------------------------------------------------------------
    // Non-command line members

    final protected SAMFileHeader samHeader;

    /**
     * an object that keeps track of the information necessary for quality score quantization
     */
    private QuantizationInfo quantizationInfo;

    /**
     * list to hold the all the covariate objects that were requested (required + standard + experimental)
     */
    private StandardCovariateList covariates;

    private RecalibrationEngine recalibrationEngine;

    private int minimumQToUse;

    protected static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
    private ReferenceDataSource referenceDataSource; // datasource for the reference. We're using a different one from the engine itself to avoid messing with its caches.
    private final static byte NO_BAQ_UNCERTAINTY = (byte) '@';


    private long accumulator;

    /**
     * This method returns a new BaseCalibratorUprooted that is configured according to the command-line
     * arguments in "toolArgs" and the input file header in "header". The command-line parser will print
     * any diagnostic messages to "out". If the command-line requested the program not be run
     * (e.g. it was "--help"), then this method returns null.
     */
    public static BaseRecalibratorWorker fromCommandLine(SAMFileHeader header, String commandLine, PrintStream out) {
        BaseRecalibratorWorker br = new BaseRecalibratorWorker(header);
        br.BRAC = new BaseRecalibrationArgumentCollection();
        boolean shouldRun = new CommandLineParser(br).parseArguments(out, Utils.escapeExpressions(commandLine));
        if (!shouldRun) return null;
        return br;
    }

    public static BaseRecalibratorWorker fromArgs(SAMFileHeader header, BaseRecalibrationArgumentCollection toolArgs) {
        return new BaseRecalibratorWorker(header, toolArgs);
    }

    private BaseRecalibratorWorker(SAMFileHeader samHeader) {
        this.samHeader = samHeader;
    }

    private BaseRecalibratorWorker(SAMFileHeader samHeader, BaseRecalibrationArgumentCollection toolArgs) {
        this.samHeader = samHeader;
        this.BRAC = toolArgs;
    }

    /**
     * Throws if knownSites is missing.
     */
    public void checkClientArguments() {
        if (null == BRAC.RAC.knownSites || BRAC.RAC.knownSites.size() == 0) {
            throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);
        }
    }

    public void onTraversalStart(File refFile) {
        accumulator = 0L;

        baq = new BAQ(BRAC.BAQGOP); // setup the BAQ object with the provided gap open penalty

        if (BRAC.RAC.FORCE_PLATFORM != null)
            BRAC.RAC.DEFAULT_PLATFORM = BRAC.RAC.FORCE_PLATFORM;

        covariates = new StandardCovariateList(BRAC.RAC, getHeaderForReads());

        initializeRecalibrationEngine();
        minimumQToUse = BRAC.PRESERVE_QSCORES_LESS_THAN;

        if (null!=refFile) {
            referenceDataSource = new ReferenceDataSource(refFile);
        }
    }

    public SAMFileHeader getHeaderForReads() {
        return samHeader;
    }

    /**
     * Initialize the recalibration engine
     */
    private void initializeRecalibrationEngine() {
        int numReadGroups = getHeaderForReads().getReadGroups().size();

        recalibrationEngine = new RecalibrationEngine(covariates, numReadGroups);
    }

    private boolean isLowQualityBase(final SAMRecord read, final int offset) {
        return read.getBaseQualities()[offset] < minimumQToUse;
    }

    public static ReadFilter readFilter() {
        return ReadFilterLibrary.WELLFORMED
                .and(MAPPING_QUALITY_NOT_ZERO)
                .and(MAPPING_QUALITY_AVAILABLE)
                .and(MAPPED)
                .and(PRIMARY_ALIGNMENT)
                .and(NOT_DUPLICATE)
                .and(PASSES_VENDOR_QUALITY_CHECK);
    }

    private static SAMRecord consolidateCigar(SAMRecord read) {
        // Always consolidate the cigar string into canonical form, collapsing zero-length / repeated cigar elements.
        // Downstream code cannot necessarily handle non-consolidated cigar strings.
        read.setCigar(AlignmentUtils.consolidateCigar(read.getCigar()));
        return read;
    }

    private SAMRecord resetOriginalBaseQualities(SAMRecord read) {
        if (!BRAC.useOriginalBaseQualities) {
            return read;
        }
        return ReadUtils.resetOriginalBaseQualities(read);
    }

    private SAMRecord setDefaultBaseQualities(SAMRecord read) {
        // if we are using default quals, check if we need them, and add if necessary.
        // 1. we need if reads are lacking or have incomplete quality scores
        // 2. we add if defaultBaseQualities has a positive value
        if (BRAC.defaultBaseQualities < 0) {
            return read;
        }
        byte reads[] = read.getReadBases();
        byte quals[] = read.getBaseQualities();
        if (quals == null || quals.length < reads.length) {
            byte new_quals[] = new byte[reads.length];
            Arrays.fill(new_quals, BRAC.defaultBaseQualities);
            read.setBaseQualities(new_quals);
        }
        return read;
    }

    public RecalibrationTables getTables() {
        return recalibrationEngine.getRecalibrationTables();
    }

    /**
     * Returns the requested covariates, possibly mutated during the execution.
     */
    public StandardCovariateList getRequestedCovariates() {
        return this.covariates;
    }

    /**
     * @param rt recalibration tables
     * @return the quantization information
     */
    public QuantizationInfo getQuantizationInfo(RecalibrationTables rt) {
        return new QuantizationInfo(rt, BRAC.RAC.QUANTIZING_LEVELS);
    }

    private ReadTransformer makeReadTransform() {
        ReadTransformer f0 = BaseRecalibratorWorker::consolidateCigar;

        ReadTransformer f =
                f0.andThen(this::setDefaultBaseQualities)
                        .andThen(this::resetOriginalBaseQualities)
                        .andThen(ReadClipper::hardClipAdaptorSequence)
                        .andThen(ReadClipper::hardClipSoftClippedBases);

        return f;
    }


    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location.
     * The "knownLocations" list doesn't need to be complete but it must includes all those that overlap with the read.
     */
    public void apply(SAMRecord originalRead, final ReferenceContext ref, List<? extends Locatable> knownLocations) {
        ReadTransformer transform = makeReadTransform();
        final SAMRecord read = transform.apply(originalRead);

        if (ReadUtils.isEmpty(read)) {
            return;
        } // the whole read was inside the adaptor so skip it

        RecalUtils.parsePlatformForRead(read, BRAC.RAC);
        if (!RecalUtils.isColorSpaceConsistent(BRAC.RAC.SOLID_NOCALL_STRATEGY, read)) { // parse the solid color space and check for color no-calls
            return; // skip this read completely
        }

        // We've checked in onTraversalStart() that we have a reference, so ref.get() is safe
        final int[] isSNP = calculateIsSNP(read, ref, originalRead);
        final int[] isInsertion = calculateIsIndel(read, EventType.BASE_INSERTION);
        final int[] isDeletion = calculateIsIndel(read, EventType.BASE_DELETION);
        final int nErrors = nEvents(isSNP, isInsertion, isDeletion);

        // note for efficiency regions we don't compute the BAQ array unless we actually have
        // some error to marginalize over.  For ILMN data ~85% of reads have no error
        final byte[] baqArray = nErrors == 0 ? flatBAQArray(read) : calculateBAQArray(read);

        if( baqArray != null ) { // some reads just can't be BAQ'ed
            final ReadCovariates covariates = RecalUtils.computeCovariates(read, this.covariates);
            final boolean[] skip = calculateSkipArray(read, knownLocations); // skip known sites of variation as well as low quality and non-regular bases
            final double[] snpErrors = calculateFractionalErrorArray(isSNP, baqArray);
            final double[] insertionErrors = calculateFractionalErrorArray(isInsertion, baqArray);
            final double[] deletionErrors = calculateFractionalErrorArray(isDeletion, baqArray);

            // aggregate all of the info into our info object, and update the data
            final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, covariates, skip, snpErrors, insertionErrors, deletionErrors);
            recalibrationEngine.updateDataForRead(info);
        }

        accumulator++;
    }

    /**
     * Compute the number of mutational events across all hasEvent vectors
     * <p>
     * Simply the sum of entries in hasEvents
     *
     * @param hasEvents a vector a vectors of 0 (no event) and 1 (has event)
     * @return the total number of events across all hasEvent arrays
     */
    protected static int nEvents(final int[]... hasEvents) {
        int n = 0;
        for (final int[] hasEvent : hasEvents) {
            n += MathUtils.sum(hasEvent);
        }
        return n;
    }

    private boolean[] calculateSkipArray(final SAMRecord read, final List<? extends Locatable> skipLocations) {
        final byte[] bases = read.getReadBases();
        final boolean[] skip = new boolean[bases.length];
        final boolean[] knownSites = calculateKnownSites(read, skipLocations);
        for (int iii = 0; iii < bases.length; iii++) {
            skip[iii] = !BaseUtils.isRegularBase(bases[iii]) || isLowQualityBase(read, iii) || knownSites[iii] || badSolidOffset(read, iii);
        }
        return skip;
    }

    protected boolean badSolidOffset(final SAMRecord read, final int offset) {
        return ReadUtils.isSOLiDRead(read) && BRAC.RAC.SOLID_RECAL_MODE != RecalUtils.SOLID_RECAL_MODE.DO_NOTHING && !RecalUtils.isColorSpaceConsistent(read, offset);
    }

    protected boolean[] calculateKnownSites(final SAMRecord read, final List<? extends Locatable> knownLocations) {
        final int readLength = read.getReadBases().length;
        final boolean[] knownSites = new boolean[readLength];
        Arrays.fill(knownSites, false);
        for (final Locatable feat : knownLocations) {
            int featureStartOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), feat.getStart(), ReadUtils.ClippingTail.LEFT_TAIL, true); // BUGBUG: should I use LEFT_TAIL here?
            if (featureStartOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED) {
                featureStartOnRead = 0;
            }

            int featureEndOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), feat.getEnd(), ReadUtils.ClippingTail.LEFT_TAIL, true);
            if (featureEndOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED) {
                featureEndOnRead = readLength;
            }

            if (featureStartOnRead > readLength) {
                featureStartOnRead = featureEndOnRead = readLength;
            }

            Arrays.fill(knownSites, Math.max(0, featureStartOnRead), Math.min(readLength, featureEndOnRead + 1), true);
        }
        return knownSites;
    }

    // TODO: can be merged with calculateIsIndel
    protected static int[] calculateIsSNP(final SAMRecord read, /*final BasesCache ref*/ final ReferenceContext ref, final SAMRecord originalRead) {
        final byte[] readBases = read.getReadBases();
        final byte[] refBases = Arrays.copyOfRange(ref.getBases(), read.getAlignmentStart() - originalRead.getAlignmentStart(), ref.getBases().length + read.getAlignmentEnd() - originalRead.getAlignmentEnd());
        int refPos = 0;


        final int[] snp = new int[readBases.length];
        int readPos = 0;
        for (final CigarElement ce : read.getCigar().getCigarElements()) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    for (int iii = 0; iii < elementLength; iii++) {
                        snp[readPos] = (BaseUtils.basesAreEqual(readBases[readPos], refBases[refPos]) ? 0 : 1);
                        readPos++;
                        refPos++;
                    }
                    break;
                case D:
                case N:
                    refPos += elementLength;
                    break;
                case I:
                case S: // ReferenceContext doesn't have the soft clipped bases!
                    readPos += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return snp;
    }

    protected static int[] calculateIsIndel(final SAMRecord read, final EventType mode) {
        final int[] indel = new int[read.getReadBases().length];
        int readPos = 0;
        for (final CigarElement ce : read.getCigar().getCigarElements()) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                case S: {
                    readPos += elementLength;
                    break;
                }
                case D: {
                    final int index = (read.getReadNegativeStrandFlag() ? readPos : readPos - 1);
                    updateIndel(indel, index, mode, EventType.BASE_DELETION);
                    break;
                }
                case I: {
                    final boolean forwardStrandRead = !read.getReadNegativeStrandFlag();
                    if (forwardStrandRead) {
                        updateIndel(indel, readPos - 1, mode, EventType.BASE_INSERTION);
                    }
                    readPos += elementLength;
                    if (!forwardStrandRead) {
                        updateIndel(indel, readPos, mode, EventType.BASE_INSERTION);
                    }
                    break;
                }
                case N:
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return indel;
    }

    private static void updateIndel(final int[] indel, final int index, final EventType mode, final EventType requiredMode) {
        if (mode == requiredMode && index >= 0 && index < indel.length)
            // protect ourselves from events at the start or end of the read (1D3M or 3M1D)
            indel[index] = 1;
    }

    protected static double[] calculateFractionalErrorArray(final int[] errorArray, final byte[] baqArray) {
        if (errorArray.length != baqArray.length) {
            throw new GATKException("Array length mismatch detected. Malformed read?");
        }

        final int BLOCK_START_UNSET = -1;

        final double[] fractionalErrors = new double[baqArray.length];
        Arrays.fill(fractionalErrors, 0.0);
        boolean inBlock = false;
        int blockStartIndex = BLOCK_START_UNSET;
        int iii;
        for (iii = 0; iii < fractionalErrors.length; iii++) {
            if (baqArray[iii] == NO_BAQ_UNCERTAINTY) {
                if (!inBlock) {
                    fractionalErrors[iii] = (double) errorArray[iii];
                } else {
                    calculateAndStoreErrorsInBlock(iii, blockStartIndex, errorArray, fractionalErrors);
                    inBlock = false; // reset state variables
                    blockStartIndex = BLOCK_START_UNSET; // reset state variables
                }
            } else {
                inBlock = true;
                if (blockStartIndex == BLOCK_START_UNSET) {
                    blockStartIndex = iii;
                }
            }
        }
        if (inBlock) {
            calculateAndStoreErrorsInBlock(iii - 1, blockStartIndex, errorArray, fractionalErrors);
        }
        if (fractionalErrors.length != errorArray.length) {
            throw new GATKException("Output array length mismatch detected. Malformed read?");
        }
        return fractionalErrors;
    }

    private static void calculateAndStoreErrorsInBlock(final int iii,
                                                       final int blockStartIndex,
                                                       final int[] errorArray,
                                                       final double[] fractionalErrors) {
        int totalErrors = 0;
        for (int jjj = Math.max(0, blockStartIndex - 1); jjj <= iii; jjj++) {
            totalErrors += errorArray[jjj];
        }
        for (int jjj = Math.max(0, blockStartIndex - 1); jjj <= iii; jjj++) {
            fractionalErrors[jjj] = ((double) totalErrors) / ((double) (iii - Math.max(0, blockStartIndex - 1) + 1));
        }
    }

    /**
     * Create a BAQ style array that indicates no alignment uncertainty
     *
     * @param read the read for which we want a BAQ array
     * @return a BAQ-style non-null byte[] counting NO_BAQ_UNCERTAINTY values
     * // TODO -- could be optimized avoiding this function entirely by using this inline if the calculation code above
     */
    protected static byte[] flatBAQArray(final SAMRecord read) {
        final byte[] baq = new byte[read.getReadLength()];
        Arrays.fill(baq, NO_BAQ_UNCERTAINTY);
        return baq;
    }

    /**
     * Compute an actual BAQ array for read, based on its quals and the reference sequence
     *
     * @param read the read to BAQ
     * @return a non-null BAQ tag array for read
     */
    private byte[] calculateBAQArray(final SAMRecord read) {
        baq.baqRead(read, referenceDataSource, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
        return BAQ.getBAQTag(read);
    }

    public Object onTraversalDone() {
        recalibrationEngine.finalizeData();

        logger.info("Calculating quantized quality scores...");
        quantizeQualityScores();

        return accumulator;
    }

    private RecalibrationTables getRecalibrationTable() {
        return recalibrationEngine.getFinalRecalibrationTables();
    }

    /**
     * go through the quality score table and use the # observations and the empirical quality score
     * to build a quality score histogram for quantization. Then use the QuantizeQual algorithm to
     * generate a quantization map (recalibrated_qual -> quantized_qual)
     */
    private void quantizeQualityScores() {
        quantizationInfo = new QuantizationInfo(getRecalibrationTable(), BRAC.RAC.QUANTIZING_LEVELS);
    }

    private void generateReport() {
        RecalUtils.outputRecalibrationReport(BRAC.RAC, quantizationInfo, getRecalibrationTable(), covariates, BRAC.RAC.SORT_BY_ALL_COLUMNS);
    }

    /**
     * Generates the report based on a finalized recalibrationTables.
     */
    public void saveReport(RecalibrationTables rt, StandardCovariateList finalRequestedCovariates) throws java.io.FileNotFoundException {
        BRAC.RAC.RECAL_TABLE = new PrintStream(BRAC.RAC.RECAL_TABLE_FILE);
        QuantizationInfo quantizationInfo = getQuantizationInfo(rt); //new QuantizationInfo(rt, BRAC.RAC.QUANTIZING_LEVELS);
        RecalUtils.outputRecalibrationReport(BRAC.RAC, quantizationInfo, rt, finalRequestedCovariates, BRAC.RAC.SORT_BY_ALL_COLUMNS);
    }
}