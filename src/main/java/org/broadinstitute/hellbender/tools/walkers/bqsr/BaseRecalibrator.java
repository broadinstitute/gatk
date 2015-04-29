package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.tribble.Feature;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.*;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import static org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.*;

/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various covariates
 * (such as read group, reported quality score, machine cycle, and nucleotide context).
 *
 * <p>
 * This walker is designed to work as the first pass in a two-pass processing step. It does a by-locus traversal operating
 * only at sites that are not in dbSNP. We assume that all reference mismatches we see are therefore errors and indicative
 * of poor base quality. This walker generates tables based on various user-specified covariates (such as read group,
 * reported quality score, cycle, and context). Since there is a large amount of data one can then calculate an empirical
 * probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.
 * The output file is a table (of the several covariate values, num observations, num mismatches, empirical quality score).
 * <p>
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and will be added for the user regardless of whether or not they were specified.
 *
 * <p>
 *
 * <h3>Input</h3>
 * <p>
 * The input read data whose base quality scores need to be assessed.
 * <p>
 * A database of known polymorphic sites to skip over.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A GATK Report file with many tables:
 * <ol>
 *     <li>The list of arguments</li>
 *     <li>The quantized qualities table</li>
 *     <li>The recalibration table by read group</li>
 *     <li>The recalibration table by quality score</li>
 *     <li>The recalibration table for all the optional covariates</li>
 * </ol>
 *
 * The GATK Report is intended to be easy to read by humans or computers. Check out the documentation of the GATKReport to learn how to manipulate this table.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T BaseRecalibrator \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -o recal_data.table
 * </pre>
 */

@CommandLineProgramProperties(
        usage = "First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        usageShort = "Generates recalibration table",
        programGroup = ReadProgramGroup.class
)
public class BaseRecalibrator extends ReadWalker {
    final protected static Logger logger = LogManager.getLogger(BaseRecalibrator.class);

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc="all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    @Argument(fullName = "bqsrBAQGapOpenPenalty", shortName="bqsrBAQGOP", doc="BQSR BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets", optional = true)
    public double BAQGOP = BAQ.DEFAULT_GOP;


    /**
     * This flag tells GATK not to modify quality scores less than this value. Instead they will be written out unmodified in the recalibrated BAM file.
     * In general it's unsafe to change qualities scores below < 6, since base callers use these values to indicate random or bad bases.
     * For example, Illumina writes Q2 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 bin can be elevated to Q8 or Q10, leading to issues downstream.
     */
    @Argument(fullName = "preserve_qscores_less_than", shortName = "preserveQ", doc = "Don't recalibrate bases with quality scores less than this threshold (with -BQSR)", optional = true)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;


    // --------------------------------------------------------------------------------------------------------------
    //
    // quality encoding checking arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * This flag tells GATK to use the original base qualities (that were in the data before BQSR/recalibration) which
     * are stored in the OQ tag, if they are present, rather than use the post-recalibration quality scores. If no OQ
     * tag is present for a read, the standard qual score will be used.
     */
    @Argument(fullName="useOriginalQualities", shortName = "OQ", doc = "Use the base quality scores from the OQ tag", optional = true)
    public Boolean useOriginalBaseQualities = false;

    /**
     * If reads are missing some or all base quality scores, this value will be used for all base quality scores.
     * By default this is set to -1 to disable default base quality assignment.
     */
    //TODO: minValue = 0, maxValue = Byte.MAX_VALUE)
    @Argument(fullName="defaultBaseQualities", shortName = "DBQ", doc = "Assign a default base quality", optional = true)
    public byte defaultBaseQualities = -1;

    /**
     * an object that keeps track of the information necessary for quality score quantization
     */
    private QuantizationInfo quantizationInfo;

    /**
     * list to hold the all the covariate objects that were requested (required + standard + experimental)
     */
    private Covariate[] requestedCovariates;

    private RecalibrationEngine recalibrationEngine;

    private int minimumQToUse;

    private static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
    private ReferenceDataSource referenceDataSource; // datasource for the reference. We're using a different one from the engine itself to avoid messing with its caches.
    private final static byte NO_BAQ_UNCERTAINTY = (byte)'@';


    private long accumulator;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    @Override
    public void onTraversalStart() {
        accumulator= 0L;

        baq = new BAQ(BAQGOP); // setup the BAQ object with the provided gap open penalty

        if (RAC.FORCE_PLATFORM != null)
            RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;

        if (RAC.knownSites.isEmpty() && !RAC.RUN_WITHOUT_DBSNP) // Warn the user if no dbSNP file or other variant mask was specified
            throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);

        requestedCovariates = getCovariatesArray();

        logger.info("The covariates being used here: ");
        for (Covariate cov : requestedCovariates) { // list all the covariates being used
            logger.info("\t" + cov.getClass().getSimpleName());
            cov.initialize(RAC); // initialize any covariate member variables using the shared argument collection
        }

        try {
            RAC.RECAL_TABLE = new PrintStream(RAC.RECAL_TABLE_FILE);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(RAC.RECAL_TABLE_FILE, e);
        }

        initializeRecalibrationEngine();
        minimumQToUse = PRESERVE_QSCORES_LESS_THAN;

        referenceDataSource = new ReferenceDataSource(referenceArguments.getReferenceFile());
    }

    private Covariate[] getCovariatesArray() {
        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalUtils.initializeCovariates(RAC); // initialize the required and optional covariates
        List<Covariate> covariatesList = new ArrayList<>(covariates.getLeft().size() + covariates.getRight().size());
        covariatesList.addAll(covariates.getLeft());
        covariatesList.addAll(covariates.getRight());
        return covariatesList.toArray(new Covariate[covariatesList.size()]);
    }

    /**
     * Initialize the recalibration engine
     */
    private void initializeRecalibrationEngine() {
        int numReadGroups = getHeaderForReads().getReadGroups().size();

        recalibrationEngine = new RecalibrationEngine(requestedCovariates, numReadGroups);
    }

    private boolean isLowQualityBase( final SAMRecord read, final int offset ) {
        return read.getBaseQualities()[offset] < minimumQToUse;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return super.makeReadFilter()
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

    private SAMRecord resetOriginalBaseQualities(SAMRecord read){
        if (! useOriginalBaseQualities) {
            return read;
        }
        return ReadUtils.resetOriginalBaseQualities(read);
    }

    private SAMRecord setDefaultBaseQualities(SAMRecord read) {
        // if we are using default quals, check if we need them, and add if necessary.
        // 1. we need if reads are lacking or have incomplete quality scores
        // 2. we add if defaultBaseQualities has a positive value
        if (defaultBaseQualities < 0) {
            return read;
        }
        byte reads[] = read.getReadBases();
        byte quals[] = read.getBaseQualities();
        if (quals == null || quals.length < reads.length) {
            byte new_quals[] = new byte[reads.length];
            Arrays.fill(new_quals, defaultBaseQualities);
            read.setBaseQualities(new_quals);
        }
        return read;
    }


    private ReadTransformer makeReadTransform(){
        ReadTransformer f0 = BaseRecalibrator::consolidateCigar;

        ReadTransformer f =
                f0.andThen(this::setDefaultBaseQualities)
                .andThen(this::resetOriginalBaseQualities)
                .andThen(ReadClipper::hardClipAdaptorSequence)
                .andThen(ReadClipper::hardClipSoftClippedBases);

        return f;
    }
    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     */
    @Override
    public void apply( SAMRecord originalRead, ReferenceContext ref, FeatureContext featureContext ) {
        ReadTransformer transform = makeReadTransform();
        final SAMRecord read = transform.apply(originalRead);

        if( ReadUtils.isEmpty(read) ) { return; } // the whole read was inside the adaptor so skip it

        RecalUtils.parsePlatformForRead(read, RAC);
        if (!RecalUtils.isColorSpaceConsistent(RAC.SOLID_NOCALL_STRATEGY, read)) { // parse the solid color space and check for color no-calls
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
            final ReadCovariates covariates = RecalUtils.computeCovariates(read, requestedCovariates);
            final boolean[] skip = calculateSkipArray(read, featureContext); // skip known sites of variation as well as low quality and non-regular bases
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
     *
     * Simply the sum of entries in hasEvents
     *
     * @param hasEvents a vector a vectors of 0 (no event) and 1 (has event)
     * @return the total number of events across all hasEvent arrays
     */
    protected static int nEvents(final int[]... hasEvents) {
        int n = 0;
        for ( final int[] hasEvent : hasEvents ) {
            n += MathUtils.sum(hasEvent);
        }
        return n;
    }

    private boolean[] calculateSkipArray( final SAMRecord read, final FeatureContext featureContext ) {
        final byte[] bases = read.getReadBases();
        final boolean[] skip = new boolean[bases.length];
        final boolean[] knownSites = calculateKnownSites(read, featureContext.getValues(RAC.knownSites));
        for( int iii = 0; iii < bases.length; iii++ ) {
            skip[iii] = !BaseUtils.isRegularBase(bases[iii]) || isLowQualityBase(read, iii) || knownSites[iii] || badSolidOffset(read, iii);
        }
        return skip;
    }

    protected boolean badSolidOffset( final SAMRecord read, final int offset ) {
        return ReadUtils.isSOLiDRead(read) && RAC.SOLID_RECAL_MODE != RecalUtils.SOLID_RECAL_MODE.DO_NOTHING && !RecalUtils.isColorSpaceConsistent(read, offset);
    }

    protected boolean[] calculateKnownSites( final SAMRecord read, final List<? extends Feature> features) {
        final int readLength = read.getReadBases().length;
        final boolean[] knownSites = new boolean[readLength];
        Arrays.fill(knownSites, false);
        for( final Feature feat : features ) {
            int featureStartOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), feat.getStart(), ReadUtils.ClippingTail.LEFT_TAIL, true); // BUGBUG: should I use LEFT_TAIL here?
            if( featureStartOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureStartOnRead = 0;
            }

            int featureEndOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), feat.getEnd(), ReadUtils.ClippingTail.LEFT_TAIL, true);
            if( featureEndOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureEndOnRead = readLength;
            }

            if( featureStartOnRead > readLength ) {
                featureStartOnRead = featureEndOnRead = readLength;
            }

            Arrays.fill(knownSites, Math.max(0, featureStartOnRead), Math.min(readLength, featureEndOnRead + 1), true);
        }
        return knownSites;
    }

    // TODO: can be merged with calculateIsIndel
    protected static int[] calculateIsSNP( final SAMRecord read, final ReferenceContext ref, final SAMRecord originalRead ) {
        final byte[] readBases = read.getReadBases();
        final byte[] refBases = Arrays.copyOfRange(ref.getBases(), read.getAlignmentStart() - originalRead.getAlignmentStart(), ref.getBases().length + read.getAlignmentEnd() - originalRead.getAlignmentEnd());
        final int[] snp = new int[readBases.length];
        int readPos = 0;
        int refPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        snp[readPos] = ( BaseUtils.basesAreEqual(readBases[readPos], refBases[refPos]) ? 0 : 1 );
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

    protected static int[] calculateIsIndel( final SAMRecord read, final EventType mode ) {
        final int[] indel = new int[read.getReadBases().length];
        int readPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                case S:
                {
                    readPos += elementLength;
                    break;
                }
                case D:
                {
                    final int index = ( read.getReadNegativeStrandFlag() ? readPos : readPos - 1 );
                    updateIndel(indel, index, mode, EventType.BASE_DELETION);
                    break;
                }
                case I:
                {
                    final boolean forwardStrandRead = !read.getReadNegativeStrandFlag();
                    if( forwardStrandRead ) {
                        updateIndel(indel, readPos - 1, mode, EventType.BASE_INSERTION);
                    }
                    readPos += elementLength;
                    if( !forwardStrandRead ) {
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
        if ( mode == requiredMode && index >= 0 && index < indel.length )
            // protect ourselves from events at the start or end of the read (1D3M or 3M1D)
            indel[index] = 1;
    }

    protected static double[] calculateFractionalErrorArray( final int[] errorArray, final byte[] baqArray ) {
        if(errorArray.length != baqArray.length ) {
            throw new GATKException("Array length mismatch detected. Malformed read?");
        }

        final int BLOCK_START_UNSET = -1;

        final double[] fractionalErrors = new double[baqArray.length];
        Arrays.fill(fractionalErrors, 0.0);
        boolean inBlock = false;
        int blockStartIndex = BLOCK_START_UNSET;
        int iii;
        for( iii = 0; iii < fractionalErrors.length; iii++ ) {
            if( baqArray[iii] == NO_BAQ_UNCERTAINTY ) {
                if( !inBlock ) {
                    fractionalErrors[iii] = (double) errorArray[iii];
                } else {
                    calculateAndStoreErrorsInBlock(iii, blockStartIndex, errorArray, fractionalErrors);
                    inBlock = false; // reset state variables
                    blockStartIndex = BLOCK_START_UNSET; // reset state variables
                }
            } else {
                inBlock = true;
                if( blockStartIndex == BLOCK_START_UNSET ) { blockStartIndex = iii; }
            }
        }
        if( inBlock ) {
            calculateAndStoreErrorsInBlock(iii-1, blockStartIndex, errorArray, fractionalErrors);
        }
        if( fractionalErrors.length != errorArray.length ) {
            throw new GATKException("Output array length mismatch detected. Malformed read?");
        }
        return fractionalErrors;
    }

    private static void calculateAndStoreErrorsInBlock( final int iii,
                                                        final int blockStartIndex,
                                                        final int[] errorArray,
                                                        final double[] fractionalErrors ) {
        int totalErrors = 0;
        for( int jjj = Math.max(0, blockStartIndex - 1); jjj <= iii; jjj++ ) {
            totalErrors += errorArray[jjj];
        }
        for( int jjj = Math.max(0, blockStartIndex - 1); jjj <= iii; jjj++ ) {
            fractionalErrors[jjj] = ((double) totalErrors) / ((double)(iii - Math.max(0, blockStartIndex - 1) + 1));
        }
    }

    /**
     * Create a BAQ style array that indicates no alignment uncertainty
     * @param read the read for which we want a BAQ array
     * @return a BAQ-style non-null byte[] counting NO_BAQ_UNCERTAINTY values
     * // TODO -- could be optimized avoiding this function entirely by using this inline if the calculation code above
     */
    protected  static byte[] flatBAQArray(final SAMRecord read) {
        final byte[] baq = new byte[read.getReadLength()];
        Arrays.fill(baq, NO_BAQ_UNCERTAINTY);
        return baq;
    }

    /**
     * Compute an actual BAQ array for read, based on its quals and the reference sequence
     * @param read the read to BAQ
     * @return a non-null BAQ tag array for read
     */
    private byte[] calculateBAQArray( final SAMRecord read ) {
        baq.baqRead(read,referenceDataSource, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
        return BAQ.getBAQTag(read);
    }

    @Override
    public Object onTraversalDone() {
        recalibrationEngine.finalizeData();

        logger.info("Calculating quantized quality scores...");
        quantizeQualityScores();

        logger.info("Writing recalibration report...");
        generateReport();
        logger.info("...done!");

        //logger.info("BaseRecalibrator was able to recalibrate " + result + " reads");
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
        quantizationInfo = new QuantizationInfo(getRecalibrationTable(), RAC.QUANTIZING_LEVELS);
    }

    private void generateReport() {
        RecalUtils.outputRecalibrationReport(RAC, quantizationInfo, getRecalibrationTable(), requestedCovariates, RAC.SORT_BY_ALL_COLUMNS);
    }
}