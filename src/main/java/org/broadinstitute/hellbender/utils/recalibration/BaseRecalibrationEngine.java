package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.Arrays;

public final class BaseRecalibrationEngine implements Serializable {
    private static final long serialVersionUID = 1L;

    protected static final Logger logger = LogManager.getLogger(BaseRecalibrationEngine.class);

    private RecalibrationArgumentCollection recalArgs;

    private RecalibrationTables recalTables;

    private SAMFileHeader readsHeader;

    /**
     * list to hold the all the covariate objects that were requested (required + standard + experimental)
     */
    private StandardCovariateList covariates;

    private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
    private static final byte NO_BAQ_UNCERTAINTY = (byte)'@';

    private long numReadsProcessed = 0L;

    /**
     * Has finalizeData() been called?
     */
    private boolean finalized = false;

    public BaseRecalibrationEngine( final RecalibrationArgumentCollection recalArgs, final SAMFileHeader readsHeader ) {
        this.recalArgs = recalArgs;
        this.readsHeader = readsHeader;

        baq = new BAQ(recalArgs.BAQGOP); // setup the BAQ object with the provided gap open penalty

        covariates = new StandardCovariateList(recalArgs, readsHeader);

        logger.info("The covariates being used here: ");
        for (final Covariate cov : covariates) { // list all the covariates being used
            logger.info('\t' + cov.getClass().getSimpleName());
        }

        final int numReadGroups = readsHeader.getReadGroups().size();
        if ( numReadGroups < 1 ) {
            throw new UserException("Number of read groups must be >= 1, but is " + numReadGroups);
        }
        recalTables = new RecalibrationTables(covariates, numReadGroups);
    }

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     */
    public void processRead( final GATKRead originalRead, final ReferenceDataSource refDS, final Iterable<? extends Locatable> knownSites ) {
        ReadTransformer transform = makeReadTransform();
        final GATKRead read = transform.apply(originalRead);

        if( read.isEmpty() ) {
            return; // the whole read was inside the adaptor so skip it
        }

        RecalUtils.parsePlatformForRead(read, readsHeader, recalArgs);
        // parse the solid color space and check for color no-calls
        if ( ! RecalUtils.isColorSpaceConsistent(recalArgs.SOLID_NOCALL_STRATEGY, read, readsHeader)) {
            return; // skip this read completely
        }

        final int[] isSNP = calculateIsSNP(read, refDS);
        final int[] isInsertion = calculateIsIndel(read, EventType.BASE_INSERTION);
        final int[] isDeletion = calculateIsIndel(read, EventType.BASE_DELETION);
        final int nErrors = nEvents(isSNP, isInsertion, isDeletion);

        // note for efficiency regions we don't compute the BAQ array unless we actually have
        // some error to marginalize over.  For ILMN data ~85% of reads have no error
        final byte[] baqArray = nErrors == 0 ? flatBAQArray(read) : calculateBAQArray(read, refDS);

        if( baqArray != null ) { // some reads just can't be BAQ'ed
            final ReadCovariates covariates = RecalUtils.computeCovariates(read, readsHeader, this.covariates);
            final boolean[] skip = calculateSkipArray(read, knownSites); // skip known sites of variation as well as low quality and non-regular bases
            final double[] snpErrors = calculateFractionalErrorArray(isSNP, baqArray);
            final double[] insertionErrors = calculateFractionalErrorArray(isInsertion, baqArray);
            final double[] deletionErrors = calculateFractionalErrorArray(isDeletion, baqArray);

            // aggregate all of the info into our info object, and update the data
            final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, covariates, skip, snpErrors, insertionErrors, deletionErrors);
            updateRecalTablesForRead(info);
        }

        numReadsProcessed++;
    }

    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to processRead have been issued.
     *
     * Assumes that all of the principal tables (by quality score) have been completely updated,
     * and walks over this data to create summary data tables like by read group table.
     */
    public void finalizeData() {
        if ( finalized ) {
            throw new IllegalStateException("FinalizeData() has already been called");
        }

        finalizeRecalibrationTables(recalTables);
        finalized = true;
    }

    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to updateDataForRead have been issued.
     *
     * Assumes that all of the principal tables (by quality score) have been completely updated,
     * and walks over this data to create summary data tables like by read group table.
     */
    public static void finalizeRecalibrationTables( final RecalibrationTables tables ) {
        Utils.nonNull(tables);
        final NestedIntegerArray<RecalDatum> byReadGroupTable = tables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> byQualTable = tables.getQualityScoreTable();

        // iterate over all values in the qual table
        for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int eventIndex = leaf.keys[2];
            final RecalDatum rgDatum = byReadGroupTable.get(rgKey, eventIndex);
            final RecalDatum qualDatum = leaf.value;

            if ( rgDatum == null ) {
                // create a copy of qualDatum, and initialize byReadGroup table with it
                byReadGroupTable.put(new RecalDatum(qualDatum), rgKey, eventIndex);
            } else {
                // combine the qual datum with the existing datum in the byReadGroup table
                rgDatum.combine(qualDatum);
            }
        }
    }

    /**
     * Get a possibly not-final recalibration table, to deal with distributed execution.
     */
    public RecalibrationTables getRecalibrationTables() {
        return recalTables;
    }

    /**
     * Get the final recalibration tables, after finalizeData() has been called
     *
     * This returns the finalized recalibration table collected by this engine.
     *
     * It is an error to call this function before finalizeData has been called
     *
     * @return the finalized recalibration table collected by this engine
     */
    public RecalibrationTables getFinalRecalibrationTables() {
        if ( ! finalized ) {
            throw new IllegalStateException("Cannot get final recalibration tables until finalizeData() has been called");
        }
        return recalTables;
    }

    public StandardCovariateList getCovariates() {
        return covariates;
    }

    public long getNumReadsProcessed() {
        return numReadsProcessed;
    }

    /**
     * Update the recalibration statistics using the information in recalInfo
     * @param recalInfo data structure holding information about the recalibration values for a single read
     */
    private void updateRecalTablesForRead( final ReadRecalibrationInfo recalInfo ) {
        if ( finalized ) {
            throw new IllegalStateException("FinalizeData() has already been called");
        }

        final GATKRead read = recalInfo.getRead();
        final ReadCovariates readCovariates = recalInfo.getCovariatesValues();
        final NestedIntegerArray<RecalDatum> qualityScoreTable = recalTables.getQualityScoreTable();

        for( int offset = 0; offset < read.getBases().length; offset++ ) {
            if( ! recalInfo.skip(offset) ) {

                for (final EventType eventType : EventType.values()) {
                    final int[] keys = readCovariates.getKeySet(offset, eventType);
                    final int eventIndex = eventType.ordinal();
                    final byte qual = recalInfo.getQual(eventType, offset);
                    final double isError = recalInfo.getErrorFraction(eventType, offset);

                    RecalUtils.incrementDatumOrPutIfNecessary(qualityScoreTable, qual, isError, keys[0], keys[1], eventIndex);

                    for (int i = 2; i < covariates.size(); i++) { //XXX the 2 is hard-wired here as the number of special covariates
                        if (keys[i] < 0) {
                            continue;
                        }

                        RecalUtils.incrementDatumOrPutIfNecessary(recalTables.getTable(i), qual, isError, keys[0], keys[1], keys[i], eventIndex);
                    }
                }
            }
        }
    }


    private ReadTransformer makeReadTransform() {
        ReadTransformer f0 = BaseRecalibrationEngine::consolidateCigar;

        ReadTransformer f = f0.andThen(this::setDefaultBaseQualities)
                .andThen(this::resetOriginalBaseQualities)
                .andThen(ReadClipper::hardClipAdaptorSequence)
                .andThen(ReadClipper::hardClipSoftClippedBases);

        return f;
    }

    private boolean isLowQualityBase( final GATKRead read, final int offset ) {
        return read.getBaseQualities()[offset] < recalArgs.PRESERVE_QSCORES_LESS_THAN;
    }

    private static GATKRead consolidateCigar( final GATKRead read ) {
        // Always consolidate the cigar string into canonical form, collapsing zero-length / repeated cigar elements.
        // Downstream code cannot necessarily handle non-consolidated cigar strings.
        read.setCigar(AlignmentUtils.consolidateCigar(read.getCigar()));
        return read;
    }

    private GATKRead resetOriginalBaseQualities( final GATKRead read ) {
        if (! recalArgs.useOriginalBaseQualities) {
            return read;
        }
        return ReadUtils.resetOriginalBaseQualities(read);
    }

    private GATKRead setDefaultBaseQualities( final GATKRead read ) {
        // if we are using default quals, check if we need them, and add if necessary.
        // 1. we need if reads are lacking or have incomplete quality scores
        // 2. we add if defaultBaseQualities has a positive value
        if (recalArgs.defaultBaseQualities < 0) {
            return read;
        }
        byte reads[] = read.getBases();
        byte quals[] = read.getBaseQualities();
        if (quals == null || quals.length < reads.length) {
            byte new_quals[] = new byte[reads.length];
            Arrays.fill(new_quals, recalArgs.defaultBaseQualities);
            read.setBaseQualities(new_quals);
        }
        return read;
    }

    /**
     * Compute the number of mutational events across all hasEvent vectors
     *
     * Simply the sum of entries in hasEvents
     *
     * @param hasEvents a vector a vectors of 0 (no event) and 1 (has event)
     * @return the total number of events across all hasEvent arrays
     */
    protected static int nEvents( final int[]... hasEvents ) {
        int n = 0;
        for ( final int[] hasEvent : hasEvents ) {
            n += MathUtils.sum(hasEvent);
        }
        return n;
    }

    private boolean[] calculateSkipArray( final GATKRead read, final Iterable<? extends Locatable> knownSites ) {
        final byte[] bases = read.getBases();
        final boolean[] skip = new boolean[bases.length];
        final boolean[] knownSitesArray = calculateKnownSites(read, knownSites);
        for( int iii = 0; iii < bases.length; iii++ ) {
            skip[iii] = !BaseUtils.isRegularBase(bases[iii]) || isLowQualityBase(read, iii) || knownSitesArray[iii] || badSolidOffset(read, iii);
        }
        return skip;
    }

    protected boolean badSolidOffset( final GATKRead read, final int offset ) {
        return ReadUtils.isSOLiDRead(read, readsHeader) && recalArgs.SOLID_RECAL_MODE != RecalUtils.SOLID_RECAL_MODE.DO_NOTHING && !RecalUtils.isColorSpaceConsistent(read, offset);
    }

    protected boolean[] calculateKnownSites( final GATKRead read, final Iterable<? extends Locatable> knownSites ) {
        final int readLength = read.getBases().length;
        final boolean[] knownSitesArray = new boolean[readLength];//initializes to all false
        for ( final Locatable knownSite : knownSites ) {
            int featureStartOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), knownSite.getStart(), ReadUtils.ClippingTail.LEFT_TAIL, true);
            if( featureStartOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureStartOnRead = 0;
            }

            int featureEndOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), knownSite.getEnd(), ReadUtils.ClippingTail.LEFT_TAIL, true);
            if( featureEndOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureEndOnRead = readLength;
            }

            if( featureStartOnRead > readLength ) {
                featureStartOnRead = featureEndOnRead = readLength;
            }

            Arrays.fill(knownSitesArray, Math.max(0, featureStartOnRead), Math.min(readLength, featureEndOnRead + 1), true);
        }
        return knownSitesArray;
    }

    // TODO: can be merged with calculateIsIndel
    protected static int[] calculateIsSNP( final GATKRead read, final ReferenceDataSource refDS ) {
        final byte[] readBases = read.getBases();
        final byte[] refBases = refDS.queryAndPrefetch(read.getContig(), read.getStart(), read.getEnd()).getBases();

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

    public static int[] calculateIsIndel( final GATKRead read, final EventType mode ) {
        final int[] indel = new int[read.getBases().length];
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
                    final int index = ( read.isReverseStrand() ? readPos : readPos - 1 );
                    updateIndel(indel, index, mode, EventType.BASE_DELETION);
                    break;
                }
                case I:
                {
                    final boolean forwardStrandRead = !read.isReverseStrand();
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
        {
            indel[index] = 1;
        }
    }

    public static double[] calculateFractionalErrorArray( final int[] errorArray, final byte[] baqArray ) {
        if ( errorArray.length != baqArray.length ) {
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
    protected  static byte[] flatBAQArray(final GATKRead read) {
        final byte[] baq = new byte[read.getLength()];
        Arrays.fill(baq, NO_BAQ_UNCERTAINTY);
        return baq;
    }

    /**
     * Compute an actual BAQ array for read, based on its quals and the reference sequence
     * @param read the read to BAQ
     * @return a non-null BAQ tag array for read
     */
    private byte[] calculateBAQArray( final GATKRead read, final ReferenceDataSource refDS ) {
        baq.baqRead(read, refDS, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
        return BAQ.getBAQTag(read);
    }


}