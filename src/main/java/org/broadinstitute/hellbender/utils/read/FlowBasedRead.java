package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedHMMEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Tail;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
/**
 * Adds flow information to the usual GATKRead. In addition to the usual read data this class keeps flowMatrix,
 * that contains probabilities for alternative hmer calls.
 *
 * Main function deals with parsing flow-specific QUAL representation readFlowMatrix.
 * Note that there is a lot of code that deals with other various formats of the representation (e.g. when the matrix
 * is coded in the tags of the BAM and is given in flow space). This code is not used in production, but was used in
 * development and testing
 *
 * A common usage pattern is to covert a GATKRead into a FlowBasedRead. Additionally
 * a SAMRecord can also be converted into a FlowBasedRead. Follows a common usage pattern:
 *
 * For a self contained example of a usage pattern, see {@link FlowBasedReadUtils#convertToFlowBasedRead(GATKRead, SAMFileHeader)}
 *
 **/

public class FlowBasedRead extends SAMRecordToGATKReadAdapter implements GATKRead, Serializable {

    final static public int     MAX_CLASS = 12;
    public static final String     DEFAULT_FLOW_ORDER = "TGCA";
    private static final long serialVersionUID = 42L;
    private final Logger logger = LogManager.getLogger(this.getClass());
    private static final OneShotLogger vestigialOneShotLogger = new OneShotLogger(FlowBasedRead.class);

    // constants
    static private final int MINIMAL_READ_LENGTH = 10; // check if this is the right number
    private final double MINIMAL_CALL_PROB = 0.1;

    // constants for clippingTagContains.
    // The tag is present when the end of the read was clipped at base calling.
    // The value of the tag is a string consisting of any one or more of the following:
    // A - adaptor clipped
    // Q - quality clipped
    // Z - "three zeros" clipped
    final  public static String  CLIPPING_TAG_NAME = "tm";

    final public static String FLOW_MATRIX_TAG_NAME = "tp";
    final public static String FLOW_MATRIX_T0_TAG_NAME = "t0";
    final public static String FLOW_MATRiX_OLD_TAG_KR = "kr";
    final public static String FLOW_MATRiX_OLD_TAG_TI = "ti";
    public static final String FLOW_MATRiX_OLD_TAG_KH = "kh";
    public static final String FLOW_MATRiX_OLD_TAG_KF = "kf";
    public static final String FLOW_MATRiX_OLD_TAG_KD = "kd";

    final public static String MAX_CLASS_READ_GROUP_TAG = "mc";

    final public static String[] HARD_CLIPPED_TAGS = new String[] {
            FLOW_MATRIX_TAG_NAME, FLOW_MATRiX_OLD_TAG_KR, FLOW_MATRiX_OLD_TAG_TI,
            FLOW_MATRiX_OLD_TAG_KH, FLOW_MATRiX_OLD_TAG_KF, FLOW_MATRiX_OLD_TAG_KD,
            FLOW_MATRIX_T0_TAG_NAME
    };

    /**
     * The sam record from which this flow based read originated
     */
    private SAMRecord samRecord;

    /**
     * The read's sequence, always in forward direction
     */
    private byte[] forwardSequence;

    /**
     * The flow key for the read - i.e. lengths of hmers in an flow order.
     *
     * For example, assuming a flow order of TGCA, and a forward sequence of GGAAT, the key will be 0,2,0,2,1
     */
    private int[] key;

    /**
     * the maping of key elements to their origin locations in the sequence. Entry n contains the offset in the sequence
     * where the hmer described by this key element starts.
     */
    private int [] flow2base;

    /**
     * The maximal length of an hmer that can be encoded (normally in the 10-12 range)
     */
    private int maxHmer;

    /**
     * The order in which flow key in encoded (See decription for key field). Flow order may be wrapped if a longer one
     * needed.
     */
    private byte[] flowOrder;

    /**
     * The probability matrix for this read. [n][m] position represents that probablity that an hmer of n length will be
     * present at the m key position. Therefore, the first dimention is in the maxHmer order, where the second dimension
     * is length(key).
     */
    private double[][] flowMatrix;

    /**
     * The validity status of the key. Certain operations may produce undefined/errornous results. This is signaled by
     * the read being marked with a validKey == false
     */
    private boolean validKey;

    /**
     * The direction of this read. After being red, the direction will also swing to be REFERENCE
     */
    private Direction direction = Direction.SYNTHESIS;

    /**
     * Was base clipping applied to this read? (normally to trim a read to the span of a haplotype)
     */
    private boolean baseClipped = false;

    /**
     * If applyBaseClipping was called, the left trimming that was actually applied to the read
     */
    private int trimLeftBase = 0 ;

    /**
     * If applyBaseClipping was called, the right trimming that was actually applied to the read
     */
    private int trimRightBase = 0 ;

    /**
     * The flow based argument collection under which this read was created
     */
    private final FlowBasedArgumentCollection fbargs;

    /**
     * This allows tools to reduce/enlarge the lower limit of read size for clipping operations.
     */
    static private int minimalReadLength = MINIMAL_READ_LENGTH;


    /**
     * Various methods for storing arrays of pre-computed qualities to be used in {@link FlowBasedHMMEngine}.
     */
    public byte[] getReadInsQuals() {
        return readInsQuals;
    }
    public void setReadInsQuals(byte[] readInsQuals) {
        this.readInsQuals = readInsQuals;
    }
    public byte[] getReadDelQuals() {
        return readDelQuals;
    }
    public void setReadDelQuals(byte[] readDelQuals) {
        this.readDelQuals = readDelQuals;
    }
    public byte[] getOverallGCP() {
        return overallGCP;
    }
    public void setOverallGCP(byte[] overallGCP) {
        this.overallGCP = overallGCP;
    }
    private byte[] readInsQuals;
    private byte[] readDelQuals;
    private byte[] overallGCP;

    public enum Direction {
        REFERENCE, SYNTHESIS
    }

    /**

     * Constructor from GATKRead. flow order, hmer and arguments
     * @param read GATK read
     * @param flowOrder flow order string (one cycle)
     * @param maxHmer maximal hmer to keep in the flow matrix
     * @param fbargs arguments that control resolution etc. of the flow matrix
     */
    public FlowBasedRead(final GATKRead read, final String flowOrder, final int maxHmer, final FlowBasedArgumentCollection fbargs) {
        this(read.convertToSAMRecord(null), flowOrder, maxHmer, fbargs);
    }

    /**
     * Same as above but constructs from SAMRecord
     * @param samRecord record from SAM file
     * @param flowOrder flow order (single cycle)
     * @param maxHmer maximal hmer to keep in the flow matrix
     * @param fbargs arguments that control resoltion of the flow matrix
     */
    public FlowBasedRead(final SAMRecord samRecord, final String flowOrder, final int maxHmer, final FlowBasedArgumentCollection fbargs) {
        super(samRecord);
        Utils.nonNull(fbargs);
        Utils.validate(FlowBasedReadUtils.isFlow(samRecord), "FlowBasedRead can only be used on flow reads. failing read: " + samRecord);
        this.fbargs = fbargs;
        this.maxHmer = maxHmer;
        this.samRecord = samRecord;
        forwardSequence = getForwardSequence();

        // read flow matrix in. note that below code contains accomodates for old formats
        if ( samRecord.hasAttribute(FLOW_MATRIX_TAG_NAME) ) {

            // this path is the production path. A flow read should contain a FLOW_MATRIX_TAG_NAME tag
            readFlowMatrix(flowOrder);

        } else {

            // NOTE: this path is vestigial and deals with old formats of the matrix
            if ( samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_KR) ) {
                readVestigialFlowMatrixFromKR(flowOrder);
            } else if ( samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_TI) ) {
                readVestigialFlowMatrixFromTI(flowOrder);
            } else {
                throw new GATKException("read missing flow matrix attribute: " + FLOW_MATRIX_TAG_NAME);
            }
        }
        implementMatrixMods(fbargs.getFlowMatrixModsInstructions());

        //Spread boundary flow probabilities when the read is unclipped
        //in this case the value of the hmer is uncertain
        if ( !fbargs.keepBoundaryFlows ) {
            if (CigarUtils.countClippedBases(samRecord.getCigar(), Tail.LEFT, CigarOperator.HARD_CLIP) == 0) {
                spreadFlowLengthProbsAcrossCountsAtFlow(findFirstNonZero(key));
            }
            if (CigarUtils.countClippedBases(samRecord.getCigar(), Tail.RIGHT, CigarOperator.HARD_CLIP) == 0) {
                spreadFlowLengthProbsAcrossCountsAtFlow(findLastNonZero(key));
            }
        }

        if ( logger.isDebugEnabled() ) {
            logger.debug("cons: name: " + samRecord.getReadName()
                    + " len: " + samRecord.getReadLength()
                    + " loc: " + samRecord.getStart() + "-" + samRecord.getEnd()
                    + " rev: " + isReverseStrand()
                    + " cigar:" + samRecord.getCigarString());
            logger.debug("     bases: " + new String(samRecord.getReadBases()));
            logger.debug("       key: " + FlowBasedKeyCodec.keyAsString(key));
        }

        validateSequence();
    }

    //since the last unclipped flow is uncertain (we give high probabilities to
    //also hmers higher than the called hmer)
    private void spreadFlowLengthProbsAcrossCountsAtFlow(final int flowToSpread) {
        if (flowToSpread<0) //boundary case when all the key is zero
            return;

        final int call = key[flowToSpread];
        if (call==0){
            throw new IllegalStateException("Boundary key value should not be zero for the spreading");
        }

        final int numberToFill = maxHmer - call+1;
        double total = 0;
        for (int i = call; i < maxHmer+1; i++)
            total += flowMatrix[i][flowToSpread];
        final double fillProb = Math.max(total / numberToFill, fbargs.fillingValue);
        for (int i = call; i < maxHmer+1; i++){
            flowMatrix[i][flowToSpread] = fillProb;
        }
    }



    // This is the code for parsing the current/production BAM format (with TP tag)
    private void readFlowMatrix(final String _flowOrder) {

        // generate key (base to flow space)
        setDirection(Direction.REFERENCE);  // base is always in reference/alignment direction

        key = FlowBasedKeyCodec.baseArrayToKey(samRecord.getReadBases(), _flowOrder);
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);
        flowOrder = FlowBasedKeyCodec.getFlowToBase(_flowOrder, key.length);

        // initialize matrix
        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = fbargs.fillingValue;
            }
        }

        // access qual, convert to flow representation
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      tp = samRecord.getSignedByteArrayAttribute(FLOW_MATRIX_TAG_NAME);
        boolean specialTreatmentForZeroCalls = false;
        final byte[]      t0 = SAMUtils.fastqToPhred(samRecord.getStringAttribute(FLOW_MATRIX_T0_TAG_NAME));
        final double[]     t0probs = new double[quals.length];
        if ((t0!=null) && fbargs.useT0Tag){
            specialTreatmentForZeroCalls = true;

            if (t0.length!=tp.length){
                throw new GATKException("Illegal read len(t0)!=len(qual): " + samRecord.getReadName());
            }
        }

        final double[]    probs = new double[quals.length];
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final double q = quals[i];
            final double p = Math.pow(10, -q/10);
            probs[i] = p;
            if (specialTreatmentForZeroCalls){
                final double qq = t0[i];
                t0probs[i] = Math.pow(10,-qq/10);
            }
        }

        // apply key and qual/tp to matrix
        int     qualOfs = 0; //converts between base -> flow
        for ( int i = 0 ; i < key.length ; i++ ) {
            final int        run = key[i];
            if (run > 0) {
                parseSingleHmer(probs, tp, i, run, qualOfs);
            }

            if ((run == 0)&&(specialTreatmentForZeroCalls)){
                parseZeroQuals(t0probs, i, qualOfs);
            }

            double totalErrorProb = 0;

            for (int k=0; k < maxHmer; k++ ){
                totalErrorProb += flowMatrix[k][i];
            }
            final double callProb = Math.max(MINIMAL_CALL_PROB, 1-totalErrorProb);
            // the probability in the recalibration is not divided by two for hmers of length 1
            flowMatrix[Math.min(run, maxHmer)][i] = callProb;
            qualOfs+=run;
        }
        applyFilteringFlowMatrix();
    }


    //convert qualities from the single hmer to a column in a flow matrix
    private void parseSingleHmer(final double[] probs, final byte[] tp, final int flowIdx,
                                 final int flowCall, final int qualOfs){
        for (int i = qualOfs ; i < qualOfs+flowCall; i++) {
            if (tp[i]!=0) {
                final int loc = Math.max(Math.min(flowCall+tp[i], maxHmer),0);
                if (flowMatrix[loc][flowIdx] == fbargs.fillingValue) {
                    flowMatrix[loc][flowIdx] = probs[i];
                } else {
                    flowMatrix[loc][flowIdx] += probs[i];
                }
            }
        }
    }

    //convert qualities from the single hmer to a column in a flow matrix
    private void parseZeroQuals(final double[] probs, final int flowIdx, final int qualOfs){
        if ((qualOfs == 0) | (qualOfs==probs.length)){ // do not report zero error probability on the edge of the read
            return;
        }
        if ((probs[qualOfs-1])==(probs[qualOfs])){
            flowMatrix[1][flowIdx] = Math.max(flowMatrix[1][flowIdx], Math.max(probs[qualOfs-1],probs[qualOfs]));
        }
    }

    public String getFlowOrder() {
        return new String(Arrays.copyOfRange(flowOrder, 0, Math.min(fbargs.flowOrderCycleLength,flowOrder.length)));
    }


    public int getMaxHmer() {
        return maxHmer;
    }

    public int getNFlows() {
        return key.length;
    }

    public Direction getDirection(){
        return direction;
    }


    private byte[] getForwardSequence(){
        if (!isReverseStrand()) {
            return samRecord.getReadBases();
        } else {
            final byte[] result = new byte[samRecord.getReadBases().length];
            System.arraycopy(samRecord.getReadBases(), 0, result, 0, result.length);
            SequenceUtil.reverseComplement(result);
            return result;
        }
    }


    private int[] getAttributeAsIntArray(final String attributeName, final boolean isSigned) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = this.samRecord.getAttribute(attributeName);

        if (attributeValue == null) {
            return null;
        } else if (attributeValue instanceof byte[]) {
            final byte[] byteArrayAttributeValue = (byte[]) attributeValue;
            final int[] ret = new int[byteArrayAttributeValue.length];
            for (int i = 0; i < ret.length; i++)
                if (!isSigned) ret[i] = byteArrayAttributeValue[i]&0xff;
                else ret[i]=byteArrayAttributeValue[i]; //converting signed byte to unsigned
            return Arrays.copyOf(ret, ret.length);
        } else if ((attributeValue instanceof int[])) {
            final int[] ret = (int[]) attributeValue;
            return Arrays.copyOf(ret, ret.length);
        } else if  (attributeValue instanceof short[]) {
            final short [] shortArrayAttributeValue = (short[]) attributeValue;
            final int[] ret = new int[shortArrayAttributeValue.length];
            for (int i = 0 ; i < shortArrayAttributeValue.length; i++ )
                ret[i] = shortArrayAttributeValue[i];
            return Arrays.copyOf(ret, ret.length);
        }else {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer array");
        }
    }


    private void validateSequence(){
        for (final int b : key) {
            if (b > maxHmer - 1) {
                validKey = false;
            }
        }
        validKey = true;
    }

    public boolean isValid() {
        return validKey;
    }

    /**
     * get a specific cell from the flow matrix. Each cell contains the probability
     * for an hmer of the given length to appear the given position in the flow key
     *
     * @param flow - position in the flow key (index into key[])
     * @param hmer - length of the hmer
     * @return
     */
    public double getProb(final int flow, final int hmer) {
        double prob = flowMatrix[hmer < maxHmer ? hmer : maxHmer][flow];
        return (prob <= 1) ? prob : 1;
    }

    /*
    * Legacy function from the time when the error probability were in flow space and when the read was clipped we had to
    * translate the clipping in the base space to the clipping in the flow space. Now does nothing (isBaseFormat() is true)
    * but is kept for R&D cases.
     */
    public void applyAlignment(){
        if ((getDirection() == Direction.SYNTHESIS) && ( isReverseStrand() )) {
            flipMatrix();
            ArrayUtils.reverse(key);
            flow2base = FlowBasedKeyCodec.getKeyToBase(key);
            SequenceUtil.reverseComplement(flowOrder);

        }

        final boolean isBase = isBaseFormat();
        final int[] basePair = {0, 0};
        final int[] clipLeftPair = !isBase ? findLeftClippingFromCigar() : basePair;
        final int[] clipRightPair = !isBase ? findRightClippingFromCigar() : basePair;
        final int clipLeft = clipLeftPair[0];
        final int leftHmerClip = clipLeftPair[1];
        final int clipRight = clipRightPair[0];
        final int rightHmerClip = clipRightPair[1];

        applyClipping(clipLeft, leftHmerClip, clipRight, rightHmerClip, false);

        setDirection(Direction.REFERENCE);

    }

    private boolean isBaseFormat() {
       return samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_TI) || samRecord.hasAttribute(FLOW_MATRIX_TAG_NAME);
    }

    private void fillFlowMatrix(final int [] kh, final int [] kf,
                                final double [] kdProbs ) {
        for ( int i = 0 ; i < kh.length; i++ ) {
            // normal matrix filling
            final int     pos = kf[i];
            final int     hmer = kh[i] & 0xff;
            if (hmer > maxHmer){
                flowMatrix[maxHmer][pos] = Math.max(flowMatrix[maxHmer][pos], kdProbs[i]);
            } else {
                flowMatrix[hmer][pos] = Math.max(flowMatrix[hmer][pos], kdProbs[i]);
            }
        }

    }

    // execute the matrix modifications
    private void implementMatrixMods(final int[] flowMatrixModsInstructions) {

        if ( flowMatrixModsInstructions != null ) {
            for (int hmer = 0; hmer < flowMatrixModsInstructions.length; hmer++) {
                final int hmer2 = flowMatrixModsInstructions[hmer];
                if (hmer2 != 0) {
                    for (int pos = 0; pos < flowMatrix[0].length; pos++) {

                        if (flowMatrix[hmer][pos] > flowMatrix[hmer2][pos]) {
                            flowMatrix[hmer2][pos] = flowMatrix[hmer][pos];
                        }

                        // if we are copying bacwards, zero out source
                        if (hmer > hmer2)
                            flowMatrix[hmer][pos] = 0;
                    }
                }
            }
        }
    }


    private void flipMatrix() {
        for ( int i = 0 ; i < flowMatrix.length; i++) {
            ArrayUtils.reverse(flowMatrix[i]);
        }
    }

    private static int findFirstNonZero(final int[] array){
        int result = -1;
        for (int i = 0 ; i < array.length; i++){
            if (array[i]!=0) {
                result = i;
                break;
            }
        }
        return result;
    }

    private static int findLastNonZero(final int[] array){
        int result = -1;
        for (int i = array.length-1 ; i >= 0; i--){
            if (array[i]!=0) {
                result = i;
                break;
            }
        }
        return result;
    }

    private static void shiftColumnUp(final double[][] matrix, final int colnum, final int shift) {
        for (int i = 0; i < matrix.length - shift; i ++ ) {
            matrix[i][colnum] = matrix[i+shift][colnum];
        }
        for (int i = matrix.length - shift; i < matrix.length; i ++ ) {
            matrix[i][colnum] = 0;
        }

    }

    public void setDirection(final Direction dir ) {
        direction = dir;
    }

    //trims base-spaced reads. Usually not needed, but kept for completeness
    public void applyBaseClipping(final int clipLeftBase, final int clipRightBase, boolean spread){
        final int[] clipLeftPair = findLeftClipping(clipLeftBase);
        final int[] clipRightPair = findRightClipping(clipRightBase);
        final int clipLeft = clipLeftPair[0];
        final int leftHmerClip = clipLeftPair[1];
        final int clipRight = clipRightPair[0];
        final int rightHmerClip = clipRightPair[1];
        if (getLength() - clipLeftBase - clipRightBase < minimalReadLength) {
            baseClipped = true;
            validKey = false;
            trimLeftBase = clipLeftBase;
            trimRightBase = clipRightBase;
        } else {
            applyClipping(clipLeft, leftHmerClip, clipRight, rightHmerClip, spread);
            baseClipped = true;
            trimLeftBase = clipLeftBase;
            trimRightBase = clipRightBase;
        }
    }

    private void applyClipping(int clipLeft, final int leftHmerClip, int clipRight, final int rightHmerClip, final boolean spread){
        if ((clipLeft < 0) || (clipRight < 0)  || (clipLeft >= getKeyLength() ) || ( clipRight >= getKeyLength())) {
            throw new IllegalStateException("Weird read clip calculated");
        }

        if ((leftHmerClip < 0) || (rightHmerClip < 0)  || (leftHmerClip >= 14 ) || ( rightHmerClip >= 14)) {
            throw new IllegalStateException("Weird read clip calculated");
        }

        final int originalLength = key.length;

        key[clipLeft]-=leftHmerClip;
        boolean shiftLeft = true;
        while (key[clipLeft] == 0) {
            clipLeft += 1 ;
            shiftLeft = false;
        }

        key[key.length - clipRight-1] -= rightHmerClip;
        boolean shiftRight = true;
        while (key[originalLength - 1- clipRight] == 0) {
            clipRight += 1 ;
            shiftRight = false;
        }

        key = Arrays.copyOfRange(key, clipLeft, originalLength - clipRight);
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);
        flowOrder = Arrays.copyOfRange(flowOrder, clipLeft, originalLength - clipRight);

        final double [][] newFlowMatrix = new double[flowMatrix.length][originalLength - clipLeft - clipRight] ;
        for ( int i = 0 ; i < newFlowMatrix.length; i++) {
            newFlowMatrix[i] = Arrays.copyOfRange(flowMatrix[i], clipLeft, originalLength - clipRight);
        }

        flowMatrix = newFlowMatrix;
        if (shiftLeft) {
            shiftColumnUp(flowMatrix, 0, leftHmerClip);
        }

        if (shiftRight) {
            shiftColumnUp(flowMatrix, flowMatrix[0].length-1, rightHmerClip);
        }

        //Spread boundary flow probabilities for the boundary hmers of the read
        //in this case the value of the genome hmer is uncertain
        if ( spread ) {
            spreadFlowLengthProbsAcrossCountsAtFlow(findFirstNonZero(key));
            spreadFlowLengthProbsAcrossCountsAtFlow(findLastNonZero(key));
        }
    }

    private int[] findLeftClippingFromCigar() {
        final List<CigarElement> cigar = getCigarElements();
        final int[] result = new int[2];
        if (cigar.size() == 0 ) {
            return result;
        }

        final CigarElement start = cigar.get(0);
        if (start.getOperator() != CigarOperator.H) {
            return result;
        }

        final int basesClipped = start.getLength();
        return findLeftClipping(basesClipped);
    }


    private int[] findLeftClipping(final int basesClipped){
        return FlowBasedReadUtils.findLeftClipping(basesClipped, flow2base, key);
    }

    private int[] findRightClippingFromCigar() {
        final List<CigarElement> cigar = getCigarElements();
        final int[] result = new int[2];
        if (cigar.size() == 0 ) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        final CigarElement end = cigar.get(cigar.size()-1);
        if (end.getOperator() != CigarOperator.H) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        final int basesClipped = end.getLength();

        return findRightClipping(basesClipped);
    }

    private int[] findRightClipping(final int basesClipped) {
        final int[] rkey = new int[key.length];
        for (int i = 0 ; i < key.length; i++ ){
            rkey[i] = key[key.length-1-i];
        }
        final int[] rflow2base = FlowBasedKeyCodec.getKeyToBase(rkey);

        return FlowBasedReadUtils.findRightClipping(basesClipped, rflow2base, rkey);
    }


    public void writeKey(final FileWriter oos)
            throws IOException {
        for (int i = 0; i < key.length; i++)
            oos.write(key[i]+"\n");

    }

    /**
     * Flow matrix logger
     *
     * This is used exclusively for testing
     *
     * @param oos
     * @throws IOException
     */
    protected void writeMatrix(final OutputStreamWriter oos)
            throws IOException {
        final DecimalFormat formatter = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

        final byte[]      bases = samRecord.getReadBases();
        int         basesOfs = 0;
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      ti = samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_TI) ? samRecord.getByteArrayAttribute(FLOW_MATRiX_OLD_TAG_TI) : (new byte[key.length*3]);
        if ( isReverseStrand() )
        {
            ArrayUtils.reverse(quals);
            ArrayUtils.reverse(ti);
        }

        for ( int col = 0 ; col < key.length ; col++ ) {
            oos.write("C,R,F,B,Bi,Q,ti\n");
            final byte base = (key[col] != 0) ? (basesOfs < bases.length ? bases[basesOfs] : (byte)'?') : (byte)'.';
            final String bi = (key[col] != 0) ? Integer.toString(basesOfs) : ".";
            final String q = (key[col] != 0) ? Integer.toString(quals[basesOfs]) : ".";
            final String Ti = (key[col] != 0) ? Integer.toString(ti[basesOfs]) : ".";
            for (int row = 0; row < flowMatrix.length; row++) {
                final String s = formatter.format(flowMatrix[row][col]);
                oos.write(""
                        + col + ","
                        + row + ","
                        + key[col] + ","
                        + (char)base + ","
                        + bi + ","
                        + q + ","
                        + Ti + ","
                        + (isReverseStrand() ? "r" : ".") + " "
                        + s + "\n");
            }
            if ( key[col] != 0 )
                basesOfs +=  key[col];
            oos.write("\n");
        }

    }


    public byte [] getFlowOrderArray() {
        return flowOrder;
    }

    public int getKeyLength() {
        return key.length;
    }

    public int[] getKey() {
        return key;
    }

    /**
     * Number of total bases that the flow based key generates
     * @return number of bases
     */
    public int totalKeyBases()  {
        int sum = 0 ;
        for (int i = 0 ; i < key.length; i++){
            sum += key[i];
        }
        return sum;
    }

    public int seqLength(){
        return forwardSequence.length;
    }
    public boolean isBaseClipped() {
        return baseClipped;
    }

    public int getTrimmedStart() {
        return trimLeftBase + getStart();
    }
    public int getTrimmedEnd() {
        return getEnd() - trimRightBase;
    }



    //functions that take care of simulating base format
    //they perform modifications on the flow matrix that are defined in applyFilteringFlowMatrix
    //this function was only applied when we tested what is the necessary information to be reported in the flow matrix
    private void applyFilteringFlowMatrix(){

        if (fbargs.disallowLargerProbs) {
            removeLargeProbs();
        }

        if (fbargs.removeLongerThanOneIndels) {
            removeLongIndels( key );
        }

        if (fbargs.removeOneToZeroProbs) {
            removeOneToZeroProbs(key);
        }

        if ((fbargs.lumpProbs)) {
            lumpProbs();
        }
        clipProbs();

        if (fbargs.symmetricIndels) {
            smoothIndels(key);
        }
        if (fbargs.onlyInsOrDel) {
            reportInsOrDel(key);
        }

        if ((fbargs.retainMaxNProbs)){
            reportMaxNProbsHmer(key);
        }

    }

    /**
     * clip probability values to fbargs.probabilityRatioThreshold
     */
    private void clipProbs() {
        for ( int i = 0 ; i < getMaxHmer(); i++ ) {
            for ( int j =0; j < getNFlows(); j++) {
                if ((flowMatrix[i][j] <= fbargs.probabilityRatioThreshold) &&
                        (key[j]!=i)) {
                    flowMatrix[i][j] = fbargs.fillingValue;
                }
            }
        }
    }

    /**
     * remove probabilities larger than 1
     */
    private void removeLargeProbs(){
        for (int i = 0; i < getNFlows(); i++){
            for (int j = 0 ; j < getMaxHmer()+1; j++) {
                if (flowMatrix[j][i] > 1) {
                    flowMatrix[j][i] = 1;
                }
            }
        }
    }

    /**
     * This is vestigial and applies only to old formats
     * @param key_kh
     */
    private void removeLongIndels(final int [] key_kh ){
        for ( int i = 0 ; i < getNFlows(); i++ ) {
            for (int j = 0; j < getMaxHmer()+1; j++){
                if (Math.abs(j-key_kh[i])>1){
                    flowMatrix[j][i] = fbargs.fillingValue;
                }
            }
        }
    }

    /**
     * This is vestigial and applies only to old formats
     * @param key_kh
     */
    private void removeOneToZeroProbs(final int [] key_kh) {
        for (int i = 0 ; i < getNFlows(); i++){
            if (key_kh[i] == 0){
                for (int j = 1; j < getMaxHmer()+1; j++){
                    flowMatrix[j][i]=fbargs.fillingValue;
                }
            }
        }
    }


    /**
     * Quantize probability values according to fbargs.probabilityQuantization and fbargs.probabilityScalingFactor
     * @param kd_probs
     */

    private void quantizeProbs(final int [] kd_probs ) {
        final int nQuants = fbargs.probabilityQuantization;
        final double bin_size = 6*fbargs.probabilityScalingFactor/(float)nQuants;
        for ( int i = 0 ; i < kd_probs.length; i++) {
            if (kd_probs[i] <=0)
                continue;
            else {
                kd_probs[i] = (byte)(bin_size * (byte)(kd_probs[i]/bin_size)+1);
            }
        }
    }

    /**
     * Smooth out probabilities by averaging with neighbours
     * @param kr
     */
    private void smoothIndels(final int [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            final int idx = kr[i];
            if (( idx > 1 ) && ( idx < maxHmer) ) {
                final double prob = (flowMatrix[idx - 1][i] + flowMatrix[idx + 1][i]) / 2;
                flowMatrix[idx - 1][i] = prob;
                flowMatrix[idx + 1][i] = prob;
            }
        }
    }

    /**
     * Only report probability of insertions or of deletions, never both
     * @param kr
     */
    private void reportInsOrDel(final int [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            final int idx = kr[i];
            if (( idx > 1 ) && ( idx < maxHmer) ) {
                if ((flowMatrix[idx-1][i] > fbargs.fillingValue) && (flowMatrix[idx+1][i] > fbargs.fillingValue)) {
                    final int fixCell = flowMatrix[idx-1][i] > flowMatrix[idx+1][i] ? idx+1 : idx-1;
                    flowMatrix[fixCell][i] = fbargs.fillingValue;
                }
            }
        }
    }

    /**
     * Combine all probabilities of insertions together and report them as probabilities of 1-mer insertion
     * Combine all probabilities of deletions together and report them as probabilities of 1-mer deletion
     */
    private void lumpProbs() {

        for (int i = 0; i < getMaxHmer(); i++) {
            for (int j = 0 ; j < getNFlows(); j ++ ) {
                final int fkey = key[j];
                if (flowMatrix[i][j]<=fbargs.fillingValue) {
                    continue;
                } else {
                    if ( (i - fkey) < -1 ){
                        flowMatrix[fkey-1][j]+=flowMatrix[i][j];
                        flowMatrix[i][j] = fbargs.fillingValue;
                    } else if ((i-fkey) > 1) {
                        flowMatrix[fkey+1][j]+=flowMatrix[i][j];
                        flowMatrix[i][j] = fbargs.fillingValue;
                    }

                }

            }
        }

    }

    /*
    Given full vector of error probabilities retain only the probabilities that can be reported in the base format
    (N+1/2 highest error probabilities)
     */
    private void reportMaxNProbsHmer(final int [] key) {
        final double [] tmpContainer = new double[maxHmer];
        for (int i = 0 ; i < key.length;i++){

            for (int j = 0 ; j < tmpContainer.length; j++) {
                tmpContainer[j] = flowMatrix[j][i];
            }
            final int k = (key[i]+1)/2;
            final double kth_highest = findKthLargest(tmpContainer, k+1);
            for (int j = 0 ; j < maxHmer; j++)
                if (flowMatrix[j][i] < kth_highest)
                    flowMatrix[j][i] = fbargs.fillingValue;
        }

    }


    private static double findKthLargest(final double[] nums, final int k) {
        final PriorityQueue<Double> q = new PriorityQueue<Double>(k);
        for(final double i: nums){
            q.offer(i);

            if(q.size()>k){
                q.poll();
            }
        }

        return q.peek();
    }

    private double[] phredToProb(final int [] kq) {
        final double [] result = new double[kq.length];
        for (int i = 0 ; i < kq.length; i++ ) {
            //disallow probabilities below fillingValue
            result[i] = Math.max(Math.pow(10, ((double)-kq[i])/fbargs.probabilityScalingFactor), fbargs.fillingValue);
        }
        return result;
    }

    /**
     * Hard clips uncertain flows (currently four first flows read that often generate noisy base calls
     * @param inputRead GATKREad
     * @param flowOrder flow order string from the SAM header
     * @param fbargs arguments
     * @return read with flowNumUncertainFlows trimmed
     */
    public static GATKRead hardClipUncertainBases(final GATKRead inputRead, final String flowOrder,
                                                  final FlowBasedArgumentCollection fbargs ){
        Utils.nonNull(fbargs);
        Utils.validateArg(fbargs.flowFirstUncertainFlowBase.length()==1, "First uncertain flow base should be of length 1");
        ReadClipper clipper = new ReadClipper(inputRead);
        final String adjustedFlowOrder = adjustFlowOrderToUncertainFlow(flowOrder, fbargs.flowFirstUncertainFlowBase.charAt(0), fbargs.flowOrderCycleLength);
        if (inputRead.isReverseStrand()) {
            final int nUncertain = nUncertainBases(inputRead, adjustedFlowOrder, fbargs.flowNumUncertainFlows, false);
            clipper.addOp(new ClippingOp(inputRead.getLength()-nUncertain, inputRead.getLength()-1));
        }
        else  {
            final int nUncertain = nUncertainBases(inputRead, adjustedFlowOrder, fbargs.flowNumUncertainFlows, true);
            clipper.addOp(new ClippingOp(0, nUncertain-1));
        }

        return clipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);

    }

    /**
     * Hard clips uncertain flows (currently four first flows read that often generate noisy base calls
     * @param inputRead GATKREad
     * @param samHeader  sam file header (to extract to flow order from
     * @param fbargs arguments
     * @return read with flowNumUncertainFlows trimmed
     */
    public static GATKRead hardClipUncertainBases(final GATKRead inputRead, final SAMFileHeader samHeader,
                                                  final FlowBasedArgumentCollection fbargs ){
        Utils.nonNull(fbargs);
        String flowOrder = samHeader.getReadGroup(inputRead.getReadGroup()).getFlowOrder();
        if (flowOrder==null){
            throw new GATKException("Unable to trim uncertain bases without flow order information");
        }
        flowOrder = flowOrder.substring(0,fbargs.flowOrderCycleLength);
        return hardClipUncertainBases(inputRead, flowOrder, fbargs);
    }

    private static String adjustFlowOrderToUncertainFlow(final String flowOrder,
                                                         final char firstUncertainFlowBase,
                                                         final int flowOrderLength){
        String result = flowOrder + flowOrder;
        final int adjustedStartPos = result.indexOf(firstUncertainFlowBase);
        return result.substring(adjustedStartPos, adjustedStartPos+flowOrderLength);
    }

    //find how many bases are output from uncertain flows
    private static int nUncertainBases(final GATKRead inputRead, final String flowOrder,
                                       final int nUncertainFlows, final boolean isForward){
        byte [] bases;
        if (isForward){
            bases = inputRead.getBases();
        } else {
            bases = ReadUtils.getBasesReverseComplement(inputRead).getBytes();
        }

        final int[] key = FlowBasedKeyCodec.baseArrayToKey(bases, flowOrder);
        final int nTrimFlows = Math.min(nUncertainFlows, key.length);
        int result = 0;
        for (int i = 0 ; i < nTrimFlows; i++){
            result += key[i];
        }
        return result;
    }

    public static void setMinimalReadLength(int minimalReadLength) {
        FlowBasedRead.minimalReadLength = minimalReadLength;
    }

    // convert qualities and ti tag to flow matrix
    private void readVestigialFlowMatrixFromTI(final String _flowOrder) {

        vestigialOneShotLogger.warn("Vestigial read format detected: " + samRecord);

        // generate key (base to flow space)
        setDirection(Direction.REFERENCE);  // base is always in reference/alignment direction
        key = FlowBasedKeyCodec.baseArrayToKey(samRecord.getReadBases(), _flowOrder);
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);
        flowOrder = FlowBasedKeyCodec.getFlowToBase(_flowOrder, key.length);

        // initialize matrix
        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = fbargs.fillingValue;;
            }
        }

        // access qual, convert to flow representation
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      ti = samRecord.getByteArrayAttribute(FLOW_MATRiX_OLD_TAG_TI);
        final double[]    probs = new double[quals.length];
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final double q = quals[i];
            final double p = QualityUtils.qualToErrorProb(q);
            probs[i] = p*2;
        }

        // apply key and qual/ti to matrix
        int     qualOfs = 0;
        for ( int i = 0 ; i < key.length ; i++ ) {
            final int        run = key[i];

            // the probability is not divided by two for hmers of length 1
            if ( run == 1 ) {
                probs[qualOfs] = probs[qualOfs]/2;
            }

            //filling the probability for the called hmer (not reported by the quals
            if ( run <= maxHmer ) {
                flowMatrix[run][i] = (run > 0) ? (1 - probs[qualOfs]) : 1;
                //require a prob. at least 0.1
                flowMatrix[run][i] = Math.max(MINIMAL_CALL_PROB, flowMatrix[run][i]);

            }

            if ( run != 0 ) {
                if ( quals[qualOfs] != 40 ) {
                    final int     run1 = (ti[qualOfs] == 0) ? (run - 1) : (run + 1);
                    if (( run1 <= maxHmer ) && (run <= maxHmer)){
                        flowMatrix[run1][i] = probs[qualOfs] / flowMatrix[run][i];
                    }
                    if (run <= maxHmer) {
                        flowMatrix[run][i] /= flowMatrix[run][i]; // for comparison to the flow space - probabilities are normalized by the key's probability
                    }
                }
                qualOfs += run;
            }

        }

        //this is just for tests of all kinds of
        applyFilteringFlowMatrix();
    }

    // code for reading BAM format where the flow matrix is stored in sparse representation in kr,kf,kh and kd tags
    // used for development of the new basecalling, but not in production code
    private void readVestigialFlowMatrixFromKR(final String _flowOrder) {

        vestigialOneShotLogger.warn("Vestigial read format detected: " + samRecord);

        key = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KR, true);

        // creates a translation from flow # to base #
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);

        // create a translation from
        flowOrder = FlowBasedKeyCodec.getFlowToBase(_flowOrder, key.length);

        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = fbargs.fillingValue;
            }
        }

        int [] kh = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KH, true);
        int [] kf = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KF, false);
        int [] kd = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KD, true);

        final int [] key_kh = key;
        final int [] key_kf = new int[key.length];
        for ( int i = 0 ; i < key_kf.length ; i++)
            key_kf[i] = i;
        final int [] key_kd = new int[key.length];

        kh = ArrayUtils.addAll(kh, key_kh);
        kf = ArrayUtils.addAll(kf, key_kf);
        kd = ArrayUtils.addAll(kd, key_kd);

        quantizeProbs(kd);

        final double [] kdProbs = phredToProb(kd);
        fillFlowMatrix( kh, kf, kdProbs);
        applyFilteringFlowMatrix();
        validateSequence();
    }

}


