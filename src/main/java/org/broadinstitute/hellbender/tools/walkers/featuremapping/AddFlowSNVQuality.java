package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.concurrent.atomic.AtomicReference;

public final class AddFlowSNVQuality  {

    public static final String BASE_QUALITY_ATTRIBUTE_NAME = "BQ";
    public static final char PHRED_ASCII_BASE = '!';

    public static final int ERROR_PROB_BAND_1LESS = 0;
    public static final int ERROR_PROB_BAND_KEY = 1;
    public static final int ERROR_PROB_BAND_1MORE = 2;
    public static final int ERROR_PROB_BANDS = 3;

    public double minErrorRate = 1e-6;
    public double maxErrorRate = 1 - minErrorRate;
    public int maxQualityScore = 60;
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    public void addBaseQuality(final GATKRead read, final SAMFileHeader hdr) {

        // convert to a flow base read
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(hdr, read);
        final FlowBasedRead   fbRead = new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, fbargs);
        final int flowOrderLength = calcFlowOrderLength(rgInfo.flowOrder);

        // generate base quality
        final AtomicReference<double[][]> snvResultRef = new AtomicReference<>();
        final double[]      baseErrorProb = generateBaseErrorProbability(fbRead, flowOrderLength, rgInfo.flowOrder.getBytes(), snvResultRef);

        // install in read
        read.setAttribute(BASE_QUALITY_ATTRIBUTE_NAME, new String(convertErrorProbToPhred(baseErrorProb, true)));
        for ( int i = 0 ; i < flowOrderLength ; i++ ) {
            final String name = ApplySNVQR.attrNameForNonCalledBase(rgInfo.flowOrder.charAt(i));
            read.setAttribute(name, new String(convertErrorProbToPhred(snvResultRef.get()[i], true)));
        }
    }

    private byte[] convertErrorProbToPhred(double[] errorProb, boolean asciiBase) {

        final byte floor = asciiBase ? (byte)PHRED_ASCII_BASE : 0;

        final byte[] phred = new byte[errorProb.length];
        for ( int i = 0 ; i < errorProb.length ; i++ ) {
            if ( errorProb[i] == 0 ) {
                phred[i] = (byte)(maxQualityScore + floor);
            } else {
                phred[i] = (byte)(Math.min((maxQualityScore), (int) (-10 * Math.log10(errorProb[i]))) + floor);
            }
        }
        return phred;
    }

    private int calcFlowOrderLength(String flowOrder) {

        final int i = flowOrder.indexOf(flowOrder.charAt(0), 1);

        return (i < 0) ? flowOrder.length() : i;
    }

    private double[] generateBaseErrorProbability(final FlowBasedRead fbRead, final int flowOrderLength, byte[] flowOrder, AtomicReference<double[][]> snvResultRef) {

        /**
         * access key and error probabilities
         * for a description of the flow probabilities see {@link FlowBasedRead#flowMatrix}
         */
        final int[]       key = fbRead.getKey();
        final double[][]  errorProbBands = extractErrorProbBands(fbRead, minErrorRate);
        final double[]    result = new double[fbRead.getBasesNoCopy().length];

        final double[][]  snvResult = new double[flowOrderLength][];
        for ( int i = 0 ; i < snvResult.length ; i++ ) {
            snvResult[i] = new double[result.length];
        }
        snvResultRef.set(snvResult);

        // loop over hmers via flow key
        int               base = 0;
        Map<Byte, Double> allBaseProb0 = new LinkedHashMap<>();
        Map<Byte, Double> allBaseProb1 = new LinkedHashMap<>();

        for ( int flow = 0 ; flow < key.length ; flow++ ) {
            if ( key[flow] != 0 ) {

                // establish initial stat
                allBaseProb0.clear();
                allBaseProb1.clear();
                int flow_i = (flow % flowOrderLength);

                // establish hmer quality
                final int         hmerLength = key[flow];
                final double[]    hmerBaseErrorProbs = generateHmerBaseErrorProbabilities(key, errorProbBands, flow, flowOrderLength, flowOrder, allBaseProb0, allBaseProb1);

                // fill in

                // first base of the read gets the original error probability of the hmer, otherwise use computed error probability
                if ( base == 0 ) {
                    result[base++] = errorProbBands[ERROR_PROB_BAND_KEY][flow];
                } else {
                    result[base++] = hmerBaseErrorProbs[0];  // first base, or only base in case of a single base hmer
                    for ( int i = 0 ; i < flowOrderLength ; i++ ) {
                        if ( allBaseProb0.containsKey(flowOrder[i]) ) {
                            snvResult[i][base - 1] = allBaseProb0.get(flowOrder[i]);
                        } else if ( i != flow_i ) {
                            snvResult[i][base - 1] = minErrorRate;
                        }
                    }
                }

                // for hmers longer than 1
                if ( hmerLength > 1 ) {

                    // skip all but last (leave with zero error probability)
                    base += (hmerLength - 2);

                    // fill last base from computed error probability
                    result[base++] = hmerBaseErrorProbs[1]; // last base, if hmer is longer than 1

                    for ( int i = 0 ; i < flowOrderLength ; i++ ) {
                        if ( allBaseProb1.containsKey(flowOrder[i]) ) {
                            final double p = allBaseProb1.get(flowOrder[i]);
                            for ( int j = 0 ; j < hmerLength - 1 ; j++ ) {
                                snvResult[i][base - 1 - j] = p;
                            }
                        } else if ( i != flow_i ) {
                            for ( int j = 0 ; j < hmerLength - 1 ; j++ ) {
                                snvResult[i][base - 1 - j] = minErrorRate;
                            }
                        }
                    }
                }

                // override result for the last base with the orignal hmer error probability
                if ( base == result.length ) {
                    result[base - 1] = errorProbBands[ERROR_PROB_BAND_KEY][flow];
                }
            }
        }

        // adjust probability of called bases so that sum will be 1
        final byte[] bases = fbRead.getBasesNoCopy();
        for ( int ofs = 0 ; ofs < bases.length ; ofs++ ) {

            // go through alt bases and accumulate p, find out index of called bases (in flow order)
            final byte calledBase = bases[ofs];
            double altP = 0;
            int calledIndex = -1;
            for (int i = 0; i < flowOrderLength; i++) {
                if ( calledBase != flowOrder[i] ) {
                    altP += snvResult[i][ofs];
                } else {
                    calledIndex = i;
                }
            }

            // install probability in called base slot
            if ( calledIndex >= 0 ) {
                snvResult[calledIndex][ofs] = Math.max(0, 1 - altP);
            }
        }

        return result;
    }

    // extract error probability bands. middle (1) band is the key prob.
    // lower (0) and high (2) are corresponding to -1 and +1 in hmer lengths
    private static double[][] extractErrorProbBands(final FlowBasedRead flowRead, final double minValue) {

        // access key
        final int[] key = flowRead.getKey();

        // allocate result
        double[][] result = new double[ERROR_PROB_BANDS][];
        for ( int i = 0 ; i < result.length ; i++ ) {
            result[i] = new double[key.length];
        }

        for ( int i = 0 ; i < key.length ; i++ ) {

            // extract key probability
            result[ERROR_PROB_BAND_KEY][i] = Math.max(flowRead.getProb(i, key[i]), minValue);

            // -1
            if ( key[i] > 0 ) {
                result[ERROR_PROB_BAND_1LESS][i] = Math.max(flowRead.getProb(i, key[i] - 1), minValue);
            } else {
                result[ERROR_PROB_BAND_1LESS][i] = minValue;
            }

            // +1
            if ( key[i] < flowRead.getMaxHmer() ) {
                result[ERROR_PROB_BAND_1MORE][i] = Math.max(flowRead.getProb(i, key[i] + 1), minValue);
            } else {
                result[ERROR_PROB_BAND_1MORE][i] = minValue;
            }
        }

        return result;
    }

    @VisibleForTesting
    protected static double[] generateHmerBaseErrorProbabilities(final int[] key, final double[][] errorProbBands, final int flow,
                                                                 final int flowOrderLength, byte[] flowOrder,
                                                                 Map<Byte, Double> allBaseProb0, Map<Byte, Double> allBaseProb1) {

        // result is left/right error probabilities
        final double[]          errorProbs = new double[2];
        final int               hmerLength = key[flow];

        errorProbs[0] = generateSidedHmerBaseErrorProbability(key, errorProbBands, flow, -1, flowOrderLength, flowOrder, allBaseProb0);
        if ( hmerLength != 1 ) {
            errorProbs[1] = generateSidedHmerBaseErrorProbability(key, errorProbBands, flow, 1, flowOrderLength, flowOrder, allBaseProb1);
        }

        return errorProbs;
    }

    private static double generateSidedHmerBaseErrorProbability(final int[] key, final double[][] errorProbBands, final int flow, final int sideIncr,
                                                                final int flowOrderLength, final byte[] flowOrder, final Map<Byte, Double> allBaseProb) {

        // create a key slice of the area around the flow/hmer.
        final int minIndex = Math.max(flow - flowOrderLength + 1, 0);
        final int maxIndex = Math.min(flow + flowOrderLength - 1, key.length - 1);
        final int[] slice = Arrays.copyOfRange(key, minIndex, maxIndex + 1);
        final int hmerLength = key[flow];

        // walk the flows towards the side until (and including) the first non-zero key
        // on hmers of length 1 we walk both sides
        final byte            flowByte = flowOrder[flow % flowOrderLength];
        final List<int[]>     slices = new LinkedList<>();
        final List<Byte>      slicesBytes = new LinkedList<>();
        final int[]           incrs = (hmerLength != 1)
                                            ? new int[] { sideIncr }
                                            : new int[] { sideIncr, -sideIncr};
        for (int incr : incrs) {
            for (int sideFlow = flow + incr; sideFlow >= 0 && sideFlow < key.length; sideFlow += incr) {

                // create a alternative key slice by incrementing sideFlow and decrementing flow
                final int[] altSlice = Arrays.copyOf(slice, slice.length);
                altSlice[sideFlow - minIndex] += 1;
                altSlice[flow - minIndex] -= 1;
                if (sliceIsValid(altSlice, flowOrderLength)) {
                    slices.add(altSlice);
                    slicesBytes.add(flowOrder[sideFlow % flowOrderLength]);
                }

                // is the sideFlow non-zero? if so, break
                if (key[sideFlow] != 0) {
                    break;
                }
            }
        }

        // at this point, we have a list of valid slices. figure out the error probability for each of them
        // and compute the base quality
        final double keyP = sliceProb(slice, minIndex, key, errorProbBands);
        double sumP = keyP;
        int sliceIndex = 0;
        for ( int[] s : slices ) {
            final double sliceP = sliceProb(s, minIndex, key, errorProbBands);
            final byte sliceByte = slicesBytes.get(sliceIndex);
            if ( allBaseProb != null ) {
                allBaseProb.put(sliceByte, sliceP);
            }
            sumP += sliceP;
            sliceIndex++;
        }
        final double ep = 1 - (keyP / sumP);

        return ep;
    }

    // compute probability for a slice
    private static double sliceProb(int[] slice, int minIndex, int[] key, double[][] errorProbBands) {

        double p = 1.0;
        for ( int i = 0 ; i < slice.length ; i++ ) {
            final int band;
            if ( slice[i] < key[i + minIndex] ) {
                band = ERROR_PROB_BAND_1LESS;
            } else if ( slice[i] > key[i + minIndex] ) {
                band = ERROR_PROB_BAND_1MORE;
            } else {
                band = ERROR_PROB_BAND_KEY;
            }
            p *= errorProbBands[band][i + minIndex];
        }
        return p;
    }

    private static boolean sliceIsValid(final int[] slice, final int flowOrderLength) {

        // look for strings of consecutive zeros in length of flowOrderLength - 1
        int     consecutiveZeros = 0;
        for ( int key : slice ) {
            if ( key != 0 ) {
                consecutiveZeros = 0;
            } else {
                consecutiveZeros++;
                if ( consecutiveZeros >= (flowOrderLength - 1) ) {
                    return false;
                }
            }
        }

        // if here, not found -> valid
        return true;
    }
}
