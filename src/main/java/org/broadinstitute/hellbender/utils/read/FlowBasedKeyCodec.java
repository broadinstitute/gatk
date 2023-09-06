package org.broadinstitute.hellbender.utils.read;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;

import java.util.ArrayList;

public class FlowBasedKeyCodec {

    /**
     * Converts base space sequence to flow space
     * @param bases base space sequence
     * @param flowOrder flow order
     * @return Array of flow values
     */

    static public int[] baseArrayToKey(final byte[] bases, final String flowOrder){

        final ArrayList<Integer> result = new ArrayList<>();
        final byte[] flowOrderBytes = flowOrder.getBytes();
        int loc = 0;
        int flowNumber = 0 ;
        final int period = flowOrderBytes.length;
        int periodGuard = 0;
        while ( loc < bases.length ) {
            final byte flowBase = flowOrderBytes[flowNumber%period];
            if ((bases[loc]!=flowBase) && ( bases[loc]!= BaseUtils.Base.N.base)) {
                result.add(0);
                if ( ++periodGuard > period )
                    throw new GATKException("baseArrayToKey periodGuard tripped, on " + new String(bases) + ", flowOrder: " + flowOrder
                    + " This probably indicates the presence of a base (value) in the sequence that is not included in the provided flow order");
            } else {
                int count = 0;
                while ( ( loc < bases.length) && ((bases[loc]==flowBase) || (bases[loc]== BaseUtils.Base.N.base)) ){
                    loc++;
                    count ++;
                }
                result.add(count);
                periodGuard = 0;
            }
            flowNumber++;
        }
        final int[] ret = new int[result.size()];
        for (int i = 0; i < result.size(); i++) {
            ret[i] = result.get(i);
        }
        return ret;
    }

    /**
     * For every flow of the key output the index of the last base that was output prior to this flow
     * @param key given key
     * @return array
     */
    static public int[] getKeyToBase(final int[] key) {
        final int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1; i < result.length; i++) {
            result[i] = result[i - 1] + key[i - 1];
        }
        return result;

    }

    /**
     * For every flow of the key output the nucleotide that is being read for this flow
     * @param flowOrder given flow order
     * @param expectedLength the length of the key (key is not provided)
     * @return array of bases
     */

    static public byte[] getFlowToBase(final String flowOrder, final int expectedLength) {
        final byte[] result = new byte[expectedLength] ;
        for ( int i = 0; i < result.length; i++ ) {
            result[i] = (byte)flowOrder.charAt(i%flowOrder.length());
        }
        return result;
    }

    /**
     * Prints the key as character-encoded string
     * @param ints (key array)
     * @return encoded string
     */
    static public String keyAsString(final int[] ints)
    {
        final StringBuilder   sb = new StringBuilder();

        for ( final int i : ints )
            sb.append((char)((i < 10) ? ('0' + i) : ('A' + i - 10)));

        return sb.toString();
    }

    /**
     * Converts a numerical array into flow space by flattening per-base elements to fill empty flows based on minimum scores.
     * This is intended to make it easier to transform per-base read statistics that are computed on the read in base-space
     * sanely into flow-space reads in a reasonable fashion that makes sense for parameters like indel-likelihoods.
     *
     * The rules are as follows:
     *  - For every run of homopolymers, the minimum score for any given base is translated to cover the whole run
     *  - For every 0 base flow, the score from the previous filled flow is copied into place.
     *  - For the beginning of the array (in the case of preceding 0-flows) the given default value is used.
     *
     * Example:
     *  A read with the following bases with a base-space calculated into flow space (ACTG) as:
     *     TTTATGC -> 0030101101
     *  With an input array like this:
     *     byte[]{1,2,3,4,5,6,7}
     *  We should expect the following result array (Which matches the key array in length)
     *     byte[]{defaultQual,defaultQual,1,1,4,4,5,6,6,7}
     *
     * @param bases base space sequence
     * @param keyLength size of the converted bases array in flow-space
     * @param baseSpacedArrayToConvert Array of base-space scores to be conformed to flow space
     * @param defaultQual default quality to use at the head of the array for non-covered flows
     * @param flowOrder flow order
     * @return Array of translated bases comparable to the keyLength
     */
    @VisibleForTesting
    static public byte[] baseArrayToKeySpace(final byte[] bases, final int keyLength, final byte[] baseSpacedArrayToConvert, final byte defaultQual, final String flowOrder){

        if (bases.length != baseSpacedArrayToConvert.length) {
            throw new IllegalArgumentException("Read and qual arrays do not match");
        }

        final byte[] result = new byte[keyLength];
        final byte[] flowOrderBytes = flowOrder.getBytes();
        int loc = 0;
        int flowNumber = 0 ;
        byte lastQual = defaultQual;
        final int period = flowOrderBytes.length;
        while ( loc < bases.length ) {
            final byte flowBase = flowOrderBytes[flowNumber%period];
            if ((bases[loc]!=flowBase) && ( bases[loc]!= BaseUtils.Base.N.base)) {
                result[flowNumber] = lastQual;
            } else {
                byte qual = Byte.MAX_VALUE;
                while ( ( loc < bases.length) && ((bases[loc]==flowBase) || (bases[loc]== BaseUtils.Base.N.base)) ){
                    qual = (byte)Math.min(baseSpacedArrayToConvert[loc], qual);
                    loc++;
                }
                result[flowNumber] = qual;
                lastQual = qual;

            }
            flowNumber++;
        }
        return result;
    }
}
