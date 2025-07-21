package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class RepeatLengthCovariate implements Covariate {
    private static final long serialVersionUID = 1L;

    protected int MAX_REPEAT_LENGTH;
    protected int MAX_STR_UNIT_LENGTH;
    private final HashMap<String, Integer> repeatLookupTable = new HashMap<>();
    private final HashMap<Integer, String> repeatReverseLookupTable = new HashMap<>();
    private int nextId = 0;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC, final List<String> readGroups) {
        MAX_STR_UNIT_LENGTH = 8;
        MAX_REPEAT_LENGTH = 20;
    }

    @Override
    public void recordValues( final GATKRead read, final SAMFileHeader header, final ReadCovariates values, final boolean recordIndelValues) {
        // store the original bases and then write Ns over low quality ones
        final byte[] originalBases = Arrays.copyOf(read.getBases(), read.getBases().length);

        final boolean negativeStrand = read.isReverseStrand();
        byte[] bases = read.getBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        // don't record reads with N's
        if (!BaseUtils.isAllRegularBases(bases))
            return;

        for (int i = 0; i < bases.length; i++) {
            final Pair<byte[], Integer> res = findTandemRepeatUnits(bases, i);
            // to merge repeat unit and repeat length to get covariate value:
            final String repeatID =  getCovariateValueFromUnitAndLength(res.getLeft(),  res.getRight());
            final int key = keyForRepeat(repeatID);

            final int readOffset = (negativeStrand ? bases.length - i - 1 : i);
            values.addCovariate(key, key, key, readOffset);
        }

        // put the original bases back in
        read.setBases(originalBases);
    }

    public Pair<byte[], Integer> findTandemRepeatUnits(byte[] readBases, int offset) {
        int maxBW = 0;
        byte[] bestBWRepeatUnit = new byte[]{readBases[offset]};
        for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
            // fix repeat unit length
            //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
            if (offset+1-str < 0)
                break;

            // get backward repeat unit and # repeats
            byte[] backwardRepeatUnit = Arrays.copyOfRange(readBases, offset - str + 1, offset + 1);
            maxBW = GATKVariantContextUtils.findNumberOfRepetitions(backwardRepeatUnit, Arrays.copyOfRange(readBases, 0, offset + 1), false);
            if (maxBW > 1) {
                bestBWRepeatUnit = Arrays.copyOf(backwardRepeatUnit, backwardRepeatUnit.length);
                break;
            }
        }
        byte[] bestRepeatUnit = bestBWRepeatUnit;
        int maxRL = maxBW;

        if (offset < readBases.length-1) {
            byte[] bestFWRepeatUnit = new byte[]{readBases[offset+1]};
            int maxFW = 0;
            for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
                // fix repeat unit length
                //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
                if (offset+str+1 > readBases.length)
                    break;

                // get forward repeat unit and # repeats
                byte[] forwardRepeatUnit = Arrays.copyOfRange(readBases, offset + 1, offset + str + 1);
                maxFW = GATKVariantContextUtils.findNumberOfRepetitions(forwardRepeatUnit, Arrays.copyOfRange(readBases, offset + 1, readBases.length), true);
                if (maxFW > 1) {
                    bestFWRepeatUnit = Arrays.copyOf(forwardRepeatUnit, forwardRepeatUnit.length);
                    break;
                }
            }
            // if FW repeat unit = BW repeat unit it means we're in the middle of a tandem repeat - add FW and BW components
            if (Arrays.equals(bestFWRepeatUnit, bestBWRepeatUnit)) {
                maxRL = maxBW + maxFW;
                bestRepeatUnit = bestFWRepeatUnit; // arbitrary
            }
            else {
                // tandem repeat starting forward from current offset.
                // It could be the case that best BW unit was differnet from FW unit, but that BW still contains FW unit.
                // For example, TTCTT(C) CCC - at (C) place, best BW unit is (TTC)2, best FW unit is (C)3.
                // but correct representation at that place might be (C)4.
                // Hence, if the FW and BW units don't match, check if BW unit can still be a part of FW unit and add
                // representations to total
                maxBW = GATKVariantContextUtils.findNumberOfRepetitions(bestFWRepeatUnit, Arrays.copyOfRange(readBases, 0, offset + 1), false);
                maxRL = maxFW + maxBW;
                bestRepeatUnit = bestFWRepeatUnit;

            }

        }
        
        if(maxRL > MAX_REPEAT_LENGTH) { maxRL = MAX_REPEAT_LENGTH; }
        return new MutablePair<>(bestRepeatUnit, maxRL);

    }

    @Override
    public String formatKey(final int key) {
        return repeatReverseLookupTable.get(key);
    }

    protected String getCovariateValueFromUnitAndLength(final byte[] repeatFromUnitAndLength, final int repeatLength) {
        return String.format("%d", repeatLength);
    }


    @Override
    public int keyFromValue(final Object value) {
        return keyForRepeat((String) value);
    }

    private int keyForRepeat(final String repeatID) {
        if ( ! repeatLookupTable.containsKey(repeatID) ) {
            repeatLookupTable.put(repeatID, nextId);
            repeatReverseLookupTable.put(nextId, repeatID);
            nextId++;
        }
        return repeatLookupTable.get(repeatID);
    }


    /**
     * Splits repeat unit and num repetitions from covariate value.
     * For example, if value if "ATG4" it returns (ATG,4)
     * @param value             Covariate value
     * @return                  Split pair
     */
    public static Pair<String,Integer> getRUandNRfromCovariate(final String value) {

        int k = 0;
        for ( k=0; k < value.length(); k++ ) {
            if (!BaseUtils.isRegularBase(value.getBytes()[k]))
                break;
        }
        Integer nr = Integer.valueOf(value.substring(k, value.length())); // will throw NumberFormatException if format illegal
        if (k == value.length() || nr <= 0)
            throw new IllegalStateException("Covariate is not of form (Repeat Unit) + Integer");

        return new MutablePair<>(value.substring(0,k), nr);
    }

    @Override
    public int maximumKeyValue() {
        // max possible values of covariate: for repeat unit, length is up to MAX_STR_UNIT_LENGTH,
        // so we have 4^MAX_STR_UNIT_LENGTH * MAX_REPEAT_LENGTH possible values
        return (1+MAX_REPEAT_LENGTH);
    }



}