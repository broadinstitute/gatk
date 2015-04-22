package org.broadinstitute.hellbender.tools.recalibration.covariates;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.NGSPlatform;
import org.broadinstitute.hellbender.utils.SequencerFlowClass;

/**
 * The Cycle covariate.
 * For Solexa the cycle is simply the position in the read (counting backwards if it is a negative strand read)
 * For 454 the cycle is the TACG flow cycle, that is, each flow grabs all the TACG's in order in a single cycle
 * For example, for the read: AAACCCCGAAATTTTTACTG
 * the cycle would be 11111111222333333344
 * For SOLiD the cycle is a more complicated mixture of ligation cycle and primer round
 */

public final class CycleCovariate implements Covariate {

    private int MAXIMUM_CYCLE_VALUE;
    public static final int CUSHION_FOR_INDELS = 4;
    private String default_platform = null;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        this.MAXIMUM_CYCLE_VALUE = RAC.MAXIMUM_CYCLE_VALUE;

        if (RAC.DEFAULT_PLATFORM != null && !NGSPlatform.isKnown(RAC.DEFAULT_PLATFORM))
            throw new UserException.CommandLineException("The requested default platform (" + RAC.DEFAULT_PLATFORM + ") is not a recognized platform.");

        if (RAC.DEFAULT_PLATFORM != null)
            default_platform = RAC.DEFAULT_PLATFORM;
    }

    // Used to pick out the covariate's value from attributes of the read
    @Override
    public void recordValues(final SAMRecord read, final ReadCovariates values) {
        final int readLength = read.getReadLength();
        final NGSPlatform ngsPlatform = default_platform == null ? NGSPlatform.fromRead(read) : NGSPlatform.fromReadGroupPL(default_platform);

        // Discrete cycle platforms
        if (ngsPlatform.getSequencerType() == SequencerFlowClass.DISCRETE) {
            final int readOrderFactor = read.getReadPairedFlag() && read.getSecondOfPairFlag() ? -1 : 1;
            final int increment;
            int cycle;
            if (read.getReadNegativeStrandFlag()) {
                cycle = readLength * readOrderFactor;
                increment = -1 * readOrderFactor;
            }
            else {
                cycle = readOrderFactor;
                increment = readOrderFactor;
            }

            final int MAX_CYCLE_FOR_INDELS = readLength - CUSHION_FOR_INDELS - 1;
            for (int i = 0; i < readLength; i++) {
                final int substitutionKey = keyFromCycle(cycle);
                final int indelKey = (i < CUSHION_FOR_INDELS || i > MAX_CYCLE_FOR_INDELS) ? -1 : substitutionKey;
                values.addCovariate(substitutionKey, indelKey, indelKey, i);
                cycle += increment;
            }
        }

        // Flow cycle platforms
        else if (ngsPlatform.getSequencerType() == SequencerFlowClass.FLOW) {
            throw new UserException("The platform (" + read.getReadGroup().getPlatform()
                    + ") associated with read group " + read.getReadGroup()
                    + " is not a supported platform.");
        }
        // Unknown platforms
        else {
            throw new UserException("The platform (" + read.getReadGroup().getPlatform()
                    + ") associated with read group " + read.getReadGroup()
                    + " is not a recognized platform. Allowable options are " + NGSPlatform.knownPlatformsString());
        }
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return Integer.parseInt(str);
    }

    @Override
    public String formatKey(final int key) {
        int cycle = key >> 1; // shift so we can remove the "sign" bit
        if ( (key & 1) != 0 ) // is the last bit set?
            cycle *= -1; // then the cycle is negative
        return String.format("%d", cycle);
    }

    @Override
    public int keyFromValue(final Object value) {
        return (value instanceof String) ? keyFromCycle(Integer.parseInt((String) value)) : keyFromCycle((Integer) value);
    }

    @Override
    public int maximumKeyValue() {
        return (MAXIMUM_CYCLE_VALUE << 1) + 1;
    }

    private int keyFromCycle(final int cycle) {
        // no negative values because values must fit into the first few bits of the long
        int result = Math.abs(cycle);
        if ( result > MAXIMUM_CYCLE_VALUE )
            throw new UserException("The maximum allowed value for the cycle is " + MAXIMUM_CYCLE_VALUE + ", but a larger cycle (" + result + ") was detected.  Please use the --maximum_cycle_value argument to increase this value (at the expense of requiring more memory to run)");

        result = result << 1; // shift so we can add the "sign" bit
        if ( cycle < 0 )
            result++; // negative cycles get the lower-most bit set
        return result;
    }
}