package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.utils.SimpleInterval;

public class CalledModeledSegment extends ModeledSegment {

    public enum Call {
        CNLOH("C"),
        AMPLIFICATION("+"),
        DELETION("-"),
        NEUTRAL("0");

        private final String outputString;

        Call(final String outputString) {
            this.outputString = outputString;
        }

        public String getOutputString() {
            return outputString;
        }
    }

    private final CalledModeledSegment.Call call;

    public CalledModeledSegment(final SimpleInterval interval, final int numPointsCopyRatio, final int numPointsAlleleFraction,
                                final SimplePosteriorSummary log2CopyRatioSimplePosteriorSummary, final SimplePosteriorSummary minorAlleleFractionSimplePosteriorSummary,
                                final CalledModeledSegment.Call call) {
        super(interval, numPointsCopyRatio, numPointsAlleleFraction, log2CopyRatioSimplePosteriorSummary, minorAlleleFractionSimplePosteriorSummary);
        this.call = call;
    }
}
