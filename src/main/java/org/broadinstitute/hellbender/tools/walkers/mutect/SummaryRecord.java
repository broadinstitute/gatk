package org.broadinstitute.hellbender.tools.walkers.mutect;

/**
 * Created by tsato on 2/8/17.
 */

public class SummaryRecord {
    public static final String TRUE_POSITIVE_COLUMN_NAME = "true-positive";
    public static final String FALSE_POSITIVE_COLUMN_NAME = "false-positive";
    public static final String FALSE_NEGATIVE_COLUMN_NAME = "false-negative";
    public static final String SENSITIVITY_COLUMN_NAME = "sensitivity";
    public static final String PRECISION_COLUMN_NAME = "precision";
    public static final String[] SUMMARY_TABLE_COLUMN_HEADER =
            {TRUE_POSITIVE_COLUMN_NAME, FALSE_POSITIVE_COLUMN_NAME, FALSE_NEGATIVE_COLUMN_NAME, SENSITIVITY_COLUMN_NAME, PRECISION_COLUMN_NAME};

    final long truePositives;
    final long falsePositives;
    final long falseNegatives;

    public SummaryRecord(final long truePositives, final long falsePositives, final long falseNegatives){
        this.truePositives = truePositives;
        this.falsePositives = falsePositives;
        this.falseNegatives = falseNegatives;
    }

    public long getTruePositives(){
        return truePositives;
    }

    public long getFalsePositives() {
        return falsePositives;
    }

    public long getFalseNegatives(){
        return falseNegatives;
    }

    public double getSensitivity(){
        return (double) truePositives / (truePositives + falseNegatives);
    }

    public double getPrecision(){
        return (double) truePositives / (truePositives + falsePositives);
    }
}

