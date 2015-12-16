package org.broadinstitute.hellbender.utils.recalibration;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Class that encapsulates the information necessary for quality score quantization for BQSR
 */
public final class QuantizationInfo implements Serializable {
    private static final long serialVersionUID = 1L;
    private List<Byte> quantizedQuals;
    private final List<Long> empiricalQualCounts;
    private int quantizationLevels;

    private QuantizationInfo(List<Byte> quantizedQuals, List<Long> empiricalQualCounts, int quantizationLevels) {
        this.quantizedQuals = quantizedQuals;
        this.empiricalQualCounts = empiricalQualCounts;
        this.quantizationLevels = quantizationLevels;
    }

    public QuantizationInfo(List<Byte> quantizedQuals, List<Long> empiricalQualCounts) {
        this(quantizedQuals, empiricalQualCounts, calculateQuantizationLevels(quantizedQuals));
    }
    
    public QuantizationInfo(final RecalibrationTables recalibrationTables, final int quantizationLevels) {
        final Long[] qualHistogram = new Long[QualityUtils.MAX_SAM_QUAL_SCORE +1]; // create a histogram with the empirical quality distribution
        for (int i = 0; i < qualHistogram.length; i++)
            qualHistogram[i] = 0L;

        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getQualityScoreTable(); // get the quality score table

        for (final RecalDatum value : qualTable.getAllValues()) {
            final RecalDatum datum = value;
            final int empiricalQual = MathUtils.fastRound(datum.getEmpiricalQuality()); // convert the empirical quality to an integer ( it is already capped by MAX_QUAL )
            qualHistogram[empiricalQual] += datum.getNumObservations(); // add the number of observations for every key
        }
        empiricalQualCounts = new ArrayList<>(Arrays.asList(qualHistogram)); // histogram with the number of observations of the empirical qualities
        quantizeQualityScores(quantizationLevels);

        this.quantizationLevels = quantizationLevels;
    }


    public void quantizeQualityScores(int nLevels) {
        QualQuantizer quantizer = new QualQuantizer(empiricalQualCounts, nLevels, QualityUtils.MIN_USABLE_Q_SCORE); // quantize the qualities to the desired number of levels
        quantizedQuals = quantizer.getOriginalToQuantizedMap(); // map with the original to quantized qual map (using the standard number of levels in the RAC)
    }

    public void noQuantization() {
        this.quantizationLevels = QualityUtils.MAX_SAM_QUAL_SCORE;
        for (int i = 0; i < this.quantizationLevels; i++)
            quantizedQuals.set(i, (byte) i);
    }

    public List<Byte> getQuantizedQuals() {
        return quantizedQuals;
    }

    public int getQuantizationLevels() {
        return quantizationLevels;
    }

    public GATKReportTable generateReportTable() {
        GATKReportTable quantizedTable;
        quantizedTable = new GATKReportTable(RecalUtils.QUANTIZED_REPORT_TABLE_TITLE, "Quality quantization map", 3, GATKReportTable.Sorting.SORT_BY_COLUMN);
        quantizedTable.addColumn(RecalUtils.QUALITY_SCORE_COLUMN_NAME, "%d");
        quantizedTable.addColumn(RecalUtils.QUANTIZED_COUNT_COLUMN_NAME, "%d");
        quantizedTable.addColumn(RecalUtils.QUANTIZED_VALUE_COLUMN_NAME, "%d");

        for (int qual = 0; qual <= QualityUtils.MAX_SAM_QUAL_SCORE; qual++) {
            quantizedTable.set(qual, RecalUtils.QUALITY_SCORE_COLUMN_NAME, qual);
            quantizedTable.set(qual, RecalUtils.QUANTIZED_COUNT_COLUMN_NAME, empiricalQualCounts.get(qual));
            quantizedTable.set(qual, RecalUtils.QUANTIZED_VALUE_COLUMN_NAME, quantizedQuals.get(qual));
        }
        return quantizedTable;
    }

    private static int calculateQuantizationLevels(List<Byte> quantizedQuals) {
        byte lastByte = -1;
        int quantizationLevels = 0;
        for (byte q : quantizedQuals) {
            if (q != lastByte) {
                quantizationLevels++;
                lastByte = q;
            }
        }
        return quantizationLevels;
    }
}
