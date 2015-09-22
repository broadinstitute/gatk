package org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr;

import org.broadinstitute.hellbender.utils.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;

import java.io.File;
import java.io.Serializable;

/**
 * Holds the recalibration information that's output from the first phase of BQSR.
 * It comes either directly
 * from the pipeline (without having to touch disk), or via a saved textual
 * report file.
 */
public final class BaseRecalOutput implements Serializable {
    private static final long serialVersionUID = 1L;
    private final RecalibrationTables recalibrationTables;
    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final StandardCovariateList covariates; // list of all covariates to be used in this calculation

    /**
     * Load the recal info from the report file.
     *
     * @param bqsrRecalFile local file with the textual recalibration information.
     */
    public BaseRecalOutput(File bqsrRecalFile) {
        final RecalibrationReport recalibrationReport = new RecalibrationReport(bqsrRecalFile);
        recalibrationTables = recalibrationReport.getRecalibrationTables();
        covariates = recalibrationReport.getCovariates();
        quantizationInfo = recalibrationReport.getQuantizationInfo();
    }

    /**
     * Set the recal info directly, e.g. after you just calculated it.
     */
    public BaseRecalOutput(RecalibrationTables recalibrationTables, QuantizationInfo quantizationInfo, StandardCovariateList covariates) {
        this.recalibrationTables = recalibrationTables;
        this.quantizationInfo = quantizationInfo;
        this.covariates = covariates;
    }

    public RecalibrationTables getRecalibrationTables() {
        return this.recalibrationTables;
    }

    public QuantizationInfo getQuantizationInfo() {
        return this.quantizationInfo;
    }

    public StandardCovariateList getCovariates() {
        return this.covariates;
    }


}
