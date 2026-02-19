package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.hadoop.fs.Path;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.ingest.bq.SamplePloidyBQWriter;
import org.broadinstitute.hellbender.tools.gvs.ingest.parquet.SamplePloidyParquetWriter;

import java.io.File;
import java.io.IOException;
import java.util.Map;

public class SamplePloidyCreator {
    private static final Logger logger = LogManager.getLogger(SamplePloidyCreator.class);

    private SamplePloidyWriter samplePloidyWriter = null;
    private static final String SAMPLE_PLOIDY_TABLE_NAME = "sample_chromosome_ploidy";
    private static final String SAMPLE_CHROMOSOME_PLOIDY_FILETYPE_PREFIX = "sample_chromosome_ploidy";
    // When we detect mixed ploidy, we allow for a small fraction due to DRAGEN bugs.  This constant determines the cutoff
    // after which we'll throw an error instead of just normalizing the ploidy
    private static final double MIXED_PLOIDY_ERROR_CUTOFF = 0.05;
    private static final String PREFIX_SEPARATOR = "_";

    private final Long sampleId;

    private final CommonCode.OutputType outputType;

    public SamplePloidyCreator(String sampleIdentifierForOutputFileName, Long sampleId, String projectId, String datasetName, File outputDirectory, MessageType schema, final CommonCode.OutputType outputType) {
        try {
            this.sampleId = sampleId;
            this.outputType = outputType;

            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project id and dataset name when writing ploidy data to BigQuery.");
                    }
                    samplePloidyWriter = new SamplePloidyBQWriter(projectId, datasetName, SAMPLE_PLOIDY_TABLE_NAME);
                    break;
                case PARQUET:
                    if (outputDirectory == null || schema == null) {
                        throw new UserException("Must specify outputDirectory and schema when writing ploidy data to Parquet.");
                    }

                    String[] sampleComponents = {sampleId.toString(), sampleIdentifierForOutputFileName};
                    String filename = SAMPLE_CHROMOSOME_PLOIDY_FILETYPE_PREFIX + String.join(PREFIX_SEPARATOR, sampleComponents) +
                            "." + outputType.toString().toLowerCase();

                    final File sampleChromosomePloidyOutputFile = new File(outputDirectory, filename);

                    samplePloidyWriter = new SamplePloidyParquetWriter(new Path(sampleChromosomePloidyOutputFile.toURI()), schema, CompressionCodecName.SNAPPY);
                    break;
                default:
                    logger.warn("Not creating sample ploidy writer for unsupported output type " + outputType);
                    break;
            }
        } catch (Exception e) {
            throw new UserException("Could not create Sample Ploidy Table Writer", e);
        }
    }


    public void apply(Map<String, Map<Integer, Long>> ploidyData, long totalRefEntries) throws IOException {
        for (final Map.Entry<String, Map<Integer, Long>> ploidyLine : ploidyData.entrySet()) {
            Map<Integer, Long> ploidiesWithCounts = ploidyLine.getValue();
            // This is the happy path we'll normally follow--no mixed ploidy detected
            if (ploidiesWithCounts.size() == 1) {
                // we know there's only one item here, so we can just send that off
                samplePloidyWriter.write(this.sampleId, SchemaUtils.encodeLocation(ploidyLine.getKey(), 0), ploidiesWithCounts.keySet().iterator().next());
                continue;
            }

            int bestPloidy = -1;
            double highestPercentage = -1;
            long highestCount = 0L;

            int secondBestPloidy = -1;
            double secondHighestPercentage = -1;
            long secondHighestCount = 0L;
            // We can detect mixed ploidy for many reasons.  Some of which come down to bugs in the GATK code.
            for (Map.Entry<Integer, Long> ploidyEntryWithCounts : ploidiesWithCounts.entrySet()) {
                double percentage = ploidyEntryWithCounts.getValue() / totalRefEntries;
                // if necessary, update our best ploidy
                if (percentage > highestPercentage) {
                    secondBestPloidy = bestPloidy;
                    secondHighestPercentage = highestPercentage;
                    secondHighestCount = highestCount;

                    bestPloidy = ploidyEntryWithCounts.getKey();
                    highestPercentage = percentage;
                    highestCount = ploidyEntryWithCounts.getValue();
                } else if (percentage > secondHighestPercentage) {
                    secondBestPloidy = ploidyEntryWithCounts.getKey();
                    secondHighestPercentage = percentage;
                    secondHighestCount = ploidyEntryWithCounts.getValue();
                }
            }

            // Decide which ploidy to keep
            // First, see if the second best ploidy is for greater than 5% of the sample (this is likely way too generous).
            // If so, there may be a deeper error going on and we should just quit
            if (secondHighestPercentage > MIXED_PLOIDY_ERROR_CUTOFF) {
                throw new UserException("Detected mixed ploidy in sample "+this.sampleId+" on chromosome "+ploidyLine.getKey()+", with second ploidy of "+secondBestPloidy+" detected in "+(secondHighestPercentage * 100)+"% ("+secondHighestCount+" total) of samples");
            }
            // It's a small enough number to just note and move on with
            logger.warn("WARNING: Detected mixed ploidy in sample "+this.sampleId+" on chromosome "+ploidyLine.getKey()+", but second ploidy of "+secondBestPloidy+" detected in only "+(secondHighestPercentage * 100)+"% ("+secondHighestCount+" total)of samples. Going with dominant ploidy of "+bestPloidy);

            if (outputType == CommonCode.OutputType.BQ) {
                samplePloidyWriter.write(this.sampleId, SchemaUtils.encodeLocation(ploidyLine.getKey(),0), bestPloidy);
            }
        }
    }

    public void commitData() {
        this.samplePloidyWriter.commitData();
    }

    public void closeTool() {
        try {
            samplePloidyWriter.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
