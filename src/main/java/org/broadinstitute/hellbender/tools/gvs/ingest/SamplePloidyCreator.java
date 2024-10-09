package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.util.BitSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;

public class SamplePloidyCreator {
    private static final Logger logger = LogManager.getLogger(SamplePloidyCreator.class);

    private PendingBQWriter samplePloidyBQJsonWriter = null;
    private static final String SAMPLE_PLOIDY_TABLE_NAME = "sample_chromosome_ploidy";

    private final Long sampleId;
    private final String projectId;
    private final String datasetName;


    public SamplePloidyCreator(Long sampleId, String projectId, String datasetName) {
        try {
            this.sampleId = sampleId;
            this.projectId = projectId;
            this.datasetName = datasetName;

            if (projectId == null || datasetName == null) {
                throw new UserException("Must specify project-id and dataset-name.");
            }
            samplePloidyBQJsonWriter = new PendingBQWriter(projectId, datasetName, SAMPLE_PLOIDY_TABLE_NAME);
        } catch (Exception e) {
            throw new UserException("Could not create VCF Header Scratch Table Writer", e);
        }
    }


    public void apply(Map<String, Map<Integer, Long>> ploidyData, long totalRefEntries) throws IOException {
        for (final Map.Entry<String, Map<Integer, Long>> ploidyLine : ploidyData.entrySet()) {
            try {
                Map<Integer, Long> ploidiesWithCounts = ploidyLine.getValue();
                // This is the happy path we'll normally follow--no mixed ploidy detected
                if (ploidiesWithCounts.size() == 1) {
                    // we know there's only one item here, so we can just sent that off
                    samplePloidyBQJsonWriter.addJsonRow(createJson(this.sampleId, SchemaUtils.encodeLocation(ploidyLine.getKey(),0), ploidiesWithCounts.keySet().iterator().next()));
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
                if (secondHighestPercentage > 0.05) {
                    throw new UserException("Detected mixed ploidy in sample "+this.sampleId+" on chromosome "+ploidyLine.getKey()+", with second ploidy of "+secondBestPloidy+" detected in "+(secondHighestPercentage * 100)+"% ("+secondHighestCount+" total) of samples");
                }
                // It's a small enough number to just note and move on with
                logger.warn("WARNING: Detected mixed ploidy in sample "+this.sampleId+" on chromosome "+ploidyLine.getKey()+", but second ploidy of "+secondBestPloidy+" detected in only "+(secondHighestPercentage * 100)+"% ("+secondHighestCount+" total)of samples. Going with dominant ploidy of "+bestPloidy);

                samplePloidyBQJsonWriter.addJsonRow(createJson(this.sampleId, SchemaUtils.encodeLocation(ploidyLine.getKey(),0), bestPloidy));
            } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
                throw new IOException("BQ exception", ex);
            }
        }
    }

    public JSONObject createJson(Long sampleId, Long chromosome, Integer ploidy) {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        record.put("chromosome", chromosome);
        record.put("ploidy", ploidy);
        return record;
    }


    public void commitData() {
        if (samplePloidyBQJsonWriter != null) {
            samplePloidyBQJsonWriter.flushBuffer();
            samplePloidyBQJsonWriter.commitWriteStreams();
        }
    }

    public void closeTool() {
        if (samplePloidyBQJsonWriter != null) {
            try {
                samplePloidyBQJsonWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close Sample Ploidy Table writer", e);
            }
        }
        if (samplePloidyBQJsonWriter != null) {
            samplePloidyBQJsonWriter.close();
        }
    }

}
