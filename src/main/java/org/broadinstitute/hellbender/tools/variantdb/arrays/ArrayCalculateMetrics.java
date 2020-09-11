package org.broadinstitute.hellbender.tools.variantdb.arrays;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.tools.walkers.annotator.ExcessHet;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;


@CommandLineProgramProperties(
        summary = "(\"CalculateMetrics\") - Calculates HWE and Call rate per site.",
        oneLineSummary = "Tool to calculate metrics from big query and upload results",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ArrayCalculateMetrics extends GATKTool {
    @Argument(
            fullName = "genotype-counts-table",
            doc = "Fully qualified name of the table where the genotype counts already exists"
    )
    private String genotypeCountsTable = null;

    @Argument(
            fullName = "output",
            doc = "TSV file that will be output with metrics per probe_id"
    )
    private GATKPath output = null;

    private SimpleXSVWriter metricsTsvWriter = null;

    public enum HeaderFieldEnum {
        probe_id,
        hwe_pval,
        call_rate,
        invariant
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        try {
            metricsTsvWriter = new SimpleXSVWriter(output.toPath(), IngestConstants.SEPARATOR);
        } catch (IOException e) {
            throw new UserException("Can't write to output file" + e);
        }
        metricsTsvWriter.setHeaderLine(Arrays.stream(HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList()));
    }

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(1000L);

        final String genotypeCountQueryString =
                "SELECT * FROM `" + genotypeCountsTable + "`";

        //Execute Query
        final TableResult result = BigQueryUtils.executeQuery(genotypeCountQueryString);

        for (final FieldValueList row : result.iterateAll()) {
            List<String> thisRow = new ArrayList<>();
            Long probeId = row.get(0).getLongValue();
            thisRow.add(String.valueOf(probeId));

            GenotypeCounts genotypeCounts = new GenotypeCounts(row.get(1).getDoubleValue(), row.get(2).getDoubleValue(), row.get(3).getDoubleValue());
            long noCalls = row.get(4).getLongValue();
            Integer sampleCount = (int) genotypeCounts.getRefs() + (int) genotypeCounts.getHets() + (int) genotypeCounts.getHoms() + (int) noCalls;
            Double excessHetPval = ExcessHet.calculateEH(genotypeCounts, sampleCount).getRight();
            thisRow.add(String.valueOf(excessHetPval));

            Double callRate = 1.0 - ((double) noCalls / sampleCount);
            thisRow.add(String.valueOf(callRate));

            Boolean invariant = genotypeCounts.getHets() + genotypeCounts.getHoms() == 0;
            thisRow.add(String.valueOf(invariant));

            metricsTsvWriter.getNewLineBuilder().setRow(thisRow);
            progressMeter.update(null);
        }
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();
        try {
            metricsTsvWriter.close();
        } catch (IOException e) {
            throw new UserException("Can't close output file" + e);
        }
    }
}
