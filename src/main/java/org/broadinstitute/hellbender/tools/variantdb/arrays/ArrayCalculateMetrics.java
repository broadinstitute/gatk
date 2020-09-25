package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.tribble.util.popgen.HardyWeinbergCalculation;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.GenotypeCountsSchema;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
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

        TableReference tableRef = new TableReference(genotypeCountsTable, GenotypeCountsSchema.GENOTYPE_COUNTS_FIELDS);

        try (final StorageAPIAvroReader reader = new StorageAPIAvroReader(tableRef)) {
            for ( final GenericRecord row : reader ) {
                List<String> thisRow = new ArrayList<>();
                // data in row should never be null
                long probeId = (Long) row.get(GenotypeCountsSchema.PROBE_ID);
                thisRow.add(String.valueOf(probeId));

                long combined_hom_var = (Long) row.get(GenotypeCountsSchema.HOM_VAR_COUNT) +
                        (Long) row.get(GenotypeCountsSchema.HET_1_2_COUNT) +
                        (Long) row.get(GenotypeCountsSchema.HOM_VAR_COUNT);

                GenotypeCounts genotypeCounts = new GenotypeCounts((Long) row.get(GenotypeCountsSchema.HOM_REF_COUNT),
                        (Long) row.get(GenotypeCountsSchema.HET_COUNT), combined_hom_var);
                long noCalls = (Long) row.get(GenotypeCountsSchema.NO_CALL_COUNT);
                int sampleCount = (int) genotypeCounts.getRefs() + (int) genotypeCounts.getHets() + (int) genotypeCounts.getHoms() + (int) noCalls;
                double hwe;
                // If there's no data set the p-value to 1.
                if (genotypeCounts.getRefs() + genotypeCounts.getHets() + genotypeCounts.getHoms() <= 0) {
                    hwe = 1;
                } else {
                    hwe = HardyWeinbergCalculation.hwCalculate((int) genotypeCounts.getRefs(), (int) genotypeCounts.getHets(), (int) genotypeCounts.getHoms());
                }
                double phredHwe = hwe == 1 ? 0 : Math.floor(-10.0 * Math.log10(hwe));
                thisRow.add(String.format("%.0f", phredHwe));

                double callRate = 1.0 - ((double) noCalls / sampleCount);
                thisRow.add(String.format("%.3f", callRate));

                boolean invariant = genotypeCounts.getHets() + genotypeCounts.getHoms() == 0;
                thisRow.add(String.valueOf(invariant));

                metricsTsvWriter.getNewLineBuilder().setRow(thisRow);
                progressMeter.update(null);
            }
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
