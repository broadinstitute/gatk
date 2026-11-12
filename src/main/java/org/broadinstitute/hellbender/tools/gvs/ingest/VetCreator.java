package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.NormalizedAllele;
import org.broadinstitute.hellbender.tools.gvs.common.RefGetAccessionProvider;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.broadinstitute.hellbender.tools.gvs.common.VrsIdComputer;
import org.broadinstitute.hellbender.tools.gvs.common.VrsNormalizer;
import org.broadinstitute.hellbender.tools.gvs.ingest.bq.VetBQWriter;
import org.broadinstitute.hellbender.tools.gvs.ingest.parquet.VetParquetFileWriter;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class VetCreator {
    private VetWriter vetWriter = null;
    private VrsAlleleCreator vrsAlleleCreator = null;
    private final Long sampleId;
    private final boolean enableVrsIds;
    private final VrsNormalizer vrsNormalizer = new VrsNormalizer();

    private static final String VET_FILETYPE_PREFIX = "vet_";
    private static final String PREFIX_SEPARATOR = "_";

    /**
     * Written to {@code vrs_allele_ids} for a spanning-deletion {@code *} ALT, which VRS cannot represent.
     * A placeholder rather than nothing, so {@code vrs_allele_ids} stays aligned entry for entry with
     * {@code vet.alt} (BigQuery arrays cannot hold NULL).
     */
    public static final String NO_VRS_ALLELE_ID = "";

    public static boolean doRowsExistFor(CommonCode.OutputType outputType, String projectId, String datasetName, String tableNumber, Long sampleId) {
        if (outputType != CommonCode.OutputType.BQ) return false;
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, VET_FILETYPE_PREFIX + tableNumber, SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
    }

    public VetCreator(String inputVcfFileName, Long sampleId, String tableNumber, final File outputDirectory, final CommonCode.OutputType outputType, final String projectId, final String datasetName, final boolean forceLoadingFromNonAlleleSpecific, final boolean skipLoadingVqsrFields, final MessageType parquetSchema, final boolean enableVrsIds, final VrsAlleleCreator vrsAlleleCreator) {
        this.sampleId = sampleId;
        this.enableVrsIds = enableVrsIds;
        this.vrsAlleleCreator = vrsAlleleCreator;

        try {
            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    vetWriter = new VetBQWriter(projectId, datasetName, VET_FILETYPE_PREFIX + tableNumber, skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific);
                    break;
                case PARQUET:
                    String filename = getOutputFileName(tableNumber, sampleId, inputVcfFileName, outputType);
                    final File parquetOutputFile = new File(outputDirectory, filename);
                    vetWriter = new VetParquetFileWriter(new Path(parquetOutputFile.toURI()), skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific, parquetSchema, CompressionCodecName.SNAPPY);
                    break;
            }
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void apply(VariantContext variant, org.broadinstitute.hellbender.engine.ReferenceContext referenceContext) throws IOException {
        final int start = variant.getStart();
        final long location = SchemaUtils.encodeLocation(variant.getContig(), start);

        final List<VrsAlleleRecord> vrsRecords = computeVrsAlleleRecords(variant, referenceContext);

        if (vetWriter != null) {
            vetWriter.write(location, variant, sampleId, vrsAlleleIds(vrsRecords));
        }
        if (vrsAlleleCreator != null) {
            vrsAlleleCreator.apply(vrsRecords.stream()
                    .filter(java.util.Objects::nonNull)
                    .collect(java.util.stream.Collectors.toList()));
        }
    }

    /**
     * The {@code vrs_allele_ids} value for a vet row: one entry per record, {@link #NO_VRS_ALLELE_ID} where
     * the record is null (a spanning-deletion ALT).
     */
    static List<String> vrsAlleleIds(final List<VrsAlleleRecord> vrsRecords) {
        return vrsRecords.stream()
                .map(r -> r == null ? NO_VRS_ALLELE_ID : r.vrsAlleleId)
                .collect(java.util.stream.Collectors.toList());
    }

    /**
     * Normalize each non-ref alternate allele and compute its full VRS record.
     * The result is aligned entry for entry with the vet {@code alt} column, which drops only
     * {@code <NON_REF>}: a spanning-deletion {@code *} ALT gets a null entry, since VRS cannot represent it.
     * Returns an empty list when VRS IDs are disabled, reference is unavailable,
     * or the contig is not a recognized chromosome.
     */
    List<VrsAlleleRecord> computeVrsAlleleRecords(final VariantContext variant, final org.broadinstitute.hellbender.engine.ReferenceContext referenceContext) {
        if (!enableVrsIds || referenceContext == null || !referenceContext.hasBackingDataSource()
                || !RefGetAccessionProvider.isKnownContig(variant.getContig())) {
            return Collections.emptyList();
        }
        final String refgetAccession = RefGetAccessionProvider.getAccession(variant.getContig());
        final List<VrsAlleleRecord> records = new ArrayList<>();
        for (final Allele altAllele : variant.getAlternateAlleles()) {
            if (altAllele.isNonRefAllele()) {
                continue; // skip <NON_REF> and <*>, as the vet alt column does
            }
            if (altAllele.equals(Allele.SPAN_DEL)) {
                records.add(null); // keeps the list aligned with vet alt, which keeps *
                continue;
            }
            final long vrsStart = (long) variant.getStart() - 1;
            final long vrsEnd   = vrsStart + variant.getReference().length();
            final NormalizedAllele normalized = vrsNormalizer.normalize(
                variant.getContig(), vrsStart, vrsEnd, altAllele.getBaseString(), referenceContext
            );
            final String locationId = VrsIdComputer.computeLocationId(
                refgetAccession, normalized.start, normalized.end);
            final VrsAlleleRecord record;
            if (normalized.state.type == NormalizedAllele.SequenceState.StateType.REFERENCE_LENGTH_EXPRESSION) {
                final String alleleId = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(
                    refgetAccession, normalized.start, normalized.end,
                    normalized.state.length, normalized.state.repeatSubunitLength);
                record = new VrsAlleleRecord(alleleId, locationId, refgetAccession,
                    normalized.start, normalized.end,
                    normalized.state.length, normalized.state.repeatSubunitLength);
            } else {
                final String alleleId = VrsIdComputer.computeAlleleId(
                    refgetAccession, normalized.start, normalized.end, normalized.state.sequence);
                record = new VrsAlleleRecord(alleleId, locationId, refgetAccession,
                    normalized.start, normalized.end, normalized.state.sequence);
            }
            records.add(record);
        }
        return records;
    }


    public static String getOutputFileName(String tableNumber, Long sampleId, String inputVcfFileName, CommonCode.OutputType outputType) {
        String[] sampleComponents = {tableNumber, sampleId.toString(), inputVcfFileName};
        return VET_FILETYPE_PREFIX + String.join(PREFIX_SEPARATOR, sampleComponents) +
                "." + outputType.toString().toLowerCase();
    }

    public void commitData() {
        if (vetWriter != null) {
            vetWriter.commitData();
        }
        if (vrsAlleleCreator != null) {
            vrsAlleleCreator.commitData();
        }
    }

    public void closeTool() {
        try {
            if (vetWriter != null) vetWriter.close();
            if (vrsAlleleCreator != null) vrsAlleleCreator.closeTool();
        } catch (final Exception e) {
            throw new IllegalArgumentException(e);
        }
    }
}
