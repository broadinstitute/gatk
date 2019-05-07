package org.broadinstitute.hellbender.tools.sv.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.sv.EventCopyNumberPosteriors;
import org.broadinstitute.hellbender.tools.sv.records.StructuralVariantDepthQuality;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

public class StructuralVariantDepthQualityCollection extends AbstractLocatableCollection<LocatableMetadata, StructuralVariantDepthQuality> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END
     */
    enum StructuralVariantDepthQualityTableColumn {
        CONTIG,
        START,
        END,
        CNVID,
        SAMPLE_IDS,
        TYPE,
        P;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final String SAMPLE_ID_DELIIMETER = ",";

    private static final Function<DataLine, StructuralVariantDepthQuality> SITE_QUALITY_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(StructuralVariantDepthQualityTableColumn.CONTIG);
        final int start = dataLine.getInt(StructuralVariantDepthQualityTableColumn.START);
        final int end = dataLine.getInt(StructuralVariantDepthQualityTableColumn.END);
        final String cnvid = dataLine.get(StructuralVariantDepthQualityTableColumn.CNVID);
        final String[] samples = dataLine.get(StructuralVariantDepthQualityTableColumn.SAMPLE_IDS).split(StructuralVariantDepthQualityCollection.SAMPLE_ID_DELIIMETER);
        final String type = dataLine.get(StructuralVariantDepthQualityTableColumn.TYPE);
        final double p = dataLine.getDouble(StructuralVariantDepthQualityTableColumn.P);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        final EventCopyNumberPosteriors.EventRecord eventRecord = new EventCopyNumberPosteriors.EventRecord(cnvid, type, interval, samples);
        return new StructuralVariantDepthQuality(eventRecord, p);
    };

    private static final BiConsumer<StructuralVariantDepthQuality, DataLine> SITE_QUALITY_RECORD_TO_DATA_LINE_ENCODER = (StructuralVariantDepthQuality, dataLine) ->
            dataLine.append(StructuralVariantDepthQuality.getContig())
                    .append(StructuralVariantDepthQuality.getStart())
                    .append(StructuralVariantDepthQuality.getEnd())
                    .append(StructuralVariantDepthQuality.getEventRecord().getId())
                    .append(String.join(SAMPLE_ID_DELIIMETER, StructuralVariantDepthQuality.getEventRecord().getSamples()))
                    .append(StructuralVariantDepthQuality.getEventRecord().getType())
                    .append(StructuralVariantDepthQuality.getQuality());

    public StructuralVariantDepthQualityCollection(final File inputFile) {
        super(inputFile, StructuralVariantDepthQualityTableColumn.COLUMNS, SITE_QUALITY_RECORD_FROM_DATA_LINE_DECODER, SITE_QUALITY_RECORD_TO_DATA_LINE_ENCODER);
    }

    public StructuralVariantDepthQualityCollection(final LocatableMetadata metadata,
                                                   final List<StructuralVariantDepthQuality> records) {
        super(metadata, records, StructuralVariantDepthQualityTableColumn.COLUMNS, SITE_QUALITY_RECORD_FROM_DATA_LINE_DECODER, SITE_QUALITY_RECORD_TO_DATA_LINE_ENCODER);
    }
}