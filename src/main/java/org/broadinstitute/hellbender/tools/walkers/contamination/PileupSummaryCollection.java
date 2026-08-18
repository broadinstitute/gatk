package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractSampleLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Simple data structure to pass and read/write a List of {@link PileupSummary} objects.
 * All {@link PileupSummary} fields must be specified if reading/writing from/to file.
 */
public final class PileupSummaryCollection extends AbstractSampleLocatableCollection<PileupSummary> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, POSITION, REF_COUNT, ALT_COUNT, OTHER_ALT_COUNT, REF_ALLELE, ALT_ALLELE, ALLELE_FREQUENCY
     */
    enum PileupSummaryTableColumn {
        CONTIG,
        POSITION,
        REF_COUNT,
        ALT_COUNT,
        OTHER_ALT_COUNT,
        REF_NUCLEOTIDE,
        ALT_NUCLEOTIDE,
        ALLELE_FREQUENCY;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, PileupSummary> PILEUP_SUMMARY_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(PileupSummaryTableColumn.CONTIG);
        final int position = dataLine.getInt(PileupSummaryTableColumn.POSITION);
        final int refReadCount = dataLine.getInt(PileupSummaryTableColumn.REF_COUNT);
        final int altReadCount = dataLine.getInt(PileupSummaryTableColumn.ALT_COUNT);
        final int otherAltReadCount = dataLine.getInt(PileupSummaryTableColumn.OTHER_ALT_COUNT);
        final byte refAllele = dataLine.get(PileupSummaryTableColumn.REF_NUCLEOTIDE).getBytes()[0];
        final byte altAllele = dataLine.get(PileupSummaryTableColumn.ALT_NUCLEOTIDE).getBytes()[0];
        final double alleleFrequency = dataLine.getDouble(PileupSummaryTableColumn.ALLELE_FREQUENCY);
        return new PileupSummary(contig, position, refReadCount, altReadCount, otherAltReadCount, refAllele, altAllele, alleleFrequency);
    };

    private static final BiConsumer<PileupSummary, DataLine> PILEUP_SUMMARY_RECORD_TO_DATA_LINE_ENCODER = (pileupSummary, dataLine) ->
            dataLine.append(pileupSummary.getContig())
                    .append(pileupSummary.getStart())
                    .append(pileupSummary.getRefCount())
                    .append(pileupSummary.getAltCount())
                    .append(pileupSummary.getOtherAltCount())
                    .append(Byte.toString(pileupSummary.getRefBase()))
                    .append(Byte.toString(pileupSummary.getAltBase()))
                    .append(pileupSummary.getAlleleFrequency());

    public PileupSummaryCollection(final File inputFile) {
        super(inputFile, PileupSummaryTableColumn.COLUMNS, PILEUP_SUMMARY_RECORD_FROM_DATA_LINE_DECODER, PILEUP_SUMMARY_RECORD_TO_DATA_LINE_ENCODER);
    }

    public PileupSummaryCollection(final SampleLocatableMetadata metadata,
                                  final List<PileupSummary> pileupSummaries) {
        super(metadata, pileupSummaries, PileupSummaryTableColumn.COLUMNS, PILEUP_SUMMARY_RECORD_FROM_DATA_LINE_DECODER, PILEUP_SUMMARY_RECORD_TO_DATA_LINE_ENCODER);
    }
}
