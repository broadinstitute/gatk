package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.bdgenomics.adam.converters.AlignmentRecordConverter;
import org.bdgenomics.adam.models.RecordGroupDictionary;
import org.bdgenomics.adam.models.SAMFileHeaderWritable;
import org.bdgenomics.formats.avro.AlignmentRecord;

/**
 * Implementation of the {@link GATKRead} interface for the {@link AlignmentRecord} class.
 *
 * The AlignmentRecord class is the Avro-based data model used in the Big Data Genomics project (aka ADAM).  It is also
 * designed to capture the main functionality of a BAM file.  The schema can be found in the bdg-formats repo:
 * https://github.com/bigdatagenomics/bdg-formats
 *
 * Utilities for converting to/from a SAMRecord can be found in the adam-core Maven artifact, built from the main ADAM
 * repo:
 * https://github.com/bigdatagenomics/adam
 *
 * Because it is an Avro model, it admits an efficient on-wire serialized format and is also compatible with the Parquet
 * columnar storage format.
 *
 * This adapter wraps an {@link AlignmentRecord} by first converting to a {@link SAMRecord} and then passing through
 * to {@link SAMRecordToGATKReadAdapter}.
 */
public final class BDGAlignmentRecordToGATKReadAdapter extends SAMRecordToGATKReadAdapter {
    private static final long serialVersionUID = 1L;

    private final AlignmentRecord alignmentRecord;

    public BDGAlignmentRecordToGATKReadAdapter(final AlignmentRecord alignmentRecord, final SAMFileHeader header) {
        super(new AlignmentRecordConverter().convert(alignmentRecord, SAMFileHeaderWritable.apply(header),
                RecordGroupDictionary.fromSAMHeader(header)));
        this.alignmentRecord = alignmentRecord;
    }

    public static GATKRead sparkReadAdapter(final AlignmentRecord record, final SAMFileHeader header) {
        return new BDGAlignmentRecordToGATKReadAdapter(record, header);
    }

    public static GATKRead sparkReadAdapter(final AlignmentRecord record) {
        return new BDGAlignmentRecordToGATKReadAdapter(record, null);
    }

    public AlignmentRecord convertToBDGAlignmentRecord() { return alignmentRecord; }

    @Override
    public String toString() {
        return commonToString();
    }
}
