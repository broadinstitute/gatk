package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.bdgenomics.adam.converters.AlignmentRecordConverter;
import org.bdgenomics.adam.models.SAMFileHeaderWritable;
import org.bdgenomics.formats.avro.AlignmentRecord;

import java.util.Random;
import java.util.UUID;
import java.util.concurrent.atomic.AtomicLong;

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
    private final static long uuidHighWord = new Random().nextLong();
    private final static AtomicLong uuidLowWord = new AtomicLong(0);

    private final AlignmentRecord alignmentRecord;

    public BDGAlignmentRecordToGATKReadAdapter( final AlignmentRecord alignmentRecord) {
        this(alignmentRecord, null);
    }

    public BDGAlignmentRecordToGATKReadAdapter( final AlignmentRecord alignmentRecord, final SAMFileHeader header) {
        // this is 100x faster than UUID.randomUUID()
        this(alignmentRecord, header, new UUID(uuidHighWord, uuidLowWord.incrementAndGet()));
    }

    /**
     * Constructor that allows an explicit UUID to be passed in -- only meant
     * for internal use and test class use, which is why it's package protected.
     */
    BDGAlignmentRecordToGATKReadAdapter( final AlignmentRecord alignmentRecord, final SAMFileHeader header, final UUID uuid ) {
        super(new AlignmentRecordConverter().convert(alignmentRecord, SAMFileHeaderWritable.apply(header)), uuid);
        this.alignmentRecord = alignmentRecord;
    }

    /**
     * Produces a BDGAlignmentRecordToGATKReadAdapter with a 0L,0L UUID. Spark doesn't need the UUIDs
     * and loading the reads twice (which can happen when caching is missing) prevents joining.
     * @param record Read to adapt
     * @param header SAMFileHeaderWritable corresponding to the underlying SAMRecord object
     * @return adapted Read
     */
    public static GATKRead sparkReadAdapter(final AlignmentRecord record, final SAMFileHeader header) {
        return new BDGAlignmentRecordToGATKReadAdapter(record, header, new UUID(0L, 0L));
    }

    public static GATKRead sparkReadAdapter(final AlignmentRecord record) {
        return new BDGAlignmentRecordToGATKReadAdapter(record, null, new UUID(0L, 0L));
    }

    public AlignmentRecord convertToBDGAlignmentRecord() { return alignmentRecord; }
}
