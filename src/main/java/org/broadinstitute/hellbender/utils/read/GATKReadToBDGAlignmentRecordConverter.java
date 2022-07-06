package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.bdgenomics.adam.converters.SAMRecordConverter;
import org.bdgenomics.adam.models.SequenceDictionary;
import org.bdgenomics.adam.models.ReadGroupDictionary;

/**
 * Converts a GATKRead to a BDG AlignmentRecord
 */
public class GATKReadToBDGAlignmentRecordConverter {
    private static final SAMRecordConverter converter = new SAMRecordConverter();

    private SAMFileHeader header;
    private SequenceDictionary dict;
    private ReadGroupDictionary readGroups;

//    public GATKReadToBDGAlignmentRecordConverter(SAMFileHeader header) {
//        this.header = header;
//        this.dict = SequenceDictionary.fromSAMSequenceDictionary(header.getSequenceDictionary());
//        this.readGroups = ReadGroupDictionary.fromSAMHeader(header);
//    }

    public static AlignmentRecord convert( final GATKRead gatkRead, final SAMFileHeader header ) {
//        SequenceDictionary dict = SequenceDictionary.fromSAMSequenceDictionary(header.getSequenceDictionary());
        //ReadGroupDictionary readGroups = ReadGroupDictionary.fromSAMHeader(header);
        return GATKReadToBDGAlignmentRecordConverter.convert(gatkRead, header, null, null);
    }

    public static AlignmentRecord convert(
            final GATKRead gatkRead, final SAMFileHeader header, final SequenceDictionary dict, final ReadGroupDictionary readGroups ) {
        //return converter.convert(gatkRead.convertToSAMRecord(header));
        return null;
    }

    public static AlignmentRecord convert(
            final SAMRecord sam, final SequenceDictionary dict, final ReadGroupDictionary readGroups ) {
        //return converter.convert(sam);
        return null;
    }
}
