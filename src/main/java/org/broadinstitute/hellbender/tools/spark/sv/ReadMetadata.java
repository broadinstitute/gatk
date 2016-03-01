package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A bag of data about reads:  contig name to id mapping, fragment length statistics by read group, mean length.
 */
@DefaultSerializer(ReadMetadata.Serializer.class)
public class ReadMetadata {
    private final Map<String, Short> contigNameToID;
    private final Map<String, ReadGroupFragmentStatistics> readGroupToFragmentStatistics;

    public ReadMetadata( final SAMFileHeader header,
                         final List<ReadGroupFragmentStatistics> statistics,
                         final ReadGroupFragmentStatistics noGroupStats ) {
        final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
        if ( contigs.size() > Short.MAX_VALUE ) throw new GATKException("Too many reference contigs.");
        contigNameToID = new HashMap<>(SVUtils.hashMapCapacity(contigs.size()));
        final int nContigs = contigs.size();
        for ( int contigID = 0; contigID < nContigs; ++contigID ) {
            contigNameToID.put(contigs.get(contigID).getSequenceName(), (short)contigID);
        }

        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if ( readGroups.size() != statistics.size() ) throw new GATKException("Wrong number of statistics for read groups.");
        readGroupToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(readGroups.size()));
        final int nReadGroups = readGroups.size();
        for ( int readGroupId = 0; readGroupId < nReadGroups; ++readGroupId ) {
            readGroupToFragmentStatistics.put(readGroups.get(readGroupId).getId(), statistics.get(readGroupId));
        }
        readGroupToFragmentStatistics.put(null, noGroupStats);
    }

    @SuppressWarnings("unchecked")
    private ReadMetadata( final Kryo kryo, final Input input ) {
        int contigMapSize = input.readInt();
        contigNameToID = new HashMap<>(SVUtils.hashMapCapacity(contigMapSize));
        while ( contigMapSize-- > 0 ) {
            final String contigName = kryo.readObject(input, String.class);
            final short contigId = input.readShort();
            contigNameToID.put(contigName, contigId);
        }
        int readGroupMapSize = input.readInt();
        readGroupToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(readGroupMapSize));
        while ( readGroupMapSize-- > 0 ) {
            final String readGroupName = kryo.readObjectOrNull(input, String.class);
            final ReadGroupFragmentStatistics groupStats = kryo.readObject(input, ReadGroupFragmentStatistics.class);
            readGroupToFragmentStatistics.put(readGroupName, groupStats);
        }
    }

    private void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(contigNameToID.size());
        for ( final Map.Entry<String, Short> entry : contigNameToID.entrySet() ) {
            kryo.writeObject(output, entry.getKey());
            output.writeShort(entry.getValue());
        }
        output.writeInt(readGroupToFragmentStatistics.size());
        for ( final Map.Entry<String, ReadGroupFragmentStatistics> entry : readGroupToFragmentStatistics.entrySet() ) {
            kryo.writeObjectOrNull(output, entry.getKey(), String.class);
            kryo.writeObject(output, entry.getValue());
        }
    }

    public short getContigID( final String contigName ) {
        final Short result = contigNameToID.get(contigName);
        if ( result == null ) throw new GATKException("No such contig name: "+contigName);
        return result;
    }

    public ReadGroupFragmentStatistics getStatistics( final String readGroupName ) {
        final ReadGroupFragmentStatistics stats = readGroupToFragmentStatistics.get(readGroupName);
        if ( stats == null ) throw new GATKException("No such read group name: "+readGroupName);
        return stats;
    }

    public int getMaxMedianFragmentSize() {
        return readGroupToFragmentStatistics.entrySet().stream()
                .mapToInt(entry -> Math.round(entry.getValue().getMedianFragmentSize()))
                .max()
                .orElse(0);
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( !(obj instanceof ReadMetadata) ) return false;
        final ReadMetadata that = (ReadMetadata)obj;
        return this.contigNameToID.equals(that.contigNameToID) &&
                this.readGroupToFragmentStatistics.equals(that.readGroupToFragmentStatistics);
    }

    @Override
    public int hashCode() {
        return 47*(47*contigNameToID.hashCode() + readGroupToFragmentStatistics.hashCode());
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ReadMetadata> {
        @Override
        public void write( final Kryo kryo, final Output output, final ReadMetadata readMetadata ) {
            readMetadata.serialize(kryo, output);
        }

        @Override
        public ReadMetadata read( final Kryo kryo, final Input input, final Class<ReadMetadata> klass ) {
            return new ReadMetadata(kryo, input);
        }
    }

    @DefaultSerializer(ReadGroupFragmentStatistics.Serializer.class)
    public static final class ReadGroupFragmentStatistics {
        private final float medianFragmentSize;
        private final float medianAbsoluteDeviationFragmentSize;

        public ReadGroupFragmentStatistics( final float medianFragmentSize,
                                            final float medianAbsoluteDeviationFragmentSize ) {
            this.medianFragmentSize = medianFragmentSize;
            this.medianAbsoluteDeviationFragmentSize = medianAbsoluteDeviationFragmentSize;
        }

        private ReadGroupFragmentStatistics( final Kryo kryo, final Input input ) {
            medianFragmentSize = input.readFloat();
            medianAbsoluteDeviationFragmentSize = input.readFloat();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeFloat(medianFragmentSize);
            output.writeFloat(medianAbsoluteDeviationFragmentSize);
        }

        public float getMedianFragmentSize() { return medianFragmentSize; }
        public float getMedianAbsoluteDeviationFragmentSize() { return medianAbsoluteDeviationFragmentSize; }

        @Override
        public boolean equals( final Object obj ) {
            if ( !(obj instanceof ReadGroupFragmentStatistics) ) return false;
            final ReadGroupFragmentStatistics that = (ReadGroupFragmentStatistics)obj;
            return this.medianFragmentSize == that.medianFragmentSize &&
                    this.medianAbsoluteDeviationFragmentSize == that.medianAbsoluteDeviationFragmentSize;
        }

        @Override
        public int hashCode() {
            return 47*(47*Float.hashCode(medianFragmentSize))+Float.hashCode(medianAbsoluteDeviationFragmentSize);
        }

        public static final class Serializer
                extends com.esotericsoftware.kryo.Serializer<ReadGroupFragmentStatistics> {
            @Override
            public void write( final Kryo kryo, final Output output,
                               final ReadGroupFragmentStatistics readGroupFragmentStatistics ) {
                readGroupFragmentStatistics.serialize(kryo, output);
            }

            @Override
            public ReadGroupFragmentStatistics read( final Kryo kryo, final Input input,
                                                     final Class<ReadGroupFragmentStatistics> klass ) {
                return new ReadGroupFragmentStatistics(kryo, input);
            }
        }
    }
}
