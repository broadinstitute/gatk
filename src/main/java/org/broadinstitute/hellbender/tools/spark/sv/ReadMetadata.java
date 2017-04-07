package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.*;

/**
 * A bag of data about reads:  contig name to id mapping, fragment length statistics by read group, mean length.
 */
@DefaultSerializer(ReadMetadata.Serializer.class)
public class ReadMetadata {
    private final Map<String, Integer> contigNameToID;
    private final long nReads;
    private final long maxReadsInPartition;
    private final int nPartitions;
    private final int coverage;
    private final Map<String, ReadGroupFragmentStatistics> readGroupToFragmentStatistics;
    private final static String NO_GROUP = "NoGroup";

    public ReadMetadata( final SAMFileHeader header,
                         final JavaRDD<GATKRead> reads ) {
        contigNameToID = buildContigNameToIDMap(header);

        final int nReadGroups = header.getReadGroups().size();
        final List<PartitionStatistics> perPartitionStatistics =
                reads.mapPartitions(readItr ->
                        Collections.singletonList(new PartitionStatistics(readItr, nReadGroups)).iterator())
                     .collect();
        nPartitions = perPartitionStatistics.size();
        nReads = perPartitionStatistics.stream().mapToLong(PartitionStatistics::getNReads).sum();
        maxReadsInPartition = perPartitionStatistics.stream().mapToLong(PartitionStatistics::getNReads).max().orElse(0L);
        final long nReadBases = perPartitionStatistics.stream().mapToLong(PartitionStatistics::getNBases).sum();
        final long nRefBases = header.getSequenceDictionary().getSequences()
                .stream().mapToLong(SAMSequenceRecord::getSequenceLength).sum();
        coverage = (int) ((nReadBases + nRefBases - 1) / nRefBases);
        readGroupToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(header.getReadGroups().size()));
        final Map<String, long[]> combinedMaps =
                perPartitionStatistics.stream()
                        .map(PartitionStatistics::getReadGroupToFragmentSizeCountMap)
                        .reduce(new HashMap<>(SVUtils.hashMapCapacity(nReadGroups)), ReadMetadata::combineMaps);
        for ( final Map.Entry<String, long[]> entry : combinedMaps.entrySet() ) {
            readGroupToFragmentStatistics.put(entry.getKey(),
                    new ReadGroupFragmentStatistics(entry.getValue()));
        }
    }

    @VisibleForTesting
    ReadMetadata( final SAMFileHeader header, final ReadGroupFragmentStatistics stats,
                  final int nPartitions, final long nReads, final long maxReadsInPartition, final int coverage ) {
        contigNameToID = buildContigNameToIDMap(header);
        this.nPartitions = nPartitions;
        this.nReads = nReads;
        this.maxReadsInPartition = maxReadsInPartition;
        this.coverage = coverage;
        readGroupToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(header.getReadGroups().size() + 1));
        readGroupToFragmentStatistics.put(null, stats);
        for ( final SAMReadGroupRecord readGroupRecord : header.getReadGroups() ) {
            readGroupToFragmentStatistics.put(readGroupRecord.getReadGroupId(), stats);
        }
    }

    private ReadMetadata( final Kryo kryo, final Input input ) {
        int contigMapSize = input.readInt();
        contigNameToID = new HashMap<>(SVUtils.hashMapCapacity(contigMapSize));
        while ( contigMapSize-- > 0 ) {
            final String contigName = kryo.readObject(input, String.class);
            final int contigId = input.readInt();
            contigNameToID.put(contigName, contigId);
        }

        nReads = input.readLong();
        maxReadsInPartition = input.readLong();
        nPartitions = input.readInt();
        coverage = input.readInt();

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
        for ( final Map.Entry<String, Integer> entry : contigNameToID.entrySet() ) {
            kryo.writeObject(output, entry.getKey());
            output.writeInt(entry.getValue());
        }

        output.writeLong(nReads);
        output.writeLong(maxReadsInPartition);
        output.writeInt(nPartitions);
        output.writeInt(coverage);

        output.writeInt(readGroupToFragmentStatistics.size());
        for ( final Map.Entry<String, ReadGroupFragmentStatistics> entry : readGroupToFragmentStatistics.entrySet() ) {
            kryo.writeObjectOrNull(output, entry.getKey(), String.class);
            kryo.writeObject(output, entry.getValue());
        }
    }

    public Map<String, Integer> getContigNameMap() {
        return Collections.unmodifiableMap(contigNameToID);
    }

    public int getContigID( final String contigName ) {
        final Integer result = contigNameToID.get(contigName);
        if ( result == null ) throw new GATKException("No such contig name: " + contigName);
        return result;
    }

    public long getNReads() {
        return nReads;
    }

    public int getNPartitions() {
        return nPartitions;
    }

    public long getMaxReadsInPartition() {
        return maxReadsInPartition;
    }

    public int getCoverage() {
        return coverage;
    }

    public Map<String, ReadGroupFragmentStatistics> getAllGroupStatistics() {
        return readGroupToFragmentStatistics;
    }

    public ReadGroupFragmentStatistics getStatistics( final String readGroupName ) {
        final ReadGroupFragmentStatistics stats = readGroupToFragmentStatistics.get(readGroupName);
        if ( stats == null ) throw new GATKException("No such read group name: " + readGroupName);
        return stats;
    }

    public int getMaxMedianFragmentSize() {
        return readGroupToFragmentStatistics.entrySet().stream()
                .mapToInt(entry -> entry.getValue().getMedianFragmentSize())
                .max()
                .orElse(0);
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( !(obj instanceof ReadMetadata) ) return false;
        final ReadMetadata that = (ReadMetadata) obj;
        return this.contigNameToID.equals(that.contigNameToID) &&
                this.readGroupToFragmentStatistics.equals(that.readGroupToFragmentStatistics);
    }

    @Override
    public int hashCode() {
        return 47 * (47 * contigNameToID.hashCode() + readGroupToFragmentStatistics.hashCode());
    }

    private static Map<String, long[]> combineMaps( final Map<String, long[]> accumulator,
                                                    final Map<String, long[]> element ) {
        for ( final Map.Entry<String, long[]> entry : element.entrySet() ) {
            final String readGroup = entry.getKey();
            final long[] accumCounts = accumulator.get(readGroup);
            if ( accumCounts == null ) accumulator.put(readGroup, entry.getValue());
            else {
                final long[] counts = entry.getValue();
                for ( int idx = 0; idx != accumCounts.length; ++idx ) {
                    accumCounts[idx] += counts[idx];
                }
            }
        }
        return accumulator;
    }

    private static Map<String, Integer> buildContigNameToIDMap( final SAMFileHeader header ) {
        final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
        final Map<String, Integer> contigNameToID = new HashMap<>(SVUtils.hashMapCapacity(contigs.size()));
        final int nContigs = contigs.size();
        for ( int contigID = 0; contigID < nContigs; ++contigID ) {
            contigNameToID.put(contigs.get(contigID).getSequenceName(), contigID);
        }
        return contigNameToID;
    }

    public static void writeMetadata( final ReadMetadata readMetadata,
                                      final String filename,
                                      final PipelineOptions pipelineOptions ) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(filename, pipelineOptions))) ) {
            writer.write("#reads:\t" + readMetadata.getNReads() + "\n");
            writer.write("#partitions:\t" + readMetadata.getNPartitions() + "\n");
            writer.write("max reads/partition:\t" + readMetadata.getMaxReadsInPartition() + "\n");
            writer.write("coverage:\t" + readMetadata.getCoverage() + "\n");
            for ( final Map.Entry<String, ReadMetadata.ReadGroupFragmentStatistics> entry :
                    readMetadata.getAllGroupStatistics().entrySet() ) {
                final ReadMetadata.ReadGroupFragmentStatistics stats = entry.getValue();
                String name = entry.getKey();
                if ( name == null ) name = NO_GROUP;
                writer.write("group " + name + ":\t" + stats.getMedianFragmentSize() +
                        "-" + stats.getMedianNegativeDeviation() + "+" + stats.getMedianPositiveDeviation() + "\n");
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write metadata file.", ioe);
        }
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

    @DefaultSerializer(PartitionStatistics.Serializer.class)
    public static final class PartitionStatistics {
        private final static int MAX_TRACKED_FRAGMENT_LENGTH = 10000;
        private final Map<String, long[]> readGroupToFragmentSizeCountMap;
        private final long nReads;
        private final long nBases;

        public PartitionStatistics( final Iterator<GATKRead> readItr, final int nReadGroups ) {
            readGroupToFragmentSizeCountMap = new HashMap<>(SVUtils.hashMapCapacity(nReadGroups));
            long reads = 0L;
            long bases = 0L;
            while ( readItr.hasNext() ) {
                final GATKRead read = readItr.next();
                reads += 1L;
                bases += read.getLength();
                if ( read.isFirstOfPair() && !read.isSecondaryAlignment() && !read.isSupplementaryAlignment() &&
                        !read.isUnmapped() && !read.mateIsUnmapped() &&
                        Objects.equals(read.getContig(), read.getMateContig()) ) {
                    int tLen = Math.abs(read.getFragmentLength());
                    if ( tLen > MAX_TRACKED_FRAGMENT_LENGTH ) tLen = MAX_TRACKED_FRAGMENT_LENGTH;
                    final String readGroup = read.getReadGroup();
                    final long[] counts = readGroupToFragmentSizeCountMap.computeIfAbsent(readGroup,
                                                              k -> new long[MAX_TRACKED_FRAGMENT_LENGTH+1]);
                    counts[tLen] += 1;
                }
            }
            nReads = reads;
            nBases = bases;
        }

        private PartitionStatistics( final Kryo kryo, final Input input ) {
            final boolean refs = kryo.getReferences();
            kryo.setReferences(false);
            int nEntries = input.readInt();
            readGroupToFragmentSizeCountMap = new HashMap<>(SVUtils.hashMapCapacity(nEntries));
            while ( nEntries-- > 0 ) {
                final String readGroup = kryo.readObjectOrNull(input, String.class);
                final long[] counts = new long[MAX_TRACKED_FRAGMENT_LENGTH + 1];
                for ( int idx = 0; idx <= MAX_TRACKED_FRAGMENT_LENGTH; ++idx ) {
                    counts[idx] = input.readLong();
                }
                readGroupToFragmentSizeCountMap.put(readGroup, counts);
            }
            nReads = input.readLong();
            nBases = input.readLong();
            kryo.setReferences(refs);
        }

        public long getNReads() {
            return nReads;
        }

        public long getNBases() {
            return nBases;
        }

        public Map<String, long[]> getReadGroupToFragmentSizeCountMap() {
            return readGroupToFragmentSizeCountMap;
        }

        private void serialize( final Kryo kryo, final Output output ) {
            final boolean refs = kryo.getReferences();
            kryo.setReferences(false);
            output.writeInt(readGroupToFragmentSizeCountMap.size());
            for ( final Map.Entry<String, long[]> entry : readGroupToFragmentSizeCountMap.entrySet() ) {
                kryo.writeObjectOrNull(output, entry.getKey(), String.class);
                for ( final long value : entry.getValue() ) {
                    output.writeLong(value);
                }
            }
            output.writeLong(nReads);
            output.writeLong(nBases);
            kryo.setReferences(refs);
        }

        public static final class Serializer
                extends com.esotericsoftware.kryo.Serializer<PartitionStatistics> {
            @Override
            public void write( final Kryo kryo, final Output output,
                               final PartitionStatistics partitionStatistics ) {
                partitionStatistics.serialize(kryo, output);
            }

            @Override
            public PartitionStatistics read( final Kryo kryo, final Input input,
                                             final Class<PartitionStatistics> klass ) {
                return new PartitionStatistics(kryo, input);
            }
        }
    }

    /**
     * class to track distribution of fragment lengths
     */
    @DefaultSerializer(ReadGroupFragmentStatistics.Serializer.class)
    public static final class ReadGroupFragmentStatistics {
        // the distribution of fragment lengths is often quite asymmetric around the median, so we'll calculate
        // the median deviation separately for negative and positive deviations.
        private final int medianFragmentSize;
        private final int medianNegativeDeviation;
        private final int medianPositiveDeviation;

        /**
         * Given an array that counts the number of reads having a fragment length (TLEN) equal to the array index,
         * figure out the statistics.
         */
        public ReadGroupFragmentStatistics( final long[] counts ) {
            // total number of reads
            final long total = Arrays.stream(counts).sum();

            // calculate the median fragment length
            long sum = 0L;
            int medianFragmentSize = 0; // increment this, summing counts as we go, until we've encountered 1/2 the reads
            while ( medianFragmentSize != counts.length ) {
                sum += counts[medianFragmentSize];
                if ( 2 * sum >= total ) break; // break if we've seen 1/2 of the reads -- we've discovered the median
                medianFragmentSize += 1;
            }
            this.medianFragmentSize = medianFragmentSize;

            // calculate the median negative deviation
            sum = counts[medianFragmentSize] / 2; // half the counts in the median bin go with the negative deviation
            int medianNegativeDeviation = 0; // increment this, summing counts as we walk down and away from the median bin
            // until we've seen 1/4 of the reads
            while ( 4 * sum < total && medianNegativeDeviation < medianFragmentSize ) {
                medianNegativeDeviation += 1;
                sum += counts[medianFragmentSize - medianNegativeDeviation];
            }
            this.medianNegativeDeviation = medianNegativeDeviation;

            // calculate the median positive deviation
            final int maxMedianPositiveDeviation = counts.length - medianFragmentSize - 1; // array boundary
            sum = counts[medianFragmentSize] / 2; // half the counts in the median bin go with the positive deviation
            int medianPositiveDeviation = 0; // increment this, summing counts as we walk up and away from the median bin
            // until we've seen 1/4 of the reads
            while ( 4 * sum < total && medianPositiveDeviation < maxMedianPositiveDeviation ) {
                medianPositiveDeviation += 1;
                sum += counts[medianFragmentSize + medianPositiveDeviation];
            }
            this.medianPositiveDeviation = medianPositiveDeviation;
        }

        @VisibleForTesting
        ReadGroupFragmentStatistics( final int medianFragmentSize,
                                     final int medianNegativeDeviation,
                                     final int medianPositiveDeviation ) {
            this.medianFragmentSize = medianFragmentSize;
            this.medianNegativeDeviation = medianNegativeDeviation;
            this.medianPositiveDeviation = medianPositiveDeviation;
        }

        private ReadGroupFragmentStatistics( final Kryo kryo, final Input input ) {
            medianFragmentSize = input.readInt();
            medianNegativeDeviation = input.readInt();
            medianPositiveDeviation = input.readInt();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(medianFragmentSize);
            output.writeInt(medianNegativeDeviation);
            output.writeInt(medianPositiveDeviation);
        }

        public int getMedianFragmentSize() {
            return medianFragmentSize;
        }

        public int getMedianNegativeDeviation() {
            return medianNegativeDeviation;
        }

        public int getMedianPositiveDeviation() {
            return medianPositiveDeviation;
        }

        public float getZIshScore( final int fragmentSize ) {
            if ( fragmentSize < 0 ) throw new GATKException("negative fragment size");
            final int diff = fragmentSize - medianFragmentSize;
            if ( diff == 0 ) return 0.0f;
            if ( diff > 0 ) return 1.0f * diff / medianPositiveDeviation;
            return 1.0f * diff / medianNegativeDeviation;
        }

        @Override
        public boolean equals( final Object obj ) {
            if ( !(obj instanceof ReadGroupFragmentStatistics) ) return false;
            final ReadGroupFragmentStatistics that = (ReadGroupFragmentStatistics) obj;
            return this.medianFragmentSize == that.medianFragmentSize &&
                    this.medianNegativeDeviation == that.medianNegativeDeviation &&
                    this.medianPositiveDeviation == that.medianPositiveDeviation;
        }

        @Override
        public int hashCode() {
            return 47 * (47 * (47 * medianFragmentSize + medianNegativeDeviation) + medianPositiveDeviation);
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
