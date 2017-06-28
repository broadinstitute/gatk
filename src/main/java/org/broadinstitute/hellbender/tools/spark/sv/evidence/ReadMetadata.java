package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.*;

/**
 * A bag of data about reads:  contig name to id mapping, fragment length statistics by read group, mean length.
 * The fragment length statistics pertain to a library, but they're accessed by read group.  (I.e., in the
 * readGroupToFragmentStatistics map, all read groups that are derived from a given library point to the same set of
 * statistics.)
 */
@DefaultSerializer(ReadMetadata.Serializer.class)
public class ReadMetadata {
    private final Set<Integer> crossContigIgnoreSet;
    private final Map<String, Integer> contigNameToID;
    private final long nReads;
    private final long maxReadsInPartition;
    private final int coverage;
    private final PartitionBounds[] partitionBounds;
    private final Map<String, LibraryFragmentStatistics> readGroupToFragmentStatistics;
    private final static String NO_GROUP = "NoGroup";

    public ReadMetadata( final Set<Integer> crossContigIgnoreSet,
                         final SAMFileHeader header,
                         final JavaRDD<GATKRead> mappedReads ) {
        this.crossContigIgnoreSet = crossContigIgnoreSet;
        contigNameToID = buildContigNameToIDMap(header);

        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        final int mapCapacity = SVUtils.hashMapCapacity(header.getReadGroups().size());
        final Map<String, String> readGroupToLibraryMap = new HashMap<>(mapCapacity);
        for ( final SAMReadGroupRecord groupRecord : readGroups ) {
            readGroupToLibraryMap.put(groupRecord.getId(), groupRecord.getLibrary());
        }
        final List<PartitionStatistics> perPartitionStatistics =
                mappedReads.mapPartitions(readItr ->
                        Collections.singletonList(new PartitionStatistics(readItr, readGroupToLibraryMap)).iterator())
                     .collect();
        nReads = perPartitionStatistics.stream().mapToLong(PartitionStatistics::getNReads).sum();
        maxReadsInPartition = perPartitionStatistics.stream().mapToLong(PartitionStatistics::getNReads).max().orElse(0L);
        final long nReadBases = perPartitionStatistics.stream().mapToLong(PartitionStatistics::getNBases).sum();
        final long nRefBases = header.getSequenceDictionary().getSequences()
                .stream().mapToLong(SAMSequenceRecord::getSequenceLength).sum();
        coverage = (int)((nReadBases + nRefBases/2) / nRefBases); // rounding the coverage number
        final int nPartitions = perPartitionStatistics.size();
        partitionBounds = new PartitionBounds[nPartitions];
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            final PartitionStatistics stats = perPartitionStatistics.get(idx);
            final Integer firstContigID = contigNameToID.get(stats.getFirstContig());
            final Integer lastContigID = contigNameToID.get(stats.getLastContig());
            partitionBounds[idx] = new PartitionBounds(
                    firstContigID==null ? PartitionBounds.UNMAPPED : firstContigID,
                    stats.getFirstLocation(),
                    lastContigID==null ? PartitionBounds.UNMAPPED : lastContigID,
                    stats.getLastLocation());
        }
        final Map<String, long[]> combinedMaps =
                perPartitionStatistics.stream()
                        .map(PartitionStatistics::getLibraryToFragmentSizeCountMap)
                        .reduce(new HashMap<>(mapCapacity), ReadMetadata::combineMaps);
        final Map<String, LibraryFragmentStatistics> statsMap = new HashMap<>(mapCapacity);
        for ( final Map.Entry<String, long[]> entry : combinedMaps.entrySet() ) {
            statsMap.put(entry.getKey(), new LibraryFragmentStatistics(entry.getValue()));
        }
        readGroupToFragmentStatistics = new HashMap<>(mapCapacity);
        for ( final SAMReadGroupRecord groupRecord : readGroups ) {
            readGroupToFragmentStatistics.put(groupRecord.getId(), statsMap.get(groupRecord.getLibrary()));
        }
    }

    @VisibleForTesting
    ReadMetadata( final Set<Integer> crossContigIgnoreSet, final SAMFileHeader header,
                  final LibraryFragmentStatistics stats, final PartitionBounds[] partitionBounds,
                  final long nReads, final long maxReadsInPartition, final int coverage ) {
        this.crossContigIgnoreSet = crossContigIgnoreSet;
        contigNameToID = buildContigNameToIDMap(header);
        this.nReads = nReads;
        this.maxReadsInPartition = maxReadsInPartition;
        this.coverage = coverage;
        this.partitionBounds = partitionBounds;
        readGroupToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(header.getReadGroups().size() + 1));
        readGroupToFragmentStatistics.put(null, stats);
        for ( final SAMReadGroupRecord readGroupRecord : header.getReadGroups() ) {
            readGroupToFragmentStatistics.put(readGroupRecord.getReadGroupId(), stats);
        }
    }

    private ReadMetadata( final Kryo kryo, final Input input ) {
        final int crossContigIgnoreSetSize = input.readInt();
        this.crossContigIgnoreSet = new HashSet<>(SVUtils.hashMapCapacity(crossContigIgnoreSetSize));
        for ( int idx = 0; idx != crossContigIgnoreSetSize; ++idx ) {
            crossContigIgnoreSet.add(input.readInt());
        }

        final int contigMapSize = input.readInt();
        contigNameToID = new HashMap<>(SVUtils.hashMapCapacity(contigMapSize));
        for ( int idx = 0; idx != contigMapSize; ++idx ) {
            final String contigName = input.readString();
            final int contigId = input.readInt();
            contigNameToID.put(contigName, contigId);
        }

        nReads = input.readLong();
        maxReadsInPartition = input.readLong();
        coverage = input.readInt();

        final int nPartitions = input.readInt();
        partitionBounds = new PartitionBounds[nPartitions];
        final PartitionBounds.Serializer boundsSerializer = new PartitionBounds.Serializer();
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            partitionBounds[idx] = boundsSerializer.read(kryo, input, PartitionBounds.class);
        }

        final int readGroupMapSize = input.readInt();
        readGroupToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(readGroupMapSize));
        final LibraryFragmentStatistics.Serializer statsSerializer =
                new LibraryFragmentStatistics.Serializer();
        for ( int idx = 0; idx != readGroupMapSize; ++idx ) {
            final String readGroupName = input.readString();
            final LibraryFragmentStatistics groupStats =
                    statsSerializer.read(kryo, input, LibraryFragmentStatistics.class);
            readGroupToFragmentStatistics.put(readGroupName, groupStats);
        }
    }

    private void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(crossContigIgnoreSet.size());
        for ( final Integer tigId : crossContigIgnoreSet ) {
            output.writeInt(tigId);
        }

        output.writeInt(contigNameToID.size());
        for ( final Map.Entry<String, Integer> entry : contigNameToID.entrySet() ) {
            output.writeString(entry.getKey());
            output.writeInt(entry.getValue());
        }

        output.writeLong(nReads);
        output.writeLong(maxReadsInPartition);
        output.writeInt(coverage);

        output.writeInt(partitionBounds.length);
        final PartitionBounds.Serializer boundsSerializer = new PartitionBounds.Serializer();
        for ( final PartitionBounds bounds : partitionBounds ) {
            boundsSerializer.write(kryo, output, bounds);
        }

        output.writeInt(readGroupToFragmentStatistics.size());
        final LibraryFragmentStatistics.Serializer statsSerializer =
                                    new LibraryFragmentStatistics.Serializer();
        for ( final Map.Entry<String, LibraryFragmentStatistics> entry : readGroupToFragmentStatistics.entrySet() ) {
            output.writeString(entry.getKey());
            statsSerializer.write(kryo, output, entry.getValue());
        }
    }

    public boolean ignoreCrossContigID( final int contigID ) { return crossContigIgnoreSet.contains(contigID); }

    public Map<String, Integer> getContigNameMap() {
        return Collections.unmodifiableMap(contigNameToID);
    }

    public int getContigID( final String contigName ) {
        final Integer result = contigNameToID.get(contigName);
        if ( result == null ) {
            throw new GATKException("No such contig name: " + contigName);
        }
        return result;
    }

    public long getNReads() {
        return nReads;
    }

    public int getNPartitions() { return partitionBounds.length; }
    public PartitionBounds getPartitionBounds( final int partitionIdx ) { return partitionBounds[partitionIdx]; }

    public long getMaxReadsInPartition() {
        return maxReadsInPartition;
    }

    public int getCoverage() {
        return coverage;
    }

    public Map<String, LibraryFragmentStatistics> getAllGroupStatistics() {
        return readGroupToFragmentStatistics;
    }

    public LibraryFragmentStatistics getStatistics( final String readGroupName ) {
        final LibraryFragmentStatistics stats = readGroupToFragmentStatistics.get(readGroupName);
        if ( stats == null ) {
            throw new GATKException("No such read group name: " + readGroupName);
        }
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
        return this.crossContigIgnoreSet.equals(that.crossContigIgnoreSet) &&
                this.contigNameToID.equals(that.contigNameToID) &&
                this.readGroupToFragmentStatistics.equals(that.readGroupToFragmentStatistics);
    }

    @Override
    public int hashCode() {
        int val = crossContigIgnoreSet.hashCode();
        val = 47 * val + contigNameToID.hashCode();
        val = 47 * val + readGroupToFragmentStatistics.hashCode();
        return 47 * val;
    }

    private static Map<String, long[]> combineMaps( final Map<String, long[]> accumulator,
                                                    final Map<String, long[]> element ) {
        for ( final Map.Entry<String, long[]> entry : element.entrySet() ) {
            final String readGroup = entry.getKey();
            final long[] accumCounts = accumulator.get(readGroup);
            if ( accumCounts == null ) {
                accumulator.put(readGroup, entry.getValue());
            } else {
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
                                      final String filename ) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(filename))) ) {
            writer.write("#reads:\t" + readMetadata.getNReads() + "\n");
            writer.write("#partitions:\t" + readMetadata.getNPartitions() + "\n");
            writer.write("max reads/partition:\t" + readMetadata.getMaxReadsInPartition() + "\n");
            writer.write("coverage:\t" + readMetadata.getCoverage() + "\n");
            for ( final Map.Entry<String, LibraryFragmentStatistics> entry :
                    readMetadata.getAllGroupStatistics().entrySet() ) {
                final LibraryFragmentStatistics stats = entry.getValue();
                String name = entry.getKey();
                if ( name == null ) {
                    name = NO_GROUP;
                }
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
        private final Map<String, long[]> libraryToFragmentSizeCountMap;
        private final long nReads;
        private final long nBases;
        private final String firstContig;
        private final int firstLocation;
        private final String lastContig;
        private final int lastLocation;

        public PartitionStatistics( final Iterator<GATKRead> mappedReadItr,
                                    final Map<String, String> readGroupToLibraryMap ) {
            if ( !mappedReadItr.hasNext() ) {
                libraryToFragmentSizeCountMap = Collections.emptyMap();
                nReads = 0L;
                nBases = 0L;
                firstContig = null;
                firstLocation = 0;
                lastContig = null;
                lastLocation = 0;
            } else {
                long reads = 0L;
                long bases = 0L;
                GATKRead read = mappedReadItr.next();
                firstContig = read.getContig();
                firstLocation = read.getUnclippedStart();
                libraryToFragmentSizeCountMap = new HashMap<>(SVUtils.hashMapCapacity(readGroupToLibraryMap.size()));
                while ( true ) {
                    reads += 1L;
                    bases += read.getLength();
                    if ( read.isFirstOfPair() && !read.isSecondaryAlignment() && !read.isSupplementaryAlignment() &&
                            !read.isUnmapped() && !read.mateIsUnmapped() &&
                            Objects.equals(read.getContig(), read.getMateContig()) ) {
                        int tLen = Math.abs(read.getFragmentLength());
                        if ( tLen > MAX_TRACKED_FRAGMENT_LENGTH ) {
                            tLen = MAX_TRACKED_FRAGMENT_LENGTH;
                        }
                        final String library = readGroupToLibraryMap.get(read.getReadGroup());
                        final long[] counts = libraryToFragmentSizeCountMap.computeIfAbsent(library,
                                k -> new long[MAX_TRACKED_FRAGMENT_LENGTH + 1]);
                        counts[tLen] += 1;
                    }
                    if ( !mappedReadItr.hasNext() ) break;
                    read = mappedReadItr.next();
                }
                lastContig = read.getContig();
                lastLocation = read.getUnclippedEnd() + 1;
                nReads = reads;
                nBases = bases;
            }
        }

        private PartitionStatistics( final Kryo kryo, final Input input ) {
            int nEntries = input.readInt();
            libraryToFragmentSizeCountMap = new HashMap<>(SVUtils.hashMapCapacity(nEntries));
            while ( nEntries-- > 0 ) {
                final String libName = kryo.readObjectOrNull(input, String.class);
                final long[] counts = new long[MAX_TRACKED_FRAGMENT_LENGTH + 1];
                for ( int idx = 0; idx <= MAX_TRACKED_FRAGMENT_LENGTH; ++idx ) {
                    counts[idx] = input.readLong();
                }
                libraryToFragmentSizeCountMap.put(libName, counts);
            }
            nReads = input.readLong();
            nBases = input.readLong();
            firstContig = input.readString();
            firstLocation = input.readInt();
            lastContig = input.readString();
            lastLocation = input.readInt();
        }

        public long getNReads() {
            return nReads;
        }

        public long getNBases() {
            return nBases;
        }

        public Map<String, long[]> getLibraryToFragmentSizeCountMap() {
            return libraryToFragmentSizeCountMap;
        }

        public String getFirstContig() { return firstContig; }
        public int getFirstLocation() { return firstLocation; }
        public String getLastContig() { return lastContig; }
        public int getLastLocation() { return lastLocation; }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(libraryToFragmentSizeCountMap.size());
            for ( final Map.Entry<String, long[]> entry : libraryToFragmentSizeCountMap.entrySet() ) {
                kryo.writeObjectOrNull(output, entry.getKey(), String.class);
                for ( final long value : entry.getValue() ) {
                    output.writeLong(value);
                }
            }
            output.writeLong(nReads);
            output.writeLong(nBases);
            output.writeString(firstContig);
            output.writeInt(firstLocation);
            output.writeString(lastContig);
            output.writeInt(lastLocation);
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

    /** A class to track the genomic location of the start of the first and last mapped reads in a partition. */
    @DefaultSerializer(PartitionBounds.Serializer.class)
    public static final class PartitionBounds {
        private final int firstContigID;
        private final int firstStart;
        private final int lastContigID;
        private final int lastStart;
        public final static int UNMAPPED = -1;

        public PartitionBounds( final int firstContigID, final int firstStart,
                                final int lastContigID, final int lastStart ) {
            this.firstContigID = firstContigID;
            this.firstStart = firstStart;
            this.lastContigID = lastContigID;
            this.lastStart = lastStart;
        }

        private PartitionBounds( final Kryo kryo, final Input input ) {
            this.firstContigID = input.readInt();
            this.firstStart = input.readInt();
            this.lastContigID = input.readInt();
            this.lastStart = input.readInt();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(firstContigID);
            output.writeInt(firstStart);
            output.writeInt(lastContigID);
            output.writeInt(lastStart);
        }

        public int getFirstContigID() { return firstContigID; }
        public int getFirstStart() { return firstStart; }
        public int getLastContigID() { return lastContigID; }
        public int getLastStart() { return lastStart; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<PartitionBounds> {
            @Override
            public void write( final Kryo kryo, final Output output, final PartitionBounds partitionBounds ) {
                partitionBounds.serialize(kryo, output);
            }

            @Override
            public PartitionBounds read( final Kryo kryo, final Input input, final Class<PartitionBounds> klass ) {
                return new PartitionBounds(kryo, input);
            }
        }
    }

    /**
     * class to track distribution of fragment lengths
     */
    @DefaultSerializer(LibraryFragmentStatistics.Serializer.class)
    public static final class LibraryFragmentStatistics {
        // the distribution of fragment lengths is often quite asymmetric around the median, so we'll calculate
        // the median deviation separately for negative and positive deviations.
        private final int medianFragmentSize;
        private final int medianNegativeDeviation;
        private final int medianPositiveDeviation;

        /**
         * Given an array that counts the number of reads having a fragment length (TLEN) equal to the array index,
         * figure out the statistics.
         */
        public LibraryFragmentStatistics( final long[] counts ) {
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
        LibraryFragmentStatistics( final int medianFragmentSize,
                                   final int medianNegativeDeviation,
                                   final int medianPositiveDeviation ) {
            this.medianFragmentSize = medianFragmentSize;
            this.medianNegativeDeviation = medianNegativeDeviation;
            this.medianPositiveDeviation = medianPositiveDeviation;
        }

        private LibraryFragmentStatistics( final Kryo kryo, final Input input ) {
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
            if ( fragmentSize < 0 ) {
                throw new GATKException("negative fragment size");
            }
            final int diff = fragmentSize - medianFragmentSize;
            if ( diff == 0 ) return 0.0f;
            if ( diff > 0 ) return 1.0f * diff / medianPositiveDeviation;
            return 1.0f * diff / medianNegativeDeviation;
        }

        @Override
        public boolean equals( final Object obj ) {
            if ( !(obj instanceof LibraryFragmentStatistics) ) return false;
            final LibraryFragmentStatistics that = (LibraryFragmentStatistics) obj;
            return this.medianFragmentSize == that.medianFragmentSize &&
                    this.medianNegativeDeviation == that.medianNegativeDeviation &&
                    this.medianPositiveDeviation == that.medianPositiveDeviation;
        }

        @Override
        public int hashCode() {
            return 47 * (47 * (47 * medianFragmentSize + medianNegativeDeviation) + medianPositiveDeviation);
        }

        public static final class Serializer
                extends com.esotericsoftware.kryo.Serializer<LibraryFragmentStatistics> {
            @Override
            public void write( final Kryo kryo, final Output output,
                               final LibraryFragmentStatistics libraryFragmentStatistics ) {
                libraryFragmentStatistics.serialize(kryo, output);
            }

            @Override
            public LibraryFragmentStatistics read( final Kryo kryo, final Input input,
                                                   final Class<LibraryFragmentStatistics> klass ) {
                return new LibraryFragmentStatistics(kryo, input);
            }
        }
    }
}
