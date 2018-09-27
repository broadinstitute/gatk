package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.spark.GATKRegistrator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
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

    public static final String CDF_PREFIX = "template size cumulative counts:";
    private final Set<Integer> crossContigIgnoreSet;
    private final Map<String, Integer> contigNameToID;
    private final String[] contigIDToName;
    private final Map<String, String> readGroupToLibrary;
    private final long nReads;
    private final int avgReadLen;
    private final long nRefBases;
    private final long maxReadsInPartition;
    private final float coverage;
    private final float meanBaseQuality;
    private final PartitionBounds[] partitionBounds;
    private final Map<String, LibraryStatistics> libraryToFragmentStatistics;
    private static final String NO_GROUP = "NoGroup";
    private static final float MIN_COVERAGE = 10.f;
    private static final float DEFAULT_MEAN_BASE_QUALITY_FOR_TESTING = 30.f;

    public ReadMetadata( final Set<Integer> crossContigIgnoreSet,
                         final SAMFileHeader header,
                         final int maxTrackedFragmentLength,
                         final JavaRDD<GATKRead> unfilteredReads,
                         final SVReadFilter filter,
                         final Logger logger ) {
        this.crossContigIgnoreSet = crossContigIgnoreSet;
        contigNameToID = buildContigNameToIDMap(header.getSequenceDictionary());
        contigIDToName = buildContigIDToNameArray(contigNameToID);
        readGroupToLibrary = buildGroupToLibMap(header);
        final Map<String, String> grpToLib = readGroupToLibrary;
        final List<PartitionStatistics> perPartitionStatistics =
                unfilteredReads
                    .mapPartitions(readItr ->
                        SVUtils.singletonIterator(
                                new PartitionStatistics(readItr, filter, maxTrackedFragmentLength, grpToLib)))
                    .collect();
        maxReadsInPartition = perPartitionStatistics.stream().mapToLong(PartitionStatistics::getNReads).max().orElse(0L);
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
                    stats.getLastLocation(),
                    stats.getSpan() );
        }
        final Map<String, LibraryRawStatistics> combinedMaps =
                perPartitionStatistics.stream()
                        .map(PartitionStatistics::getLibraryNameToStatisticsMap)
                        .reduce(new HashMap<>(), ReadMetadata::combineMaps);
        nReads = combinedMaps.values().stream().mapToLong(LibraryRawStatistics::getNReads).sum();
        final long nReadBases = combinedMaps.values().stream().mapToLong(LibraryRawStatistics::getNBases).sum();
        avgReadLen = (int)(nReadBases / nReads);
        nRefBases = header.getSequenceDictionary().getSequences()
                .stream().mapToLong(SAMSequenceRecord::getSequenceLength).sum();
        final float cov = (float)nReadBases / nRefBases;
        if ( cov >= MIN_COVERAGE ) coverage = cov;
        else {
            logger.warn("Apparent coverage (" + cov + ") too low.  Pretending it's 10x.");
            coverage = MIN_COVERAGE;
        }
        final long totalBaseQuals =
                combinedMaps.values().stream().mapToLong(LibraryRawStatistics::getTotalBaseQuality).sum();
        meanBaseQuality = (float)totalBaseQuals / nReadBases;
        libraryToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(combinedMaps.size()));
        combinedMaps.forEach( (libName, rawStats) ->
                libraryToFragmentStatistics.put(libName,rawStats.createLibraryStatistics(nRefBases)));
    }

    /** This constructor is for testing only.  It applies a single LibraryStatistics object to all libraries. */
    @VisibleForTesting
    ReadMetadata( final Set<Integer> crossContigIgnoreSet, final SAMFileHeader header,
                  final LibraryStatistics stats, final PartitionBounds[] partitionBounds,
                  final long nReads, final long maxReadsInPartition, final float coverage ) {
        this.crossContigIgnoreSet = crossContigIgnoreSet;
        contigNameToID = buildContigNameToIDMap(header.getSequenceDictionary());
        contigIDToName = buildContigIDToNameArray(contigNameToID);
        readGroupToLibrary = buildGroupToLibMap(header);
        this.nReads = nReads;
        nRefBases = header.getSequenceDictionary().getSequences()
                .stream().mapToLong(SAMSequenceRecord::getSequenceLength).sum();
        avgReadLen = (int)(coverage * nRefBases / nReads);
        this.maxReadsInPartition = maxReadsInPartition;
        this.coverage = coverage;
        meanBaseQuality = DEFAULT_MEAN_BASE_QUALITY_FOR_TESTING;
        this.partitionBounds = partitionBounds;
        libraryToFragmentStatistics = new HashMap<>(6);
        libraryToFragmentStatistics.put(null, stats);
        for ( final SAMReadGroupRecord readGroupRecord : header.getReadGroups() ) {
            libraryToFragmentStatistics.put(readGroupRecord.getLibrary(), stats);
        }
    }

    private ReadMetadata( final Kryo kryo, final Input input ) {
        final int crossContigIgnoreSetSize = input.readInt();
        this.crossContigIgnoreSet = new HashSet<>(SVUtils.hashMapCapacity(crossContigIgnoreSetSize));
        for ( int idx = 0; idx != crossContigIgnoreSetSize; ++idx ) {
            crossContigIgnoreSet.add(input.readInt());
        }

        final int groupMapSize = input.readInt();
        readGroupToLibrary = new HashMap<>(SVUtils.hashMapCapacity(groupMapSize));
        for ( int idx = 0; idx != groupMapSize; ++idx ) {
            final String groupName = input.readString();
            final String libName = input.readString();
            readGroupToLibrary.put(groupName, libName);
        }

        final int contigMapSize = input.readInt();
        contigNameToID = new HashMap<>(SVUtils.hashMapCapacity(contigMapSize));
        for ( int idx = 0; idx != contigMapSize; ++idx ) {
            final String contigName = input.readString();
            final int contigId = input.readInt();
            contigNameToID.put(contigName, contigId);
        }
        contigIDToName = buildContigIDToNameArray(contigNameToID);

        nReads = input.readLong();
        avgReadLen = input.readInt();
        nRefBases = input.readLong();
        maxReadsInPartition = input.readLong();
        coverage = input.readFloat();
        meanBaseQuality = input.readFloat();

        final int nPartitions = input.readInt();
        partitionBounds = new PartitionBounds[nPartitions];
        final PartitionBounds.Serializer boundsSerializer = new PartitionBounds.Serializer();
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            partitionBounds[idx] = boundsSerializer.read(kryo, input, PartitionBounds.class);
        }

        final int libMapSize = input.readInt();
        final LibraryStatistics.Serializer statsSerializer = new LibraryStatistics.Serializer();
        libraryToFragmentStatistics = new HashMap<>(SVUtils.hashMapCapacity(libMapSize));
        for ( int idx = 0; idx != libMapSize; ++idx ) {
            final String libraryName = input.readString();
            final LibraryStatistics stats = statsSerializer.read(kryo, input, LibraryStatistics.class);
            libraryToFragmentStatistics.put(libraryName, stats);
        }

    }

    private void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(crossContigIgnoreSet.size());
        for ( final Integer tigId : crossContigIgnoreSet ) {
            output.writeInt(tigId);
        }

        output.writeInt(readGroupToLibrary.size());
        for ( final Map.Entry<String, String> entry : readGroupToLibrary.entrySet() ) {
            output.writeString(entry.getKey());
            output.writeString(entry.getValue());
        }

        output.writeInt(contigNameToID.size());
        for ( final Map.Entry<String, Integer> entry : contigNameToID.entrySet() ) {
            output.writeString(entry.getKey());
            output.writeInt(entry.getValue());
        }

        output.writeLong(nReads);
        output.writeInt(avgReadLen);
        output.writeLong(nRefBases);
        output.writeLong(maxReadsInPartition);
        output.writeFloat(coverage);
        output.writeFloat(meanBaseQuality);

        output.writeInt(partitionBounds.length);
        final PartitionBounds.Serializer boundsSerializer = new PartitionBounds.Serializer();
        for ( final PartitionBounds bounds : partitionBounds ) {
            boundsSerializer.write(kryo, output, bounds);
        }

        output.writeInt(libraryToFragmentStatistics.size());
        final LibraryStatistics.Serializer statsSerializer = new LibraryStatistics.Serializer();
        for ( final Map.Entry<String, LibraryStatistics> entry : libraryToFragmentStatistics.entrySet() ) {
            output.writeString(entry.getKey());
            statsSerializer.write(kryo, output, entry.getValue());
        }

    }

    public boolean ignoreCrossContigID( final int contigID ) { return crossContigIgnoreSet.contains(contigID); }
    @VisibleForTesting Set<Integer> getCrossContigIgnoreSet() { return crossContigIgnoreSet; }

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

    public String getContigName( final int contigID ) { return contigIDToName[contigID]; }

    public String getLibraryName( final String readGroupName ) {
        if ( readGroupName == null ) return null;
        if ( !readGroupToLibrary.containsKey(readGroupName) ) {
            throw new GATKException("No such read group in header: "+readGroupName);
        }
        return readGroupToLibrary.get(readGroupName);
    }

    @VisibleForTesting Map<String, String> getReadGroupToLibraryMap() { return readGroupToLibrary; }

    public float getZishScore( final String readGroup, final int fragmentSize ) {
        return getFragmentLengthStatistics(readGroup).getZishScore(fragmentSize);
    }

    @VisibleForTesting
    LibraryStatistics getFragmentLengthStatistics( final String readGroup ) {
        return libraryToFragmentStatistics.get(getLibraryName(readGroup));
    }

    public long getNReads() { return nReads; }
    public int getAvgReadLen() { return avgReadLen; }
    public long getNRefBases() { return nRefBases; }
    public int getNPartitions() { return partitionBounds.length; }
    public PartitionBounds getPartitionBounds( final int partitionIdx ) { return partitionBounds[partitionIdx]; }
    @VisibleForTesting PartitionBounds[] getAllPartitionBounds() { return partitionBounds; }

    public long getMaxReadsInPartition() { return maxReadsInPartition; }
    public float getCoverage() {
        return coverage;
    }
    public float getMeanBaseQuality() { return meanBaseQuality; }
    public float getAccurateKmerCoverage( final int kSize ) {
        // i.e., p = (1 - 10**(-Q/10))**K
        final float probabilityOfAnAccurateKmer =
                (float)Math.exp(kSize * Math.log1p(-Math.pow(10.,meanBaseQuality/-10.)));
        return coverage * probabilityOfAnAccurateKmer;
    }

    public int getMedianPartitionSpan() {
        final int[] spans = new int[partitionBounds.length];
        for ( int idx = 0; idx != partitionBounds.length; ++idx ) {
            spans[idx] = partitionBounds[idx].getSpan();
        }
        Arrays.sort(spans);
        return spans[partitionBounds.length/2];
    }

    public Map<String, LibraryStatistics> getAllLibraryStatistics() { return libraryToFragmentStatistics; }

    public LibraryStatistics getLibraryStatistics( final String libraryName ) {
        final LibraryStatistics stats = libraryToFragmentStatistics.get(libraryName);
        if ( stats == null ) {
            throw new GATKException("No such library: " + libraryName);
        }
        return stats;
    }

    public int getMaxMedianFragmentSize() {
        return libraryToFragmentStatistics.entrySet().stream()
                .mapToInt(entry -> entry.getValue().getMedian())
                .max()
                .orElse(0);
    }

    private static Map<String, LibraryRawStatistics> combineMaps( final Map<String, LibraryRawStatistics> accumulator,
                                                    final Map<String, LibraryRawStatistics> element ) {
        for ( final Map.Entry<String, LibraryRawStatistics> entry : element.entrySet() ) {
            final String libraryName = entry.getKey();
            final LibraryRawStatistics accumulatorStats = accumulator.get(libraryName);
            if ( accumulatorStats == null ) {
                accumulator.put(libraryName, entry.getValue());
            } else {
                LibraryRawStatistics.reduce(accumulatorStats, entry.getValue());
            }
        }
        return accumulator;
    }

    public static Map<String, Integer> buildContigNameToIDMap( final SAMSequenceDictionary dictionary ) {
        final List<SAMSequenceRecord> contigs = dictionary.getSequences();
        final Map<String, Integer> contigNameToID = new HashMap<>(SVUtils.hashMapCapacity(contigs.size()));
        final int nContigs = contigs.size();
        for ( int contigID = 0; contigID < nContigs; ++contigID ) {
            contigNameToID.put(contigs.get(contigID).getSequenceName(), contigID);
        }
        return contigNameToID;
    }

    public static String[] buildContigIDToNameArray( final Map<String, Integer> nameToIDMap ) {
        final String[] result = new String[nameToIDMap.size()];
        for ( final Map.Entry<String, Integer> entry : nameToIDMap.entrySet() ) {
            result[entry.getValue()] = entry.getKey();
        }
        return result;
    }

    public static Map<String, String> buildGroupToLibMap( final SAMFileHeader header ) {
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        final int mapCapacity = SVUtils.hashMapCapacity(header.getReadGroups().size());
        final Map<String, String> readGroupToLibraryMap = new HashMap<>(mapCapacity);
        for ( final SAMReadGroupRecord groupRecord : readGroups ) {
            readGroupToLibraryMap.put(groupRecord.getId(), groupRecord.getLibrary());
        }
        return readGroupToLibraryMap;
    }

    public static void writeMetadata(final ReadMetadata readMetadata,
                                     final String filename ) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(filename))) ) {
            writer.write("#reads:\t" + readMetadata.getNReads() + "\n");
            writer.write("#partitions:\t" + readMetadata.getNPartitions() + "\n");
            writer.write("max reads/partition:\t" + readMetadata.getMaxReadsInPartition() + "\n");
            writer.write("coverage:\t" + readMetadata.getCoverage() + "\n");
            writer.write( "meanQ:\t" + readMetadata.getMeanBaseQuality() + "\n");
            writer.write("\nLibrary Statistics\n");
            for ( final Map.Entry<String, LibraryStatistics> entry :
                    readMetadata.getAllLibraryStatistics().entrySet() ) {
                final LibraryStatistics stats = entry.getValue();
                String name = entry.getKey();
                if ( name == null ) {
                    name = NO_GROUP;
                }
                final int median = stats.getMedian();
                writer.write(name + ":\t" + median + "-" + stats.getNegativeMAD() +
                                "+" + stats.getPositiveMAD() + "\t" + stats.getCoverage() + "\t" +
                                stats.getMeanBaseQuality() + "\t" +
                                stats.getNReads() + "\t" + stats.getReadStartFrequency() + "\n");
                final IntHistogram.CDF templateSizeCDF = stats.getCDF();
                final int cdfSize = templateSizeCDF.size();
                final long totalObservations = templateSizeCDF.getTotalObservations();
                writer.write(CDF_PREFIX);
                for(int idx = 0; idx < cdfSize; ++idx) {
                    final long cumulativeCounts = Math.round(templateSizeCDF.getFraction(idx) * totalObservations);
                    writer.write("\t" + cumulativeCounts);
                }
                writer.write("\n");
            }
            final PartitionBounds[] partitionBounds = readMetadata.partitionBounds;
            writer.write("\nPartition Boundaries\n");
            for ( int idx = 0; idx != partitionBounds.length; ++idx ) {
                final PartitionBounds bounds = partitionBounds[idx];
                writer.write(idx + "\t" + bounds.firstContigID + "\t" + bounds.getFirstStart() + "\t" +
                        bounds.getLastContigID() + "\t" + bounds.getLastEnd() + "\n");
            }
            writer.write("contigs map:\n");
            try {
                for (int i = 0; i < readMetadata.contigIDToName.length; ++i) {
                    writer.write(i + ":" + readMetadata.contigIDToName[i] + "\n");
                }
            } catch ( final IOException ex ) {
                throw new GATKException("Can't write metadata contig entry", ex);
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write metadata file.", ioe);
        }
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ReadMetadata> {

        /**
         * Placed at the beginning of each standalone serialization of the meta-data
         * identifies file that most likely are valid serializations.
         * <p>
         * Is not a 100% guarantee as one construct an invalid file that contains,
         * but is extremely-unlikely to be accidental.
         * </p>
         */
        public static final String MAGIC_STRING = "9wdgy2yEbw0jg";

        /**
         * Serialization version.
         * <p>
         *     This follows the {@link #MAGIC_STRING} in the serialization file indicating
         *     the format version.
         * </p>
         * <p>
         *     This number/string should be change every time the serialization format changes
         *     so that we can catch incompatible deserialization early.
         * </p>
         */
        public static final String VERSION_STRING = "0.1";

        @Override
        public void write( final Kryo kryo, final Output output, final ReadMetadata readMetadata ) {
            readMetadata.serialize(kryo, output);
        }

        @Override
        public ReadMetadata read( final Kryo kryo, final Input input, final Class<ReadMetadata> klass ) {
            return new ReadMetadata(kryo, input);
        }

        /**
         * Serializes a read-metadata by itself into a file.
         * @param meta the read-metadata to serialize.
         * @param whereTo the name of the file or resource where it will go to.
         * @throws IllegalArgumentException if either {@code meta} or {@code whereTo}
         * is {@code null}.
         * @throws UserException if there was a problem during serialization.
         */
        public static void writeStandalone(final ReadMetadata meta, final String whereTo) {
            try {
                final OutputStream outputStream = BucketUtils.createFile(whereTo);
                final OutputStream actualStream = IOUtil.hasBlockCompressedExtension(whereTo)
                        ? new GzipCompressorOutputStream(outputStream) : outputStream;
                final Output output = new Output(actualStream);
                final Kryo kryo = new Kryo();
                final Serializer serializer = new Serializer();
                output.writeString(MAGIC_STRING);
                output.writeString(VERSION_STRING);
                serializer.write(kryo, output, meta);
                output.close();
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(whereTo, ex);
            }
        }

        /**
         * Reads a read-metadata from a file or resource.
         * @param whereFrom the file or resource containing the read-metadata.
         * @throws IllegalArgumentException if {@code whereFrom} is {@code null}.
         * @throws UserException if any problem occurred where deserializing.
         * @return never {@code null}.
         */
        public static ReadMetadata readStandalone(final String whereFrom) {
            try (final InputStream inputStream = BucketUtils.openFile(whereFrom);
                 final Input input = new Input(inputStream)) {
                final Kryo kryo = new Kryo();
                final Serializer serializer = new Serializer();
                final String magicString = input.readString();
                if (!Objects.equals(MAGIC_STRING, magicString)) {
                    throw new UserException.BadInput("Bad file format in " + whereFrom +
                            "; it does not seem to be a valid read-metadata serialization");
                }
                final String versionString = input.readString();
                if (!Objects.equals(VERSION_STRING, versionString)) {
                    throw new UserException.BadInput("Bad file format in " + whereFrom +
                            "; it contains an incompatible version " + versionString + "(expected: " + VERSION_STRING + ")");
                }
                final ReadMetadata result = serializer.read(kryo, input, ReadMetadata.class);
                if (result == null) {
                    throw new UserException.BadInput("Missing read-metadata in " + whereFrom);
                }
                return result;

            } catch (final Exception ex) {
                throw new UserException.CouldNotCreateOutputFile(whereFrom, ex);
            }
        }
    }

    @DefaultSerializer(LibraryRawStatistics.Serializer.class)
    public static final class LibraryRawStatistics {
        private final IntHistogram fragmentSizes;
        private long nReads;
        private long nBases;
        private long totalBaseQuality;

        public LibraryRawStatistics( final int maxTrackedValue ) {
            fragmentSizes = new IntHistogram(maxTrackedValue);
            nReads = nBases = 0;
        }

        private LibraryRawStatistics( final Kryo kryo, final Input input ) {
            fragmentSizes = new IntHistogram.Serializer().read(kryo, input, IntHistogram.class);
            nReads = input.readLong();
            nBases = input.readLong();
            totalBaseQuality = input.readLong();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            new IntHistogram.Serializer().write(kryo, output, fragmentSizes);
            output.writeLong(nReads);
            output.writeLong(nBases);
            output.writeLong(totalBaseQuality);
        }

        public void addRead( final int readLength, final int summedQuals, final int templateLength, final boolean isTemplateLengthTestable ) {
            nReads += 1;
            nBases += readLength;
            totalBaseQuality += summedQuals;
            if ( isTemplateLengthTestable ) {
                fragmentSizes.addObservation(Math.abs(templateLength));
            }
        }

        public long getNReads() { return nReads; }
        public long getNBases() { return nBases; }
        public long getTotalBaseQuality() { return totalBaseQuality; }

        public LibraryStatistics createLibraryStatistics( final long nRefBases ) {
            return new LibraryStatistics(fragmentSizes.getCDF(), nBases, nReads, totalBaseQuality, nRefBases);
        }

        // assumes that the BAM partitioning has no overlap --
        // we'll see each primary line exactly once across all partitions
        public static LibraryRawStatistics reduce( final LibraryRawStatistics stats1,
                                                   final LibraryRawStatistics stats2 ) {
            stats1.fragmentSizes.addObservations(stats2.fragmentSizes);
            stats1.nBases += stats2.nBases;
            stats1.nReads += stats2.nReads;
            stats1.totalBaseQuality += stats2.totalBaseQuality;
            return stats1;
        }

        public static final class Serializer
                extends com.esotericsoftware.kryo.Serializer<LibraryRawStatistics> {
            @Override
            public void write( final Kryo kryo, final Output output,
                               final LibraryRawStatistics libraryRawStatistics ) {
                libraryRawStatistics.serialize(kryo, output);
            }

            @Override
            public LibraryRawStatistics read( final Kryo kryo, final Input input,
                                              final Class<LibraryRawStatistics> klass ) {
                return new LibraryRawStatistics(kryo, input);
            }
        }
    }

    @DefaultSerializer(PartitionStatistics.Serializer.class)
    public static final class PartitionStatistics {
        private final Map<String, LibraryRawStatistics> libraryNameToStatisticsMap;
        private final String firstContig;
        private final int firstLocation;
        private final String lastContig;
        private final int lastLocation;
        private final int span;

        public PartitionStatistics( final Iterator<GATKRead> unfilteredReadItr,
                                    final SVReadFilter filter,
                                    final int maxTrackedFragmentLength,
                                    final Map<String, String> readGroupToLibraryMap ) {
            final Iterator<GATKRead> mappedReadItr = filter.applyFilter(unfilteredReadItr, SVReadFilter::isMappedPrimary);
            libraryNameToStatisticsMap = new HashMap<>();
            if ( !mappedReadItr.hasNext() ) {
                firstContig = lastContig = null;
                firstLocation = lastLocation = -1;
                span = 0;
                return;
            }
            GATKRead mappedRead = mappedReadItr.next();
            firstContig = mappedRead.getContig();
            firstLocation = mappedRead.getStart();
            String currentContig = firstContig;
            int currentSpan = -firstLocation;
            while ( true ) {
                final String libraryName = readGroupToLibraryMap.get(mappedRead.getReadGroup());
                final boolean isTestable = filter.isTemplateLenTestable(mappedRead);
                int summedQuals = 0;
                for ( final byte qual : mappedRead.getBaseQualities() ) {
                    summedQuals += qual;
                }
                libraryNameToStatisticsMap
                        .computeIfAbsent(libraryName, key -> new LibraryRawStatistics(maxTrackedFragmentLength))
                    .addRead(CigarUtils.countAlignedBases(mappedRead.getCigar()), summedQuals, mappedRead.getFragmentLength(), isTestable);
                if ( !mappedReadItr.hasNext() ) break;
                final int endPos = mappedRead.getEnd() + 1;
                mappedRead = mappedReadItr.next();
                if ( !mappedRead.getContig().equals(currentContig) ) {
                    currentSpan += endPos;
                    currentContig = mappedRead.getContig();
                    currentSpan -= mappedRead.getStart();
                }
            }

            lastContig = mappedRead.getContig();
            lastLocation = mappedRead.getEnd() + 1;
            span = currentSpan + lastLocation;
        }

        private PartitionStatistics( final Kryo kryo, final Input input ) {
            final LibraryRawStatistics.Serializer rawStatsSerializer = new LibraryRawStatistics.Serializer();
            int nEntries = input.readInt();
            libraryNameToStatisticsMap = new HashMap<>(SVUtils.hashMapCapacity(nEntries));
            while ( nEntries-- > 0 ) {
                final String libName = input.readString();
                final LibraryRawStatistics rawStats = rawStatsSerializer.read(kryo, input, LibraryRawStatistics.class);
                libraryNameToStatisticsMap.put(libName, rawStats);
            }
            firstContig = input.readString();
            firstLocation = input.readInt();
            lastContig = input.readString();
            lastLocation = input.readInt();
            span = input.readInt();
        }

        public long getNReads() {
            return libraryNameToStatisticsMap.values().stream().mapToLong(LibraryRawStatistics::getNReads).sum();
        }

        public Map<String, LibraryRawStatistics> getLibraryNameToStatisticsMap() {
            return libraryNameToStatisticsMap;
        }

        public String getFirstContig() { return firstContig; }
        public int getFirstLocation() { return firstLocation; }
        public String getLastContig() { return lastContig; }
        public int getLastLocation() { return lastLocation; }
        public int getSpan() { return span; }

        private void serialize( final Kryo kryo, final Output output ) {
            final LibraryRawStatistics.Serializer rawStatsSerializer = new LibraryRawStatistics.Serializer();
            output.writeInt(libraryNameToStatisticsMap.size());
            for ( final Map.Entry<String, LibraryRawStatistics> entry : libraryNameToStatisticsMap.entrySet() ) {
                output.writeString(entry.getKey());
                rawStatsSerializer.write(kryo, output, entry.getValue());
            }
            output.writeString(firstContig);
            output.writeInt(firstLocation);
            output.writeString(lastContig);
            output.writeInt(lastLocation);
            output.writeInt(span);
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
        private final int lastEnd;
        private final int span;
        public final static int UNMAPPED = Integer.MAX_VALUE;

        public PartitionBounds( final int firstContigID, final int firstStart,
                                final int lastContigID, final int lastEnd, final int span ) {
            this.firstContigID = firstContigID;
            this.firstStart = firstStart;
            this.lastContigID = lastContigID;
            this.lastEnd = lastEnd;
            this.span = span;
        }

        private PartitionBounds( final Kryo kryo, final Input input ) {
            this.firstContigID = input.readInt();
            this.firstStart = input.readInt();
            this.lastContigID = input.readInt();
            this.lastEnd = input.readInt();
            this.span = input.readInt();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(firstContigID);
            output.writeInt(firstStart);
            output.writeInt(lastContigID);
            output.writeInt(lastEnd);
            output.writeInt(span);
        }

        public int getFirstContigID() { return firstContigID; }
        public int getFirstStart() { return firstStart; }
        public int getLastContigID() { return lastContigID; }
        public int getLastEnd() { return lastEnd; }
        public int getSpan() { return span; }

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
}
