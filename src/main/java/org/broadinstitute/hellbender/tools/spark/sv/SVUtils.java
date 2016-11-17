package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.util.*;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * Useful scraps of this and that.
 */
public final class SVUtils {

    /**
     * Read a file of kmers.
     * Each line must be exactly SVConstants.KMER_SIZE characters long, and must match [ACGT]*.
     */
    public static Set<SVKmer> readKmersFile(final int kSize, final String kmersFile, final PipelineOptions popts ) {
        final Set<SVKmer> kmers;

        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(kmersFile, popts))) ) {
            final long fileLength = BucketUtils.fileSize(kmersFile, popts);
            kmers = new HopscotchSet<>((int)(fileLength/(kSize+1)));
            String line;
            while ( (line = rdr.readLine()) != null ) {
                if ( line.length() != kSize ) {
                    throw new GATKException("SVKmerLong kill set contains a line of length " + line.length() +
                            " but we were expecting K=" + kSize);
                }

                final SVKmerizer kmerizer = new SVKmerizer(line, kSize, new SVKmerLong(kSize));
                if ( !kmerizer.hasNext() ) {
                    throw new GATKException("Unable to kmerize the kmer kill set string '" + line + "'.");
                }

                kmers.add(kmerizer.next());
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read kmers from "+kmersFile, ioe);
        }

        return kmers;
    }

    /** Write kmers to file. */
    public static <KType extends SVKmer> void writeKmersFile(final int kSize, final String kmersFile, final PipelineOptions popts,
                                                             final Collection<KType> kmers ) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(kmersFile, popts))) ) {
            for ( final KType kmer : kmers ) {
                writer.write(kmer.toString(kSize));
                writer.write('\n');
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to write kmers to "+kmersFile, ioe);
        }
    }

    /** Read intervals from file. */
    public static List<SVInterval> readIntervalsFile( final String intervalsFile, final PipelineOptions popts,
                                                      final Map<String, Integer> contigNameMap ) {
        final List<SVInterval> intervals;
        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(intervalsFile, popts))) ) {
            final int INTERVAL_FILE_LINE_LENGTH_GUESS = 25;
            final long sizeGuess = BucketUtils.fileSize(intervalsFile, popts)/INTERVAL_FILE_LINE_LENGTH_GUESS;
            intervals = new ArrayList<>((int)sizeGuess);
            String line;
            int lineNo = 0;
            while ( (line = rdr.readLine()) != null ) {
                ++lineNo;
                final String[] tokens = line.split("\t");
                if ( tokens.length != 3 ) {
                    throw new GATKException("Interval file "+intervalsFile+" line "+
                                            lineNo+" did not contain 3 columns: "+line);
                }
                try {
                    final Integer contigId = contigNameMap.get(tokens[0]);
                    if ( contigId == null ) throw new GATKException("contig name "+tokens[0]+" not in dictionary");
                    final int start = Integer.valueOf(tokens[1])-1;
                    final int end = Integer.valueOf(tokens[2]);
                    intervals.add(new SVInterval(contigId, start, end));
                }
                catch ( final Exception e ) {
                    throw new GATKException("Unable to parse interval file "+intervalsFile+" line "+lineNo+": "+line, e);
                }
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read intervals from "+intervalsFile, ioe);
        }
        return intervals;
    }

    /** Write intervals to a file. */
    public static void writeIntervalsFile( final String intervalsFile, final PipelineOptions popts,
                                           final Collection<SVInterval> intervals, final List<String> contigNames ) {
        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(intervalsFile, popts)))) {
            for (final SVInterval interval : intervals) {
                final String seqName = contigNames.get(interval.getContig());
                writer.write(seqName + "\t" + (interval.getStart()+1) + "\t" + interval.getEnd() + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write intervals file " + intervalsFile, ioe);
        }
    }

    /** return a good initialCapacity for a HashMap that will hold a given number of elements */
    public static int hashMapCapacity( final int nElements )
    {
        return (int)((nElements*4L)/3) + 1;
    }

    /** count the number of items available from an iterator */
    public static <T> int iteratorSize( final Iterator<T> itr ) {
        int result = 0;
        while ( itr.hasNext() ) { result += 1; itr.next(); }
        return result;
    }

    /**
     * Provides a stream collector that will collect items into an array list with a given initial capacity.
     */
    public static <T> Collector<T, ?, ArrayList<T>> arrayListCollector(final int size) {
        return Collectors.toCollection( () -> new ArrayList<>(size));
    }

    /**
     * Create an RDD from the reference sequences.
     * The reference sequences are transformed into a single, large collection of byte arrays. The collection is then
     * parallelized into an RDD.
     * Each contig that exceeds a size given by REF_RECORD_LEN is broken into a series of REF_RECORD_LEN chunks with a
     * K-1 base overlap between successive chunks. (I.e., for K=63, the last 62 bases in chunk n match the first 62
     * bases in chunk n+1) so that we don't miss any kmers due to the chunking -- we can just kmerize each record
     * independently.
     */
    public static JavaRDD<byte[]> getRefRDD(final JavaSparkContext ctx,
                                             final int kSize,
                                             final ReferenceMultiSource ref,
                                             final PipelineOptions options,
                                             final SAMSequenceDictionary readsDict,
                                             final int ref_record_len,
                                             final int ref_records_per_partition) {
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
        if ( dict == null ) throw new GATKException("No reference dictionary available");

        final int effectiveRecLen = ref_record_len - kSize + 1;
        final List<byte[]> sequenceChunks = new ArrayList<>();
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String seqName = rec.getSequenceName();
            final int seqLen = rec.getSequenceLength();
            final SimpleInterval interval = new SimpleInterval(seqName, 1, seqLen);
            try {
                final byte[] bases = ref.getReferenceBases(options, interval).getBases();
                for ( int start = 0; start < seqLen; start += effectiveRecLen ) {
                    sequenceChunks.add(Arrays.copyOfRange(bases, start, Math.min(start+ref_record_len, seqLen)));
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't get reference sequence bases for " + interval, ioe);
            }
        }

        return ctx.parallelize(sequenceChunks, sequenceChunks.size()/ref_records_per_partition+1);
    }
}
