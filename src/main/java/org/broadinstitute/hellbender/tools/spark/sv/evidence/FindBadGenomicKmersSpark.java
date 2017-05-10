package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVDUSTFilteredKmerizer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerLong;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

/**
 * SparkTool to identify 63-mers in the reference that occur more than 3 times.
 */
@CommandLineProgramProperties(summary="Find the set of high copy number kmers in a reference.",
        oneLineSummary="find ref kmers with high copy number",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class FindBadGenomicKmersSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    @VisibleForTesting static final int MAX_KMER_FREQ = 3;
    private static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    private static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

    @Argument(doc = "file for ubiquitous kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "kmer size", fullName = "kSize", optional = true)
    private int kSize = SVConstants.KMER_SIZE;

    @Argument(doc = "maximum kmer DUST score", fullName = "kmerMaxDUSTScore")
    private int maxDUSTScore = SVConstants.MAX_DUST_SCORE;

    @Argument(doc = "additional high copy kmers (mitochondrion, e.g.) fasta file name",
            fullName = "highCopyFasta", optional = true)
    private String highCopyFastaFilename;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /** Get the list of high copy number kmers in the reference, and write them to a file. */
    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader hdr = getHeaderForReads();
        SAMSequenceDictionary dict = null;
        if ( hdr != null ) dict = hdr.getSequenceDictionary();
        final PipelineOptions options = getAuthenticatedGCSOptions();
        final ReferenceMultiSource referenceMultiSource = getReference();
        Collection<SVKmer> killList = findBadGenomicKmers(ctx, kSize, maxDUSTScore, referenceMultiSource, options, dict);
        if ( highCopyFastaFilename != null ) {
            killList = uniquify(killList, processFasta(kSize, maxDUSTScore, highCopyFastaFilename, options));
        }

        SVUtils.writeKmersFile(kSize, outputFile, killList);
    }

    /** Find high copy number kmers in the reference sequence */
    @VisibleForTesting
    static List<SVKmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                                    final int kSize,
                                                    final int maxDUSTScore,
                                                    final ReferenceMultiSource ref,
                                                    final PipelineOptions options,
                                                    final SAMSequenceDictionary readsDict ) {
        // Generate reference sequence RDD.
        final JavaRDD<byte[]> refRDD = SVUtils.getRefRDD(ctx, kSize, ref, options, readsDict,
                                                                REF_RECORD_LEN, REF_RECORDS_PER_PARTITION);

        // Find the high copy number kmers
        return processRefRDD(kSize, maxDUSTScore, MAX_KMER_FREQ, refRDD);
    }

    /**
     * A <Kmer,count> pair.
     */
    @DefaultSerializer(KmerAndCount.Serializer.class)
    @VisibleForTesting
    final static class KmerAndCount extends SVKmerLong implements Map.Entry<SVKmer, Integer> {
        private int count;

        KmerAndCount( final SVKmerLong kmer ) { this(kmer,1); }

        KmerAndCount( final SVKmerLong kmer, final int count ) {
            super(kmer);
            this.count = count;
        }

        private KmerAndCount( final Kryo kryo, final Input input ) {
            super(kryo, input);
            count = input.readInt();
        }

        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(count);
        }

        @Override
        public SVKmer getKey() { return new SVKmerLong(this); }
        @Override
        public Integer getValue() { return count; }
        @Override
        public Integer setValue( final Integer count ) {
            final Integer result = this.count;
            this.count = count;
            return result;
        }
        public int grabCount() { return count; }
        public void bumpCount() { count += 1; }
        public void bumpCount( final int extra ) { count += extra; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<KmerAndCount> {
            @Override
            public void write( final Kryo kryo, final Output output, final KmerAndCount kmerAndInterval) {
                kmerAndInterval.serialize(kryo, output);
            }

            @Override
            public KmerAndCount read(final Kryo kryo, final Input input,
                                     final Class<KmerAndCount> klass ) {
                return new KmerAndCount(kryo, input);
            }
        }
    }

    /**
     * Do a map/reduce on an RDD of genomic sequences:
     * Kmerize, mapping to a pair <kmer,1>, reduce by summing values by key, filter out <kmer,N> where
     * N <= MAX_KMER_FREQ, and collect the high frequency kmers back in the driver.
     */
    @VisibleForTesting static List<SVKmer> processRefRDD( final int kSize,
                                                          final int maxDUSTScore,
                                                          final int maxKmerFreq,
                                                          final JavaRDD<byte[]> refRDD ) {
        final int nPartitions = refRDD.getNumPartitions();
        final int hashSize = 2*REF_RECORDS_PER_PARTITION;
        final int arrayCap = REF_RECORDS_PER_PARTITION/100;
        return refRDD
                .mapPartitions(seqItr -> {
                    final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerCounts = new HopscotchMap<>(hashSize);
                    while ( seqItr.hasNext() ) {
                        final byte[] seq = seqItr.next();
                        SVDUSTFilteredKmerizer.stream(seq, kSize, maxDUSTScore, new SVKmerLong())
                                .map(kmer -> kmer.canonical(kSize))
                                .forEach(kmer -> {
                                    final KmerAndCount entry = kmerCounts.find(kmer);
                                    if ( entry == null ) kmerCounts.add(new KmerAndCount((SVKmerLong)kmer));
                                    else entry.bumpCount();
                                });
                    }
                    return kmerCounts.iterator();
                })
                .mapToPair(entry -> new Tuple2<>(entry.getKey(), entry.getValue()))
                .partitionBy(new HashPartitioner(nPartitions))
                .mapPartitions(pairItr -> {
                    final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerCounts =
                            new HopscotchMap<>(hashSize);
                    while ( pairItr.hasNext() ) {
                        final Tuple2<SVKmer, Integer> pair = pairItr.next();
                        final SVKmer kmer = pair._1();
                        final int count = pair._2();
                        KmerAndCount entry = kmerCounts.find(kmer);
                        if ( entry == null ) kmerCounts.add(new KmerAndCount((SVKmerLong)kmer, count));
                        else entry.bumpCount(count);
                    }
                    final List<SVKmer> highFreqKmers = new ArrayList<>(arrayCap);
                    for ( KmerAndCount kmerAndCount : kmerCounts ) {
                        if ( kmerAndCount.grabCount() > maxKmerFreq ) highFreqKmers.add(kmerAndCount.getKey());
                    }
                    return highFreqKmers.iterator();
                })
                .collect();
    }

    @VisibleForTesting static List<SVKmer> processFasta( final int kSize,
                                                         final int maxDUSTScore,
                                                         final String fastaFilename,
                                                         final PipelineOptions options ) {
        try ( BufferedReader rdr = new BufferedReader(new InputStreamReader(BucketUtils.openFile(fastaFilename))) ) {
            final List<SVKmer> kmers = new ArrayList<>((int) BucketUtils.fileSize(fastaFilename));
            String line;
            final StringBuilder sb = new StringBuilder();
            final SVKmer kmerSeed = new SVKmerLong();
            while ( (line = rdr.readLine()) != null ) {
                if ( line.charAt(0) != '>' ) sb.append(line);
                else if ( sb.length() > 0 ) {
                    SVDUSTFilteredKmerizer.stream(sb,kSize,maxDUSTScore,kmerSeed)
                            .map(kmer -> kmer.canonical(kSize)).forEach(kmers::add);
                    sb.setLength(0);
                }
            }
            if ( sb.length() > 0 ) {
                SVDUSTFilteredKmerizer.stream(sb,kSize,maxDUSTScore,kmerSeed)
                        .map(kmer -> kmer.canonical(kSize)).forEach(kmers::add);
            }
            return kmers;
        }
        catch ( IOException ioe ) {
            throw new GATKException("Can't read high copy kmers fasta file "+fastaFilename, ioe);
        }
    }

    private static Collection<SVKmer> uniquify(final Collection<SVKmer> coll1, final Collection<SVKmer> coll2 ) {
        final HopscotchSet<SVKmer> kmers = new HopscotchSet<>(coll1.size() + coll2.size());
        kmers.addAll(coll1);
        kmers.addAll(coll2);
        return kmers;
    }
}
