package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Identifies sequences that occur at high frequency in a reference
 *
 * <p>Search the reference for kmers (fixed-length substrings) that occur more than a specified number of times,
 * and list them to an output file.  The resulting output file is appropriate for use as the --kmers-to-ignore
 * input file by the StructuralVariationDiscoveryPipelineSpark tool, which will ignore these kmers when trying
 * to produce candidate reads for local assemblies.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A reference.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A text file describing the ubiquitous kmers.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk FindBadGenomicKmersSpark \
 *     -R reference.fasta \
 *     -O kmers_to_ignore.txt
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Identifies sequences that occur at high frequency in a reference",
        summary =
        "Search the reference for kmers (fixed-length substrings) that occur more than a specified number of times," +
        " and list them to an output file.  The resulting output file is appropriate for use as the --kmers-to-ignore" +
        " input file by the StructuralVariationDiscoveryPipelineSpark tool, which will ignore these kmers when trying" +
        " to produce candidate reads for local assemblies.",
        programGroup = ReferenceProgramGroup.class)
public final class FindBadGenomicKmersSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    public static final int REF_RECORD_LEN = 10000;

    public static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

    public static final int MAX_KMER_FREQ = 3;

    @Argument(doc = "file for ubiquitous kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "kmer size", fullName = "k-size")
    private int kSize = StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.KMER_SIZE;

    @Argument(doc = "maximum kmer DUST score", fullName = "kmer-max-dust-score")
    private int maxDUSTScore = StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.MAX_DUST_SCORE;

    @Argument(doc = "additional high copy kmers (mitochondrion, e.g.) fasta file name",
            fullName = "high-copy-fasta", optional = true)
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
        final ReferenceMultiSparkSource referenceMultiSource = getReference();
        Collection<SVKmer> killList = findBadGenomicKmers(ctx, kSize, maxDUSTScore, referenceMultiSource, dict);
        if ( highCopyFastaFilename != null ) {
            killList = SVUtils.uniquify(killList, processFasta(kSize, maxDUSTScore, highCopyFastaFilename));
        }

        SVFileUtils.writeKmersFile(outputFile, kSize, killList);
    }

    /** Find high copy number kmers in the reference sequence */
    @VisibleForTesting
    static List<SVKmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                             final int kSize,
                                             final int maxDUSTScore,
                                             final ReferenceMultiSparkSource ref,
                                             final SAMSequenceDictionary readsDict ) {
        // Generate reference sequence RDD.
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
        if ( dict == null ) throw new GATKException("No reference dictionary available");
        final JavaRDD<byte[]> refRDD = SVReferenceUtils.getReferenceBasesRDD(ctx, kSize, ref, dict,
                                                        REF_RECORD_LEN, REF_RECORDS_PER_PARTITION);

        // Find the high copy number kmers
        return collectUbiquitousKmersInReference(kSize, maxDUSTScore, MAX_KMER_FREQ, refRDD);
    }

    /**
     * Do a map/reduce on an RDD of genomic sequences:
     * Kmerize, mapping to a pair <kmer,1>, reduce by summing values by key, filter out <kmer,N> where
     * N <= MAX_KMER_FREQ, and collect the high frequency kmers back in the driver.
     */
    @VisibleForTesting
    static List<SVKmer> collectUbiquitousKmersInReference(final int kSize,
                                                          final int maxDUSTScore,
                                                          final int maxKmerFreq,
                                                          final JavaRDD<byte[]> refRDD) {
        Utils.nonNull(refRDD, "reference bases RDD is null");
        Utils.validateArg(kSize > 0, "provided kmer size is non positive");
        Utils.validateArg(maxDUSTScore > 0, "provided DUST filter score is non positive");
        Utils.validateArg(maxKmerFreq > 0, "provided kmer frequency is non positive");

        final int nPartitions = refRDD.getNumPartitions();
        final int hashSize = 2*REF_RECORDS_PER_PARTITION;
        return refRDD
                .mapPartitions(seqItr -> {
                    final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerCounts = new HopscotchMap<>(hashSize);
                    while ( seqItr.hasNext() ) {
                        final byte[] seq = seqItr.next();
                        SVDUSTFilteredKmerizer.canonicalStream(seq, kSize, maxDUSTScore, new SVKmerLong())
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
                    final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerCounts = new HopscotchMap<>(hashSize);
                    while ( pairItr.hasNext() ) {
                        final Tuple2<SVKmer, Integer> pair = pairItr.next();
                        final SVKmer kmer = pair._1();
                        final int count = pair._2();
                        KmerAndCount entry = kmerCounts.find(kmer);
                        if ( entry == null ) kmerCounts.add(new KmerAndCount((SVKmerLong)kmer, count));
                        else entry.bumpCount(count);
                    }
                    return kmerCounts.stream()
                            .filter(kmerAndCount -> kmerAndCount.grabCount() > maxKmerFreq)
                            .map(KmerAndCount::getKey).iterator();
                })
                .collect();
    }

    @VisibleForTesting static List<SVKmer> processFasta( final int kSize,
                                                         final int maxDUSTScore,
                                                         final String fastaFilename) {
        try ( BufferedReader rdr = new BufferedReader(new InputStreamReader(BucketUtils.openFile(fastaFilename))) ) {
            final List<SVKmer> kmers = new ArrayList<>((int) BucketUtils.fileSize(fastaFilename));
            String line;
            final StringBuilder sb = new StringBuilder();
            final SVKmer kmerSeed = new SVKmerLong();
            while ( (line = rdr.readLine()) != null ) {
                if ( line.charAt(0) != '>' ) sb.append(line);
                else if ( sb.length() > 0 ) {
                    SVDUSTFilteredKmerizer.canonicalStream(sb,kSize,maxDUSTScore,kmerSeed).forEach(kmers::add);
                    sb.setLength(0);
                }
            }
            if ( sb.length() > 0 ) {
                SVDUSTFilteredKmerizer.canonicalStream(sb,kSize,maxDUSTScore,kmerSeed).forEach(kmers::add);
            }
            return kmers;
        }
        catch ( IOException ioe ) {
            throw new GATKException("Can't read high copy kmers fasta file "+fastaFilename, ioe);
        }
    }

}
