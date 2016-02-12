package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * SparkTool to identify 63-mers in the reference that occur more than 3 times.
 */
@CommandLineProgramProperties(summary="Find the set of high copy number kmers in a reference.",
        oneLineSummary="find ref kmers with high copy number",
        programGroup = SparkProgramGroup.class)
public final class FindBadGenomicKmersSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    @VisibleForTesting static final Long MAX_KMER_FREQ = 3L;
    private static final int REF_RECORD_LEN = 10000;
    private static final int CHUNK_SIZE = 3000000;

    @Argument(doc = "file for ubiquitous kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputFile;

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
        final List<SVKmer> killList = findBadGenomicKmers(ctx, getReference(), getAuthenticatedGCSOptions(), dict);
        writeKmersToFile(new File(outputFile), killList);
    }

    /** Find high copy number kmers in the reference sequence */
    public static List<SVKmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                                    final ReferenceMultiSource ref,
                                                    final GCSOptions gcsOptions,
                                                    final SAMSequenceDictionary readsDict ) {
        // Generate reference text file.
        final File refSequenceFile = writeRefAsOverlappingRecords(ref, gcsOptions, readsDict);

        // Find the high copy number kmers
        final List<SVKmer> killList = processRefFile(ctx, refSequenceFile);

        if ( !refSequenceFile.delete() ) throw new GATKException("Failed to delete "+refSequenceFile.getPath());

        return killList;
    }

    /** Write kmer list to file. */
    public static void writeKmersToFile( final File kmersFile, final List<SVKmer> kmerList ) {
        try ( final Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(kmersFile))) ) {
            for ( final SVKmer kmer : kmerList ) {
                writer.write(kmer.toString(SVConstants.KMER_SIZE));
                writer.write('\n');
            }
        }
        catch ( final IOException e ) {
            throw new GATKException("Unable to write ubiquitous kmer kill file.", e);
        }
    }

    /**
     * Turn a text file of overlapping records from a reference sequence into an RDD, and do a classic map/reduce:
     * Kmerize, mapping to a pair <kmer,1>, reduce by summing values by key, filter out <kmer,N> where
     * N <= MAX_KMER_FREQ, and collect the high frequency kmers back in the driver.
     */
    @VisibleForTesting static List<SVKmer> processRefFile( final JavaSparkContext ctx, final File refSequenceFile ) {
        final int nParts = (int)(refSequenceFile.length()/CHUNK_SIZE + 1);
        final JavaRDD<String> genomicBases = ctx.textFile("file://" + refSequenceFile.getAbsolutePath(), nParts);
        return genomicBases
                        .flatMapToPair(seq ->
                                SVKmerizer.stream(seq, SVConstants.KMER_SIZE)
                                        .map(kmer -> new Tuple2<>(kmer.canonical(SVConstants.KMER_SIZE),1))
                                        .collect(Collectors.toCollection(() -> new ArrayList<>(CHUNK_SIZE))))
                        .reduceByKey(( v1, v2 ) -> v1 + v2)
                        .filter(kv -> kv._2 > MAX_KMER_FREQ)
                        .map(kv -> kv._1)
                        .collect();
    }

    /**
     * Write all the reference sequences to a text file as a series of REF_RECORD_LEN length lines.
     * For each contig, records are written with a K-1 base overlap
     * (i.e., for K=63, the last 62 bases in line n match the first 62 bases in line n+1) so that we don't miss any
     * kmers due to line breaks -- we can just kmerize each record.
     */
    private static File writeRefAsOverlappingRecords( final ReferenceMultiSource ref,
                                                      final GCSOptions gcsOptions,
                                                      final SAMSequenceDictionary readsDict ) {
        final File refSequenceFile;
        try {
            refSequenceFile = File.createTempFile("genomic_sequence", "txt");
            refSequenceFile.deleteOnExit();
            try ( final OutputStream os = new BufferedOutputStream(new FileOutputStream(refSequenceFile)) ) {
                final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
                if ( dict == null ) throw new GATKException("No reference dictionary available.");

                final int effectiveRecLen = REF_RECORD_LEN - SVConstants.KMER_SIZE + 1;
                for ( final SAMSequenceRecord rec : dict.getSequences() ) {
                    final String seqName = rec.getSequenceName();
                    final int seqLen = rec.getSequenceLength();
                    final SimpleInterval interval = new SimpleInterval(seqName, 1, seqLen);
                    final byte[] bases = ref.getReferenceBases(gcsOptions, interval).getBases();
                    for ( int start = 0; start < seqLen; start += effectiveRecLen ) {
                        os.write(bases, start, Math.min(REF_RECORD_LEN, seqLen - start));
                        os.write('\n');
                    }
                }
            }
        }
        catch ( final IOException e ) {
            throw new GATKException("Unable to write reference text file for kmerizing.", e);
        }

        return refSequenceFile;
    }
}
