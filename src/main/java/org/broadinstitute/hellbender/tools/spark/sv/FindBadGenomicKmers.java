package org.broadinstitute.hellbender.tools.spark.sv;

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

/**
 * SparkTool to identify 63-mers in the reference that occur more than 3 times.
 *
 * Created by tsharpe on 12/4/15.
 */
@CommandLineProgramProperties(summary="Use spark to find the set of high copy number kmers in a reference.",
        oneLineSummary="find ref kmers with high copy number",
        programGroup = SparkProgramGroup.class)
public final class FindBadGenomicKmers extends GATKSparkTool
{
    private static final long serialVersionUID = 1L;
    private static final Long MAX_KMER_FREQ = 3L;
    private static final int KLEN = 63;
    private static final int REF_REC_LEN = 10000;
    private static final long CHUNK_SIZE = 3000000;

    @Argument(doc = "file for ubiquitous kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputFile;

    @Override
    public boolean requiresReference()
    {
        return true;
    }

    // Get the list of high copy number kmers in the reference, and write them to a file.
    @Override
    protected void runTool( final JavaSparkContext ctx )
    {
        SAMFileHeader hdr = getHeaderForReads();
        SAMSequenceDictionary dict = null;
        if ( hdr != null ) dict = hdr.getSequenceDictionary();
        List<Kmer> killList = findBadGenomicKmers(ctx,getReference(),dict);
        final File killFile = new File(outputFile);
        try
        {
            final Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(killFile)));

            for ( Kmer kmer : killList )
            {
                writer.write(kmer.toString(KLEN));
                writer.write('\n');
            }
            writer.close();
        }
        catch ( IOException e )
        { throw new GATKException("Unable to write ubiquitous kmer kill file.", e); }
    }

    // Find high copy number kmers in the reference sequence
    public static List<Kmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                                  final ReferenceMultiSource ref,
                                                  final SAMSequenceDictionary readsDict )
    {
        // Write all the reference sequences to a text file as a series of 10,000-base lines.
        // For each contig, records are written with a K-1 base overlap
        // (i.e., the last 62 bases in line n match the first 62 bases in line n+1) so that we don't miss any kmers
        // due to line breaks -- we can just kmerize each record.
        File refSequenceFile;
        try
        {
            refSequenceFile = File.createTempFile("genomic_sequence","txt");
            refSequenceFile.deleteOnExit();
            final OutputStream os = new BufferedOutputStream(new FileOutputStream(refSequenceFile));

            final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
            if ( dict == null )
                throw new GATKException("No reference dictionary available.");

            final int effectiveRecLen = REF_REC_LEN - KLEN + 1;
            for ( final SAMSequenceRecord rec : dict.getSequences() )
            {
                final String seqName = rec.getSequenceName();
                final int seqLen = rec.getSequenceLength();
                final byte[] bases = ref.getReferenceBases(null,new SimpleInterval(seqName,1,seqLen)).getBases();
                for ( int start = 0; start < seqLen; start += effectiveRecLen )
                {
                    os.write(bases,start,Math.min(REF_REC_LEN,seqLen-start));
                    os.write('\n');
                }
            }

            os.close();
        }
        catch ( IOException e )
        { throw new GATKException("Unable to write reference text file for kmerizing.", e); }

        // Turn that text file into an RDD, and do a classic map/reduce:
        // Kmerize, mapping to a pair <kmer,1>, reduce by summing values by key, filter out <kmer,N> where N < 3,
        // and collect the high frequency kmers back in the driver.
        final int nParts = (int)(refSequenceFile.length()/CHUNK_SIZE + 1);
        final JavaRDD<String> genomicBases = ctx.textFile(refSequenceFile.getPath(),nParts);
        List<Kmer> killList =
                genomicBases
                        .flatMapToPair(seq ->
                        {
                            if ( seq.length() < KLEN )
                                return Collections.emptyList();

                            final List<Tuple2<Kmer,Integer>> kmers = new ArrayList<>(seq.length() - KLEN + 1);
                            StringKmerizer itr = new StringKmerizer(seq,KLEN);
                            while ( itr.hasNext() )
                                kmers.add(new Tuple2<>(itr.next().canonical(KLEN),1));

                            return kmers;
                        })
                        .reduceByKey(( v1, v2 ) -> v1 + v2)
                        .flatMap(kv -> kv._2 > MAX_KMER_FREQ ? Collections.singletonList(kv._1) : Collections.emptyList())
                        .collect();
        if ( !refSequenceFile.delete() )
            throw new GATKException("Failed to delete "+refSequenceFile.getPath());
        return killList;
    }
}
