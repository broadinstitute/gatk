package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.*;
import java.util.*;

/**
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

    public static List<Kmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                           final ReferenceMultiSource ref,
                                           final SAMSequenceDictionary readsDict )
    {
        final File hadoopFile = new File("genomicbases.txt");
        hadoopFile.deleteOnExit();
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
        if ( dict == null )
            throw new GATKException("No reference dictionary available.");
        final int effectiveRecLen = REF_REC_LEN - KLEN + 1;
        try
        {
            final OutputStream os = new BufferedOutputStream(new FileOutputStream(hadoopFile));
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

        final JavaRDD<String> genomicBases = ctx.textFile(hadoopFile.getPath(),300);
        List<Kmer> killList =
            genomicBases
                .flatMapToPair(seq ->
                {
                    if ( seq.length() < KLEN )
                        return Collections.emptyList();

                    final List<Tuple2<Kmer,Integer>> kmers = new ArrayList<>(seq.length() - KLEN + 1);
                    Kmerizer itr = new Kmerizer(seq,KLEN);
                    while ( itr.hasNext() )
                        kmers.add(new Tuple2<>(itr.next().canonical(KLEN),1));

                    return kmers;
                })
                .reduceByKey(( v1, v2 ) -> v1 + v2)
                .flatMap(kv -> kv._2 > MAX_KMER_FREQ ? Collections.singletonList(kv._1) : Collections.emptyList())
                .collect();
        if ( !hadoopFile.delete() )
            throw new GATKException("Failed to delete "+hadoopFile.getPath());
        return killList;
    }

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    public boolean requiresReference()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx )
    {
        List<Kmer> killList = findBadGenomicKmers(ctx,getReference(),getHeaderForReads().getSequenceDictionary());
        final File killFile = new File("ubiquitous_kmers.txt");
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
        { throw new GATKException("Unable to write kmer kill file.", e); }
    }
}
