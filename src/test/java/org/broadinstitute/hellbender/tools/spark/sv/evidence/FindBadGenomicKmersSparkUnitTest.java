package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerLong;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;

public class FindBadGenomicKmersSparkUnitTest extends GATKBaseTest {

    private static final int KMER_SIZE = StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.KMER_SIZE;
    private static final String REFERENCE_FILE_NAME = hg19MiniReference;

    @Test(groups = "sv")
    public void badKmersTest() {

        final byte[] polyA = new byte[KMER_SIZE]; Arrays.fill(polyA, (byte)'A');
        final byte[] polyC = new byte[KMER_SIZE]; Arrays.fill(polyC, (byte)'C');
        final byte[] polyT = new byte[KMER_SIZE]; Arrays.fill(polyT, (byte)'T');

        int nTimes = FindBadGenomicKmersSpark.MAX_KMER_FREQ;
        final List<byte[]> sequenceChunks = new ArrayList<>(nTimes*2+1);

        // add polyA and polyC the max number of times possible while evading the trigger
        while ( nTimes-- > 0L ) {
            sequenceChunks.add(polyA);
            sequenceChunks.add(polyC);
        }
        // tip polyA over the edge by writing its reverse complement
        sequenceChunks.add(polyT);

        final JavaRDD<byte[]> refRDD = SparkContextFactory.getTestSparkContext().parallelize(sequenceChunks);
        final List<SVKmer> badKmers = FindBadGenomicKmersSpark.collectUbiquitousKmersInReference(KMER_SIZE,
                                                                             Integer.MAX_VALUE,
                                                                             FindBadGenomicKmersSpark.MAX_KMER_FREQ,
                                                                             refRDD);

        // should have just one bad kmer:  polyA
        Assert.assertEquals(badKmers.size(), 1);
        Assert.assertEquals(badKmers.get(0), SVKmerizer.toKmer(polyA,new SVKmerLong()));
    }

    @Test(groups = "sv")
    public void miniRefTest() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReferenceMultiSparkSource ref = new ReferenceMultiSparkSource(
                REFERENCE_FILE_NAME, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(null);
        if ( dict == null ) throw new GATKException("No reference dictionary available.");

        final Map<SVKmer, Long> kmerMap = new LinkedHashMap<>();
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final SimpleInterval interval = new SimpleInterval(rec.getSequenceName(), 1, rec.getSequenceLength());
            final byte[] bases = ref.getReferenceBases(interval).getBases();
            SVKmerizer.canonicalStream(bases, KMER_SIZE, new SVKmerLong())
                    .forEach(kmer -> kmerMap.put(kmer, kmerMap.getOrDefault(kmer, 0L) + 1));
        }
        kmerMap.entrySet().removeIf( x -> x.getValue() <= FindBadGenomicKmersSpark.MAX_KMER_FREQ);

        final List<SVKmer> badKmers =
                FindBadGenomicKmersSpark.findBadGenomicKmers(ctx, KMER_SIZE, Integer.MAX_VALUE, ref, null);
        final Set<SVKmer> badKmerSet = new HashSet<>(badKmers);
        Assert.assertEquals(badKmers.size(), badKmerSet.size());
        Assert.assertEquals(badKmerSet, kmerMap.keySet());
    }
}
