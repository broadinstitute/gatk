package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.BroadcastJoinReadsWithRefBases;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.util.*;

//TODO: docs
// Notes:  GATKSparkTool.getReferenceSequenceDictionary()
//   SparkCommandLineProgram.getAuthenticatedGCSOptions()

public class OxoQScorer {
    /** This scorer assumes that we are only interested in a context that is symmetric on both sides and a fixed number of bases. */
    public final static int EXPECTED_CONTEXT_PADDING = 1;

    public static OxoQBinKey createOxoQBinKey(final byte artifactRef, final byte artifactAlt,
                                   final byte contextPre, final byte contextPost,
                                   final boolean isOriginalMolecule) {
        return new OxoQBinKey(new ArtifactMode(artifactRef, artifactAlt),
                            new Tuple2<>(contextPre, contextPost), isOriginalMolecule);
    }

    /**
     * Reference window function for creating OxoQBinKeys. For each read, returns an interval representing the span of
     * reference bases matching that read with a padding of one on either end.
     *
     * Implemented as a static class rather than an anonymous class or lambda due to serialization issues in spark.
     */
    public final static class OxoQBinReferenceWindowFunction implements SerializableFunction<GATKRead, SimpleInterval> {
        private static final long serialVersionUID = 111L;

        SAMSequenceDictionary samSequenceDictionary;

        public OxoQBinReferenceWindowFunction(final SAMSequenceDictionary samSequenceDictionary) {
            this.samSequenceDictionary = samSequenceDictionary;
        }

        @Override
        public SimpleInterval apply( GATKRead read ) {
            return new SimpleInterval(read).expandWithinContig(EXPECTED_CONTEXT_PADDING, samSequenceDictionary);
        }
    }


    /**
     * Will automatically filter out unmapped reads.
     *
     * Reads MUST have unmapped reads removed.
     * TODO: Finish docs
     * @param reads
     */
    public static double scoreReads(final JavaRDD<GATKRead> reads, final ReferenceMultiSource reference, final JavaSparkContext ctx) {

        // For each read.  Determine if original molecule:
        //  (Orientation==F) xor (Read# == 1)

        // Create 2 x 2 matrix isOriginal x ref alt
        // F_o is "is original molecule"
        //  error of F_o = F_o_alt/(F_o_ref + F_o_alt)
        //  error of R_o = R_o_alt/(R_o_ref + R_o_alt)
        //
        // This is being represented as a Map
        final Map<String, OxoQBinKey> stringOxoQBinKeyMap = createStringOxoQBinKeyMap();

        // Make the mapping act as a cache and broadcast it, since it will be read only.
        Broadcast<Map<String, OxoQBinKey>> bStringOxoQBinKeyMap = ctx.broadcast(stringOxoQBinKeyMap);
        JavaPairRDD<GATKRead, ReferenceBases> readsWithReferenceBases = BroadcastJoinReadsWithRefBases.addBases(reference, reads);
        final Map<OxoQBinKey, Long> binCountsByKey = readsWithReferenceBases
                .flatMap(r -> OxoQScorer.createOxoQBinKeys(r, bStringOxoQBinKeyMap.getValue())).countByValue();
        final Set<OxoQBinKey> allKeys = binCountsByKey.keySet();
        final long fORef = allKeys.stream()
                .filter(k -> k.isOriginalMolecule() && (k.getArtifactMode().getAlt() == k.getArtifactMode().getRef()))
                .mapToLong(k -> binCountsByKey.get(k)).sum();
        final long fOAlt = allKeys.stream()
                .filter(k -> k.isOriginalMolecule() && (k.getArtifactMode().getAlt() != k.getArtifactMode().getRef()))
                .mapToLong(k -> binCountsByKey.get(k)).sum();
        final long rORef = allKeys.stream()
                .filter(k -> !k.isOriginalMolecule() && (k.getArtifactMode().getAlt() == k.getArtifactMode().getRef()))
                .mapToLong(k -> binCountsByKey.get(k)).sum();
        final long rOAlt = allKeys.stream()
                .filter(k -> !k.isOriginalMolecule() && (k.getArtifactMode().getAlt() != k.getArtifactMode().getRef()))
                .mapToLong(k -> binCountsByKey.get(k)).sum();

        final double errorFO = (double)fOAlt/ (double)(fOAlt + fORef);
        final double errorRO = (double)rOAlt/ (double)(rOAlt + rORef);

        return -10 * Math.log10(Math.max(errorFO-errorRO, Math.pow(10, -10)));
    }

    @VisibleForTesting
    static Map<String, OxoQBinKey> createStringOxoQBinKeyMap() {
        // Create a mapping of strings of all possible OxoQBinKeys to a single instance of an OxoQBinKey
        final Map<String, OxoQBinKey> stringOxoQBinKeyMap = new HashMap<>();

        // Go through all possible combinations of context, artifact mode, and isOriginalMolecule
        final String preBases = "TCGA";
        final String postBases = "TCGA";
        final String artifactRefs = "TCGA";
        final String artifactAlts = "TCGA";
        final boolean[] isOriginalMolecules={true, false};

        for (byte preBase: preBases.getBytes()) {
            for (byte postBase: postBases.getBytes()) {
                for (byte artifactRef: artifactRefs.getBytes()) {
                    for (byte artifactAlt: artifactAlts.getBytes()) {
                        for (boolean isOriginalMolecule: isOriginalMolecules) {
                            final String k = createOxoQBinCacheKey(preBase, postBase, artifactRef, artifactAlt, isOriginalMolecule);
                            final OxoQBinKey v = OxoQScorer.createOxoQBinKey(artifactRef, artifactAlt, preBase, postBase, isOriginalMolecule);
                            stringOxoQBinKeyMap.put(k,v);
                        }
                    }
                }
            }
        }
        return stringOxoQBinKeyMap;
    }

    private static String createOxoQBinCacheKey(byte preBase, byte postBase, byte artifactRef, byte artifactAlt, boolean isOriginalMolecule) {
        return String.valueOf(preBase) + String.valueOf(postBase) + String.valueOf(artifactRef) + String.valueOf(artifactAlt) + String.valueOf(isOriginalMolecule);
    }

    @VisibleForTesting
    static List<OxoQBinKey> createOxoQBinKeys(final Tuple2<GATKRead, ReferenceBases> readWithRef, final Map<String, OxoQBinKey> stringOxoQBinKeyMap) {

        final GATKRead read = readWithRef._1();

        final ReferenceBases ref = readWithRef._2();
        final byte[] refBases = ref.getBases();

        List<OxoQBinKey> result = new LinkedList<>();
        final boolean isOriginalMolecule = read.isFirstOfPair() ^ !read.isReverseStrand();

        final AlignmentStateMachine genomeIterator = new AlignmentStateMachine(read);
        genomeIterator.stepForwardOnGenome();

        // Since we need to derive context from the read, we (effectively) skip the first base.
        int readOffset = genomeIterator.getReadOffset(); // position of corresponding reference base in the read
        int genomeOffset = genomeIterator.getGenomeOffset(); // position of corresponding reference base in the read

        while (genomeIterator.stepForwardOnGenome() != null) {

            final int prevReadOffset = readOffset;
            final int prevGenomeOffset = genomeOffset;

            readOffset = genomeIterator.getReadOffset();
            genomeOffset = genomeIterator.getGenomeOffset();

            if (!( ((readOffset - prevReadOffset) > 0) && ((genomeOffset - prevGenomeOffset) > 0) )) {
                continue;
            }

            if (readOffset >= (read.getLength()-1)) {
                break;
            }

            // Corresponding index in the reference bases.
            final int indexForRefBase = genomeIterator.getGenomePosition()-ref.getInterval().getStart();
            if (indexForRefBase >= refBases.length) {
                break;
            }
            final byte modeRef = refBases[indexForRefBase];
            final byte modeAlt = read.getBase(readOffset);
            final String k = OxoQScorer.createOxoQBinCacheKey(read.getBase(readOffset - 1), read.getBase(readOffset + 1), modeRef, modeAlt, isOriginalMolecule);
            final OxoQBinKey v = stringOxoQBinKeyMap.getOrDefault(k, null);
            if (v != null) {
                result.add(v);
            }
        }

        return result;
    }

    // Do not allow instantiation of this class
    private OxoQScorer() {};
}
