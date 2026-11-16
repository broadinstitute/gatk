package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.IOException;
import java.util.Arrays;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * One sample's counts of the sites of a block of the genome by reference context and depth bin (see
 * {@link F1R2FilterUtils#binDepth}): of sites without alt reads, and of sites with one alt read by its base and orientation
 */
final class BlockSiteCounts {
    private static final Map<String, Integer> CONTEXT_INDICES = IntStream.range(0, F1R2FilterConstants.NUM_KMERS).boxed()
            .collect(Collectors.toMap(F1R2FilterConstants.ALL_KMERS::get, Function.identity()));

    private final int maxDepth;

    // the depth of each bin, and the bin of each depth up to the maximum
    private final int[] binDepths;
    private final int[] binsOfDepths;

    private final int[][] refSites;             // [context][depth bin]
    private final int[][][][] oneAltReadSites;  // [context][alt base][orientation][depth bin]

    BlockSiteCounts(final int maxDepth) {
        this.maxDepth = maxDepth;
        binDepths = F1R2FilterUtils.getDepthBins(maxDepth);
        binsOfDepths = IntStream.rangeClosed(0, maxDepth)
                .map(depth -> depth == 0 ? 0 : Arrays.binarySearch(binDepths, F1R2FilterUtils.binDepth(depth, maxDepth))).toArray();
        refSites = new int[F1R2FilterConstants.NUM_KMERS][binDepths.length];
        oneAltReadSites = new int[F1R2FilterConstants.NUM_KMERS][Nucleotide.STANDARD_BASES.size()][ReadOrientation.values().length][binDepths.length];
    }

    // the depths of the bins
    int[] getBinDepths() {
        return binDepths;
    }

    void addRefSite(final String context, final int depth) {
        refSites[CONTEXT_INDICES.get(context)][bin(depth)]++;
    }

    void addOneAltReadSite(final String context, final Nucleotide altBase, final ReadOrientation orientation, final int depth) {
        oneAltReadSites[CONTEXT_INDICES.get(context)][altBase.ordinal()][orientation.ordinal()][bin(depth)]++;
    }

    private int bin(final int depth) {
        return binsOfDepths[Math.min(depth, maxDepth)];
    }

    /**
     * Write a row to the ref table for each reference context with sites without alt reads, and a row to the alt table
     * for each context, alt base and orientation with sites with one alt read, of the block that starts at the given
     * position of the contig.  Then reset the counts to zero for the next block.
     */
    void writeRowsAndReset(final TableWriter<F1R2HistogramTable.Row> refTable, final TableWriter<F1R2HistogramTable.Row> altTable,
                           final String contig, final int blockStart) throws IOException {
        for (int n = 0; n < F1R2FilterConstants.NUM_KMERS; n++) {
            final String context = F1R2FilterConstants.ALL_KMERS.get(n);
            writeRowAndReset(refTable, new F1R2HistogramTable.Row(contig, blockStart, context, null, null, binDepths, refSites[n]));
            for (final Nucleotide altBase : Nucleotide.STANDARD_BASES) {
                for (final ReadOrientation orientation : ReadOrientation.values()) {
                    final int[] counts = oneAltReadSites[n][altBase.ordinal()][orientation.ordinal()];
                    writeRowAndReset(altTable, new F1R2HistogramTable.Row(contig, blockStart, context, altBase, orientation, binDepths, counts));
                }
            }
        }
    }

    private static void writeRowAndReset(final TableWriter<F1R2HistogramTable.Row> table, final F1R2HistogramTable.Row row) throws IOException {
        if (Arrays.stream(row.counts()).anyMatch(count -> count > 0)) {
            table.writeRecord(row);
            Arrays.fill(row.counts(), 0);
        }
    }
}
