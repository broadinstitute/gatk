package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.EdgeFactory;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Like the ReadThreadingGraph but without multiplke copies of non-unique kmers
 */
public final class PlainDeBruijnGraph extends AbstractReadThreadingGraph {

    protected static final Logger logger = LogManager.getLogger(PlainDeBruijnGraph.class);

    private static final long serialVersionUID = 1l;

    PlainDeBruijnGraph(final int kmerSize, final byte minBaseQualityToUseInAssembly) {
        super(kmerSize, false, minBaseQualityToUseInAssembly, 1, -1);
    }

    @Override
    protected void preprocessReads() { }

    @Override
    protected boolean shouldRemoveReadsAfterGraphConstruction() {
        return true;
    }

    @Override
    public boolean isLowQualityGraph() {
        return false;
    }

    @Override
    protected void trackKmer(final Kmer kmer, final MultiDeBruijnVertex newVertex) {
        kmerToVertexMap.putIfAbsent(kmer, newVertex);
    }

    @Override
    protected boolean isThreadingStart(final Kmer kmer, final boolean startThreadingOnlyAtExistingVertex) { return true; }

    @Override
    public SeqGraph toSequenceGraph() {
        buildGraphIfNecessary();
        return super.toSequenceGraph();
    }

    @Override
    protected MultiDeBruijnVertex getNextKmerVertexForChainExtension(final Kmer kmer, final boolean isRef, final MultiDeBruijnVertex prevVertex) {
        return getKmerVertex(kmer, true);
    }

    @Override
    public void postProcessForHaplotypeFinding(final File debugGraphOutputPath, final Locatable refHaplotype) { }

    @Override
    public String toString() {
        return "PlainDeBruijnGraph{kmerSize=" + kmerSize + '}';
    }

}