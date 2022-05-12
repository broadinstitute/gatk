package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.io.File;

public class SharedSequenceMergerUnitTest extends GATKBaseTest {

    @Test
    public void testMergeComplex(){
        final SeqGraph original = new SeqGraph(11);
        final SeqVertex v1 = new SeqVertex("TOP");
        final SeqVertex v2 = new SeqVertex("T");
        final SeqVertex v3 = new SeqVertex("T");
        final SeqVertex v4 = new SeqVertex("C");
        final SeqVertex v5 = new SeqVertex("BOTTOM");

        original.addVertices(v1, v2, v3, v4, v5);
        original.addEdges(v1, v2, v5);
        original.addEdges(v1, v3, v5);
        original.addEdges(v1, v4, v3);

        original.getEdge(v1,v3).setIsRef(true);
        original.getEdge(v3,v5).setIsRef(true);

        original.printGraph(new File("test.mergeComplet.pre.dot" ), 0);
        final MergeCommonSuffices originalMerged = new MergeCommonSuffices(original.clone());
        originalMerged.transformUntilComplete();
        originalMerged.getGraph().printGraph(new File("test.mergeComplet.post.dot" ), 0);
    }
}