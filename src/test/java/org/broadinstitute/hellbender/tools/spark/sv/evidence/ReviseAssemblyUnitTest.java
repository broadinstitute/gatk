package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly.Contig;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly.Connection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ReviseAssemblyUnitTest extends GATKBaseTest {
    @Test(groups = "sv")
    void testShadowedContigs() {
        final String[] testSeqs = {
                // contig that will do some shadowing
                "GTATAATTTAGGCTTTCTGTGTCCTGTGTGACCTCCGTCCCATACACGTTTTTGTTTATTTTTACAACCCTTCATACACGTAAAAGCCGCCCATCCTGCCA",
                // not shadowed -- there's some extra 5' sequence
            "TCATGTATAATTTAGGCTTTCTGTGTCCTGTGTGACCTCCGTCCCATACACGTTTTTGTTTATTTTTACAACCCTTCATACACGTAAAAGCCG",
                // not shadowed -- there's some extra 3' sequence
                    "AATTTAGGCTTTCTGTGTCCTGTGTGACCTCCGTCCCATACACGTTTTTGTTTATTTTTACAACCCTTCATACACGTAAAAGCCGCCCATCCTGCCATCAC",
                // not shadowed -- too many mismatches (9 of 99)
                 "CCCCCCCCCGGCTTTCTGTGTCCTGTGTGACCTCCGTCCCATACACGTTTTTGTTTATTTTTACAACCCTTCATACACGTAAAAGCCGCCCATCCTGCC",
                // shadowed, despite 1 mismatch (in 100 bases)
                 "CATAATTTAGGCTTTCTGTGTCCTGTGTGACCTCCGTCCCATACACGTTTTTGTTTATTTTTACAACCCTTCATACACGTAAAAGCCGCCCATCCTGCCA",
                // RC also shadowed
                 "TGGCAGGATGGGCGGCTTTTACGTGTATGAAGGGTTGTAAAAATAAACAAAAACGTGTATGGGACGGAGGTCACACAGGACACAGAAAGCCTAAATTATG"
        };
        final FermiLiteAssembly assembly = new FermiLiteAssembly(
            Arrays.stream(testSeqs)
                .map(seq -> new Contig(seq.getBytes(), null, 1))
                .collect(SVUtils.arrayListCollector(testSeqs.length)));
        final FermiLiteAssembly unshadowedAssembly = FermiLiteAssemblyHandler.removeShadowedContigs(assembly);
        Assert.assertEquals(unshadowedAssembly.getContigs(), assembly.getContigs().subList(0, 4));
    }

    @Test(groups = "sv")
    void testRemoveUnbranched() {
        // test assembly has graph structure A->C, B->C, C->D.  C and D can be replaced with CD.
        final Contig tigA = new Contig("ACACTTTT".getBytes(), null, 1);
        final Contig tigB = new Contig("TGTGTTTT".getBytes(), null, 1);
        final Contig tigC = new Contig("TTTTAAAA".getBytes(), null, 1);
        final Contig tigD = new Contig("AAAAGCGC".getBytes(), null, 1);
        tigA.setConnections(Collections.singletonList(new Connection(tigC, 4, false, false)));
        tigB.setConnections(Collections.singletonList(new Connection(tigC, 4, false, false)));
        final List<Connection> tigCConnections = new ArrayList<>(3);
        tigCConnections.add(new Connection(tigA, 4, true, true));
        tigCConnections.add(new Connection(tigB, 4, true, true));
        tigCConnections.add(new Connection(tigD, 4, false, false));
        tigC.setConnections(tigCConnections);
        tigD.setConnections(Collections.singletonList(new Connection(tigC, 4, true, true)));
        final List<Contig> tigList = new ArrayList<>(5);
        tigList.add(tigA);
        tigList.add(tigB);
        tigList.add(tigC);
        tigList.add(tigD);
        final FermiLiteAssembly assembly = new FermiLiteAssembly(tigList);
        final FermiLiteAssembly noUnbranchedAssembly = FermiLiteAssemblyHandler.removeUnbranchedConnections(assembly);
        Assert.assertEquals(noUnbranchedAssembly.getNContigs(), 3);
        Assert.assertEquals(noUnbranchedAssembly.getContig(0), tigA);
        Assert.assertEquals(noUnbranchedAssembly.getContig(1), tigB);
        Assert.assertEquals(new String(noUnbranchedAssembly.getContig(2).getSequence()), "TTTTAAAAGCGC");

        // now add another contig E->D and make sure we don't create CD (since D now has multiple predecessors)
        final Contig tigE = new Contig("GAGAAAAA".getBytes(), null, 1);
        tigE.setConnections(Collections.singletonList(new Connection(tigD, 4, false, false)));
        tigList.add(tigE);
        final List<Connection> tigDConnections = new ArrayList<>(2);
        tigDConnections.add(tigD.getConnections().get(0));
        tigDConnections.add(new Connection(tigE, 4, true, true));
        tigD.setConnections(tigDConnections);
        final FermiLiteAssembly assembly2 = new FermiLiteAssembly(tigList);
        final FermiLiteAssembly noUnbranchedAssembly2 = FermiLiteAssemblyHandler.removeUnbranchedConnections(assembly2);
        Assert.assertEquals(noUnbranchedAssembly2.getContigs(), assembly2.getContigs());
    }

    @Test(groups = "sv")
    void testNoCrossingUnphasedContigs() {
        // test assembly has the structure A->C, B->C, C->D, C->E.  expanded contigs should be AC, BC, CD, and CE.
        final Contig tigA = new Contig("ACACTTTT".getBytes(), null, 1);
        final Contig tigB = new Contig("TGTGTTTT".getBytes(), null, 1);
        final Contig tigC = new Contig("TTTTAAAA".getBytes(), null, 1);
        final Contig tigD = new Contig("AAAAGCGC".getBytes(), null, 1);
        final Contig tigE = new Contig("AAAATATA".getBytes(), null, 1);
        tigA.setConnections(Collections.singletonList(new Connection(tigC, 4, false, false)));
        tigB.setConnections(Collections.singletonList(new Connection(tigC, 4, false, false)));
        tigD.setConnections(Collections.singletonList(new Connection(tigC, 4, true, true)));
        tigE.setConnections(Collections.singletonList(new Connection(tigC, 4, true, true)));
        final List<Connection> tigCConnections = new ArrayList<>(4);
        tigCConnections.add(new Connection(tigA, 4, true, true));
        tigCConnections.add(new Connection(tigB, 4, true, true));
        tigCConnections.add(new Connection(tigD, 4, false, false));
        tigCConnections.add(new Connection(tigE, 4, false, false));
        tigC.setConnections(tigCConnections);
        final List<Contig> tigList = new ArrayList<>(5);
        tigList.add(tigA);
        tigList.add(tigB);
        tigList.add(tigC);
        tigList.add(tigD);
        tigList.add(tigE);
        final FermiLiteAssembly assembly = new FermiLiteAssembly(tigList);
        final FermiLiteAssembly expandedAssembly = FermiLiteAssemblyHandler.expandAssemblyGraph(assembly);
        Assert.assertEquals(expandedAssembly.getNContigs(), 4);
        Assert.assertEquals(new String(expandedAssembly.getContig(2).getSequence()), "ACACTTTTAAAA");
        Assert.assertEquals(new String(expandedAssembly.getContig(3).getSequence()), "TGTGTTTTAAAA");
        Assert.assertEquals(new String(expandedAssembly.getContig(0).getSequence()), "TTTTAAAAGCGC");
        Assert.assertEquals(new String(expandedAssembly.getContig(1).getSequence()), "TTTTAAAATATA");
    }

    @Test(groups = "sv")
    void testCycleResolution() {
        // test assembly has the structure A->B, B->B.  expanded contig should be ABB.
        final Contig tigA = new Contig("ACACTTTT".getBytes(), null, 1);
        final Contig tigB = new Contig("TTTTTTTT".getBytes(), null, 1);
        tigA.setConnections(Collections.singletonList(new Connection(tigB, 4, false, false)));
        final List<Connection> tigBConnections = new ArrayList<>(3);
        tigBConnections.add(new Connection(tigA, 4, true, true));
        tigBConnections.add(new Connection(tigB, 4, false, false));
        tigBConnections.add(new Connection(tigB, 4, true, true));
        tigB.setConnections(tigBConnections);
        final List<Contig> tigList = new ArrayList<>(2);
        tigList.add(tigA);
        tigList.add(tigB);
        final FermiLiteAssembly assembly = new FermiLiteAssembly(tigList);
        final FermiLiteAssembly expandedAssembly = FermiLiteAssemblyHandler.expandAssemblyGraph(assembly);
        Assert.assertEquals(expandedAssembly.getNContigs(), 1);
        Assert.assertEquals(new String(expandedAssembly.getContig(0).getSequence()), "ACACTTTTTTTTTTTT");
    }

    @Test(groups = "sv")
    void testCycleResolutionAcrossPhasing() {
        // test assembly has the structure A->B, B->B, B->C.  expanded contigs should be AB, BB and BC.
        final Contig tigA = new Contig("ACACTTTT".getBytes(), null, 1);
        final Contig tigB = new Contig("TTTTTTTT".getBytes(), null, 1);
        final Contig tigC = new Contig("TTTTGTGT".getBytes(), null, 1);
        tigA.setConnections(Collections.singletonList(new Connection(tigB, 4, false, false)));
        tigC.setConnections(Collections.singletonList(new Connection(tigB, 4, true, true)));
        final List<Connection> tigBConnections = new ArrayList<>(4);
        tigBConnections.add(new Connection(tigA, 4, true, true));
        tigBConnections.add(new Connection(tigB, 4, false, false));
        tigBConnections.add(new Connection(tigB, 4, true, true));
        tigBConnections.add(new Connection(tigC, 4, false, false));
        tigB.setConnections(tigBConnections);
        final List<Contig> tigList = new ArrayList<>(3);
        tigList.add(tigA);
        tigList.add(tigB);
        tigList.add(tigC);
        final FermiLiteAssembly assembly = new FermiLiteAssembly(tigList);
        final FermiLiteAssembly expandedAssembly = FermiLiteAssemblyHandler.expandAssemblyGraph(assembly);
        Assert.assertEquals(expandedAssembly.getNContigs(), 3);
        Assert.assertEquals(new String(expandedAssembly.getContig(0).getSequence()), "TTTTTTTTTTTT");
        Assert.assertEquals(new String(expandedAssembly.getContig(1).getSequence()), "TTTTTTTTGTGT");
        Assert.assertEquals(new String(expandedAssembly.getContig(2).getSequence()), "ACACTTTTTTTT");
    }

    @Test(groups = "sv")
    void testReviseAssembly() {
        final String fastqFile = getToolTestDataDir() + "/test.fastq";
        final List<SVFastqUtils.FastqRead> readList = SVFastqUtils.readFastqFile(fastqFile);
        final FermiLiteAssembly initialAssembly = new FermiLiteAssembler().createAssembly(readList);
        Assert.assertEquals(initialAssembly.getNContigs(), 9);
        final FermiLiteAssembly revisedAssembly =
                FermiLiteAssemblyHandler.reviseAssembly(initialAssembly, true, true);
        Assert.assertEquals(revisedAssembly.getNContigs(), 6);
    }
}
