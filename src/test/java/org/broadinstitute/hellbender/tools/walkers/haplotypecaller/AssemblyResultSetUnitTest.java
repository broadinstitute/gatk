package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.testutils.TestingReadThreadingGraph;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests for {@link AssemblyResultSet}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class AssemblyResultSetUnitTest extends GATKBaseTest {
    private GenomeLocParser genomeLocParser;
    private SAMFileHeader header;

    @BeforeClass
    public void init() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }


    @Test
    public void testEmptyResultSet() {
        final AssemblyResultSet subject = new AssemblyResultSet();

        Assert.assertEquals(subject.getHaplotypeList().size(), 0);
        Assert.assertEquals(subject.getHaplotypeCount(), 0);
        Assert.assertEquals(subject.getReferenceHaplotype(), null);
        Assert.assertEquals(subject.getFullReferenceWithPadding(), null);
        Assert.assertEquals(subject.getPaddedReferenceLoc(), null);
        Assert.assertEquals(subject.getRegionForGenotyping(), null);
        Assert.assertEquals(subject.getUniqueReadThreadingGraph(10), null);
        Assert.assertFalse(subject.hasMultipleKmerSizes());
    }

    @Test
    public void testAddReferenceHaplotype() {

        final Haplotype ref = new Haplotype("ACGT".getBytes(),true);
        ref.setGenomeLocation(genomeLocParser.createGenomeLoc("1",1,ref.length() + 1 ));
        final AssemblyResultSet subject = new AssemblyResultSet();

        Assert.assertTrue(subject.add(ref));
        Assert.assertFalse(subject.add(ref));

        Assert.assertEquals(subject.getReferenceHaplotype(), ref);
        Assert.assertEquals(subject.getHaplotypeCount(), 1);
        Assert.assertEquals(subject.getHaplotypeList().size(), 1);
    }

    @Test(dataProvider="assemblyResults")
    public void testAddManyHaplotypes(final List<AssemblyResult> assemblyResults,
                                      final List<List<Haplotype>> haplotypes) {
        final AssemblyResultSet subject = new AssemblyResultSet();
        for (int i = 0; i < haplotypes.size(); i++) {
            final int haplotypeCountBefore = subject.getHaplotypeCount();
            final List<Haplotype> haplos = haplotypes.get(i);
            final AssemblyResult ar = assemblyResults.get(i);
            for (final Haplotype h : haplos) {
                Assert.assertTrue(subject.add(h, ar));
                Assert.assertFalse(subject.add(h, ar));
                if (h.isReference())
                    Assert.assertEquals(subject.getReferenceHaplotype(), h);
            }
            final int haplotypeCountAfter = subject.getHaplotypeCount();
            Assert.assertEquals(haplos.size(), haplotypeCountAfter - haplotypeCountBefore);
            Assert.assertTrue(subject.getMaximumKmerSize() >= ar.getKmerSize());
            Assert.assertTrue(subject.getMinimumKmerSize() <= ar.getKmerSize());
            Assert.assertEquals(subject.getUniqueReadThreadingGraph(ar.getKmerSize()), ar.getThreadingGraph());
        }
    }

    @Test(dataProvider="trimmingData")
    public void testTrimTo(final Map<Haplotype,AssemblyResult> haplotypesAndResultSets, final AssemblyRegion original) {
        final AssemblyResultSet subject = new AssemblyResultSet();
        for (final Map.Entry<Haplotype,AssemblyResult> entry : haplotypesAndResultSets.entrySet())
            subject.add(entry.getKey(),entry.getValue());
        subject.setRegionForGenotyping(original);
        final SimpleInterval originalLocation = original.getPaddedSpan();
        final int length = originalLocation.size();
        final SimpleInterval newLocation = new SimpleInterval(originalLocation.getContig(),
                                                  originalLocation.getStart() + length / 2,
                                                  originalLocation.getEnd() - length / 2);
        final AssemblyRegion newRegion = original.trim(newLocation, newLocation);

        final Map<Haplotype,Haplotype> originalHaplotypesByTrimmed = new HashMap<>(haplotypesAndResultSets.size());
        for (final Haplotype h : haplotypesAndResultSets.keySet())
            originalHaplotypesByTrimmed.put(h.trim(newRegion.getPaddedSpan()), h);

        final AssemblyResultSet trimmed = subject.trimTo(newRegion);

        for (final Haplotype h : trimmed.getHaplotypeList()) {
            Assert.assertEquals(h.getGenomeLocation(), newLocation);
            Assert.assertEquals(h.getBases().length, newLocation.size());
        }
    }

    @DataProvider(name = "swParameters")
    public Object[][] swParameters() {
        return new Object[][] {
                {SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS, true},
                {SmithWatermanAlignmentConstants.STANDARD_NGS, true},
                {SmithWatermanAlignmentConstants.ORIGINAL_DEFAULT, true},
                {new SWParameters(200, -600, -260, -11), false},
                // a mismatch ties with a one-base deletion plus a one-base insertion
                {new SWParameters(200, -520, -260, -11), false}
        };
    }

    @Test(dataProvider = "swParameters")
    public void testMismatchBeatsDeletionPlusInsertion(final SWParameters parameters, final boolean expected) {
        Assert.assertEquals(AssemblyResultSet.mismatchBeatsDeletionPlusInsertion(parameters), expected);
    }

    // A given SNV substituted into the reference haplotype gets the cigar that Smith-Waterman would find, which is all
    // matches whenever a mismatch beats a deletion plus an insertion
    @Test(dataProvider = "swParameters")
    public void testGivenSnvCigar(final SWParameters parameters, final boolean mismatchBeatsDeletionPlusInsertion) {
        final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);
        final byte[] refBases = new RandomDNA(13).nextBases(100);
        for (int offset = 0; offset < refBases.length; offset++) {
            for (final byte altBase : BaseUtils.BASES) {
                if (altBase == refBases[offset]) {
                    continue;
                }
                final AssemblyResultSet resultSet = referenceOnlyResultSet(refBases);
                final Event snv = new Event("1", REFERENCE_HAPLOTYPE_START + offset, Allele.create(refBases[offset], true), Allele.create(altBase));
                resultSet.addGivenAlleles(Collections.singletonList(snv), 0, aligner, parameters);

                Assert.assertEquals(resultSet.getHaplotypeCount(), 2);
                final Haplotype altHaplotype = resultSet.getHaplotypeList().stream().filter(h -> !h.isReference()).findFirst().get();
                Assert.assertEquals(altHaplotype.getBases()[offset], altBase);
                Assert.assertEquals(altHaplotype.getCigar(), CigarUtils.calculateCigar(refBases, altHaplotype.getBases(), aligner, parameters, SWOverhangStrategy.INDEL));
                if (mismatchBeatsDeletionPlusInsertion) {
                    Assert.assertEquals(altHaplotype.getCigar().toString(), refBases.length + "M");
                }
            }
        }
    }

    private static final int REFERENCE_HAPLOTYPE_START = 2;

    // an assembly result set with only a reference haplotype, as when we skip assembly, following AssemblyBasedCallerUtilsUnitTest
    private static AssemblyResultSet referenceOnlyResultSet(final byte[] refBases) {
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        refHaplotype.setAlignmentStartHapwrtRef(1);
        refHaplotype.setCigar(new Cigar(Collections.singletonList(new CigarElement(refBases.length, CigarOperator.M))));
        refHaplotype.setGenomeLocation(new SimpleInterval("1", REFERENCE_HAPLOTYPE_START, REFERENCE_HAPLOTYPE_START + refBases.length - 1));

        final AssemblyResultSet resultSet = new AssemblyResultSet();
        resultSet.setPaddedReferenceLoc(new SimpleInterval("1", 1, REFERENCE_HAPLOTYPE_START + refBases.length - 1));
        resultSet.setFullReferenceWithPadding(ArrayUtils.addAll(new byte[] {'A'}, refBases));
        resultSet.add(refHaplotype);
        return resultSet;
    }

    @DataProvider(name="trimmingData")
    public Iterator<Object[]> trimmingData() {
        final AssemblyRegion activeRegion = new AssemblyRegion(new SimpleInterval("1",1000,1100), 25, header);
        final int length = activeRegion.getPaddedSpan().size();
        final RandomDNA rnd = new RandomDNA(13); // keep it prepoducible by fixing the seed to lucky 13.
        final AssemblyRegionTestDataSet actd = new AssemblyRegionTestDataSet(10,new String(rnd.nextBases(length)),new String[] {
                "Civar:*1T*" }, new String[0], new byte[0], new byte[0], new byte[0]);

        final List<Haplotype> haplotypes = actd.haplotypeList();
        for (final Haplotype h : haplotypes)
            h.setGenomeLocation(activeRegion.getPaddedSpan());

        final ReadThreadingGraph rtg = new ReadThreadingGraph(10);
        for (final Haplotype h : haplotypes)
            rtg.addSequence("seq-" + Math.abs(h.hashCode()), h.getBases(), h.isReference());
        final SeqGraph seqGraph = rtg.toSequenceGraph();
        final AssemblyResult ar = new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION,seqGraph, rtg);
        final Map<Haplotype,AssemblyResult> result =
                new HashMap<>();
        for (final Haplotype h : haplotypes)
            result.put(h,ar);
        return Collections.singleton(new Object[]{result, activeRegion}).iterator();

    }




    @DataProvider(name="assemblyResults")
    public Iterator<Object[]> assemblyResults() {
        final int size = THREE_KS_GRAPH_AND_HAPLOTYPES.length * (1 + TEN_KS_GRAPH_AND_HAPLOTYPES.length);
        final Object[][] result = new Object[size][];

        for (int i = 0; i < THREE_KS_GRAPH_AND_HAPLOTYPES.length; i++) {
            final ReadThreadingGraph rtg = new TestingReadThreadingGraph((String) THREE_KS_GRAPH_AND_HAPLOTYPES[i][0]);
            final AssemblyResult ar = new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION,rtg.toSequenceGraph(), rtg);
            final Object[] haplotypeStrings = (Object[]) THREE_KS_GRAPH_AND_HAPLOTYPES[i][1];
            final Haplotype[] haplotypes = new Haplotype[haplotypeStrings.length];
            for (int j = 0; j < haplotypeStrings.length; j++) {
                haplotypes[j] = new Haplotype(((String)haplotypeStrings[j]).getBytes(),j == 0);
                haplotypes[j].setGenomeLocation(genomeLocParser.createGenomeLoc("1",1,haplotypes[j].length() + 1));
            }
            result[i] = new Object[] { Collections.singletonList(ar), Arrays.asList(Arrays.asList(haplotypes))};
            for (int j = 0; j < TEN_KS_GRAPH_AND_HAPLOTYPES.length; j++) {
                final ReadThreadingGraph rtg10 = new TestingReadThreadingGraph((String) TEN_KS_GRAPH_AND_HAPLOTYPES[j][0]);
                final AssemblyResult ar10 = new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION,rtg10.toSequenceGraph(),rtg10);
                final Object[] haplotypeStrings10 = (Object[]) TEN_KS_GRAPH_AND_HAPLOTYPES[j][1];
                final Haplotype[] haplotype10 = new Haplotype[haplotypeStrings10.length];
                for (int k = 0; k < haplotypeStrings10.length; k++) {
                    haplotype10[k] = new Haplotype(((String)haplotypeStrings10[k]).getBytes(),false);
                    haplotype10[k].setGenomeLocation(genomeLocParser.createGenomeLoc("1", 1, haplotype10[k].length() + 1));
                }

                result[THREE_KS_GRAPH_AND_HAPLOTYPES.length + i * TEN_KS_GRAPH_AND_HAPLOTYPES.length + j] = new Object[] { Arrays.asList(ar, ar10),
                        Arrays.asList(Arrays.asList(haplotypes), Arrays.asList(haplotype10)) };
            }
        }
        return Arrays.asList(result).iterator();
    }


    private static final Object[][] THREE_KS_GRAPH_AND_HAPLOTYPES = new Object[][] {
            {"[ks=3]{REF: ACT}",new Object[] {"ACT"}},
            {"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                    "{ (3) -> A -> G ->          (2) }" +
                    "{  (1) -> A -> G ->  (2) }",new Object[] {"ACTTGA","ACTAGGA","ACTTAGGA"}},
            {"[ks=3]{REF: ACT -> C(1) -> G}{ACT -> C(1) -> G}{ACT -> C(1) -> G}", new Object[] {"ACTCG"}} ,
            {"[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> C -> G -> T }" +
                    "{A(1) -> T -> A(2) }", new Object[] {"ACTAGACGT","ACTATACGT"}}  ,
            {"[ks=3]{REF: ACT -> A -> T(2) -> C -> A -> G -> T -> A -> C -> G -> T -> A(1) -> T}" +
                    "{ ACT -> A -> T(2) -> C -> T -> A -> C -> G -> T -> A(1) -> T}",
                           new Object[] {"ACTATCAGTACGTAT","ACTATCTACGTAT"}} ,
            {"[ks=3]{REF: ACT -> A -> T    -> C -> A -> G -> T -> A -> C -> G -> T -> A    -> T}",
                           new Object[] {"ACTATCAGTACGTAT"}},
            {"[ks=3]{REF: ACT -> A -> T(1) }" +
                    "{ ACT -> A -> T(1) }", new Object[] {"ACTAT"}},
            {"[ks=3]{REF: TTT -> A(1) -> C -> T(2)}{ A(1) -> T(2) } ", new Object[] {"TTTACT","TTTAT"}}
    };

    private static final Object[][] TEN_KS_GRAPH_AND_HAPLOTYPES = new Object[][] {
            {"[ks=10]{ACTAGTAAAT -> A -> T -> A -> A -> T -> A", new Object[] {"ACTAGTAAATATAATA"}},
            {"[ks=10]{ATAGTAATAA(1) -> A -> C -> T -> A(2) -> C}{ (1) -> C -> C -> C -> A(2) -> C}",
                    new Object[] {"ATAGTAATAAACTAC","ATAGTAATAACCCAC"}},

    };

}
