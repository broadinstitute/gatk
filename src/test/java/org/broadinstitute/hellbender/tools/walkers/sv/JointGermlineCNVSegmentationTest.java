package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class JointGermlineCNVSegmentationTest extends CommandLineProgramTest {

    private final JointGermlineCNVSegmentation walker = new JointGermlineCNVSegmentation();
    private SampleDB pedigree;
    private Set<String> allosomalContigs;
    private Set<String> badAllosomes;
    private final int refCopyNumber = 2;

    @Test
    public void testResolveVariantContexts() {
        VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("sample1");

        //make first variant sample 1
        genotypeBuilder.alleles(Arrays.asList(Allele.REF_N, GATKSVVCFConstants.DEL_ALLELE)).attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1);
        builder.chr("1").start(1000).stop(2000).alleles(Arrays.asList(Allele.REF_N, GATKSVVCFConstants.DEL_ALLELE)).genotypes(genotypeBuilder.make());
        final VariantContext sample1_var1 = builder.make();

        //then second variant in sample 1, to update the sample-copy number map
        builder.start(3000).stop(5000).alleles(Arrays.asList(Allele.REF_N, GATKSVVCFConstants.DEL_ALLELE)).genotypes(genotypeBuilder.make());
        final VariantContext sample1_var2 = builder.make();

        //then an overlapping variant in sample 2
        final GenotypeBuilder gb_sample2 = new GenotypeBuilder("sample2");
        gb_sample2.alleles(Arrays.asList(Allele.REF_N, GATKSVVCFConstants.DEL_ALLELE)).attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1);
        builder.start(4000).stop(5000).alleles(Arrays.asList(Allele.REF_N, GATKSVVCFConstants.DEL_ALLELE)).genotypes(genotypeBuilder.make());
        final VariantContext sample2 = builder.make();

        final SortedSet<String> sampleSet = new TreeSet<>();
        sampleSet.addAll(Arrays.asList("sample1","sample2"));
        final List<VariantContext> resolvedVCs = walker.resolveVariantContexts(Collections.emptySet(), 2, null, sampleSet,
                Arrays.asList(sample1_var1, sample1_var2, sample2));

        Assert.assertEquals(Integer.parseInt(resolvedVCs.get(2).getGenotype("sample1").getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 1);
    }

    @Test
    public void testUpdateGenotypes() {
        final SortedSet<String> allSamples = new TreeSet<>();
        allSamples.addAll(pedigree.getSamples().stream().map(Sample::getID).collect(Collectors.toSet()));
        final Map<String, JointGermlineCNVSegmentation.CopyNumberAndEndRecord> sampleCopyNumbers = new LinkedHashMap<>();
        final VariantContext vcOut = JointGermlineCNVSegmentation.updateGenotypes(allosomalContigs, refCopyNumber, pedigree, allSamples,
                JointGermlineCNVSegmentation.buildVariantContext(SVTestUtils.call1, ReferenceUtils.createReferenceReader(new GATKPath(GATKBaseTest.hg38Reference))),
                sampleCopyNumbers);
        Assert.assertEquals(vcOut.getStart(), SVTestUtils.call1.getPositionA());
    }

    @BeforeMethod
    private void initializePedigree() {
        final SampleDBBuilder pedBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        pedBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(new GATKPath(getToolTestDataDir() + "overlapping.ped")));
        pedigree = pedBuilder.getFinalSampleDB();

        allosomalContigs = new LinkedHashSet<>();
        allosomalContigs.addAll(Arrays.asList("X","Y","chrX","chrY"));

        badAllosomes = new LinkedHashSet<>();
        badAllosomes.addAll(Arrays.asList("X","Y","chrX","chrY","Z","chrZ"));
    }

    @DataProvider
    public Object[][] samplesForPloidyQuery() {
        return new Object[][]{
                {"NA00000", "X", null, 2},
                {"NA00000", "Y", null, 0},
                {"NA00000", "1", null, 2},
                {"sample1", "X", SVTestUtils.sample1Ploidy1, 1}, //sample1 isn't in ped, so use GT ploidy
                {"sample1", "1", SVTestUtils.sample1, 2}, //sample1 isn't in ped, but contig is autosome
        };
    }

    @Test(dataProvider = "samplesForPloidyQuery")
    public void testSamplePloidy(final String sampleName, final String contig, final Genotype g, final int expected) {
        Assert.assertEquals(JointGermlineCNVSegmentation.getSamplePloidy(allosomalContigs, refCopyNumber, pedigree, sampleName, contig, g), expected);
    }

    //allosomalContigs contains something other than X and Y -- we only support mammals for the time being
    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testBadAllosomes() {
        JointGermlineCNVSegmentation.getSamplePloidy(badAllosomes, refCopyNumber, pedigree, "NA00000", "Z", null);
    }

    //sample not in pedigree and no genotype provided
    @Test(expectedExceptions = {IllegalStateException.class})
    public void testBadInputs() {
        JointGermlineCNVSegmentation.getSamplePloidy(allosomalContigs, refCopyNumber, pedigree, "NA99999", "X", null);
    }

}