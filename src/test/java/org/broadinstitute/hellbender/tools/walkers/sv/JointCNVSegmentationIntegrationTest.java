package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.GATKBaseTest.toolsTestDir;
import static org.broadinstitute.hellbender.testutils.BaseTest.createTempFile;
import static org.testng.Assert.*;

public class JointCNVSegmentationIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/gcnv-postprocess");

    private static final List<File> SEGMENTS_VCF_CORRECT_OUTPUTS = Arrays.asList(
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_000.vcf"),
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_001.vcf"),
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_002.vcf"));

    @DataProvider
    public Object[][] postprocessOutputs() {
        return new Object[][] {
               new Object[]{SEGMENTS_VCF_CORRECT_OUTPUTS, 5, 6},
        };
    }

    @DataProvider
    public Object[][] overlappingSamples() {
        return new Object[][] {
            new Object[] {
                    Arrays.asList(new File(getToolTestDataDir() + "HG00365.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG01623.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG01789.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG02165.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG02221.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA07357.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA11829.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12005.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12046.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12814.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12873.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA18946.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA18997.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA19428.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA19456.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA20502.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA21120.overlaps.vcf.gz"))
            }
        };
    }

    @Test(dataProvider = "postprocessOutputs")
    public void testThreeGCNVSamples(final List<File> inputVcfs, final int expectedCountDefault, final int expectedCountNoFilter) {
        final File output = createTempFile("threeSamples", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "threeSamples.interval_list")
                .addIntervals(new File(getToolTestDataDir() + "threeSamples.interval_list"));

        inputVcfs.forEach(vcf -> args.addVCF(vcf));

        runCommandLine(args, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> withQStreshold = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(withQStreshold.getRight().size(), expectedCountDefault);
        Assert.assertTrue(withQStreshold.getRight().stream().noneMatch(vc -> vc.getContig().equals("1") || vc.getContig().equals("Y")));  //1 and Y are reference on all samples
        Assert.assertTrue(withQStreshold.getRight().stream().noneMatch(vc -> vc.getReference().equals(Allele.REF_N))); //by supplying a reference we can fill in ref bases

        final File output2 = createTempFile("threeSamples.noQSthreshold",".vcf");

        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addOutput(output2)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointCNVSegmentation.MIN_QUALITY_LONG_NAME, 0)
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "threeSamples.interval_list")
                .addIntervals(new File(getToolTestDataDir() + "threeSamples.interval_list"));
        inputVcfs.forEach(vcf -> args2.addVCF(vcf));

        runCommandLine(args2, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> withoutQStreshold = VariantContextTestUtils.readEntireVCFIntoMemory(output2.getAbsolutePath());
        Assert.assertEquals(withoutQStreshold.getRight().size(), expectedCountNoFilter); //extra variant at X:227988

        //another test to make sure adjacent events with different copy numbers don't get merged/defragmented?
    }

    @Test
    public void testDefragmentation() {
        final File output = createTempFile("defragmented",".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .addVCF(getToolTestDataDir() + "NA20533.fragmented.segments.vcf.gz")
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "intervals.chr13.interval_list")
                .addInterval("13:52951204-115064572");

        runCommandLine(args, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> defragmentedEvents = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        //events get combined into a single event of 62113369 bp
        Assert.assertEquals(defragmentedEvents.getRight().size(), 1);
        Assert.assertEquals(defragmentedEvents.getRight().get(0).getAttributeAsInt(GATKSVVCFConstants.SVLEN,0), 62113369);

        final File output2 = createTempFile("notDefragmented",".vcf");
        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addOutput(output2)
                .addReference(GATKBaseTest.b37Reference)
                .addVCF(getToolTestDataDir() + "adjacentDifferentCN.vcf")
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "intervals.chr8snippet.interval_list")
                .addInterval("8:190726-666104");

        runCommandLine(args2, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> notDefragmentedEvents = VariantContextTestUtils.readEntireVCFIntoMemory(output2.getAbsolutePath());
        //events have different copy numbers and don't get merged
        Assert.assertEquals(notDefragmentedEvents.getRight().size(), 3);
        Assert.assertEquals(Integer.parseInt(notDefragmentedEvents.getRight().get(0).getGenotype("NA20520").getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 3);
        Assert.assertEquals(Integer.parseInt(notDefragmentedEvents.getRight().get(1).getGenotype("NA20520").getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 4);
        Assert.assertEquals(Integer.parseInt(notDefragmentedEvents.getRight().get(2).getGenotype("NA20520").getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 3);
    }

    @Test(dataProvider = "overlappingSamples")
    public void testOverlappingEvents(final List<File> inputVcfs) {
        final File output = createTempFile("overlaps", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "intervals.chr22.interval_list")
                .addInterval("22:22,538,114-23,538,437");

        inputVcfs.forEach(vcf -> args.addVCF(vcf));

        runCommandLine(args, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> overlappingEvents = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(overlappingEvents.getRight().size(), 6);
        //do copy number checks on genotypes
        //at the start of the contig, all homRef genotypes should be CN2
        final VariantContext vc0 = overlappingEvents.getRight().get(0);
        for (final Genotype g : vc0.getGenotypes()) {
            if (g.isHomRef()) {
                Assert.assertEquals(Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 2);
            } else {
                Assert.assertNotEquals(Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 2);
            }
        }

        //all the events in this example are deletions
        for (final VariantContext vc : overlappingEvents.getRight()) {
            Assert.assertTrue(vc.getAlternateAlleles().size() == 1 && vc.getAlternateAllele(0).equals(GATKSVVCFConstants.DEL_ALLELE));
            Assert.assertTrue(vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY));
            Assert.assertFalse(vc.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY, "").contains(",")); //no zero ACs for uncalled alts
            Assert.assertTrue(vc.hasAttribute(GATKSVVCFConstants.SVTYPE) && vc.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals("DEL"));
            Assert.assertTrue(vc.hasAttribute(GATKSVVCFConstants.SVLEN));
        }

        //in NA11829 variant events are not overlapping, so there should be a CN2 homRef in between
        final List<String> samplesWithOverlaps = Arrays.asList("HG00365", "HG01789", "HG02221", "NA07357", "NA12005", "NA12873", "NA18997", "NA19428", "NA21120");
        final List<String> samplesWithGaps = Arrays.asList("NA11829");

        //all of these samples have an event that overlaps the next event, which is not called in that sample
        boolean sawVariant;
        for (final String sample : samplesWithOverlaps) {
            sawVariant = false;
            for (final VariantContext vc : overlappingEvents.getRight()) {
                if (!sawVariant && !vc.getGenotype(sample).isHomRef()) {
                    sawVariant = true;
                    continue;
                }
                if (sawVariant) {
                    Assert.assertTrue(vc.getGenotype(sample).isHomRef()
                            && (Integer.parseInt(vc.getGenotype(sample).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()) != 2));
                    break;
                }
            }
        }

        //these samples have a variant that doesn't overlap the next call
        for (final String sample : samplesWithGaps) {
            sawVariant = false;
            for (final VariantContext vc : overlappingEvents.getRight()) {
                if (!sawVariant && !vc.getGenotype(sample).isHomRef()) {
                    sawVariant = true;
                    continue;
                }
                if (sawVariant) {
                    Assert.assertTrue(vc.getGenotype(sample).isHomRef()
                            && (Integer.parseInt(vc.getGenotype(sample).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()) == 2));
                    break;
                }
            }
        }

        //all the hom var deletions should be CN0
        for (final VariantContext vc : overlappingEvents.getRight()) {
            for (final Genotype g : vc.getGenotypes()) {
                if (g.isHomVar()) {
                    Assert.assertTrue(g.hasAnyAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT));
                    Assert.assertEquals(Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 0);
                }
            }
        }
    }
}