package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class JointGermlineCNVSegmentationIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/gcnv-postprocess");
    private static final double DEFAULT_DEFRAG_SAMPLE_OVERLAP = 0.5;
    private static final double DEFAULT_DEFRAG_PADDING_FRACTION = 0.8;

    private static final List<File> SEGMENTS_VCF_CORRECT_OUTPUTS = Arrays.asList(
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_000.vcf"),
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_001.vcf"),
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_002.vcf"));

    @DataProvider
    public Object[][] postprocessOutputs() {
        return new Object[][] {
               new Object[]{SEGMENTS_VCF_CORRECT_OUTPUTS, 6, 7},
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
                    new File(getToolTestDataDir() + "NA21120.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA00000.overlaps.vcf"))
            }
        };
    }

    @Test(dataProvider = "postprocessOutputs")
    public void testThreeGCNVSamples(final List<File> inputVcfs, final int expectedCountDefault, final int expectedCountNoFilter) {
        final File output = createTempFile("threeSamples", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointGermlineCNVSegmentation.MODEL_CALL_INTERVALS_LONG_NAME, getToolTestDataDir() + "threeSamples.interval_list")
                .addIntervals(new File(getToolTestDataDir() + "threeSamples.interval_list"))
                .add(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, getToolTestDataDir() + "threeSamples.ped");

        inputVcfs.forEach(vcf -> args.addVCF(vcf));

        runCommandLine(args, JointGermlineCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> withQStreshold = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final List<VariantContext> withThresholdVariants = withQStreshold.getRight();
        Assert.assertEquals(withQStreshold.getRight().size(), expectedCountDefault);
        Assert.assertTrue(withQStreshold.getRight().stream().noneMatch(vc -> vc.getContig().equals("1")));  //1 is reference on all samples
        //Note that many of these are the same in PostprocessGermlineCNVCallsIntegrationTest::testQualScoreCalculationWithBreakpoints

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        final List<String> variantKeys = variants.stream().map(VariantContextTestUtils::keyForVariant).collect(Collectors.toList());

        final List<String> expectedKeys = Arrays.asList(
                "2:230925-231288 G*, [<DEL>]",
                "2:233003-234369 A*, [<DUP>]",
                "3:1415190-1415854 T*, [<DEL>]",
                "X:223929-224644 A*, [<DUP>]",
                "X:230719-230984 C*, [<DUP>]",
                "Y:1521543-1684258 N*, [<DEL>]");

        Assert.assertTrue(variantKeys.containsAll(expectedKeys));

        VariantContext variant = withThresholdVariants.get(0);
        validateAlleleCount(variant, 2);
        Assert.assertTrue(variant.getGenotype("SAMPLE_001").isHomVar());
        validateCopyNumber(variant, "SAMPLE_001", 0);

        variant = withThresholdVariants.get(1);
        validateAlleleCount(variant, 1);
        validateCopyNumber(variant,"SAMPLE_001", 4);

        variant = withThresholdVariants.get(2);
        validateAlleleCount(variant, 2);
        validateCopyNumber(variant, "SAMPLE_002", 0);

        variant = withThresholdVariants.get(3);
        validateAlleleCount(variant, 1);
        validateCopyNumber(variant, "SAMPLE_001", 3);
        //this is a male sample, so a dupe on X is phased and doesn't have to be a no-call
        Assert.assertTrue(variant.getGenotype("SAMPLE_001").getPloidy() == 1 && variant.getGenotype("SAMPLE_001").isHomVar());

        variant = withThresholdVariants.get(4);
        validateAlleleCount(variant, 1);
        validateCopyNumber(variant, "SAMPLE_001", 2);
        //this is a male sample, so a dupe on X is phased and doesn't have to be a no-call
        Assert.assertTrue(variant.getGenotype("SAMPLE_001").getPloidy() == 1 && variant.getGenotype("SAMPLE_001").isHomVar());

        //The chrY entry starts in PAR1, so it does get a legit N
        variant = withThresholdVariants.get(5);
        validateAlleleCount(variant, 1);
        validateCopyNumber(variant, "SAMPLE_000", 0);
        validateCopyNumber(variant, "SAMPLE_001", 0);
        validateCopyNumber(variant, "SAMPLE_002", 0);
        Assert.assertTrue(variant.getGenotype("SAMPLE_000").getPloidy() == 1 && variant.getGenotype("SAMPLE_000").isNoCall());
        Assert.assertTrue(variant.getGenotype("SAMPLE_001").getPloidy() == 1 && variant.getGenotype("SAMPLE_001").isHomVar());
        Assert.assertTrue(variant.getGenotype("SAMPLE_002").getPloidy() == 1 && variant.getGenotype("SAMPLE_002").isNoCall());

        final File output2 = createTempFile("threeSamples.noQSthreshold",".vcf");

        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addOutput(output2)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointGermlineCNVSegmentation.MIN_QUALITY_LONG_NAME, 0)
                .add(JointGermlineCNVSegmentation.MODEL_CALL_INTERVALS_LONG_NAME, getToolTestDataDir() + "threeSamples.interval_list")
                .addIntervals(new File(getToolTestDataDir() + "threeSamples.interval_list"))
                .add(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, getToolTestDataDir() + "threeSamples.ped");
        inputVcfs.forEach(vcf -> args2.addVCF(vcf));

        runCommandLine(args2, JointGermlineCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> withoutQStreshold = VariantContextTestUtils.readEntireVCFIntoMemory(output2.getAbsolutePath());
        Assert.assertEquals(withoutQStreshold.getRight().size(), expectedCountNoFilter); //extra variant at X:227988
        //check extra variant, make sure other variants are same
    }

    @Test
    public void testDefragmentation() {
        final File output = createTempFile("defragmented",".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .addVCF(getToolTestDataDir() + "NA20533.fragmented.segments.vcf.gz")
                .add(JointGermlineCNVSegmentation.MODEL_CALL_INTERVALS_LONG_NAME, getToolTestDataDir() + "intervals.chr13.interval_list")
                .addInterval("13:52951204-115064572")
                .add(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, getToolTestDataDir() + "NA20533.ped")  //this sample actually appears Turner (X0), but doesn't matter for chr13 sample actually appears Turner (X0), but doesn't matter for chr13
                .add(JointGermlineCNVSegmentation.MIN_SAMPLE_NUM_OVERLAP_LONG_NAME, DEFAULT_DEFRAG_SAMPLE_OVERLAP)
                .add(JointGermlineCNVSegmentation.DEFRAGMENTATION_PADDING_LONG_NAME, DEFAULT_DEFRAG_PADDING_FRACTION);

        runCommandLine(args, JointGermlineCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> defragmentedEvents = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        //events get combined into a single event of 62113369 bp
        Assert.assertEquals(defragmentedEvents.getRight().size(), 1);
        Assert.assertEquals(defragmentedEvents.getRight().get(0).getAttributeAsInt(GATKSVVCFConstants.SVLEN,0), 62113369);

        final File output2 = createTempFile("notDefragmented",".vcf");
        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addOutput(output2)
                .addReference(GATKBaseTest.b37Reference)
                .addVCF(getToolTestDataDir() + "adjacentDifferentCN.vcf")
                .add(JointGermlineCNVSegmentation.MODEL_CALL_INTERVALS_LONG_NAME, getToolTestDataDir() + "intervals.chr8snippet.interval_list")
                .addInterval("8:190726-666104")
                .add(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, getToolTestDataDir() + "NA20520.ped") //this sample actually appears Turner (X0), but doesn't matter for chr13 sample actually appears Turner (X0), but doesn't matter for chr13
                .add(JointGermlineCNVSegmentation.MIN_SAMPLE_NUM_OVERLAP_LONG_NAME, DEFAULT_DEFRAG_SAMPLE_OVERLAP)
                .add(JointGermlineCNVSegmentation.DEFRAGMENTATION_PADDING_LONG_NAME, DEFAULT_DEFRAG_PADDING_FRACTION);

        runCommandLine(args2, JointGermlineCNVSegmentation.class.getSimpleName());

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
                .add(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, getToolTestDataDir() + "overlapping.ped")
                .add(JointGermlineCNVSegmentation.MODEL_CALL_INTERVALS_LONG_NAME, getToolTestDataDir() + "intervals.chr22.interval_list")
                .addInterval("22:22,538,114-23,538,437")
                .add(JointGermlineCNVSegmentation.CLUSTERING_INTERVAL_OVERLAP_LONG_NAME, 0.8)
                .add(JointGermlineCNVSegmentation.CLUSTERING_BREAKEND_WINDOW_LONG_NAME, 0);

        inputVcfs.forEach(vcf -> args.addVCF(vcf));

        runCommandLine(args, JointGermlineCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> overlappingEvents = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(overlappingEvents.getRight().size(), 7);
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

        //spot check some AFs and make sure integer division isn't borked
        Assert.assertEquals(overlappingEvents.getRight().get(0).getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0), 0.056);
        Assert.assertEquals(overlappingEvents.getRight().get(1).getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0), 0.139);
        Assert.assertEquals(overlappingEvents.getRight().get(6).getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0), 0.028);
    }

    private void validateCopyNumber(final VariantContext variant, final String sampleName, final int expectedCopyNumber) {
        Assert.assertEquals(Integer.parseInt(variant.getGenotype(sampleName).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), expectedCopyNumber);
    }

    private void validateAlleleCount(final VariantContext variant, final int expectedAlleleCount) {
        Assert.assertEquals(Integer.parseInt(variant.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY,"")), expectedAlleleCount);
    }
}