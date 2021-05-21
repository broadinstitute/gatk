package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableSet;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.variant.utils.VCFHeaderReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class VariantAnnotatorIntegrationTest extends CommandLineProgramTest {
    // small VCF with a few IDs filled in.  Contains multiallelic sites.
    private final File BASIC_INPUT = getTestFile("input.vcf");

    // the above subsetted to indel sites with some INFO and FILTER fields added.
    private final File INDELS = getTestFile("indels.vcf");

    private static final String FOO = "foo";
    private static final String FOO2 = FOO + 2;

    private static final String FOO_FILTER = FOO + ".FILTER";
    private static final String FOO_ID = FOO + ".ID";
    private static final String AC = VCFConstants.ALLELE_COUNT_KEY;
    private static final String FOO_AC = FOO + "." + AC;

    private static final double GENEROUS_RELATIVE_TOLERANCE = 0.2;

    // the purpose of VariantAnnotator is to emulate the annotations that HaplotypeCaller and Mutect2 would emit after
    // local assembly and calculating likelihoods with Pair-HMM.  Thus our basic test is to compare the annotations of
    // VariantAnnotator to those of Mutect2
    @Test
    public void testAgainstMutect2() throws Exception {
        final File bam = new File(largeFileTestDir + "mutect/dream_synthetic_bams/", "tumor_4.bam");
        final File mutect2AnnotatedVcf = createTempFile("mutect2-annotated", ".vcf");
        final File mutect2UnannotatedVcf = createTempFile("mutect2-unannotated", ".vcf");
        final File reannotatedVcf = createTempFile("reannotated", ".vcf");
        final File reference = new File(b37Reference);

        // run Mutect2 on an interval that has indels and multiallelics
        final ArgumentsBuilder mutect2Args = new ArgumentsBuilder()
                .addInput(bam)
                .addOutput(mutect2AnnotatedVcf)
                .addReference(reference)
                .addInterval(new SimpleInterval("20", 50_000_000, 52_000_000))
                .add(StandardArgumentDefinitions.ENABLE_ALL_ANNOTATIONS, true);

        final ArgumentsBuilder mutect2ArgsWithoutAnnotations = new ArgumentsBuilder()
                .addInput(bam)
                .addOutput(mutect2UnannotatedVcf)
                .addReference(reference)
                .addInterval(new SimpleInterval("20", 50_000_000, 52_000_000))
                .add(StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, true);

        new Main().instanceMain(makeCommandLineArgs(mutect2Args.getArgsList(), Mutect2.class.getSimpleName()));

        new Main().instanceMain(makeCommandLineArgs(mutect2ArgsWithoutAnnotations.getArgsList(), Mutect2.class.getSimpleName()));

        // annotate the unannotated Mutect2 output
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(bam)
                .addVCF(mutect2UnannotatedVcf)
                .addReference(reference)
                .addOutput(reannotatedVcf)
                .add(StandardArgumentDefinitions.ENABLE_ALL_ANNOTATIONS, true);

        runCommandLine(args.getArgsArray());

        // check a representative sampling of annotations
        final Set<VariantAnnotation> representativeInfoIntegerAnnotations = ImmutableSet.of(new BaseQuality(), new MappingQuality(), new FragmentLength(), new ChromosomeCounts());
        final Set<VariantAnnotation> representativeInfoStringAnnotations = ImmutableSet.of(new ReferenceBases());
        final Set<VariantAnnotation> representativeInfoBooleanAnnotations = ImmutableSet.of(new TandemRepeat());

        final Set<GenotypeAnnotation> representativeGenotypeAnnotations = ImmutableSet.of(new DepthPerAlleleBySample(), new DepthPerSampleHC(), new AlleleFraction());

        // check the header
        final Set<String> infoHeaderKeys = getHeaderFromFile(reannotatedVcf).getInfoHeaderLines().stream().map(VCFInfoHeaderLine::getID).collect(Collectors.toSet());
        Stream.of(representativeInfoBooleanAnnotations, representativeInfoIntegerAnnotations, representativeInfoStringAnnotations)
                .flatMap(Set::stream)
                .flatMap(annotation -> annotation.getKeyNames().stream())
                .forEach(key -> Assert.assertTrue(infoHeaderKeys.contains(key)));

        final Set<String> formatHeaderKeys = getHeaderFromFile(reannotatedVcf).getFormatHeaderLines().stream().map(VCFFormatHeaderLine::getID).collect(Collectors.toSet());
        representativeGenotypeAnnotations.stream().flatMap(a -> a.getKeyNames().stream()).forEach(key -> Assert.assertTrue(formatHeaderKeys.contains(key)));

        // check the annotated values
        final List<VariantContext> mutect2AnnotatedVariants = VariantContextTestUtils.getVariantContexts(mutect2AnnotatedVcf);
        final List<VariantContext> reannotatedVariants = VariantContextTestUtils.getVariantContexts(reannotatedVcf);

        Assert.assertEquals(mutect2AnnotatedVariants.size(), reannotatedVariants.size());

        final MutableInt matchingInfoFieldValues = new MutableInt(0);
        final MutableInt nonMatchingInfoFieldValues = new MutableInt(0);

        final MutableInt matchingFormatFieldValues = new MutableInt(0);
        final MutableInt nonMatchingFormatFieldValues = new MutableInt(0);

        for (int n = 0; n < reannotatedVariants.size(); n++) {
            final VariantContext mutect2VC = mutect2AnnotatedVariants.get(n);
            final VariantContext reannotatedVC = reannotatedVariants.get(n);
            Assert.assertEquals(keyForVariant(mutect2VC), keyForVariant(reannotatedVC));

            // check INFO field annotations
            representativeInfoBooleanAnnotations.stream()
                    .flatMap(a -> a.getKeyNames().stream())
                    .forEach(key -> Assert.assertEquals(mutect2VC.hasAttribute(key), reannotatedVC.hasAttribute(key)));

            representativeInfoIntegerAnnotations.stream()
                    .flatMap(a -> a.getKeyNames().stream())
                    .forEach(key -> {
                                final double[] mutect2Values = Doubles.toArray(mutect2VC.getAttributeAsDoubleList(key, 1.0));
                                final double[] reannotatedValues = Doubles.toArray(reannotatedVC.getAttributeAsDoubleList(key, 1.0));
                        (approximatelyEqual(mutect2Values, reannotatedValues) ? matchingInfoFieldValues : nonMatchingInfoFieldValues).increment();
                            });

            representativeInfoStringAnnotations.stream()
                    .flatMap(a -> a.getKeyNames().stream())
                    .forEach(key -> Assert.assertEquals(mutect2VC.getAttribute(key), reannotatedVC.getAttribute(key)));

            // check FORMAT annotations
            final Genotype mutect2Genotype = mutect2VC.getGenotype(0);
            final Genotype reannotatedGenotype = reannotatedVC.getGenotype(0);

            representativeGenotypeAnnotations.stream()
                    .flatMap(a -> a.getKeyNames().stream())
                    .forEach(key -> {
                        final double[] mutect2Values = VariantContextGetters.getAttributeAsDoubleArray(mutect2Genotype, key, () -> null, 1.0);
                        final double[] reannotatedValues = VariantContextGetters.getAttributeAsDoubleArray(reannotatedGenotype, key, () -> null, 1.0);
                        (approximatelyEqual(mutect2Values, reannotatedValues) ? matchingFormatFieldValues : nonMatchingFormatFieldValues).increment();
                    });
        }

        // VariantAnnotator can't hope to match M2 and HC perfectly since it doesn't go throught the whole process of assembly and realignment.
        // The best we can expect is that annotations usually match.  Performance as of this writing (GATK 4.1.5.0) is that less than 1 in 40 INFO
        // and less than 1 in 15 FORMAT annotations are off
        Assert.assertTrue(matchingInfoFieldValues.intValue() > 40 * nonMatchingInfoFieldValues.intValue());
        Assert.assertTrue(matchingFormatFieldValues.intValue() > 15 * nonMatchingFormatFieldValues.intValue());
    }

    // test the comp annotation with one and two files
    @Test
    public void testComp() {
        final File inputVCF = BASIC_INPUT;
        final File outputVCF = createTempFile("output", ".vcf");
        final File compVcf = INDELS;

        final Set<String> unfilteredCompVariants = VariantContextTestUtils.streamVcf(compVcf)
                .filter(VariantContext::isNotFiltered)
                .map(VariantAnnotatorIntegrationTest::keyForVariant)
                .collect(Collectors.toSet());

        final ArgumentsBuilder argsForOneComp = new ArgumentsBuilder()
                .addVCF(inputVCF)
                .addOutput(outputVCF)
                .add(StandardArgumentDefinitions.COMPARISON_LONG_NAME + ":" + FOO, compVcf);

        runCommandLine(argsForOneComp.getArgsList());

        for (final VariantContext outputVC : VariantContextTestUtils.getVariantContexts(outputVCF)) {
            Assert.assertEquals(outputVC.hasAttribute(FOO), unfilteredCompVariants.contains(keyForVariant(outputVC)));
        }

        // add the input as a comp -- every site should get this annotation
        final ArgumentsBuilder argsForTwoComps = new ArgumentsBuilder(argsForOneComp.getArgsArray())
                .add(StandardArgumentDefinitions.COMPARISON_LONG_NAME + ":" + FOO2, inputVCF);

        runCommandLine(argsForTwoComps.getArgsList());

        for (final VariantContext outputVC : VariantContextTestUtils.getVariantContexts(outputVCF)) {
            Assert.assertEquals(outputVC.hasAttribute(FOO), unfilteredCompVariants.contains(keyForVariant(outputVC)));
            Assert.assertTrue(outputVC.hasAttribute(FOO2));
        }

    }

    // test dbSNP annotation
    // the input VCF -- up to the ID field -- is as follows:
    //  20	10002353	rs372143558     in dbSNP with same tag
    //  20	10002374	.               not in dbSNP
    //  20	10002443	rs1;rs2         in dbSNP with tag rs188831105
    //  20	10002458	.               in dbSNP with tag rs34527371
    //  20	10002460	rstestnot       not in dbSNP
    //  20	10002470	rstest          in dbSNP with tag rs2327260
    //  20	10002470	rstest          in dbSNP with tag rs3062214
    //  20	10002478	.               in dbSNP with tag rs33961276
    //  we expect the correct tag to be added to the ID field if necessary and the "DB" INFO tag to be added iff the variant is in dbSNP
    @Test
    public void testDBSNP() {
        final File inputVCF = BASIC_INPUT;
        final File outputVCF = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(inputVCF)
                .addOutput(outputVCF)
                .add(DbsnpArgumentCollection.DBSNP_LONG_NAME, dbsnp_138_b37_20_21_vcf);

        runCommandLine(args.getArgsList());

        final List<VariantContext> inputVCs = VariantContextTestUtils.getVariantContexts(inputVCF);
        final List<VariantContext> outputVCs = VariantContextTestUtils.getVariantContexts(outputVCF);

        Assert.assertTrue(inputVCs.size() == outputVCs.size());
        for (int n = 0; n < inputVCs.size(); n++) {
            final VariantContext inputVC = inputVCs.get(n);
            final VariantContext outputVC = outputVCs.get(n);

            // all but two variants are in dbSNP, one of which has the right tag
            final boolean inDbSNP = !(n == 1 || n == 4);
            final boolean inDbSNPWithSameKey = n == 0;
            final boolean inDbSNPWithDifferentKey = inDbSNP && !inDbSNPWithSameKey;

            // check the INFO tag
            Assert.assertEquals(outputVC.hasAttribute(VCFConstants.DBSNP_KEY), inDbSNP);

            // check the ID field has the correct ID appended if necessary
            if (inDbSNP) {
                Assert.assertTrue(!inputVC.hasID() || outputVC.getID().startsWith(inputVC.getID()));
                Assert.assertEquals(outputVC.getID().length() > inputVC.getID().length(), inDbSNPWithDifferentKey);
            }
        }
    }

    // test the -resource argument and associated ID, FILTER and INFO field annotation
    @Test
    public void testResource() {
        final File inputVCF = BASIC_INPUT;
        final File outputVCF = createTempFile("output", ".vcf");
        final File resourceVcf = INDELS;

        final Map<Integer, VariantContext> resourceVariantsByStart = VariantContextTestUtils.streamVcf(resourceVcf)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc));

        final Map<String, VariantContext> resourceVariantsByStartAndAlleles = VariantContextTestUtils.streamVcf(resourceVcf)
                .collect(Collectors.toMap(VariantAnnotatorIntegrationTest::keyForVariant, vc -> vc));

        final String[] baseArgs = new ArgumentsBuilder()
                .addVCF(inputVCF)
                .addOutput(outputVCF)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":" + FOO, resourceVcf)
                .getArgsArray();

        // test the --expression foo.FILTER annotation
        final ArgumentsBuilder filterArgs = new ArgumentsBuilder(baseArgs).add(VariantAnnotator.EXPRESSION_LONG_NAME, FOO_FILTER);
        runCommandLine(filterArgs.getArgsList());
        for (final VariantContext outputVC : VariantContextTestUtils.getVariantContexts(outputVCF)) {
            final VariantContext resourceVariant = resourceVariantsByStart.get(outputVC.getStart());
            if (resourceVariant != null) {
                Assert.assertTrue(outputVC.hasAttribute(FOO_FILTER));
                Assert.assertEquals(outputVC.getAttribute(FOO_FILTER), resourceVariant.isNotFiltered() ? "PASS" : String.join(";", new ArrayList<>(resourceVariant.getFilters())));
            } else {
                Assert.assertFalse(outputVC.hasAttribute(FOO_FILTER));
            }
        }

        // test the --expression foo.ID annotation
        final ArgumentsBuilder idArgs = new ArgumentsBuilder(baseArgs).add(VariantAnnotator.EXPRESSION_LONG_NAME, FOO_ID);
        runCommandLine(idArgs.getArgsList());
        for (final VariantContext outputVC : VariantContextTestUtils.getVariantContexts(outputVCF)) {
            final VariantContext resourceVariant = resourceVariantsByStart.get(outputVC.getStart());
            if (resourceVariant != null && resourceVariant.hasID()) {
                Assert.assertTrue(outputVC.hasAttribute(FOO_ID));
                Assert.assertEquals(outputVC.getAttribute(FOO_ID), resourceVariant.getID());
            } else {
                Assert.assertFalse(outputVC.hasAttribute(FOO_ID));
            }
        }

        // test the --expression foo.AC annotation as an example of a resource INFO field that can be array-valued
        final ArgumentsBuilder infoArgs = new ArgumentsBuilder(baseArgs).add(VariantAnnotator.EXPRESSION_LONG_NAME, FOO_AC);
        runCommandLine(infoArgs.getArgsList());
        for (final VariantContext outputVC : VariantContextTestUtils.getVariantContexts(outputVCF)) {
            final VariantContext resourceVariant = resourceVariantsByStartAndAlleles.get(keyForVariant(outputVC));
            if (resourceVariant != null && resourceVariant.hasAttribute(AC)) {
                Assert.assertTrue(outputVC.hasAttribute(FOO_AC));
                Assert.assertEquals(outputVC.getAttribute(FOO_AC), resourceVariant.getAttribute(AC));
            } else {
                Assert.assertFalse(outputVC.hasAttribute(FOO_AC));
            }
        }
    }

    @Test
    public void testNoReads() {
        final File inputVCF = BASIC_INPUT;
        final File outputVCF = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(inputVCF)
                .addOutput(outputVCF)
                .add(StandardArgumentDefinitions.ENABLE_ALL_ANNOTATIONS, true)
                .addInterval(new SimpleInterval("20", 1, 1));

        runCommandLine(args.getArgsList());
    }

    @Test(expectedExceptions = UserException.class)
    public void testValidationReadsDontMatch() {
        final File inputVCF = BASIC_INPUT;
        final File outputVCF = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(inputVCF)
                .addOutput(outputVCF)
                .addInput(new File(largeFileTestDir + "mutect/dream_synthetic_bams/", "tumor_4.bam")) // This bam corresponds to a sample that is not in the input VCF
                .add(StandardArgumentDefinitions.ENABLE_ALL_ANNOTATIONS, true)
                .addInterval(new SimpleInterval("20", 1, 1));

        runCommandLine(args.getArgsList());
    }

    @Test
    public void testDeNovo() {
        final File inputVCF = getTestFile("trioGGVCF.vcf.gz");
        final File outputVCF = createTempFile("output", ".vcf");
        final File pedigree = new File(getToolTestDataDir(), "trio.ped");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(inputVCF)
                .addOutput(outputVCF)
                .add(StandardArgumentDefinitions.ANNOTATION_LONG_NAME, PossibleDeNovo.class.getSimpleName())
                .add(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, pedigree);

        runCommandLine(args.getArgsList());

        final int[] lociWithLowConfidenceDeNovo = VariantContextTestUtils.streamVcf(outputVCF)
                .filter(vc -> vc.hasAttribute(GATKVCFConstants.LO_CONF_DENOVO_KEY))
                .mapToInt(VariantContext::getStart)
                .toArray();

        final int[] lociWithHighConfidenceDeNovo = VariantContextTestUtils.streamVcf(outputVCF)
                .filter(vc -> vc.hasAttribute(GATKVCFConstants.HI_CONF_DENOVO_KEY))
                .mapToInt(VariantContext::getStart)
                .toArray();

        // known possible de novo sites
        Assert.assertEquals(lociWithLowConfidenceDeNovo, new int[] {10088967});
        Assert.assertEquals(lociWithHighConfidenceDeNovo, new int[] {10130767, 10197999});
    }

    private static String keyForVariant( final VariantContext variant ) {
        return String.format("%s:%d-%d %s", variant.getContig(), variant.getStart(), variant.getEnd(), variant.getAlleles());
    }

    private static VCFHeader getHeaderFromFile(final File vcfFile) throws IOException {
        try (SeekablePathStream stream = new SeekablePathStream(vcfFile.toPath())) {
            return VCFHeaderReader.readHeaderFrom(stream);
        }
    }

    private boolean approximatelyEqual(final double[] array1, final double[] array2) {
        if (array1 == null || array2 == null) {
            return array1 == null && array2 == null;
        }
        if (array1.length != array2.length) {
            return false;
        }

        final double diffSum = MathUtils.sum(MathUtils.applyToArray(MathArrays.ebeSubtract(array1, array2), Math::abs));
        final double absSum = MathUtils.sum(MathArrays.ebeAdd(MathUtils.applyToArray(array1, Math::abs), MathUtils.applyToArray(array2, Math::abs)));
        return diffSum / absSum < GENEROUS_RELATIVE_TOLERANCE;
    }
}

