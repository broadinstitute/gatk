package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.cluster.*;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.TestException;

import java.util.*;
import java.util.stream.Collectors;

public class SVTestUtils {

    public static final ReferenceSequenceFile hg38Reference = ReferenceUtils.createReferenceReader(new GATKPath(GATKBaseTest.hg38Reference));
    public final static SAMSequenceDictionary hg38Dict = hg38Reference.getSequenceDictionary();
    public final static int chr1Length = hg38Dict.getSequence("chr1").getSequenceLength();

    private final static GenomeLocParser glParser = new GenomeLocParser(SVTestUtils.hg38Dict);
    public static final SVCollapser<SVCallRecord> defaultCollapser =
            new CanonicalSVCollapser(
                    hg38Reference,
                    CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
                    CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
                    CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN);

    public static final String PESR_ALGORITHM = "pesr";
    public static final List<String> DEPTH_ONLY_ALGORITHM_LIST = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);
    public static final List<String> PESR_ONLY_ALGORITHM_LIST = Collections.singletonList(PESR_ALGORITHM);

    public static CanonicalSVLinkage<SVCallRecord> getNewDefaultLinkage() {
        final CanonicalSVLinkage<SVCallRecord> linkage = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, false);
        linkage.setDepthOnlyParams(defaultDepthOnlyParameters);
        linkage.setMixedParams(defaultMixedParameters);
        linkage.setEvidenceParams(defaultEvidenceParameters);
        return linkage;
    }

    public static SVClusterEngine<SVCallRecord> getNewDefaultSingleLinkageEngine() {
        return new SVClusterEngine<>(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, defaultCollapser, getNewDefaultLinkage());
    }

    public static SVClusterEngine<SVCallRecord> getNewDefaultMaxCliqueEngine() {
        return new SVClusterEngine<>(SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE, defaultCollapser, getNewDefaultLinkage());
    }

    public static final ClusteringParameters defaultDepthOnlyParameters = ClusteringParameters.createDepthParameters(0.8, 0, 0);
    public static final ClusteringParameters defaultMixedParameters = ClusteringParameters.createMixedParameters(0.8, 1000, 0);
    public static final ClusteringParameters defaultEvidenceParameters = ClusteringParameters.createPesrParameters(0.5, 500, 0);

    public static final SVClusterEngine<SVCallRecord> defaultSingleLinkageEngine = getNewDefaultSingleLinkageEngine();
    public static final SVClusterEngine<SVCallRecord> defaultMaxCliqueEngine = getNewDefaultMaxCliqueEngine();

    public final static int start = 10001;

    public final static int length = 10000;

    public final static int length1b = (int)Math.round(defaultDepthOnlyParameters.getReciprocalOverlap() * length);

    public final static int start1b = start + length - length1b;

    //separated from end of call1 by defragmenter padding (in bin space, according to bins defined by targetIntervals below)
    public final static int start2 = start + length*14/9;

    public final static int start3 = start + (int)Math.round(Math.floor((1- defaultDepthOnlyParameters.getReciprocalOverlap())*length));

    //make intervals like xxxx----xxxx----xxxx----xxxx----xxxx
    public final static List<GenomeLoc> targetIntervals = new ArrayList<>(
            Arrays.asList(glParser.createGenomeLoc("chr1", 1, length/2),  //left edge call
                    glParser.createGenomeLoc("chr1", start, start + length/9),  //start of call1
                    glParser.createGenomeLoc("chr1", start + length*2/9, start + length*3/9),
                    glParser.createGenomeLoc("chr1", start + length*4/9, start + length*5/9),
                    glParser.createGenomeLoc("chr1", start + length*6/9, start + length*7/9),
                    glParser.createGenomeLoc("chr1", start + length*8/9, start + length),  //end of call1
                    glParser.createGenomeLoc("chr1", start + length*10/9, start + length*11/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start + length*12/9, start + length*13/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start2, start2 + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length*2/9, start2 + length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length*4/9, start2 + length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length*6/9, start2 + length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length*8/9, start2 + length),
                    glParser.createGenomeLoc("chr1", start2 + length, start2 + length + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*2/9, start2 + length+ length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*4/9, start2 + length+ length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*6/9, start2 + length+ length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*8/9, start2 + length+ length),
                    glParser.createGenomeLoc("chr1", chr1Length - 99, chr1Length),
                    glParser.createGenomeLoc("chr10", start, start + length - 1)))
            ;

    public final static GenotypeBuilder sample1 = new GenotypeBuilder("sample1",
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1);

    public final static GenotypeBuilder sample1Ploidy1 = new GenotypeBuilder("sample1",
            Lists.newArrayList(Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1);

    public final static GenotypeBuilder sample1_CN0 = new GenotypeBuilder("sample1",
            Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);

    public final static GenotypeBuilder sample2 = new GenotypeBuilder("sample2",
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 3);

    public final static GenotypeBuilder sample3 = new GenotypeBuilder("sample3",
            Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);

    public final static SVCallRecord makeRecord(final String id,
                                                 final String contigA,
                                                 final int positionA,
                                                 final Boolean strandA,
                                                 final String contigB,
                                                 final int positionB,
                                                 final Boolean strandB,
                                                 final StructuralVariantType type,
                                                 final Integer length,
                                                 final List<String> algorithms,
                                                 final List<Allele> alleles,
                                                 final List<GenotypeBuilder> genotypeBuilders) {
        final Allele refAllele = Allele.create(ReferenceUtils.getRefBaseAtPosition(hg38Reference, contigA, positionA), true);
        final List<Allele> newAlleles = replaceRefAlleles(alleles, refAllele);
        final List<Genotype> genotypes = new ArrayList<>(genotypeBuilders.size());
        for (final GenotypeBuilder builder : genotypeBuilders) {
            genotypes.add(makeGenotypeWithRefAllele(builder, refAllele));
        }
        return new SVCallRecord(id, contigA, positionA, strandA, contigB, positionB, strandB, type, length, algorithms,
                newAlleles, genotypes);
    }

    public static final Genotype makeGenotypeWithRefAllele(final GenotypeBuilder builder, final Allele refAllele) {
        final List<Allele> alleles = replaceRefAlleles(builder.make().getAlleles(), refAllele);
        builder.alleles(alleles);
        builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, alleles.size());
        return builder.make();
    }

    public final static SVCallRecord rightEdgeCall = makeRecord("rightEdgeCall", "chr1", chr1Length - 99, null,
            "chr1", chr1Length, null,
            StructuralVariantType.CNV, 100,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord leftEdgeCall = makeRecord("leftEdgeCall", "chr1", 1, null,
            "chr1", length/2, null,
            StructuralVariantType.CNV, length/2,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord nonDepthOnly = makeRecord("nonDepthOnly", "chr1", start, null,
            "chr1", start + length -1, null,
            StructuralVariantType.CNV, length,
            Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, "PE"),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1 = makeRecord("call1", "chr1", start, null,
            "chr1", start + length -1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1b = makeRecord("call1b", "chr1", start1b, null,
            "chr1", start1b + length1b - 1, null,
            StructuralVariantType.CNV, length1b,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call2 = makeRecord("call2", "chr1", start2, null,
            "chr1", start2 + length - 1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call3 = makeRecord("call3", "chr1", start2+length, null,
            "chr1", start2 + length + length - 1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call4_chr10 = makeRecord("call4_chr10", "chr10", start, null,
            "chr10", start + length - 1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1_CN1 = makeRecord("call1_CN1", "chr1", start, null,
            "chr1", start + length - 1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1));

    public final static SVCallRecord call2_CN0 = makeRecord("call2_CN0", "chr1", start2, null,
            "chr1", start2 + length - 1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1_CN0));

    public final static SVCallRecord sameBoundsSampleMismatch = makeRecord("sameBoundsSampleMismatch", "chr1", start, null,
            "chr1", start + length - 1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1, sample3));

    public final static List<Genotype> threeGenotypes = Arrays.asList(sample1.make(), sample2.make(), sample3.make());

    public final static SVCallRecord inversion = makeRecord("inversion", "chr1", 10000, true, "chr1", 20000, true,
            StructuralVariantType.INV, 10001, Arrays.asList("SR", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthOnly = makeRecord("depthOnly", "chr1", 10000, null, "chr1", 20000, null,
            StructuralVariantType.CNV, 10001, Collections.singletonList("depth"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL), Collections.singletonList(sample1));

    public static final SVCallRecord depthAndStuff = makeRecord("depthAndStuff", "chr1", 10000, null, "chr1", 20000, null,
                    StructuralVariantType.CNV, 10001, Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthAndStuff2 = makeRecord("depthAndStuff2", "chr1", 10000 - defaultEvidenceParameters.getWindow(), null, "chr1", 20000, null,
            StructuralVariantType.CNV, 10001 + defaultEvidenceParameters.getWindow(), Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthAndStuff3 = makeRecord("depthAndStuff3", "chr1", 10000 + defaultEvidenceParameters.getWindow(), null, "chr1", 20000, null,
            StructuralVariantType.CNV, 10001 - defaultEvidenceParameters.getWindow(), Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public final static SVCallRecord overlapsCall1 = makeRecord("overlapsCall1", "chr1", start3, null,
            "chr1", start3 + length -1, null,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public static void assertEqualsExceptMembership(final SVCallRecord one, final SVCallRecord two) {
        assertEqualsExceptExcludedAttributes(one, two, Collections.singletonList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
    }

    public static void assertEqualsExceptMembershipAndGT(final SVCallRecord one, final SVCallRecord two) {
        assertEqualsExceptExcludedAttributes(one, two, Collections.singletonList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY), true);
    }

    private static Genotype nullGT(final Genotype gt) {
        return new GenotypeBuilder(gt).alleles(null).make();
    }

    public static void assertEqualsExceptExcludedAttributes(final SVCallRecord one, final SVCallRecord two,
                                                            final Collection<String> excludedAttributeKeys) {
        assertEqualsExceptExcludedAttributes(one, two, excludedAttributeKeys, false);
    }

    public static void assertEqualsExceptExcludedAttributes(final SVCallRecord one, final SVCallRecord two,
                                                            final Collection<String> excludedAttributeKeys,
                                                            final boolean ignoreGT) {
        if (one == two) return;
        Assert.assertEquals(one.getAlgorithms(), two.getAlgorithms());
        Assert.assertEquals(one.getAlleles(), two.getAlleles());

        // Exclude MEMBERS field
        final Map<String, Object> attributesOne = new HashMap<>(one.getAttributes());
        final Map<String, Object> attributesTwo = new HashMap<>(two.getAttributes());
        for (final String key : excludedAttributeKeys) {
            attributesOne.remove(key);
            attributesTwo.remove(key);
        }
        Assert.assertEquals(attributesOne, attributesTwo);

        Assert.assertEquals(one.getPositionAInterval(), two.getPositionAInterval());
        Assert.assertEquals(one.getPositionBInterval(), two.getPositionBInterval());
        Assert.assertEquals(one.getStrandA(), two.getStrandA());
        Assert.assertEquals(one.getStrandB(), two.getStrandB());
        Assert.assertEquals(one.getId(), two.getId());
        Assert.assertEquals(one.getLength(), two.getLength());
        Assert.assertEquals(one.getType(), two.getType());
        Assert.assertEquals(one.getGenotypes().size(), two.getGenotypes().size());
        for (int i = 0; i < one.getGenotypes().size(); i++) {
            final Genotype g1 = ignoreGT ? nullGT(one.getGenotypes().get(i)) : one.getGenotypes().get(i);
            final Genotype g2 = ignoreGT ? nullGT(two.getGenotypes().get(i)) : two.getGenotypes().get(i);
            VariantContextTestUtils.assertGenotypesAreEqual(g1, g2);
        }
    }

    public static void assertContainsAllIgnoreRefAlleleBase(final GenotypesContext one, final GenotypesContext two,
                                                            final boolean ignoreGT) {
        Assert.assertTrue(one.size() >= two.size());
        if ( !one.isEmpty()) {
            Assert.assertTrue(one.getSampleNames().containsAll(two.getSampleNames()), "sample names set");
            final Set<String> samples = two.getSampleNames();
            for ( final String sample : samples ) {
                assertSameGenotypeExceptRefBase(one.get(sample), two.get(sample), ignoreGT);
            }
        }
    }

    public static void assertSameGenotypeExceptRefBase(final Genotype g1, final Genotype g2, final boolean ignoreGT) {
        final GenotypeBuilder builderOne = new GenotypeBuilder(g1);
        final GenotypeBuilder builderTwo = new GenotypeBuilder(g2);
        if (ignoreGT) {
            builderOne.alleles(null);
            builderTwo.alleles(null);
        } else {
            builderOne.alleles(replaceRefAlleles(g1.getAlleles(), Allele.REF_N));
            builderTwo.alleles(replaceRefAlleles(g2.getAlleles(), Allele.REF_N));
        }
        VariantContextTestUtils.assertGenotypesAreEqual(builderOne.make(), builderTwo.make());
    }

    public static List<Allele> replaceRefAlleles(final List<Allele> alleles, final Allele replace) {
        return alleles.stream().map(a -> a.isReference() ? replace : a).collect(Collectors.toList());
    }

    public static Genotype buildHomGenotypeWithPloidy(final Allele allele, final int ploidy) {
        return new GenotypeBuilder()
                .alleles(buildHomAlleleListWithPloidy(allele, ploidy))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy)
                .make();
    }

    public static List<Allele> buildHomAlleleListWithPloidy(final Allele allele, final int ploidy) {
        return Collections.nCopies(ploidy, allele);
    }

    public static List<Allele> buildBiallelicListWithPloidy(final Allele altAllele, final Allele refAllele, final int copyNumber, final int ploidy) {
        final int numAlt;
        if (altAllele.equals(Allele.SV_SIMPLE_DEL)) {
            numAlt = Math.max(ploidy - copyNumber, 0);
        } else if (altAllele.equals(Allele.SV_SIMPLE_DUP)) {
            numAlt = Math.max(copyNumber - ploidy, 0);
        } else {
            throw new TestException("Unsupported alt allele: " + altAllele.getDisplayString());
        }
        final List<Allele> alleles = new ArrayList<>(ploidy);
        final int numRef = ploidy - numAlt;
        for (int i = 0; i < numRef; i++) {
            alleles.add(refAllele);
        }
        for (int i = 0; i < numAlt; i++) {
            alleles.add(altAllele);
        }
        return alleles;
    }

    public static SVCallRecord newCallRecordWithAlleles(final List<Allele> genotypeAlleles, final List<Allele> variantAlleles,
                                                        final StructuralVariantType svtype, final Integer expectedCopyNumber,
                                                        final Integer copyNumber) {
        GenotypeBuilder builder = new GenotypeBuilder("sample").alleles(genotypeAlleles);
        if (expectedCopyNumber != null) {
            builder = builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
        }
        if (copyNumber != null) {
            builder = builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
        }
        return new SVCallRecord("", "chr1", 100, getValidTestStrandA(svtype), "chr1", 199, getValidTestStrandB(svtype),
                svtype, 100, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                variantAlleles,
                Collections.singletonList(builder.make()),
                Collections.emptyMap());
    }

    public static SVCallRecord newNamedDeletionRecordWithAttributes(final String id, final Map<String, Object> attributes) {
        return new SVCallRecord(id, "chr1", 100, true, "chr1", 199, false,
                StructuralVariantType.DEL,
                100, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                attributes);
    }

    public static final Map<String, Object> keyValueArraysToMap(final String[] keys, final Object[] values) {
        final Map<String, Object> map = new HashMap<>();
        for (int i = 0; i < keys.length; i++) {
            map.put(keys[i], values[i]);
        }
        return map;
    }

    // Note that strands may not be valid
    public static SVCallRecord newCallRecordWithLengthAndTypeAndChrom2(final Integer length, final StructuralVariantType svtype, final String chrom2) {
        final int positionB = length == null ? 1 : length;
        return new SVCallRecord("", "chr1", 1, getValidTestStrandA(svtype), chrom2, positionB, getValidTestStrandB(svtype),
                svtype, length, PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(), Collections.emptyList(),
                Collections.emptyMap());
    }

    public static SVCallRecord newDeletionCallRecordWithId(final String id) {
        return new SVCallRecord(id, "chr1", 1, true, "chr1", 100, false,
                StructuralVariantType.DEL, 100, PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap());
    }

    public static SVCallRecord newDeletionCallRecordWithIdAndAlgorithms(final String id, final List<String> algorithms) {
        return new SVCallRecord(id, "chr1", 1, true, "chr1", 100, false,
                StructuralVariantType.DEL, 100, algorithms, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap());
    }

    // Note strands and length may not be set properly
    public static SVCallRecord newCallRecordWithIntervalAndType(final int start, final int end, final StructuralVariantType svtype) {
        return new SVCallRecord("", "chr1", start, getValidTestStrandA(svtype), "chr1", end, getValidTestStrandB(svtype),
                svtype, getLength(start, end, svtype), PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap());
    }

    public static Integer getLength(final int start, final int end, final StructuralVariantType type) {
        if (type.equals(StructuralVariantType.BND) || type.equals(StructuralVariantType.INS)) {
            return null;
        }
        return end - start + 1;
    }

    public static SVCallRecord newBndCallRecordWithStrands(final boolean strandA, final boolean strandB) {
        return new SVCallRecord("", "chr1", 1000, strandA, "chr1", 1999, strandB, StructuralVariantType.BND, null,
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap());
    }

    public static SVCallRecord newCnvCallRecordWithStrands(final Boolean strandA, final Boolean strandB) {
        return new SVCallRecord("", "chr1", 1000, strandA, "chr1", 1999, strandB, StructuralVariantType.CNV, 1000,
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap());
    }

    public static SVCallRecord newCallRecordWithCoordinates(final String id, final String chrA, final int posA, final String chrB, final int posB) {
        return new SVCallRecord(id, chrA, posA, true, chrB, posB, false, StructuralVariantType.BND, null,
                Collections.singletonList("peser"),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap());
    }

    public static SVCallRecord newCallRecordWithCoordinatesAndType(final String id, final String chrA, final int posA, final String chrB, final int posB, final StructuralVariantType type) {
        return new SVCallRecord(id, chrA, posA, true, chrB, posB, false, type, getLength(posA, posB, type),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap());
    }

    public static SVCallRecord newCallRecordWithAlgorithms(final List<String> algorithms) {
        return new SVCallRecord("", "chr1", 1000, true, "chr1", 1000, false, StructuralVariantType.INS, length,
                algorithms,
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap());
    }

    public static SVCallRecord newCallRecordInsertionWithLength(final Integer length) {
        return new SVCallRecord("", "chr1", 1000, true, "chr1", 1000, false, StructuralVariantType.INS, length,
                PESR_ONLY_ALGORITHM_LIST,
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap());
    }

    public static VariantContext newVariantContext(final String id, final String chrA, final int posA,
                                                   final int end, final List<Allele> alleles,
                                                   final List<Genotype> genotypes, final Integer svlen,
                                                   final String strands, final StructuralVariantType svtype,
                                                   final List<String> algorithms, final String chr2, final Integer end2,
                                                   final Map<String, Object> extendedAttributes) {
        final Map<String, Object> attributes = new HashMap<>();
        attributes.put(GATKSVVCFConstants.SVLEN, svlen);
        attributes.put(GATKSVVCFConstants.SVTYPE, svtype);
        attributes.put(GATKSVVCFConstants.STRANDS_ATTRIBUTE, strands);
        attributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, algorithms);
        if (chr2 != null && end2 != null) {
            attributes.put(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, chr2);
            attributes.put(GATKSVVCFConstants.END2_ATTRIBUTE, end2);
        }
        attributes.putAll(extendedAttributes);
        return new VariantContextBuilder()
                .id(id)
                .start(posA)
                .chr(chrA)
                .stop(end)
                .alleles(alleles)
                .genotypes(genotypes)
                .attributes(attributes)
                .make();
    }

    public static Boolean getValidTestStrandA(final StructuralVariantType type) {
        if (type == StructuralVariantType.DUP) {
            return Boolean.FALSE;
        } else if (type == StructuralVariantType.CNV) {
            return null;
        } else if (type == StructuralVariantType.INV) {
            return Boolean.TRUE;
        }
        return Boolean.TRUE;
    }

    public static Boolean getValidTestStrandB(final StructuralVariantType type) {
        if (type == StructuralVariantType.DUP) {
            return Boolean.TRUE;
        } else if (type == StructuralVariantType.CNV) {
            return null;
        } else if (type == StructuralVariantType.INV) {
            return Boolean.TRUE;
        }
        return Boolean.FALSE;
    }
}
