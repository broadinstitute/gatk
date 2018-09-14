package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.mockito.Mockito;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static org.mockito.Mockito.when;

public class ImpreciseVariantDetectorUnitTest extends GATKBaseTest {


    private final Logger localLogger = LogManager.getLogger(ImpreciseVariantDetectorUnitTest.class);
    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @DataProvider(name = "evidenceTargetLinksAndVariants")
    public Object[][] getEvidenceTargetLinksAndVariants() {

        final VariantContext impreciseDeletion = new VariantContextBuilder()
                .id("DEL_IMPRECISE_20_950_1050_1975_2025")
                .chr("20").start(1000).stop(2000)
                .alleles("N", SimpleSVType.ImpreciseDeletion.createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL))
                .attribute(VCFConstants.END_KEY, 2000)
                .attribute(GATKSVVCFConstants.SVTYPE, SimpleSVType.SupportedType.DEL.toString())
                .attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, 7)
                .attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, 5)
                .attribute(GATKSVVCFConstants.IMPRECISE, true)
                .attribute(GATKSVVCFConstants.CIPOS, "-50,50")
                .attribute(GATKSVVCFConstants.CIEND, "-25,25")
                .attribute(GATKSVVCFConstants.SVLEN, -1000)
                .make();

        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {
                Arrays.asList(
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 950, 1050), true),
                                new StrandedInterval(new SVInterval(0, 1975, 2025), false),
                                5, 7, new HashSet<>(), new HashSet<>())),
                Arrays.asList( impreciseDeletion ) }
        );

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "evidenceTargetLinksAndVariants", groups = "sv")
    public void testProcessEvidenceTargetLinks(final List<EvidenceTargetLink> etls,
                                               final List<VariantContext> expectedVariants) {
        final int impreciseEvidenceVariantCallingThreshold =
                new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection().impreciseVariantEvidenceThreshold;

        final int maxCallableImpreciseVariantDeletionSize =
                new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection().maxCallableImpreciseVariantDeletionSize;

        final ReferenceMultiSparkSource referenceMultiSource = new ReferenceMultiSparkSource(twoBitRefURL, ReferenceWindowFunctions.IDENTITY_FUNCTION);

        ReadMetadata metadata = Mockito.mock(ReadMetadata.class);
        when(metadata.getMaxMedianFragmentSize()).thenReturn(300);
        when(metadata.getContigName(0)).thenReturn("20");

        PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTree = new PairedStrandedIntervalTree<>();
        etls.forEach(e -> evidenceTree.put(e.getPairedStrandedIntervals(), e));

        final List<VariantContext> impreciseVariants =
                ImpreciseVariantDetector.callImpreciseDeletionFromEvidenceLinks(evidenceTree,
                        metadata,
                        referenceMultiSource,
                        impreciseEvidenceVariantCallingThreshold,
                        maxCallableImpreciseVariantDeletionSize,
                        localLogger);

        VariantContextTestUtils.assertEqualVariants(impreciseVariants, expectedVariants);
    }

}