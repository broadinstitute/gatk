package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class MergeMutect2CallsWithMC3IntegrationTest extends CommandLineProgramTest {

    @Test
    public void test() {
        final String mc3Vcf = toolsTestDir + "/validation/mc3/MC3.vcf";
        final String m2Vcf = toolsTestDir + "/validation/mc3/M2.vcf";
        final File outputVcf = createTempFile("output", ".vcf");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, m2Vcf,
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, mc3Vcf,
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputVcf.toString(),
        };
        runCommandLine(args);

        final Map<String, VariantContext> m2Variants =
                StreamSupport.stream(new FeatureDataSource<VariantContext>(m2Vcf).spliterator(), false)
                .collect(Collectors.toMap(vc -> keyForVariant(vc), vc -> vc));

        final Map<String, VariantContext> mc3Variants =
                StreamSupport.stream(new FeatureDataSource<VariantContext>(mc3Vcf).spliterator(), false)
                        .collect(Collectors.toMap(vc -> keyForVariant(vc), vc -> vc));

        final Map<String, VariantContext> mergedVariants =
                StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false)
                        .collect(Collectors.toMap(vc -> keyForVariant(vc), vc -> vc));

        // all m2 variants should be in output except perhaps filtered variants not in mc3
        Assert.assertTrue(m2Variants.entrySet().stream()
                .allMatch(entry -> (entry.getValue().isFiltered() && !mc3Variants.containsKey(entry.getKey())) || mergedVariants.containsKey(entry.getKey())));

        // all mc3 variants should be in merged vcf
        Assert.assertTrue(mc3Variants.keySet().stream().allMatch(mergedVariants::containsKey));

        // merged variants should contain M2 as a center iff they were in the m2 vcf
        // they should have an M2 filters field only iff it was filtered in M2
        mergedVariants.entrySet().forEach(entry -> {
            final boolean hasM2Center = entry.getValue().getAttributeAsString(MergeMutect2CallsWithMC3.CENTERS_KEY, "")
                    .contains(MergeMutect2CallsWithMC3.M2_CENTER_NAME);
            final boolean inM2Vcf = m2Variants.containsKey(entry.getKey());
            final boolean hasM2Filters = entry.getValue().hasAttribute(MergeMutect2CallsWithMC3.M2_FILTERS_KEY);
            final boolean filteredInM2Vcf = inM2Vcf && m2Variants.get(entry.getKey()).isFiltered();
            Assert.assertEquals(hasM2Center, inM2Vcf);
            Assert.assertEquals(hasM2Filters, filteredInM2Vcf);
        });
    }


    private static String keyForVariant( final VariantContext variant ) {
        return String.format("%s:%d %s", variant.getContig(), variant.getStart(), variant.getEnd(),
                variant.getAlleles());
    }

}