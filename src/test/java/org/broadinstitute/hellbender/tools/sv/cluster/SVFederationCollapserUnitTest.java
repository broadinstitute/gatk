package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class SVFederationCollapserUnitTest {

    private static final SVFederationCollapser COLLAPSER = new SVFederationCollapser(
            SVTestUtils.hg38Reference,
            CanonicalSVCollapser.AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE,
            CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE,
            CanonicalSVCollapser.FlagFieldLogic.OR);

    @DataProvider(name = "collapseFiltersData")
    public Object[][] collapseFiltersData() {
        return new Object[][]{
                {List.of(Set.of("LOW_Q"), Set.of("BOTHSIDES_SUPPORT")), Set.of("LOW_Q", "BOTHSIDES_SUPPORT")},
                {List.of(Collections.emptySet(), Set.of("LOW_Q")), Collections.emptySet()},
                {List.of(Collections.emptySet(), Collections.emptySet()), Collections.emptySet()}
        };
    }

    @Test(dataProvider = "collapseFiltersData")
    public void testCollapseFilters(final List<Set<String>> filterSets, final Set<String> expected) {
        final List<SVCallRecord> records = filterSets.stream()
                .map(filters -> new SVCallRecord("record" + filters.hashCode(), "chr1", 1000, true, "chr1", 1999, false,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                        1000, Collections.emptyList(), SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                        List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL), Collections.emptyList(),
                        Collections.emptyMap(), filters, null, SVTestUtils.hg38Dict))
                .collect(Collectors.toList());

        Assert.assertEquals(COLLAPSER.collapseFilters(records), expected);
    }
}
