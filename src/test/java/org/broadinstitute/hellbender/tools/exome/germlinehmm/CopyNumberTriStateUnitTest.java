package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Unit tests for {@link CopyNumberTriState}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CopyNumberTriStateUnitTest {
    @Test
    public void testAlleleListsLogic() {
        Assert.assertEquals(CopyNumberTriState.ALL_ALLELES.size(), CopyNumberTriState.values().length);
        Assert.assertEquals(CopyNumberTriState.ALL_ALLELES.subList(1, CopyNumberTriState.ALL_ALLELES.size()),
                CopyNumberTriState.ALTERNATIVE_ALLELES);
        Assert.assertFalse(CopyNumberTriState.ALTERNATIVE_ALLELES.stream().anyMatch(Allele::isReference));
        Assert.assertTrue(CopyNumberTriState.ALL_ALLELES.get(0).isReference());
        Assert.assertEquals(CopyNumberTriState.ALL_ALLELES.get(0), CopyNumberTriState.NEUTRAL.allele);
        final List<Allele> expectedOrder = new ArrayList<>();
        for (final CopyNumberTriState state : CopyNumberTriState.values()) {
            if (state == CopyNumberTriState.NEUTRAL) {
                continue;
            }
            expectedOrder.add(state.allele);
        }
        Assert.assertEquals(CopyNumberTriState.ALTERNATIVE_ALLELES, expectedOrder);
    }
}
