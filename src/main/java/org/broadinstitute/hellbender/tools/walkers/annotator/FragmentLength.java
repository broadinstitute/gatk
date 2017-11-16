package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.OptionalInt;

/**
 * Median fragment length of reads supporting each allele.
 *
 * Created by David Benjamin on 3/20/17.
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Median fragment length of reads supporting each allele (MFRL)")
public class FragmentLength extends PerAlleleAnnotation implements StandardMutectAnnotation {
    public static final String KEY = "MFRL";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : GATKProtectedMathUtils.median(Ints.toArray(values));
    }

    @Override
    protected boolean includeRefAllele() { return true; }

    @Override
    protected String getVcfKey() { return KEY; }

    @Override
    protected String getDescription() { return "median fragment length"; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);
        //abs because fragment lengths are negative if mate comes first
        return OptionalInt.of(Math.abs(read.getFragmentLength()));
    }
}
