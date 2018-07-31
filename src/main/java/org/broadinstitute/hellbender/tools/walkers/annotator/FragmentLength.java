package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.OptionalInt;

/**
 * Median fragment length of reads supporting each allele.
 *
 * <p>The output is an array containing, for each allele, the median fragment length (this is synonymous with "insert size") over all reads that best match that allele.  Fragment length here is defined to be non-negative and is zero for unpaired reads.  Thus is annotation is only useful with paired-end sequencing data.</p>
 *
 * <p>For example, if the left read in a pair starts at position 100 and the right ends at position 300 (inclusive) then the fragment length is 201 regardless of which is read 1 and which is read 2, i.e. it is <i>not</i> -201.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Median fragment length of reads supporting each allele (MFRL)")
public class FragmentLength extends PerAlleleAnnotation implements StandardMutectAnnotation {
    public static final String KEY = "MFRL";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : MathUtils.median(Ints.toArray(values));
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
