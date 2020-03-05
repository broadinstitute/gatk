package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.IntStream;


/**
 * Apply a read-based annotation that reports the number of Ns seen at a given site. This is intended for use on consensus called data.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of Ns at the pileup")
public class CountNs extends InfoFieldAnnotation {
    /**
     * Calculate annotations for each allele based on given VariantContext and likelihoods for a given genotype's sample
     * and add the annotations to the GenotypeBuilder.  By default annotations are only calculated for alt alleles but
     * implementations may override the {@code includeRefAllele()} method.  See parent class docs in {@link GenotypeAnnotation}.
     */

    public Map<String, Object> annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( likelihoods == null ) {
            return Collections.emptyMap();
        }
        long Count = IntStream.range(0, likelihoods.numberOfSamples()).boxed()
                .flatMap(n -> likelihoods.sampleEvidence(n).stream())
                .filter(read -> doesReadHaveN(read, vc)).count();

        return ImmutableMap.of(GATKVCFConstants.N_COUNT_KEY, Count);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.N_COUNT_KEY));
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.N_COUNT_KEY); }

    @VisibleForTesting
    static Boolean doesReadHaveN(final GATKRead read, final VariantContext vc) {
        final Optional<Byte> readBase = ReadUtils.getReadBaseAtReferenceCoordinate(read, vc.getStart());
        return readBase.isPresent() && readBase.get() == 'N';
    }
}
