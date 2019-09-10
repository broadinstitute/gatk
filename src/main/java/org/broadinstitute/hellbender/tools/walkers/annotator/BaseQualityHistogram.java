package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.TreeMultiset;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;

public class BaseQualityHistogram extends InfoFieldAnnotation {

    public static final String KEY = "BQHIST";

    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( likelihoods == null ) {
            return Collections.emptyMap();
        }

        final Map<Allele, TreeMultiset<Integer>> values = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> TreeMultiset.create()));

        Utils.stream(likelihoods.bestAllelesBreakingTies())
                .filter(ba -> ba.isInformative() && isUsableRead(ba.evidence))
                .forEach(ba -> BaseQuality.getBaseQuality(ba.evidence, vc).ifPresent(v -> values.get(ba.allele).add(v)));


        final List<Integer> distinctBaseQualities = likelihoods.alleles().stream()
                .flatMap(a -> values.get(a).stream())
                .distinct()
                .sorted()
                .collect(Collectors.toList());

        final List<Integer> output = new ArrayList<>();

        for (final int qual : distinctBaseQualities) {
            output.add(qual);
            likelihoods.alleles().forEach(allele -> output.add(values.get(allele).count(qual)));
        }


        return ImmutableMap.of(KEY, output);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer,
                "Base quality counts for each allele represented sparsely as alternating entries of qualities and counts for each allele." +
                "For example [10,1,0,20,0,1] means one ref base with quality 10 and one alt base with quality 20."));
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(KEY); }

    private static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }
}
