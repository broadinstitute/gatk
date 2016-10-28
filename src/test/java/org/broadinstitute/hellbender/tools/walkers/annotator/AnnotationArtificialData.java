package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Created by davidben on 10/16/16.
 */
public class AnnotationArtificialData {
    private static final double MATCH_LIKELIHOOD = -1.0;

    public static GATKRead makeRead(final int qual, final int mappingQuality) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
        read.setMappingQuality(mappingQuality);
        read.setBaseQualities(Utils.dupBytes((byte) qual, 10));
        return read;
    }

    public static ReadLikelihoods<Allele> makeLikelihoods(final String sample,
                                                          final List<GATKRead> refReads,
                                                          final double refReadAltLikelihood,
                                                          final Allele refAllele,
                                                          final Allele altAllele) {
        return makeLikelihoods(sample, refReads, Collections.emptyList(), refReadAltLikelihood, 0, refAllele, altAllele);
    }

    public static ReadLikelihoods<Allele> makeLikelihoods(final String sample,
                                                          final List<GATKRead> refReads,
                                                          final List<GATKRead> altReads,
                                                          final double refReadAltLikelihood,
                                                          final double altReadRefLikelihood,
                                                          final Allele refAllele,
                                                          final Allele altAllele) {
        return makeLikelihoods(sample, refReads, altReads, Collections.emptyList(), refReadAltLikelihood, altReadRefLikelihood,
                0, refAllele, altAllele);
    }

    public static ReadLikelihoods<Allele> makeLikelihoods(final String sample,
                                                          final List<GATKRead> refReads,
                                                          final List<GATKRead> altReads,
                                                          final List<GATKRead> uninformativeReads,
                                                          final double refReadAltLikelihood,
                                                          final double altReadRefLikelihood,
                                                          final double badReadAltLikelihood,
                                                          final Allele refAllele,
                                                          final Allele altAllele) {
        final List<GATKRead> reads = ListUtils.union(ListUtils.union(refReads, altReads), uninformativeReads);
        final ReadLikelihoods<Allele> likelihoods = initializeReadLikelihoods(sample, refAllele, altAllele, reads);

        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        int readIndex = 0;
        for (int i = 0; i < refReads.size(); i++) {
            matrix.set(0, readIndex, MATCH_LIKELIHOOD);
            matrix.set(1, readIndex, refReadAltLikelihood);
            readIndex++;
        }

        for (int i = 0; i < altReads.size(); i++) {
            matrix.set(0, readIndex, altReadRefLikelihood);
            matrix.set(1, readIndex, MATCH_LIKELIHOOD);
            readIndex++;
        }

        for (int i = 0; i < uninformativeReads.size(); i++) {
            matrix.set(0, readIndex, MATCH_LIKELIHOOD);
            matrix.set(1, readIndex, badReadAltLikelihood);
            readIndex++;
        }

        return likelihoods;
    }

    private static ReadLikelihoods<Allele> initializeReadLikelihoods(String sample, Allele refAllele, Allele altAllele, List<GATKRead> reads) {
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(refAllele, altAllele));
        return new ReadLikelihoods<>(sampleList, alleleList, readsBySample);
    }
}
