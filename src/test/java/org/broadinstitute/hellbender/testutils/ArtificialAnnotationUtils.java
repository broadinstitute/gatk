package org.broadinstitute.hellbender.testutils;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 10/16/16.
 */
public class ArtificialAnnotationUtils {
    private static final double MATCH_LIKELIHOOD = -1.0;
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final String SAMPLE = "sample1";

    public static GATKRead makeRead(final int qual, final int mappingQuality) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
        read.setMappingQuality(mappingQuality);
        read.setBaseQualities(Utils.dupBytes((byte) qual, 10));
        return read;
    }

    public static VariantContext makeVC() {
        final GenotypesContext testGC = GenotypesContext.create(2);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
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
        final ReadLikelihoods<Allele> likelihoods = initializeReadLikelihoods(sample, new IndexedAlleleList<>(Arrays.asList(refAllele, altAllele)), reads);

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

    public static ReadLikelihoods<Allele> makeTriAllelicLikelihoods(final String sample,
                                                                    final List<GATKRead> refReads,
                                                                    final List<GATKRead> alt1Reads,
                                                                    final List<GATKRead> alt2Reads,
                                                                    final List<GATKRead> uninformativeReads,
                                                                    final double refReadAltLikelihood,
                                                                    final double alt1ReadRefLikelihood,
                                                                    final double alt2ReadRefLikelihood,
                                                                    final double badReadAltLikelihood,
                                                                    final Allele refAllele,
                                                                    final Allele alt1Allele,
                                                                    final Allele alt2Allele) {
        final List<GATKRead> reads = ListUtils.union(ListUtils.union(refReads, ListUtils.union(alt1Reads,alt2Reads)), uninformativeReads);
        final ReadLikelihoods<Allele> likelihoods = initializeReadLikelihoods(sample, new IndexedAlleleList<>(Arrays.asList(refAllele, alt1Allele, alt2Allele)), reads);

        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        int readIndex = 0;
        for (int i = 0; i < refReads.size(); i++) {
            matrix.set(0, readIndex, MATCH_LIKELIHOOD);
            matrix.set(1, readIndex, refReadAltLikelihood);
            matrix.set(2, readIndex, refReadAltLikelihood);
            readIndex++;
        }

        for (int i = 0; i < alt1Reads.size(); i++) {
            matrix.set(0, readIndex, alt1ReadRefLikelihood);
            matrix.set(1, readIndex, MATCH_LIKELIHOOD);
            matrix.set(2, readIndex, alt1ReadRefLikelihood);
            readIndex++;
        }

        for (int i = 0; i < alt2Reads.size(); i++) {
            matrix.set(0, readIndex, alt2ReadRefLikelihood);
            matrix.set(1, readIndex, alt2ReadRefLikelihood);
            matrix.set(2, readIndex, MATCH_LIKELIHOOD);
            readIndex++;
        }

        for (int i = 0; i < uninformativeReads.size(); i++) {
            matrix.set(0, readIndex, MATCH_LIKELIHOOD);
            matrix.set(1, readIndex, badReadAltLikelihood);
            matrix.set(2, readIndex, badReadAltLikelihood);
            readIndex++;
        }

        return likelihoods;
    }

    private static ReadLikelihoods<Allele> initializeReadLikelihoods(String sample, AlleleList<Allele> alleleList, List<GATKRead> reads) {
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample));
        return new ReadLikelihoods<>(sampleList, alleleList, readsBySample);
    }

    // TODO Add a way to get the likelihoods from this path
    public static VariantContext makeHetAlleleVariantContext(final int refDepth, final int altDepth) {
        final int[] expectedAD = {refDepth, altDepth};

        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> refReads = IntStream.range(0, refDepth).mapToObj(i -> makeRead(30, 5)).collect(Collectors.toList());
        final List<GATKRead> altReads = IntStream.range(0, altDepth).mapToObj(i -> makeRead(30, 5)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, -100.0, -100.0, REF, ALT);

        return new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();
    }
}
