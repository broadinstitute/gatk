package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.apache.commons.lang.mutable.MutableInt;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 *  Count of read pairs in the F1R2 and F2R1 configurations supporting the reference and alternate alleles
 *
 *  <p>This is an annotation that gathers information about the read pair configuration for the reads supporting each
 *  allele. It can be used along with downstream filtering steps to identify and filter out erroneous variants that occur
 *  with higher frequency in one read pair orientation.</p>
 *
 *  <h3>References</h3>
 *  <p>For more details about the mechanism of oxoG artifact generation, see <a href='http://www.ncbi.nlm.nih.gov/pubmed/23303777' target='_blank'>
 *      "Discovery and characterization of artefactual mutations in deep coverage targeted capture sequencing data due to oxidative DNA damage during sample preparation."
 *  by Costello et al.</a></p>
 */
public final class OxoGReadCounts extends GenotypeAnnotation implements StandardMutectAnnotation {

    private static final Logger logger = LogManager.getLogger(OxoGReadCounts.class);

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.F1R2_KEY, GATKVCFConstants.F2R1_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F1R2_KEY),
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F2R1_KEY));
    }

    @Override
    public void annotate(final ReferenceContext refContext,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final ReadLikelihoods<Allele> likelihoods){
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        if (g == null || likelihoods == null) {
            return;
        }

        final Map<Allele, MutableInt> f1r2Counts = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        final Map<Allele, MutableInt> f2r1Counts = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        Utils.stream(likelihoods.bestAlleles(g.getSampleName()))
                .filter(ba -> ba.isInformative() && isUsableRead(ba.read))
                .forEach(ba -> (isF2R1(ba.read) ? f2r1Counts : f1r2Counts).get(ba.allele).increment());

        final int[] f1r2 = vc.getAlleles().stream().mapToInt(a -> f1r2Counts.get(a).intValue()).toArray();

        final int[] f2r1 = vc.getAlleles().stream().mapToInt(a -> f2r1Counts.get(a).intValue()).toArray();

        gb.attribute(GATKVCFConstants.F1R2_KEY, f1r2);
        gb.attribute(GATKVCFConstants.F2R1_KEY, f2r1);
    }

    /**
     *  Annotate the given variant context with the OxoG read count attributes, directly from the read pileup.
     *
     *  This method may be slow and should be considered EXPERIMENTAL, especially with regard to indels and complex/mixed
     *    variants.
     *
     * @param vc variant context for the genotype.  Necessary so that we can see all alleles.
     * @param gb genotype builder to put the annotations into.
     * @param readPileup pileup of the reads at this vc.  Note that this pileup does not have to match the
     *                   genotype.  In other words, this tool does not check that the pileup was generated from the
     *                   genotype sample.
     */
    public static void annotateSingleVariant(final VariantContext vc, final GenotypeBuilder gb,
                                             final ReadPileup readPileup, int meanBaseQualityCutoff) {
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        // Create a list of unique alleles
        final List<Allele> variantAllelesWithDupes = vc.getAlleles();
        final Set<Allele> alleleSet = new LinkedHashSet<>(variantAllelesWithDupes);
        final List<Allele> variantAlleles = new ArrayList<>(alleleSet);

        // Initialize the mappings
        final Map<Allele, MutableInt> f1r2Counts = variantAlleles.stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        final Map<Allele, MutableInt> f2r1Counts = variantAlleles.stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        final List<Allele> referenceAlleles = variantAlleles.stream().filter(a -> a.isReference() && !a.isSymbolic()).collect(Collectors.toList());
        final List<Allele> altAlleles = variantAlleles.stream().filter(a -> a.isNonReference() && !a.isSymbolic()).collect(Collectors.toList());

        if (referenceAlleles.size() != 1) {
            logger.warn("Number of reference alleles does not equal  for VC: " + vc);
        }

        // We MUST have exactly 1 non-symbolic reference allele and a read pileup,
        if ((referenceAlleles.size() == 1) && (readPileup != null) && !referenceAlleles.get(0).isSymbolic()) {
            final Allele referenceAllele = referenceAlleles.get(0);
            Utils.stream(readPileup)
                    .filter(pe -> isUsableRead(pe.getRead()))
                    .forEach(pe -> incrementCounts(pe, f1r2Counts, f2r1Counts, referenceAllele, altAlleles, meanBaseQualityCutoff));
        }

        final int[] f1r2 = variantAlleles.stream().mapToInt(a -> f1r2Counts.get(a).intValue()).toArray();

        final int[] f2r1 = variantAlleles.stream().mapToInt(a -> f2r1Counts.get(a).intValue()).toArray();

        gb.attribute(GATKVCFConstants.F1R2_KEY, f1r2);
        gb.attribute(GATKVCFConstants.F2R1_KEY, f2r1);
    }

    /**
     *  If the allele is not in the count mappings, then it is not counted.  No exception will be thrown
     *  Modifies count variables in place.
     *
     * @param pileupElement pileup overlapping the alleles
     * @param f1r2Counts a mapping of allele to f1r2 counts
     * @param f2r1Counts a mapping of allele to f2r1 counts
     */
    private static void incrementCounts(final PileupElement pileupElement, final Map<Allele, MutableInt> f1r2Counts,
                                    final Map<Allele, MutableInt> f2r1Counts, final Allele referenceAllele,
                                        final List<Allele> altAlleles, int minBaseQualityCutoff) {

        final Map<Allele, MutableInt> countMap = isF2R1(pileupElement.getRead()) ? f2r1Counts : f1r2Counts;

        final boolean isRef = referenceAllele.basesMatch(getBasesForAlleleInRead(pileupElement, referenceAllele))
                && !pileupElement.isBeforeDeletionStart() && !pileupElement.isBeforeInsertion();

        Allele pileupAllele = null;
        if (!isRef) {

            for (Allele altAllele : altAlleles) {
                final VariantContext.Type variantType = GATKProtectedVariantContextUtils.typeOfVariant(referenceAllele, altAllele);

                if (variantType == VariantContext.Type.INDEL) {
                    if (isIndelInThePileupElement(pileupElement, referenceAllele, altAllele)) {
                        pileupAllele = altAllele;
                    }

                } else if (variantType == VariantContext.Type.MNP || variantType == VariantContext.Type.SNP) {
                    if (altAllele.basesMatch(getBasesForAlleleInRead(pileupElement, altAllele))) {
                        pileupAllele = altAllele;
                    }
                }

            }

        } else {
            pileupAllele = referenceAllele;
        }

        if (pileupAllele == null) {
            return;
        }

        if (getMinBaseQualityForAlleleInRead(pileupElement, pileupAllele) < minBaseQualityCutoff) {
            return;
        }

        if (countMap.containsKey(pileupAllele)) {
            countMap.get(pileupAllele).increment();
        }
    }

    private static boolean isIndelInThePileupElement(final PileupElement pileupElement, final Allele referenceAllele, final Allele altAllele) {
        boolean isAltAlleleInThePileup = false;

        // Check insertion
        if (pileupElement.isBeforeInsertion()) {
            final int insertionLength = pileupElement.getLengthOfImmediatelyFollowingIndel();
            if (insertionLength == pileupElement.getLengthOfImmediatelyFollowingIndel()) {
                final String insertionBases = pileupElement.getBasesOfImmediatelyFollowingInsertion();
                // edge case: ignore a deletion immediately preceding an insertion as p.getBasesOfImmediatelyFollowingInsertion() returns null [EB]
                if (insertionBases != null) {
                    final boolean isMatch = Allele.extend(referenceAllele, insertionBases.getBytes()).basesMatch(altAllele);
                    if (isMatch) {
                        isAltAlleleInThePileup = true;
                    }
                }
            }
        }

        // Check deletion
        if (pileupElement.isBeforeDeletionStart()) {
            final int deletionLength = pileupElement.getLengthOfImmediatelyFollowingIndel();
            if ((referenceAllele.getBases().length - altAllele.getBases().length) == deletionLength) {
                isAltAlleleInThePileup = true;
            }
        }
        return isAltAlleleInThePileup;
    }

    private static byte[] getBasesForAlleleInRead(final PileupElement pileupElement, final Allele allele) {
        return ArrayUtils.subarray(pileupElement.getRead().getBases(), pileupElement.getOffset(), pileupElement.getOffset() + allele.getBases().length);
    }

    private static int getMinBaseQualityForAlleleInRead(final PileupElement pileupElement, final Allele allele) {
        final byte[] alleleBases = allele.getBases();
        final byte[] pileupBaseQualities = ArrayUtils.subarray(pileupElement.getRead().getBaseQualities(), pileupElement.getOffset(), pileupElement.getOffset() + alleleBases.length);
        final OptionalInt minQuality = IntStream.range(0, pileupBaseQualities.length).map(i -> Byte.toUnsignedInt(pileupBaseQualities[i])).min();
        if (!minQuality.isPresent()) {
            return -1;
        } else {
            return minQuality.getAsInt();
        }
    }


    protected static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    protected static boolean isF2R1(final GATKRead read) {
        return read.isReverseStrand() == read.isFirstOfPair();
    }
}