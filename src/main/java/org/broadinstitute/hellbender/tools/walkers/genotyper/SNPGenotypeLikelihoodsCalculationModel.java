package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.*;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.genotyper.DiploidGenotype;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class SNPGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {

    private final boolean useAlleleFromVCF;

    private final double[] likelihoodSums = new double[4];

    private final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap;

    protected SNPGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
        useAlleleFromVCF = UAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES;
        perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
    }

    public VariantContext getLikelihoods(final FeatureContext features,
                                         final ReferenceContext ref,
                                         final Map<String, AlignmentContext> contexts,
                                         final AlignmentContext.ReadOrientation contextType,
                                         final List<Allele> allAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final Map<String, PerReadAlleleLikelihoodMap> sampleLikelihoodMap) {

        sampleLikelihoodMap.clear(); // not used in SNP model, sanity check to delete any older data

        final byte refBase = ref.getBase();
        final int indexOfRefBase = BaseUtils.simpleBaseToBaseIndex(refBase);
        // handle non-standard reference bases
        if ( indexOfRefBase == -1 ) {
            return null;
        }
        final Allele refAllele = Allele.create(refBase, true);

        // calculate the GLs
        ArrayList<SampleGenotypeData> GLs = new ArrayList<>(contexts.size());
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            // Down-sample with bias according to the contamination level (global or per file)
            ReadPileup pileup = sample.getValue().stratify(contextType).getBasePileup();
            final Double contamination =  UAC.getSampleContamination().get(sample.getKey());
            if( contamination > 0.0 ) {//no need to enter if no contamination reduction
                pileup = perReadAlleleLikelihoodMap.createPerAlleleDownsampledBasePileup(pileup, contamination);
            }
            if ( useBAQedPileup ) {
                pileup = createBAQedPileup(pileup);
            }

            // create the GenotypeLikelihoods object
            final DiploidSNPGenotypeLikelihoods GL = new DiploidSNPGenotypeLikelihoods(UAC.PCR_error);
            final int nGoodBases = GL.add(pileup, true, true, UAC.MIN_BASE_QUALTY_SCORE);
            if ( nGoodBases > 0 ) {
                GLs.add(new SampleGenotypeData(sample.getKey(), GL, getFilteredDepth(pileup)));
            }
        }

        // start making the VariantContext
        final SimpleInterval loc = ref.getInterval();
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(refAllele);


        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), loc.getEnd(), alleles);
        // find the alternate allele(s) that we should be using
        if ( allAllelesToUse != null ) {
            alleles.addAll(allAllelesToUse.subList(1,allAllelesToUse.size()));   // this includes ref allele
        } else if ( useAlleleFromVCF ) {
            final VariantContext vc = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(features, ref.getInterval(), true, logger, UAC.alleles);

            // ignore places where we don't have a SNP
            if ( vc == null || !vc.isSNP() ) {
                return null;
            }

            // make sure a user isn't passing the REF base in as an ALT
            if ( vc.hasAlternateAllele(refAllele, true) ) {
                throw new UserException.BadInput("Alternate allele '" + (char) refBase + "' passed in is the same as the reference at location " + vc.getContig() + ":" + vc.getStart());
            }

            alleles.addAll(vc.getAlternateAlleles());
        } else {

            alleles.addAll(determineAlternateAlleles(refBase, GLs));

            // if there are no non-ref alleles...
            if ( alleles.size() == 1 ) {
                // if we only want variants, then we don't need to calculate genotype likelihoods
                if ( UAC.outputMode == OutputMode.EMIT_VARIANTS_ONLY ) {
                    return builder.make();
                }
                else {
                    // otherwise, choose any alternate allele (it doesn't really matter)
                    alleles.add(Allele.create(BaseUtils.baseIndexToSimpleBase(indexOfRefBase == 0 ? 1 : 0)));
                }
            }
        }

        // create the alternate alleles and the allele ordering (the ordering is crucial for the GLs)
        final int numAlleles = alleles.size();
        final int numAltAlleles = numAlleles - 1;

        final int[] alleleOrdering = new int[numAlleles];
        int alleleOrderingIndex = 0;
        int numLikelihoods = 0;
        for ( Allele allele : alleles ) {
            alleleOrdering[alleleOrderingIndex++] = BaseUtils.simpleBaseToBaseIndex(allele.getBases()[0]);
            numLikelihoods += alleleOrderingIndex;
        }
        builder.alleles(alleles);

        // create the PL ordering to use based on the allele ordering.
        final int[] PLordering = new int[numLikelihoods];
        for ( int i = 0; i <= numAltAlleles; i++ ) {
            for ( int j = i; j <= numAltAlleles; j++ ) {
                // As per the VCF spec: "the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.
                // In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
                PLordering[(j * (j+1) / 2) + i] = DiploidGenotype.createDiploidGenotype(alleleOrdering[i], alleleOrdering[j]).ordinal();
            }
        }

        // create the genotypes; no-call everyone for now
        final GenotypesContext genotypes = GenotypesContext.create();
        final int ploidy = UAC.genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy);
        for ( SampleGenotypeData sampleData : GLs ) {
            final double[] allLikelihoods = sampleData.GL.getLikelihoods();
            final double[] myLikelihoods = new double[numLikelihoods];

            for ( int i = 0; i < numLikelihoods; i++ ) {
                myLikelihoods[i] = allLikelihoods[PLordering[i]];
            }

            // normalize in log space so that max element is zero.
            final GenotypeBuilder gb = new GenotypeBuilder(sampleData.name);
            final double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(myLikelihoods, false, true);
            gb.PL(genotypeLikelihoods);
            gb.DP(sampleData.depth);
            gb.alleles(noCall);
            if (UAC.annotateAllSitesWithPLs) {
                gb.attribute(GATKVCFConstants.PL_FOR_ALL_SNP_ALLELES_KEY, GenotypeLikelihoods.fromLog10Likelihoods(MathUtils.normalizeFromLog10(allLikelihoods, false, true)));
            }
            genotypes.add(gb.make());
        }

        return builder.genotypes(genotypes).make();
    }

    // determines the alleles to use
    protected List<Allele> determineAlternateAlleles(final byte ref, final List<SampleGenotypeData> sampleDataList) {

        final int baseIndexOfRef = BaseUtils.simpleBaseToBaseIndex(ref);
        final int PLindexOfRef = DiploidGenotype.createDiploidGenotype(ref, ref).ordinal();
        for ( int i = 0; i < 4; i++ ) {
            likelihoodSums[i] = 0.0;
        }

        // based on the GLs, find the alternate alleles with enough probability
        for ( SampleGenotypeData sampleData : sampleDataList ) {
            final double[] likelihoods = sampleData.GL.getLikelihoods();
            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            if ( PLindexOfBestGL != PLindexOfRef ) {
                GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindexOfBestGL);
                if ( alleles.alleleIndex1 != baseIndexOfRef ) {
                    likelihoodSums[alleles.alleleIndex1] += likelihoods[PLindexOfBestGL] - likelihoods[PLindexOfRef];
                }
                // don't double-count it
                if ( alleles.alleleIndex2 != baseIndexOfRef && alleles.alleleIndex2 != alleles.alleleIndex1 ) {
                    likelihoodSums[alleles.alleleIndex2] += likelihoods[PLindexOfBestGL] - likelihoods[PLindexOfRef];
                }
            }
        }

        final List<Allele> allelesToUse = new ArrayList<Allele>(3);
        for ( int i = 0; i < 4; i++ ) {
            if ( likelihoodSums[i] > 0.0 ) {
                allelesToUse.add(Allele.create(BaseUtils.baseIndexToSimpleBase(i), false));
            }
        }

        return allelesToUse;
    }

    public ReadPileup createBAQedPileup( final ReadPileup pileup ) {
        final List<PileupElement> BAQedElements = new ArrayList<PileupElement>();
        for( final PileupElement PE : pileup ) {
            final PileupElement newPE = new BAQedPileupElement( PE );
            BAQedElements.add( newPE );
        }
        return new ReadPileup(pileup.getLocation(), BAQedElements);
    }

    public static class BAQedPileupElement extends PileupElement {
        public BAQedPileupElement( final PileupElement PE ) {
            super(PE);
        }

        @Override
        public byte getQual() {
            if ( isDeletion() ) {
                return super.getQual();
            }
            else {
                return BAQ.calcBAQFromTag(getRead(), offset, true);
            }
        }
    }

    private static class SampleGenotypeData {

        public final String name;
        public final DiploidSNPGenotypeLikelihoods GL;
        public final int depth;

        public SampleGenotypeData(final String name, final DiploidSNPGenotypeLikelihoods GL, final int depth) {
            this.name = name;
            this.GL = GL;
            this.depth = depth;
        }
    }
}
