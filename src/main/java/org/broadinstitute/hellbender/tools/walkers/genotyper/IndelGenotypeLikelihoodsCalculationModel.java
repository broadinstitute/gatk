package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AlignmentContextUtils;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

public class IndelGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {
    private static final int HAPLOTYPE_SIZE = 80;

    private boolean DEBUG = false;
    private boolean ignoreSNPAllelesWhenGenotypingIndels = false;
    private PairHMMIndelErrorModel pairModel;


    private Map<Allele, Haplotype> haplotypeMap;

    private List<Allele> alleleList = new ArrayList<>();


    protected IndelGenotypeLikelihoodsCalculationModel(final UnifiedArgumentCollection UAC,
                                                       final Logger logger) {
        super(UAC, logger);
        pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                UAC.OUTPUT_DEBUG_INDEL_INFO, UAC.pairHMM);
        DEBUG = UAC.OUTPUT_DEBUG_INDEL_INFO;
        haplotypeMap = new LinkedHashMap<>();
        ignoreSNPAllelesWhenGenotypingIndels = UAC.IGNORE_SNP_ALLELES;
    }

    protected static List<Allele> computeConsensusAlleles(final ReferenceContext ref,
                                                 final Map<String, AlignmentContext> contexts,
                                                 final AlignmentContextUtils.ReadOrientation contextType,
                                                 final UnifiedArgumentCollection UAC) {
        final ConsensusAlleleCounter counter = new ConsensusAlleleCounter(true, UAC.MIN_INDEL_COUNT_FOR_GENOTYPING, UAC.MIN_INDEL_FRACTION_PER_SAMPLE);
        return counter.computeConsensusAlleles(ref, contexts, contextType);
    }

    private final static EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.INDEL, VariantContext.Type.MIXED);


    public VariantContext getLikelihoods(final FeatureContext tracker,
                                         final ReferenceContext ref,
                                         final Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final List<Allele> allAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser,
                                         final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        final Locatable loc = ref.getInterval();
        if (contextType == AlignmentContextUtils.ReadOrientation.COMPLETE) {
            // starting a new site: clear allele list
            haplotypeMap.clear();
            perReadAlleleLikelihoodMap.clear(); // clean mapping sample-> per read, per allele likelihoods
            alleleList = getInitialAlleleList(tracker, ref, contexts, contextType, UAC, ignoreSNPAllelesWhenGenotypingIndels);
            if (alleleList.isEmpty()) {
                return null;
            }
        }

        getHaplotypeMapFromAlleles(alleleList, ref, loc, haplotypeMap); // will update haplotypeMap adding elements
        if (haplotypeMap == null || haplotypeMap.isEmpty()) {
            return null;
        }

        // start making the VariantContext
        // For all non-snp VC types, VC end location is just startLocation + length of ref allele including padding base.
        final int endLoc = loc.getStart() + alleleList.get(0).length() - 1;
        final int eventLength = getEventLength(alleleList);

        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), endLoc, alleleList);

        // create the genotypes; no-call everyone for now
        final GenotypesContext genotypes = GenotypesContext.create();
        final int ploidy = UAC.genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy);

        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.

        for (final Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            final AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            if (!perReadAlleleLikelihoodMap.containsKey(sample.getKey())){
                // no likelihoods have been computed for this sample at this site
                perReadAlleleLikelihoodMap.put(sample.getKey(), new PerReadAlleleLikelihoodMap());
            }
            final ReadPileup pileup = context.getBasePileup();
            if (pileup != null) {
                final GenotypeBuilder b = new GenotypeBuilder(sample.getKey());
                final double[] genotypeLikelihoods = pairModel.computeDiploidReadHaplotypeLikelihoods(pileup, haplotypeMap, ref, eventLength, perReadAlleleLikelihoodMap.get(sample.getKey()), UAC.getSampleContamination(sample.getKey()));
                b.PL(genotypeLikelihoods);
                b.alleles(noCall);
                b.DP(getFilteredDepth(pileup));
                genotypes.add(b.make());

                if (DEBUG) {
                    System.out.format("Sample:%s Alleles:%s GL:", sample.getKey(), alleleList.toString());
                    for (int k = 0; k < genotypeLikelihoods.length; k++) {
                        System.out.format("%1.4f ", genotypeLikelihoods[k]);
                    }
                    System.out.println();
                }
            }
        }

        return builder.genotypes(genotypes).make();
    }

    public static void getHaplotypeMapFromAlleles(final List<Allele> alleleList,
                                                 final ReferenceContext ref,
                                                 final Locatable loc,
                                                 final Map<Allele, Haplotype> haplotypeMap) {
        // protect against having an indel too close to the edge of a contig
        if (loc.getStart() <= HAPLOTYPE_SIZE) {
            haplotypeMap.clear();
        }// check if there is enough reference window to create haplotypes (can be an issue at end of contigs)
        else if (ref.getWindow().getEnd() < loc.getEnd() + HAPLOTYPE_SIZE) {
            haplotypeMap.clear();
        } else if (alleleList.isEmpty()) {
            haplotypeMap.clear();
        } else {
            final int eventLength = getEventLength(alleleList);
            final int hsize = ref.getWindow().size() - Math.abs(eventLength) - 1;
            final int numPrefBases = ref.getInterval().getStart() - ref.getWindow().getStart() + 1;

            if (hsize <= 0)  // protect against event lengths larger than ref window sizes
            {
                haplotypeMap.clear();
            } else {
                haplotypeMap.putAll(Haplotype.makeHaplotypeListFromAlleles(alleleList, loc.getStart(),
                        ref, hsize, numPrefBases));
            }
        }
    }

    public static int getEventLength(final List<Allele> alleleList) {
        final Allele refAllele = alleleList.get(0);
        Allele altAllele = alleleList.get(1);
        // look for alt allele that has biggest length distance to ref allele
        int maxLenDiff = 0;
        for (final Allele a : alleleList) {
            if (a.isNonReference()) {
                final int lenDiff = Math.abs(a.getBaseString().length() - refAllele.getBaseString().length());
                if (lenDiff > maxLenDiff) {
                    maxLenDiff = lenDiff;
                    altAllele = a;
                }
            }
        }

        return altAllele.getBaseString().length() - refAllele.getBaseString().length();

    }
    
    public static List<Allele> getInitialAlleleList(final FeatureContext tracker,
                                                    final ReferenceContext ref,
                                                    final Map<String, AlignmentContext> contexts,
                                                    final AlignmentContextUtils.ReadOrientation contextType,
                                                    final UnifiedArgumentCollection UAC,
                                                    final boolean ignoreSNPAllelesWhenGenotypingIndels) {

        List<Allele> alleles = new ArrayList<>();
        if (UAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES) {
            VariantContext vc = null;
            for (final VariantContext vc_input : tracker.getValues(UAC.alleles, new SimpleInterval(ref.getInterval()))) {
                if (vc_input != null &&
                        allowableTypes.contains(vc_input.getType()) &&
                        ref.getWindow().getStart() == vc_input.getStart()) {
                    vc = vc_input;
                    break;
                }
            }
           // ignore places where we don't have a variant
            if (vc == null) {
                return alleles;
            }

            if (ignoreSNPAllelesWhenGenotypingIndels) {
                // if there's an allele that has same length as the reference (i.e. a SNP or MNP), ignore it and don't genotype it
                for (final Allele a : vc.getAlleles()) {
                    if (a.isNonReference() && a.getBases().length == vc.getReference().getBases().length)
                        continue;
                    else
                        alleles.add(a);
                }

            } else {
                alleles.addAll(vc.getAlleles());
            }

        } else {
            alleles = computeConsensusAlleles(ref, contexts, contextType, UAC);
        }
        return alleles;
    }

    // Overload function in GenotypeLikelihoodsCalculationModel so that, for an indel case, we consider a deletion as part of the pileup,
    // so that per-sample DP will include deletions covering the event.
    protected int getFilteredDepth(final ReadPileup pileup) {
        int count = 0;
        for (final PileupElement p : pileup) {
            if (p.isDeletion() || BaseUtils.isRegularBase(p.getBase())) {
                count++;
            }
        }

        return count;
    }

}
