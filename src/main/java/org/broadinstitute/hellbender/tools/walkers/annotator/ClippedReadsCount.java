package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Counts number of soft-clipped reads that support each allele.  *
 * <h3> Caveats </h3>
 * This annotation can only be calculated by Mutect2 and HaplotypeCaller. In case
 * the soft clipping is reverted during the flow
 */

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of soft clipped reads per allele")

public class ClippedReadsCount implements GenotypeAnnotation {
    private final static Logger logger = LogManager.getLogger(ClippedReadsCount.class);

    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.SOFT_CLIP_LEFT_COUNT_KEY, GATKVCFConstants.SOFT_CLIP_RIGHT_COUNT_KEY);
    }

    @Override
    public VCFCompoundHeaderLine.SupportedHeaderLineType annotationType() {
        return GenotypeAnnotation.super.annotationType();
    }

    @Override
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return GenotypeAnnotation.super.getDescriptions();
    }

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        List<GATKRead> allReads = likelihoods.sampleEvidence(likelihoods.indexOfSample(g.getSampleName())).stream().collect(Collectors.toList());
        allReads.addAll(likelihoods.filteredSampleEvidence(likelihoods.indexOfSample(g.getSampleName())).stream().collect(Collectors.toList()));
        List<GATKRead> leftClippedReads = allReads.stream().filter( rd -> (rd.getStart() <= rd.getEnd()) && rd.overlaps(vc) && wasReadClipped(rd, false )).collect(Collectors.toList());
        List<GATKRead> rightClippedReads = allReads.stream().filter( rd -> (rd.getStart() <= rd.getEnd()) && rd.overlaps(vc) && wasReadClipped(rd, true )).collect(Collectors.toList());
        ReadPileup leftClippedPileup = new ReadPileup(ref.getInterval(),leftClippedReads);
        Map<Allele, Integer> leftClippedCounts = PileupBasedAlleles.getPileupAlleleCounts(vc, leftClippedPileup);
        final int[] counts = new int[vc.getNAlleles()];
        counts[0] = leftClippedCounts.get(vc.getReference()); //first one in AD is always ref
        for (int i = 0; i < vc.getNAlleles() -1; i++) {
            counts[i + 1] = leftClippedCounts.get(vc.getAlternateAllele(i));
        }
        gb.attribute(getKeyNames().get(0), counts.clone());

        ReadPileup rightClippedPileup = new ReadPileup(ref.getInterval(),rightClippedReads);
        Map<Allele, Integer> rightClippedCounts = PileupBasedAlleles.getPileupAlleleCounts(vc, rightClippedPileup);
        for (int i = 0 ; i < counts.length; i++){
            counts[i] = 0;
        }
        counts[0] = rightClippedCounts.get(vc.getReference()); //first one in AD is always ref
        for (int i = 0; i < vc.getNAlleles() -1; i++) {
            counts[i + 1] = rightClippedCounts.get(vc.getAlternateAllele(i));
        }

        gb.attribute(getKeyNames().get(1), counts);
    }

    //collect reads that are now softclipped or were softclipped before reversion of the softclipping or were softclipped before hardclipping
    //softclipped bases
    private boolean wasReadClipped(final GATKRead read, boolean rightClipping){

        if ((!rightClipping) &&
                read.hasAttribute(ReadClipper.ORIGINAL_SOFTCLIP_TAG) &&
                read.getAttributeAsString(ReadClipper.ORIGINAL_SOFTCLIP_TAG).contains(ReadClipper.LEFT_SOFTCLIPPING_MARK)) {
            return true;
        }
        if ((!rightClipping) && (!read.isUnmapped()) && (read.getCigar().getFirstCigarElement().getOperator()== CigarOperator.SOFT_CLIP)){
            return true;
        }

        if ((rightClipping) &&
                read.hasAttribute(ReadClipper.ORIGINAL_SOFTCLIP_TAG) &&
                read.getAttributeAsString(ReadClipper.ORIGINAL_SOFTCLIP_TAG).contains(ReadClipper.RIGHT_SOFTCLIPPING_MARK)) {
            return true;
        }
        if ((rightClipping) && (!read.isUnmapped()) && (read.getCigar().getLastCigarElement().getOperator()== CigarOperator.SOFT_CLIP)){
            return true;
        }
        return false;
    }
}
