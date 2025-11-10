package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of MQ0 reads per allele")
public class MappingQuality0Count implements GenotypeAnnotation{
    private final static Logger logger = LogManager.getLogger(ClippedReadsCount.class);

    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.MQ0_COUNT_KEY);
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
        List<GATKRead> mq0Reads = allReads.stream().filter( rd -> (rd.getStart() <= rd.getEnd()) && rd.overlaps(vc) && (rd.getMappingQuality() == 0)).collect(Collectors.toList());
        ReadPileup mq0Pileup = new ReadPileup(ref.getInterval(),mq0Reads);
        Map<Allele, Integer> mq0Counts = PileupBasedAlleles.getPileupAlleleCounts(vc, mq0Pileup);
        final int[] counts = new int[vc.getNAlleles()];
        counts[0] = mq0Counts.get(vc.getReference()); //first one in AD is always ref
        for (int i = 0; i < vc.getNAlleles() -1; i++) {
            counts[i + 1] = mq0Counts.get(vc.getAlternateAllele(i));
        }
        gb.attribute(getKeyNames().get(0), counts);
    }
}
