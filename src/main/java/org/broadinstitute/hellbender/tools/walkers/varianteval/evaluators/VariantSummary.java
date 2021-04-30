package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

@Analysis(description = "1000 Genomes Phase I summary of variants table")
public class VariantSummary extends VariantEvaluator implements StandardEval {
    final protected static Logger logger = LogManager.getLogger(VariantSummary.class);

    /** Indels with size greater than this value are tallied in the CNV column */
    private final static int MAX_INDEL_LENGTH = 50;
    private final static double MIN_CNV_OVERLAP = 0.5;

    public enum Type {
        SNP, INDEL, CNV
    }

    public VariantSummary(VariantEvalEngine engine) {
        super(engine);

        this.knownCNVs = getEngine().getVariantEvalArgs().getKnownCNVsFile();

        nSamples = getEngine().getSampleNamesForEvaluation().size();
        countsPerSample = new TypeSampleMap(getEngine().getSampleNamesForEvaluation());
        transitionsPerSample = new TypeSampleMap(getEngine().getSampleNamesForEvaluation());
        transversionsPerSample = new TypeSampleMap(getEngine().getSampleNamesForEvaluation());
        allVariantCounts = new TypeSampleMap(getEngine().getSampleNamesForEvaluation());
        knownVariantCounts = new TypeSampleMap(getEngine().getSampleNamesForEvaluation());
        depthPerSample = new TypeSampleMap(getEngine().getSampleNamesForEvaluation());
    }

    private FeatureInput<Feature> knownCNVs = null;

    // basic counts on various rates found
    @DataPoint(description = "Number of samples", format = "%d")
    public long nSamples = 0;

    @DataPoint(description = "Number of processed loci", format = "%d")
    public long nProcessedLoci = 0;

    @DataPoint(description = "Number of SNPs", format = "%d")
    public long nSNPs = 0;
    @DataPoint(description = "Overall TiTv ratio", format = "%.2f")
    public double TiTvRatio = 0;
    @DataPoint(description = "SNP Novelty Rate", format = "%s")
    public String SNPNoveltyRate = "NA";
    @DataPoint(description = "Mean number of SNPs per individual", format = "%d")
    public long nSNPsPerSample = 0;
    @DataPoint(description = "Mean TiTv ratio per individual", format = "%.2f")
    public double TiTvRatioPerSample = 0;
    @DataPoint(description = "Mean depth of coverage per sample at SNPs", format = "%.1f")
    public double SNPDPPerSample = 0;

    @DataPoint(description = "Number of Indels", format = "%d")
    public long nIndels = 0;
    @DataPoint(description = "Indel Novelty Rate", format = "%s")
    public String IndelNoveltyRate = "NA";
    @DataPoint(description = "Mean number of Indels per individual", format = "%d")
    public long nIndelsPerSample = 0;
    @DataPoint(description = "Mean depth of coverage per sample at Indels", format = "%.1f")
    public double IndelDPPerSample = 0;

    @DataPoint(description = "Number of SVs", format = "%d")
    public long nSVs = 0;
    @DataPoint(description = "SV Novelty Rate", format = "%s")
    public String SVNoveltyRate = "NA";
    @DataPoint(description = "Mean number of SVs per individual", format = "%d")
    public long nSVsPerSample = 0;

    private TypeSampleMap allVariantCounts, knownVariantCounts;
    private TypeSampleMap countsPerSample;
    private TypeSampleMap transitionsPerSample, transversionsPerSample;
    private TypeSampleMap depthPerSample;

    private final static String ALL = "ALL";

    private class TypeSampleMap extends EnumMap<Type, Map<String, Integer>> {
        private static final long serialVersionUID = 1L;

        public TypeSampleMap(final Collection<String> samples) {
            super(Type.class);
            for ( Type type : Type.values() ) {
                Map<String, Integer> bySample = new HashMap<String, Integer>(samples.size());
                for ( final String sample : samples ) {
                    bySample.put(sample, 0);
                }
                bySample.put(ALL, 0);
                this.put(type, bySample);
            }
        }

        public final void inc(final Type type, final String sample) {
            final int count = this.get(type).get(sample);
            get(type).put(sample, count + 1);
        }

        public final int all(Type type) {
            return get(type).get(ALL);
        }

        public final int meanValue(Type type) {
            long sum = 0;
            int n = 0;
            for ( final Map.Entry<String, Integer> pair : get(type).entrySet() ) {
                if ( pair.getKey() != ALL)  { // truly must be string ==
                    n++;
                    sum += pair.getValue();
                }
            }
            return (int)(Math.round(sum / (1.0 * n)));
        }

        public final double ratioValue(Type type, TypeSampleMap denoms, boolean allP) {
            double sum = 0;
            int n = 0;
            for ( final String sample : get(type).keySet() ) {
                if ( (allP && sample == ALL) || (!allP && sample != ALL) ) { // truly must be string ==
                    final long num = get(type).get(sample);
                    final long denom = denoms.get(type).get(sample);
                    sum += ratio(num, denom);
                    n++;
                }
            }

            return n > 0 ? sum / (1.0 * n) : 0.0;
        }
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    private Type getType(VariantContext vc) {
        switch (vc.getType()) {
            case SNP:
                return Type.SNP;
            case INDEL:
                for ( int l : vc.getIndelLengths() )
                    if ( Math.abs(l) > MAX_INDEL_LENGTH )
                        return Type.CNV;
                return Type.INDEL;
            case SYMBOLIC:
                return Type.CNV;
            default:
                //throw new UserException.BadInput("Unexpected variant context type: " + vc);
                return null;
        }
    }

    private boolean overlapsKnownCNV(VariantContext cnv, VariantEvalContext context) {
        if ( knownCNVs != null ) {
            List<Feature> overlaps = context.queryFeaturesIncludingOverlapping(knownCNVs, new SimpleInterval(cnv.getContig(), cnv.getStart(), cnv.getEnd()));
            GenomeLocParser parser = new GenomeLocParser(context.getSequenceDictionaryForDrivingVariants());
            GenomeLoc loc1 = parser.createGenomeLoc(cnv);
            for (Feature vc : overlaps) {
                GenomeLoc loc2 = parser.createGenomeLoc(vc);
                final double overlapP = loc1.reciprocialOverlapFraction(loc2);
                if ( overlapP > MIN_CNV_OVERLAP )
                    return true;
            }
        }

        return false;
    }

    @Override
    public void update2(final VariantContext eval, final VariantContext comp, final VariantEvalContext context) {
        if ( eval == null || (getEngine().getVariantEvalArgs().ignoreAC0Sites() && eval.isMonomorphicInSamples()) )
            return;

        final Type type = getType(eval);
        if ( type == null )
            return;

        TypeSampleMap titvTable = null;

        // update DP, if possible
        if ( eval.hasAttribute(VCFConstants.DEPTH_KEY) )
            depthPerSample.inc(type, ALL);

        // update counts
        allVariantCounts.inc(type, ALL);

        // type specific calculations
        if ( type == Type.SNP && eval.isBiallelic() ) {
            titvTable = GATKVariantContextUtils.isTransition(eval) ? transitionsPerSample : transversionsPerSample;
            titvTable.inc(type, ALL);
        }

        // novelty calculation
        if ( comp != null || (type == Type.CNV && overlapsKnownCNV(eval, context)))
            knownVariantCounts.inc(type, ALL);

        // per sample metrics
        for (final Genotype g : eval.getGenotypes()) {
            if ( ! g.isNoCall() && ! g.isHomRef() ) {
                countsPerSample.inc(type, g.getSampleName());

                // update transition / transversion ratio
                if ( titvTable != null ) titvTable.inc(type, g.getSampleName());

                if ( g.hasDP() )
                    depthPerSample.inc(type, g.getSampleName());
            }
        }
    }

    private String noveltyRate(Type type) {
        final int all = allVariantCounts.all(type);
        final int known = knownVariantCounts.all(type);
        return Utils.formattedPercent(all - known, all);
    }

    @Override
    public void finalizeEvaluation() {
        nProcessedLoci = getEngine().getnProcessedLoci();
        nSNPs = allVariantCounts.all(Type.SNP);
        nIndels = allVariantCounts.all(Type.INDEL);
        nSVs = allVariantCounts.all(Type.CNV);

        TiTvRatio = transitionsPerSample.ratioValue(Type.SNP, transversionsPerSample, true);
        TiTvRatioPerSample = transitionsPerSample.ratioValue(Type.SNP, transversionsPerSample, false);

        nSNPsPerSample = countsPerSample.meanValue(Type.SNP);
        nIndelsPerSample = countsPerSample.meanValue(Type.INDEL);
        nSVsPerSample = countsPerSample.meanValue(Type.CNV);

        SNPNoveltyRate = noveltyRate(Type.SNP);
        IndelNoveltyRate = noveltyRate(Type.INDEL);
        SVNoveltyRate = noveltyRate(Type.CNV);

        SNPDPPerSample = depthPerSample.meanValue(Type.SNP);
        IndelDPPerSample = depthPerSample.meanValue(Type.INDEL);
    }

    @Override
    public boolean requiresTerritoryToBeSpecified() {
        return true;
    }
}