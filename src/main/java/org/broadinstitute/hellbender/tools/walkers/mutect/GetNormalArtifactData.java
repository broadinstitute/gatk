package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;
import org.broadinstitute.hellbender.tools.walkers.contamination.PileupSummary;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk GetNormalArtifactData \
 *   -I tumor.bam \
 *   -I normal.bam \
 *   -normal normal_sample \
 *   -L intervals.list \
 *   -O normal-artifact.table
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Collects data for training normal artifact filter",
        oneLineSummary = "Collects data for training normal artifact filter",
        programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
public class GetNormalArtifactData extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional=false)
    private File outputTable;

    @Argument(fullName = M2ArgumentCollection.NORMAL_SAMPLE_LONG_NAME, shortName = M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME,
            doc = "BAM sample name of normal.  May be URL-encoded as output by GetSampleName with -encode argument.")
    protected List<String> normalSamples = new ArrayList<>();

    public static final String ERROR_PROB_NAME = "error-prob";
    private static double DEFAULT_ERROR_PROB = 0.001;
    @Argument(fullName = ERROR_PROB_NAME, doc = "Error probability for p-values", optional = true)
    protected double errorProb = DEFAULT_ERROR_PROB;

    private List<NormalArtifactRecord> data = new ArrayList<>();

    private SAMFileHeader header;

    private final Random rng = Utils.getRandomGenerator();

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresIntervals() {
        return false;
    }

    @Override
    public boolean requiresFeatures() {
        return false;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public void onTraversalStart() {
        header = getHeaderForReads();
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final ReadPileup pileup = alignmentContext.getBasePileup();
        final ReadPileup normalPileup = pileup.makeFilteredPileup(pe -> normalSamples.contains(ReadUtils.getSampleName(pe.getRead(), header)));
        final byte refBase = referenceContext.getBase();
        final int[] normalCounts = getBaseCounts(normalPileup, refBase);

        final int bestNormalAllele = MathUtils.maxElementIndex(normalCounts);
        final int normalAltCount = normalCounts[bestNormalAllele];

        // skip cases of no evidence in normal or likely germline
        if (normalAltCount == 0 || normalAltCount > 0.2 * normalPileup.size()) {
            return;
        }

        final ReadPileup tumorPileup = pileup.makeFilteredPileup(pe -> !normalSamples.contains(ReadUtils.getSampleName(pe.getRead(), header)));
        final int tumorAltCount = getBaseCounts(tumorPileup, refBase)[bestNormalAllele];

        // p value for this tumor alt count or greater
        // we don't want to bloat our data with a lot of sites with a sequencing error in the normal and little or nothing in the tumor
        // at the same time, we must include them.  Thus we downsample and record the downsampling in order to upsample later
        // when the p value is not significant
        final double tumorPValue = 1 - new BinomialDistribution(tumorPileup.size(), errorProb).cumulativeProbability(tumorAltCount - 1);
        final double downsampleProb = Math.max(1 - tumorPValue, 0.05);

        if (rng.nextDouble() > downsampleProb) {
            return;
        } else if (tumorAltCount > 0.5 * tumorPileup.size()) {
            return;
        }

        final String type = bestNormalAllele < 4 ? "SNV" : "INDEL";
        data.add(new NormalArtifactRecord(normalAltCount, normalPileup.size(), tumorAltCount, tumorPileup.size(), downsampleProb, type));

    }

    @Override
    public Object onTraversalSuccess() {
        NormalArtifactRecord.writeToFile(data, outputTable);
        return "SUCCESS";
    }

    /**
     * Get counts of A, C, G, T, before insertion start, before deletion start in order, which returns a int[6] vector with counts according
     * to BaseUtils.simpleBaseToBaseIndex for each base, and with indices 4 and 5 for pileup elements preceding insertions
     * and deletions.  Insertions and deletions themselves are not counted to avoid overcounting.
     */
    private static int[] getBaseCounts(final ReadPileup pileup, final byte refBase) {
        final int[] counts = new int[6];

        for (final PileupElement pe : pileup) {
            if (pe.isDeletion()) {
                continue;
            } else if (pe.isBeforeInsertion()) {
                counts[4]++;
            } else if (pe.isBeforeDeletionStart()) {
                counts[5]++;
            } else if (pe.getBase() != refBase) {
                final int index = BaseUtils.simpleBaseToBaseIndex(pe.getBase());
                if (index != -1) {
                    counts[index]++;
                }
            }
        }

        return counts;
    }
}
