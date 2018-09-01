package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.primitives.Ints;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.tools.walkers.readorientation.AltSiteRecord.AltSiteRecordTableWriter;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientation.F1R2;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientation.F2R1;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2FilterConstants.*;

/**
 * At each genomic locus, count the number of F1R2/F2R1 alt reads.
 * {@link LearnReadOrientationModel} uses the tsv output of this tool
 *
 * <h3>Usage Example</h3>
 *
 * gatk CollectF1R2Counts \
 *   -R GRCh38.fasta \
 *   -I tumor.bam \
 *   -O tumor-artifact-prior-table.tsv \
 *   -alt-table tumor-alt.tsv \
 *   -ref-hist tumor-ref.metrics \
 *   -alt-hist tumor-alt.metrics
 */

@CommandLineProgramProperties(
        summary = "Collect F1R2 read counts for the Mutect2 orientation bias mixture model filter",
        oneLineSummary = "Collect F1R2 read counts for the Mutect2 orientation bias mixture model filter",
        programGroup = CoverageAnalysisProgramGroup.class
)

public class CollectF1R2Counts extends LocusWalker {
    public static final String ALT_DATA_TABLE_LONG_NAME = "alt-table";
    public static final String ALT_DEPTH1_HISTOGRAM_LONG_NAME = "alt-hist";
    public static final String REF_SITE_METRICS_LONG_NAME = "ref-hist";
    public static final String MIN_MEDIAN_MQ_LONG_NAME = "median-mq";
    public static final String MIN_BASE_QUALITY_LONG_NAME = "min-bq";
    public static final String MAX_DEPTH_LONG_NAME = "max-depth";

    @Argument(fullName = MIN_MEDIAN_MQ_LONG_NAME, doc = "skip sites with median mapping quality below this value", optional = true)
    private int MINIMUM_MEDIAN_MQ = 30;

    @Argument(fullName = MIN_BASE_QUALITY_LONG_NAME, doc = "exclude bases below this quality from pileup", optional = true)
    private int MINIMUM_BASE_QUALITY = 20;

    @Argument(fullName = ALT_DATA_TABLE_LONG_NAME, doc = "a tab-separated output table of pileup data over alt sites")
    private File altDataTable = null;

    @Argument(fullName = REF_SITE_METRICS_LONG_NAME, doc = "a metrics file with overall summary metrics and reference context-specific depth histograms")
    private File refMetricsOutput = null;

    @Argument(fullName = ALT_DEPTH1_HISTOGRAM_LONG_NAME, doc = "a histogram of alt sites with alt depth = 1")
    private File altMetricsOutput = null;

    @Argument(fullName = MAX_DEPTH_LONG_NAME, doc = "sites with depth higher than this value will be grouped", optional = true)
    private int maxDepth = F1R2FilterConstants.DEFAULT_MAX_DEPTH;

    // For each reference context, count ref sites in a histogram keyed by depth
    private Map<String, Histogram<Integer>> refSiteHistograms = new HashMap<>(ALL_KMERS.size());

    // Store the total depths of alt sites with alt depth = 1 separately to save memory
    private DepthOneHistograms depthOneAltHistograms;

    private AltSiteRecordTableWriter altTableWriter;

    private final MetricsFile<?, Integer> refMetricsFile = getMetricsFile();

    private final MetricsFile<?, Integer> altMetricsFile = getMetricsFile();

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public void onTraversalStart() {
        // Initialize for each reference the histogram of the counts of reference sites by depth
        ALL_KMERS.forEach(context -> {
            Histogram<Integer> emptyRefHistogram = F1R2FilterUtils.createRefHistogram(context, maxDepth);
            refSiteHistograms.put(context, emptyRefHistogram);
        });

        depthOneAltHistograms = new DepthOneHistograms(maxDepth);
        // Intentionally not use try-with-resources so that the writer stays open outside of the try block
        try {
            altTableWriter = new AltSiteRecordTableWriter(altDataTable);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception creating a writer for %s", altDataTable), e);
        }

    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final int position = referenceContext.getInterval().getStart();
        final String refContext = referenceContext.getKmerAround(position, F1R2FilterConstants.REF_CONTEXT_PADDING);
        if (refContext == null){
            return;
        }
        final Nucleotide refBase = F1R2FilterUtils.getMiddleBase(refContext);

        if (refContext.contains("N") || refContext.length() != F1R2FilterConstants.REFERENCE_CONTEXT_SIZE) {
            return;
        }

        if (refContext == null) {
            logger.warn(String.format("Skipped a site with null reference at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), refContext));
            return;
        }

        final ReadPileup pileup = alignmentContext.getBasePileup().makeFilteredPileup(pe -> pe.getQual() > MINIMUM_BASE_QUALITY);
        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        if (!isPileupGood(pileup)) {
            return;
        }

        // Make a copy of base counts and update the counts of ref to -1. Now the maxElementIndex of the array gives us
        // the alt base.
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = -1;
        final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        // If the site is ref, we simply update the coverage histogram
        if (referenceSite) {
            refSiteHistograms.get(refContext).increment(Math.min(depth, maxDepth));
            return;
        }

        // If we got here, we have an alt site with a single alt base
        final Nucleotide altBase = Nucleotide.decode(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        final int refCount = baseCounts[refBase.ordinal()];
        final int altCount = baseCounts[altBaseIndex];
        Utils.validate(altCount > 0, "We must have a nonzero alt read but got " + altCount);

        final int refF1R2 = pileup.getNumberOfElements(pe -> Nucleotide.decode(pe.getBase()) == refBase && ReadUtils.isF1R2(pe.getRead()));
        final int altF1R2 = pileup.getNumberOfElements(pe -> Nucleotide.decode(pe.getBase()) == altBase && ReadUtils.isF1R2(pe.getRead()));

        if (altCount == 1) {
            final ReadOrientation type = altF1R2 == 1 ? F1R2 : F2R1;
            depthOneAltHistograms.increment(refContext, altBase, type, depth);
            return;
        }

        try {
            altTableWriter.writeRecord(new AltSiteRecord(refContext, refCount, altCount, refF1R2, altF1R2, altBase));
        } catch (IOException e) {
            throw new UserException("Encountered an IO Exception writing to the alt data table", e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        refSiteHistograms.values().forEach(h -> refMetricsFile.addHistogram(h));
        refMetricsFile.write(refMetricsOutput);

        depthOneAltHistograms.getHistograms().forEach(h -> altMetricsFile.addHistogram(h));
        altMetricsFile.write(altMetricsOutput);

        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (altTableWriter != null) {
            try {
                altTableWriter.close();
            } catch (IOException e) {
                throw new UserException("Encountered an IO exception while closing the alt table writer", e);
            }
        }
    }

    /**
     * Use a series of heuristics to detect a bad pileup.
     */
    private boolean isPileupGood(final ReadPileup pileup){
        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        List<Integer> mappingQualities = Ints.asList(pileup.getMappingQuals());

        // If more than 1% of the reads is indel then consider this site an indel
        final int indelThreshold = depth/100;
        boolean isIndel = pileup.getNumberOfElements(pe -> pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart()) > indelThreshold;

        // If depth (the sum of base counts) is 0 but the pileup is non-empty, that means all the reads
        // have deleted bases at this particular locus
        isIndel = isIndel || depth == 0 && pileup.size() > 0;

        return depth > 0 && ! isIndel && MathUtils.median(mappingQualities) >= MINIMUM_MEDIAN_MQ;

    }
}
