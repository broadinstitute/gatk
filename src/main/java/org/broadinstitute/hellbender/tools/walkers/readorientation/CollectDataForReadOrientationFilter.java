package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.primitives.Ints;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.apache.commons.math3.util.Pair;
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

import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactType.F1R2;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactType.F2R1;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientationFilterConstants.*;

/**
 * Created by Takuto Sato on 7/26/17.
 */

@CommandLineProgramProperties(
        summary = "Collect data from a tumor bam for Mutect2 Read Orientation Filter",
        oneLineSummary = "Data collection for Mutect2 Read Orientation Filter",
        programGroup = CoverageAnalysisProgramGroup.class
)

public class CollectDataForReadOrientationFilter extends LocusWalker {
    public static final String ALT_DATA_TABLE_SHORT_NAME = "alt-table";
    public static final String ALT_DATA_TABLE_LONG_NAME = "alt-data-table";

    public static final String REF_SITE_METRICS_SHORT_NAME = "ref-table";
    public static final String REF_SITE_METRICS_LONG_NAME = "ref-histogram-table";

    public static final String MIN_MEDIAN_MQ_SHORT_NAME = "median-mq";

    public static final String MIN_BASE_QUALITY_SHORT_NAME = "bq";
    public static final String MIN_BASE_QUALITY_LONG_NAME = "min-bq";

    @Argument(fullName = MIN_MEDIAN_MQ_SHORT_NAME,
            shortName = MIN_MEDIAN_MQ_SHORT_NAME,
            doc = "skip sites with median mapping quality below this value", optional = true)
    private static int MINIMUM_MEDIAN_MQ = 20;

    @Argument(fullName = MIN_BASE_QUALITY_LONG_NAME,
            shortName = MIN_BASE_QUALITY_SHORT_NAME,
            doc = "exclude bases below this quality from pileup", optional = true)
    private static int MINIMUM_BASE_QUALITY = 20;

    @Argument(fullName = ALT_DATA_TABLE_LONG_NAME,
            shortName = ALT_DATA_TABLE_SHORT_NAME,
            doc = "a tab-separated output table of pileup data over alt sites")
    private static File altDataTable = null;

    @Argument(fullName = REF_SITE_METRICS_LONG_NAME,
            shortName = REF_SITE_METRICS_SHORT_NAME,
            doc = "a metrics file with overall summary metrics and reference context-specific depth histograms")
    private static File refMetricsOutput = null;

    @Argument(fullName = "alt-histograms",
            shortName = "alt-hist",
            doc = "")
    private static File altMetricsOutput = null;

    // For computational efficiency, for each reference context, we build a depth histogram over ref sites
    private static Map<String, Histogram<Integer>> refSiteHistograms = new HashMap<>(ALL_KMERS.size());

    // Maps Context -> (Alt Allele, F1R2) Pair -> Histogram
    private static Map<String, Map<Pair<Nucleotide, ArtifactType>, Histogram<Integer>>> singleAltTransitionHistograms = new HashMap<>(ALL_KMERS.size());

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
        ALL_KMERS.forEach(context -> {
            Histogram<Integer> emptyRefHistogram = new Histogram<>("depth", context);
            emptyRefHistogram.prefillBins(bins);
            refSiteHistograms.put(context, emptyRefHistogram);


            // Initialize for each context the (Alt Allele, Artifact Type) -> Histogram map
            singleAltTransitionHistograms.put(context, new HashMap<>((REGULAR_BASES.size() - 1)*ArtifactType.values().length));

            for (Nucleotide altAllele : REGULAR_BASES){
                // Skip e.g. AGT -> AGT because G is not an alt allele
                if (altAllele == Nucleotide.valueOf(context.substring(MIDDLE_INDEX, MIDDLE_INDEX+1))){
                    continue;
                }

                for (ArtifactType artifactType : ArtifactType.values()){
                    singleAltTransitionHistograms.get(context).put(new Pair<>(altAllele, artifactType),
                            initializeAltHistogram(context, altAllele, artifactType));
                }
            }
        });

        // Intentionally not use try-with-resources so that the writer stays open outside of the try block
        try {
            altTableWriter = new AltSiteRecordTableWriter(altDataTable);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception creating a writer for %s", altDataTable), e);
        }

    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext){
        // referenceContext always comes with a window of a single base, so
        // manually expand the window and get the 3-mer for now.
        // TODO: implement getBasesInInterval() in referenceContext. Maybe simplify to getKmer(int k)?
        // TODO: this is still relevant (10/2). I shouldn't mess with the internal state of the ref context object
        referenceContext.setWindow(ReadOrientationFilterConstants.REF_CONTEXT_PADDING_ON_EACH_SIDE, ReadOrientationFilterConstants.REF_CONTEXT_PADDING_ON_EACH_SIDE);
        final String refContext = new String(referenceContext.getBases());
        if (refContext.contains("N") || refContext.length() != ReadOrientationFilterConstants.REFERENCE_CONTEXT_SIZE) {
            return;
        }

        if (refContext == null){
            logger.warn(String.format("Skipped a site with null reference at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), refContext));
            return;
        }

        final ReadPileup pileup = alignmentContext.getBasePileup().makeFilteredPileup(pe -> pe.getQual() > MINIMUM_BASE_QUALITY);
        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        if (! isPileupGood(pileup, referenceContext)){
            return;
        }

        // Make a copy of base counts and update the counts of ref to -1. Now the maxElementIndex of the array gives us
        // the alt base.
        final Nucleotide refBase = Nucleotide.valueOf(refContext.getBytes()[ReadOrientationFilterConstants.MIDDLE_INDEX]);
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = -1;
        final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        // if the site is ref, we simply update the coverage histogram
        if (referenceSite){
            refSiteHistograms.get(refContext).increment(depth <= ReadOrientationFilterConstants.maxDepthForHistograms ? depth : ReadOrientationFilterConstants.maxDepthForHistograms);
            return;
        }

        // if we got here, we have an alt site
        final Nucleotide altBase = Nucleotide.valueOf(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        final int[] f1R2Counts = ReadOrientationFilterConstants.REGULAR_BASES.stream().mapToInt(base -> pileup.getNumberOfElements(
                pe -> Nucleotide.valueOf(pe.getBase()) == base && ReadUtils.isF1R2(pe.getRead()))).toArray();

        final int sumNonRefDepths = ReadOrientationFilterConstants.REGULAR_BASES.stream().filter(base -> base != refBase)
                .mapToInt(base -> baseCounts[base.ordinal()])
                .sum();


        // If there's only one alt read at depth > minimumDepthForOptimization, we store the data differently to save space
        if (sumNonRefDepths <= 1){
            if (sumNonRefDepths == 0){
                logger.warn(String.format("Encountered an alt site with 0 alt depth at %s", pileup.getLocation().toString()));
            }

            final int sumNonRefF1R2Depths = ReadOrientationFilterConstants.REGULAR_BASES.stream().filter(base -> base != refBase)
                    .mapToInt(base -> f1R2Counts[base.ordinal()])
                    .sum();
            Utils.validate(Arrays.asList(0, 1).contains(sumNonRefF1R2Depths), "sum of non-ref F1R2 depths must be 0 or 1");
            final ArtifactType type = sumNonRefF1R2Depths == 1 ? F1R2 : F2R1;

            singleAltTransitionHistograms.get(refContext).get(new Pair<>(altBase, type))
                    .increment(depth <= ReadOrientationFilterConstants.maxDepthForHistograms ? depth : ReadOrientationFilterConstants.maxDepthForHistograms);
            return;
        }

        try {
            altTableWriter.writeRecord(new AltSiteRecord(alignmentContext.getContig(), alignmentContext.getStart(),
                    refContext, baseCounts, f1R2Counts, depth, altBase));
        } catch (IOException e) {
            throw new UserException("Encountered an IO Exception writing to the alt data table", e);
        }

        return;
    }

    @Override
    public Object onTraversalSuccess() {
        refSiteHistograms.values().forEach(h -> refMetricsFile.addHistogram(h));
        refMetricsFile.write(refMetricsOutput);

        singleAltTransitionHistograms.values().forEach(hs -> hs.values().forEach(h -> altMetricsFile.addHistogram(h)));
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
    private boolean isPileupGood(final ReadPileup pileup, final ReferenceContext referenceContext){
        // This case should not happen, as AlignmentContext should come filtered, but it does happen once in a while

        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        List<Integer> mappingQualities = Ints.asList(pileup.getMappingQuals());

        boolean isIndel = pileup.getNumberOfElements(pe -> pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart()) > 0;

        // If depth (the sum of base counts) is 0 but the pileup is non-empty, that means all the reads
        // have deleted bases at this particular locus
        isIndel = isIndel || depth == 0 && pileup.size() > 0;

        return depth > 0 && ! isIndel && MathUtils.median(mappingQualities) >= MINIMUM_MEDIAN_MQ;

    }
}
