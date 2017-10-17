package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.primitives.Ints;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasUtils;
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
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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

    public static final int REF_CONTEXT_PADDING_ON_EACH_SIDE = 1;

    public static final int REFERENCE_CONTEXT_SIZE = 2 * REF_CONTEXT_PADDING_ON_EACH_SIDE + 1; // aka 3

    public static final int MIDDLE_INDEX = REF_CONTEXT_PADDING_ON_EACH_SIDE;

    // we put reference site depths above this value in the last bin of the histogram
    public static final int MAX_REF_DEPTH = 200;

    // the list of all possible kmers, where k = REFERENCE_CONTEXT_SIZE
    static final List<String> ALL_KMERS = SequenceUtil.generateAllKmers(REFERENCE_CONTEXT_SIZE).stream()
            .map(String::new).collect(Collectors.toList());
    static final List<String> ALL_KMERS_MODULO_REVERSE_COMPLEMENT = ALL_KMERS.stream()
            .map(context -> new TreeSet<>(Arrays.asList(context, SequenceUtil.reverseComplement(context))))
            .distinct()
            .map(s -> s.first().compareTo(s.last()) < 0 ? s.first() : s.last())
            .collect(Collectors.toList());


    // for computational efficiency, for each reference context, we build a depth histogram over ref sites
    private static Map<String, Histogram<Integer>> refSiteHistograms = new HashMap<>(ALL_KMERS.size());

    private AltSiteRecordTableWriter altTableWriter;

    private final MetricsFile<?, Integer> refMetricsFile = getMetricsFile();

    public static final List<Nucleotide> REGULAR_BASES = Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T);

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
        final Integer[] allBins = IntStream.rangeClosed(1, MAX_REF_DEPTH).boxed().toArray( Integer[]::new );
        ALL_KMERS.forEach(context -> {
            Histogram<Integer> emptyHistogram = new Histogram<>("depth", context);
            emptyHistogram.prefillBins(allBins);
            refSiteHistograms.put(context, emptyHistogram);
        });

        // intentionally not use try-with-resources so that the writer stays open outside of the try block
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
        referenceContext.setWindow(REF_CONTEXT_PADDING_ON_EACH_SIDE, REF_CONTEXT_PADDING_ON_EACH_SIDE);
        final String refContext = new String(referenceContext.getBases());
        if (refContext.contains("N") || refContext.length() != REFERENCE_CONTEXT_SIZE) {
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
        final Nucleotide refBase = Nucleotide.valueOf(refContext.getBytes()[MIDDLE_INDEX]);
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = -1;
        final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        // if the site is ref, we simply update the coverage histogram
        if (referenceSite){
            refSiteHistograms.get(refContext).increment(depth <= MAX_REF_DEPTH ? depth : MAX_REF_DEPTH);
            return;
        }

        // if we got here, we have an alt site
        final Nucleotide altBase = Nucleotide.valueOf(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        final int[] altF1R2Counts = REGULAR_BASES.stream().mapToInt(base -> pileup.getNumberOfElements(
                pe -> Nucleotide.valueOf(pe.getBase()) == base && ReadUtils.isF1R2(pe.getRead()))).toArray();

        try {
            altTableWriter.writeRecord(new AltSiteRecord(alignmentContext.getContig(), alignmentContext.getStart(),
                    refContext, baseCounts, altF1R2Counts, depth, altBase));
        } catch (IOException e) {
            throw new UserException("Encountered an IO Exception writing to the alt data table", e);
        }

        return;
    }

    @Override
    public Object onTraversalSuccess() {
        refSiteHistograms.values().forEach(h -> refMetricsFile.addHistogram(h));
        refMetricsFile.write(refMetricsOutput);
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
