package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.collections4.IteratorUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Annotates intervals with GC content, and optionally, mappability and segmental-duplication content.
 * The output may optionally be used as input to {@link CreateReadCountPanelOfNormals}, {@link DenoiseReadCounts},
 * and {@link GermlineCNVCaller}.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Reference FASTA file
 *     </li>
 *     <li>
 *         Intervals to be annotated. Supported formats are described in
 *         <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=1319">Article#1319</a>.
 *         The argument {@code interval-merging-rule} must be set to {@link IntervalMergingRule#OVERLAPPING_ONLY}
 *         and all other common arguments for interval padding or merging must be set to their defaults.
 *     </li>
 *     <li>
 *         (Optional) Umap single-read mappability track.
 *         This is a BED file in .bed or .bed.gz format that identifies uniquely mappable regions of the genome.
 *         The track should correspond to the appropriate read length and overlapping intervals must be merged.
 *         See <a href ="https://bismap.hoffmanlab.org/">https://bismap.hoffmanlab.org/</a>.  If scores are provided,
 *         intervals will be annotated with the length-weighted average; note that NaN scores will be taken as unity.
 *         Otherwise, scores for covered and uncovered intervals will be taken as unity and zero, respectively.
 *     </li>
 *     <li>
 *         (Optional) Segmental-duplication track.
 *         This is a BED file in .bed or .bed.gz format that identifies segmental-duplication regions of the genome.
 *         Overlapping intervals must be merged.  If scores are provided, intervals will be annotated with the
 *         length-weighted average; note that NaN scores will be taken as unity. Otherwise, scores for covered and
 *         uncovered intervals will be taken as unity and zero, respectively.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Annotated-intervals file.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a sequence dictionary,
 *         a row specifying the column headers for the contained annotations (see {@link CopyNumberAnnotations}
 *         for possible annotations), and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk AnnotateIntervals \
 *          -R reference.fa \
 *          -L intervals.interval_list \
 *          --interval-merging-rule OVERLAPPING_ONLY \
 *          -O annotated_intervals.tsv
 * </pre>
 *
 * <pre>
 *     gatk AnnotateIntervals \
 *          -R reference.fa \
 *          -L intervals.interval_list \
 *          --mappability-track mappability.bed.gz \
 *          --segmental-duplication-track segmental_duplication.bed.gz \
 *          --interval-merging-rule OVERLAPPING_ONLY \
 *          -O annotated_intervals.tsv
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Annotates intervals with GC content, mappability, and segmental-duplication content",
        oneLineSummary = "Annotates intervals with GC content, mappability, and segmental-duplication content",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class AnnotateIntervals extends GATKTool {
    private static final int DEFAULT_FEATURE_QUERY_LOOKAHEAD_IN_BP = 1_000_000;

    public static final String MAPPABILITY_TRACK_PATH_LONG_NAME = "mappability-track";
    public static final String SEGMENTAL_DUPLICATION_TRACK_PATH_LONG_NAME = "segmental-duplication-track";
    public static final String FEATURE_QUERY_LOOKAHEAD = "feature-query-lookahead";

    @Argument(
            doc = "Output file for annotated intervals.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputAnnotatedIntervalsFile;

    @Argument(
            doc = "Path to Umap single-read mappability track in .bed or .bed.gz format (see https://bismap.hoffmanlab.org/).  " +
                    "Overlapping intervals must be merged.",
            fullName = MAPPABILITY_TRACK_PATH_LONG_NAME,
            optional = true
    )
    private FeatureInput<BEDFeature> mappabilityTrackPath;

    @Argument(
            doc = "Path to segmental-duplication track in .bed or .bed.gz format.  " +
                    "Overlapping intervals must be merged.",
            fullName = SEGMENTAL_DUPLICATION_TRACK_PATH_LONG_NAME,
            optional = true
    )
    private FeatureInput<BEDFeature> segmentalDuplicationTrackPath;

    @Argument(
            doc = "Number of bases to cache when querying feature tracks.",
            fullName = FEATURE_QUERY_LOOKAHEAD,
            optional = true
    )
    private int featureQueryLookahead = DEFAULT_FEATURE_QUERY_LOOKAHEAD_IN_BP;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    private List<SimpleInterval> intervals;
    private SAMSequenceDictionary sequenceDictionary;
    private ReferenceDataSource reference;
    private FeatureManager features;
    private List<IntervalAnnotator<?>> annotators = new ArrayList<>();
    private AnnotatedIntervalCollection annotatedIntervals;

    @Override
    public void onTraversalStart() {
        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);

        logger.info("Loading intervals for annotation...");
        sequenceDictionary = getBestAvailableSequenceDictionary();
        intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);

        logger.info("Loading resources for annotation...");
        reference = ReferenceDataSource.of(referenceArguments.getReferencePath());  //the GATKTool ReferenceDataSource is package-protected, so we cannot access it directly
        features = new FeatureManager(                                              //the GATKTool FeatureManager is package-protected, so we cannot access it directly
                this,
                featureQueryLookahead,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer,
                referenceArguments.getReferencePath());

        // always perform GC-content annotation
        logger.info("Adding GC-content annotator...");
        annotators.add(new GCContentAnnotator());

        // add optional annotators
        if (mappabilityTrackPath != null) {
            checkForOverlaps(features, mappabilityTrackPath, sequenceDictionary);
            logger.info("Adding mappability annotator...");
            annotators.add(new MappabilityAnnotator(mappabilityTrackPath));
        }
        if (segmentalDuplicationTrackPath != null) {
            checkForOverlaps(features, segmentalDuplicationTrackPath, sequenceDictionary);
            logger.info("Adding segmental-duplication-content annotator...");
            annotators.add(new SegmentalDuplicationContentAnnotator(segmentalDuplicationTrackPath));
        }

        logger.info("Annotating intervals...");
    }

    private static void checkForOverlaps(final FeatureManager featureManager,
                                         final FeatureInput<BEDFeature> featureTrackPath,
                                         final SAMSequenceDictionary sequenceDictionary) {
        final List<BEDFeature> features = IteratorUtils.toList(featureManager.getFeatureIterator(featureTrackPath));
        final GenomeLocParser parser = new GenomeLocParser(sequenceDictionary);
        final List<GenomeLoc> genomeLocs = IntervalUtils.genomeLocsFromLocatables(parser, features);
        final List<GenomeLoc> mergedGenomeLocs = IntervalUtils.mergeIntervalLocations(
                genomeLocs, IntervalMergingRule.OVERLAPPING_ONLY);
        if (genomeLocs.size() != mergedGenomeLocs.size()) {
            throw new UserException.BadInput(String.format("Feature track %s contains overlapping intervals; " +
                    "these should be merged.", featureTrackPath));
        }
    }
    @Override
    public void traverse() {
        final List<AnnotatedInterval> annotatedIntervalList = new ArrayList<>(intervals.size());
        for (final SimpleInterval interval : intervals) {
            if (interval.getLengthOnReference() == 0) {
                throw new UserException.BadInput(String.format("Interval cannot have zero length: %s", interval));
            }
            final ReferenceContext referenceContext = new ReferenceContext(reference, interval);
            final FeatureContext featureContext = new FeatureContext(features, interval);
            final AnnotationMap annotations = new AnnotationMap(annotators.stream()
                    .collect(Collectors.mapping(
                            a -> Pair.of(
                                    a.getAnnotationKey(),
                                    a.applyAndValidate(interval, referenceContext, featureContext)),
                            Collectors.toList())));
            annotatedIntervalList.add(new AnnotatedInterval(interval, annotations));
            progressMeter.update(interval);
        }
        annotatedIntervals = new AnnotatedIntervalCollection(new SimpleLocatableMetadata(sequenceDictionary), annotatedIntervalList);
    }

    @Override
    public Object onTraversalSuccess() {
        reference.close();
        features.close();
        logger.info(String.format("Writing annotated intervals to %s...", outputAnnotatedIntervalsFile));
        annotatedIntervals.write(outputAnnotatedIntervalsFile);
        return super.onTraversalSuccess();
    }

    /**
     * If additional annotators are added to this tool, they should follow this interface.
     * Validation that the required resources are available should be performed before
     * calling {@link IntervalAnnotator#apply}.
     */
    abstract static class IntervalAnnotator<T> {
        public abstract AnnotationKey<T> getAnnotationKey();

        /**
         * @param interval  assumed to have non-zero length
         */
        abstract T apply(final Locatable interval,
                         final ReferenceContext referenceContext,
                         final FeatureContext featureContext);

        /**
         * @param interval  assumed to have non-zero length
         */
        T applyAndValidate(final Locatable interval,
                           final ReferenceContext referenceContext,
                           final FeatureContext featureContext) {
            return getAnnotationKey().validate(apply(interval, referenceContext, featureContext));
        }
    }

    private static class GCContentAnnotator extends IntervalAnnotator<Double> {
        @Override
        public AnnotationKey<Double> getAnnotationKey() {
            return CopyNumberAnnotations.GC_CONTENT;
        }

        /**
         * @param interval  assumed to have non-zero length
         */
        @Override
        Double apply(final Locatable interval,
                     final ReferenceContext referenceContext,
                     final FeatureContext featureContext) {
            final Nucleotide.Counter counter = new Nucleotide.Counter();
            counter.addAll(referenceContext.getBases());
            final long gcCount = counter.get(Nucleotide.C) + counter.get(Nucleotide.G);
            final long atCount = counter.get(Nucleotide.A) + counter.get(Nucleotide.T);
            final long totalCount = gcCount + atCount;
            return totalCount == 0 ? Double.NaN : gcCount / (double) totalCount;
        }
    }

    /**
     * If scores are provided, intervals will be annotated with the length-weighted average; note that NaN scores will
     * be taken as unity.  Otherwise, scores for covered and uncovered intervals will be taken as unity and zero,
     * respectively.
     */
    abstract static class BEDLengthWeightedAnnotator extends IntervalAnnotator<Double> {
        private final FeatureInput<BEDFeature> trackPath;

        BEDLengthWeightedAnnotator(final FeatureInput<BEDFeature> trackPath) {
            this.trackPath = trackPath;
        }

        /**
         * @param interval  assumed to have non-zero length
         */
        @Override
        Double apply(final Locatable interval,
                     final ReferenceContext referenceContext,
                     final FeatureContext featureContext) {
            double lengthWeightedSum = 0.;
            final List<BEDFeature> features = featureContext.getValues(trackPath);
            for (final BEDFeature feature : features) {
                final double scoreOrNaN = (double) feature.getScore();
                final double score = Double.isNaN(scoreOrNaN) ? 1. : scoreOrNaN;    // missing or NaN score -> score = 1
                lengthWeightedSum += score *
                        CoordMath.getOverlap(
                                feature.getStart(), feature.getEnd() - 1,       // zero-based
                                interval.getStart(), interval.getEnd());        // one-based
            }
            return lengthWeightedSum / interval.getLengthOnReference();
        }
    }

    private static class MappabilityAnnotator extends BEDLengthWeightedAnnotator {
        MappabilityAnnotator(final FeatureInput<BEDFeature> mappabilityTrackPath) {
            super(mappabilityTrackPath);
        }

        @Override
        public AnnotationKey<Double> getAnnotationKey() {
            return CopyNumberAnnotations.MAPPABILITY;
        }
    }

    private static class SegmentalDuplicationContentAnnotator extends BEDLengthWeightedAnnotator {
        SegmentalDuplicationContentAnnotator(final FeatureInput<BEDFeature> segmentalDuplicationTrackPath) {
            super(segmentalDuplicationTrackPath);
        }

        @Override
        public AnnotationKey<Double> getAnnotationKey() {
            return CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT;
        }
    }
}
