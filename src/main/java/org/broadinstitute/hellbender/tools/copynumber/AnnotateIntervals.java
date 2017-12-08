package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotationSet;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Annotate intervals with GC content.  The output may optionally be used as input to
 * {@link CreateReadCountPanelOfNormals} or {@link DenoiseReadCounts}.  In the former case,
 * using the resulting panel as input to {@link DenoiseReadCounts} will perform explicit GC-bias correction.
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>
 *         Reference file.
 *     </li>
 *     <li>
 *         Intervals to be annotated.
 *         The argument {@code interval-merging-rule} must be set to {@link IntervalMergingRule#OVERLAPPING_ONLY}
 *         and all other common arguments for interval padding or merging must be set to their defaults.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         GC-content annotated-intervals file.
 *         This is a TSV with a SAM-style header containing a sequence dictionary,
 *         a row specifying the column headers contained in {@link AnnotatedIntervalCollection.AnnotatedIntervalTableColumn},
 *         and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Examples</h3>
 *
 * <pre>
 *     gatk AnnotateIntervals \
 *          -R reference.fa \
 *          -L intervals.interval_list \
 *          --interval-merging-rule OVERLAPPING_ONLY \
 *          -O annotated_intervals.tsv
 * </pre>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Annotate intervals with GC content.",
        oneLineSummary = "Annotate intervals with GC content.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class AnnotateIntervals extends GATKTool {
    @Argument(
            doc = "Output file for annotated intervals.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected File outputAnnotatedIntervalsFile;

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
    private final GCContentAnnotator gcContentAnnotator = new GCContentAnnotator();
    private AnnotatedIntervalCollection annotatedIntervals;

    @Override
    public void onTraversalStart() {
        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);

        logger.info("Loading intervals for annotation...");
        sequenceDictionary = getBestAvailableSequenceDictionary();
        intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
        reference = ReferenceDataSource.of(referenceArguments.getReferenceFile());  //the GATKTool ReferenceDataSource is package-protected, so we cannot access it directly
        logger.info("Annotating intervals...");
    }

    @Override
    public void traverse() {
        final List<AnnotatedInterval> annotatedIntervalList = new ArrayList<>(intervals.size());
        intervals.forEach(interval -> {
            annotatedIntervalList.add(new AnnotatedInterval(
                    interval,
                    new AnnotationSet(gcContentAnnotator.apply(
                            interval, null, new ReferenceContext(reference, interval), null))));
            progressMeter.update(interval);
        });
        annotatedIntervals = new AnnotatedIntervalCollection(new SimpleLocatableMetadata(sequenceDictionary), annotatedIntervalList);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info(String.format("Writing annotated intervals to %s...", outputAnnotatedIntervalsFile));
        annotatedIntervals.write(outputAnnotatedIntervalsFile);
        return super.onTraversalSuccess();
    }

    //if additional annotators are added to this tool, they should follow this interface
    //(and validation that the required resources are available should be performed)
    private interface IntervalAnnotator<T> {
        T apply(final Locatable interval,
                final ReadsContext readsContext,
                final ReferenceContext referenceContext,
                final FeatureContext featureContext);
    }

    private class GCContentAnnotator implements IntervalAnnotator<Double> {
        @Override
        public Double apply(final Locatable interval,
                            final ReadsContext readsContext,
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
}
