package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.copynumber.annotation.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.annotation.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.annotation.AnnotationSet;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.utils.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Annotates intervals with GC content.  The output may optionally be used by
 * {@link CreateReadCountPanelOfNormals} and {@link DenoiseReadCounts} to perform explicit GC-bias correction.
 *
 * <h3>Examples</h3>

 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" AnnotateIntervals \
 *   -L intervals.interval_list \
 *   --reference ref_fasta.fa \
 *   --output annotated_intervals.tsv
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Annotate intervals with GC content.  " +
                "Overlapping intervals will be merged, but no other padding or merging specified by command-line arguments is allowed.",
        oneLineSummary = "Annotate intervals with GC content.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class AnnotateIntervals extends GATKTool {
    @Argument(
            doc = "Output annotated-intervals file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected File outputFile;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    private List<SimpleInterval> intervals;
    private ReferenceDataSource reference;
    private final GCContentAnnotator gcContentAnnotator = new GCContentAnnotator();
    private AnnotatedIntervalCollection annotatedIntervals;

    @Override
    public void onTraversalStart() {
        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);

        logger.info("Loading intervals for annotation...");
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
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
        annotatedIntervals = new AnnotatedIntervalCollection(annotatedIntervalList);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info(String.format("Writing annotated intervals to %s...", outputFile));
        annotatedIntervals.write(outputFile);
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
