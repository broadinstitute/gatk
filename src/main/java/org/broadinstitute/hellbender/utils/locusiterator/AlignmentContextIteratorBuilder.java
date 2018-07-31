package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.iterators.IntervalOverlappingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Create an iterator for traversing alignment contexts in a specified manner.  This class should be able to support both spark and non-spark LocusWalker
 */
public class AlignmentContextIteratorBuilder {

    protected static final Logger logger = LogManager.getLogger(AlignmentContextIteratorBuilder.class);

    private boolean isEmitEmptyLoci;
    private boolean isKeepUniqueReadListInLibs;
    private boolean isIncludeDeletions;
    private boolean isIncludeNs;
    private LIBSDownsamplingInfo downsamplingInfo;

    public boolean isEmitEmptyLoci() {
        return isEmitEmptyLoci;
    }

    public void setEmitEmptyLoci(boolean emitEmptyLoci) {
        isEmitEmptyLoci = emitEmptyLoci;
    }

    public boolean isKeepUniqueReadListInLibs() {
        return isKeepUniqueReadListInLibs;
    }

    public void setKeepUniqueReadListInLibs(boolean keepUniqueReadListInLibs) {
        isKeepUniqueReadListInLibs = keepUniqueReadListInLibs;
    }

    public boolean isIncludeDeletions() {
        return isIncludeDeletions;
    }

    public void setIncludeDeletions(boolean includeDeletions) {
        isIncludeDeletions = includeDeletions;
    }

    public boolean isIncludeNs() {
        return isIncludeNs;
    }

    public void setIncludeNs(boolean includeNs) {
        isIncludeNs = includeNs;
    }

    public LIBSDownsamplingInfo getDownsamplingInfo() {
        return downsamplingInfo;
    }

    public void setDownsamplingInfo(LIBSDownsamplingInfo downsamplingInfo) {
        this.downsamplingInfo = downsamplingInfo;
    }

    public AlignmentContextIteratorBuilder() {
        isEmitEmptyLoci = false;
        isKeepUniqueReadListInLibs = false;
        isIncludeDeletions = true;
        isIncludeNs = false;
        downsamplingInfo = LocusIteratorByState.NO_DOWNSAMPLING;
    }

    /**
     *  Have this builder return an instance of the iterator as specified by the caller.
     *
     * @param readIterator iterator of sorted GATK reads.  Not {@code null}
     * @param header SAM file header to use.  Not {@code null}
     * @param intervalsForTraversal the intervals to generate alignment contexts over.
     * @param dictionary the SAMSequenceDictionary being used for this traversal.  This can be the same as the reference.  {@code null} is supported, but will often lead to invalid parameter combinations.
     * @param isReference {@code true} if the specified dictionary came from a reference.  {@code false} otherwise.  If dictionary is {@code null}, this parameter is ignored.
     * @return iterator that produces AlignmentContexts ready for consumption (e.g. by a {@link org.broadinstitute.hellbender.engine.LocusWalker})
     */
    public Iterator<AlignmentContext> build(final Iterator<GATKRead> readIterator, final SAMFileHeader header, final List<SimpleInterval> intervalsForTraversal, final SAMSequenceDictionary dictionary, final boolean isReference) {
        Utils.nonNull(header, "Header cannot be null");
        Utils.nonNull(readIterator, "Read iterator cannot be null");
        final boolean isDefinitelyReference = (dictionary != null) && isReference ;
        return createAlignmentContextIterator(intervalsForTraversal, header, readIterator, dictionary, downsamplingInfo,
                isDefinitelyReference, isEmitEmptyLoci, isKeepUniqueReadListInLibs, isIncludeDeletions, isIncludeNs);
    }

    /**
     *  Create the appropriate instance of an alignment context spliterator based on the input parameters.
     *
     *  Please note that this wrapper is still tied to {@link LocusIteratorByState} and some parameters are being passed directly to that class.
     *
     * @param intervalsForTraversal the intervals to generate alignment contexts over.
     * @param header SAM file header to use
     * @param readIterator iterator of sorted GATK reads
     * @param dictionary the SAMSequenceDictionary being used for this traversal.  This can be the same as the reference.  {@code null} is supported, but will often lead to invalid parameter combinations.
     * @param downsamplingInfo how to downsample (for {@link LocusIteratorByState})
     * @param isReference the dictionary specified above is a reference, {@code false} if no reference being used or it is not a reference.
     * @param emitEmptyLoci whether loci with no coverage should be emitted.  In this case, the AlignmentContext will be empty (not null).
     * @param isKeepUniqueReadListInLibs if true, we will keep the unique reads from the samIterator and make them
     *                                       available via the transferReadsFromAllPreviousPileups interface (this parameter is specific to {@link LocusIteratorByState})
     * @param isIncludeDeletions include reads with deletion on the loci in question
     * @param isIncludeNs include reads with N on the loci in question
     * @return iterator that produces AlignmentContexts ready for consumption (e.g. by a {@link org.broadinstitute.hellbender.engine.LocusWalker})
     */
    private static Iterator<AlignmentContext> createAlignmentContextIterator(final List<SimpleInterval> intervalsForTraversal,
                                                                            final SAMFileHeader header,
                                                                               final Iterator<GATKRead> readIterator,
                                                                               final SAMSequenceDictionary dictionary,
                                                                               final LIBSDownsamplingInfo downsamplingInfo,
                                                                               final boolean isReference,
                                                                               boolean emitEmptyLoci,
                                                                               boolean isKeepUniqueReadListInLibs,
                                                                               boolean isIncludeDeletions,
                                                                               boolean isIncludeNs) {

        // get the samples from the read groups
        final Set<String> samples = header.getReadGroups().stream()
                .map(SAMReadGroupRecord::getSample)
                .collect(Collectors.toSet());

        // get the LIBS
        final LocusIteratorByState libs = new LocusIteratorByState(readIterator, downsamplingInfo, isKeepUniqueReadListInLibs, samples, header, isIncludeDeletions, isIncludeNs);

        List<SimpleInterval> finalIntervals = intervalsForTraversal;
        validateEmitEmptyLociParameters(emitEmptyLoci, dictionary, intervalsForTraversal, isReference);
        if (emitEmptyLoci) {

            // If no intervals were specified, then use the entire reference (or best available sequence dictionary).
            if (!areIntervalsSpecified(finalIntervals)) {
                finalIntervals = IntervalUtils.getAllIntervalsForReference(dictionary);
            }
            final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(finalIntervals.iterator());
            return new IntervalAlignmentContextIterator(libs, intervalLocusIterator, header.getSequenceDictionary());
        } else if (areIntervalsSpecified(finalIntervals)) {
            return new IntervalOverlappingIterator<>(libs, finalIntervals, header.getSequenceDictionary());
        } else {
            // prepare the iterator
            return libs;
        }
    }

    private static boolean areIntervalsSpecified(final List<SimpleInterval> finalIntervals) {
        return finalIntervals != null;
    }

    /**
     *  The emit empty loci parameter comes with several pitfalls when used incorrectly.  Here we check and either give
     *   warnings or errors.
     */
    static private void validateEmitEmptyLociParameters(boolean emitEmptyLoci, final SAMSequenceDictionary dictionary, final List<SimpleInterval> intervals, final boolean isReference) {
        if (emitEmptyLoci) {
            if ((dictionary == null) && !isReference) {
                throw new UserException.MissingReference("No sequence dictionary nor reference specified.  Therefore, emitting empty loci is impossible and this tool cannot be run.  The easiest fix here is to specify a reference dictionary.");
            }
            if (!isReference && !areIntervalsSpecified(intervals)) {
                logger.warn("****************************************");
                logger.warn("* Running this tool without a reference nor intervals can yield unexpected results, since it will emit results for loci with no reads.  A sequence dictionary has been found and the intervals will be derived from this.  The easiest way avoid this message is to specify a reference.");
                logger.warn("****************************************");
            }
        }
    }
}
