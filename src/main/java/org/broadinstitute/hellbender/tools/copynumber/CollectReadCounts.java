package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multiset;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.HDF5SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Collects read counts at specified intervals.  The count for each interval is calculated by counting
 * the number of read starts that lie in the interval.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SAM format read data
 *     </li>
 *     <li>
 *         Intervals at which counts will be collected.
 *         The argument {@code interval-merging-rule} must be set to {@link IntervalMergingRule#OVERLAPPING_ONLY}
 *         and all other common arguments for interval padding or merging must be set to their defaults.
 *     </li>
 *     <li>
 *         Output file format.  This can be used to select TSV or HDF5 output.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Counts file.
 *         By default, the tool produces HDF5 format results. This can be changed with the {@code format} option
 *         to TSV format. Using HDF5 files with {@link CreateReadCountPanelOfNormals}
 *         can decrease runtime, by reducing time spent on IO, so this is the default output format.
 *         The HDF5 format contains information in the paths defined in {@link HDF5SimpleCountCollection}. HDF5 files may be viewed using
 *         <a href="https://support.hdfgroup.org/products/java/hdfview/">hdfview</a> or loaded in python using
 *         <a href="http://www.pytables.org/">PyTables</a> or <a href="http://www.h5py.org/">h5py</a>.
 *         The TSV format has a SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers contained in
 *         {@link SimpleCountCollection.SimpleCountTableColumn}, and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk CollectReadCounts \
 *          -I sample.bam \
 *          -L intervals.interval_list \
 *          --interval-merging-rule OVERLAPPING_ONLY \
 *          -O sample.counts.hdf5
 * </pre>
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collects read counts at specified intervals",
        oneLineSummary = "Collects read counts at specified intervals",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
public final class CollectReadCounts extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(CollectReadCounts.class);

    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;

    enum Format {
        TSV, HDF5
    }

    public static final String FORMAT_LONG_NAME = "format";

    @Argument(
            doc = "Output file for read counts.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputCountsFile = null;

    @Argument(
            doc = "Output file format.",
            fullName = FORMAT_LONG_NAME,
            optional = true
    )
    private Format format = Format.HDF5;

    /**
     * Metadata contained in the BAM file.
     */
    private SampleLocatableMetadata metadata;

    private List<SimpleInterval> intervals;

    private String currentContig = null;

    /**
     * Overlap detector used to determine when read starts overlap with input intervals.
     */
    private CachedOverlapDetector<SimpleInterval> intervalCachedOverlapDetector;

    private Multiset<SimpleInterval> intervalMultiset;

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.MAPPED);
        readFilters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        return readFilters;
    }

    @Override
    public void onTraversalStart() {
        metadata = MetadataUtils.fromHeader(getHeaderForReads(), Metadata.Type.SAMPLE_LOCATABLE);
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        //this check is currently redundant, since the master dictionary is taken from the reads;
        //however, if any other dictionary is added in the future, such a check should be performed
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(metadata.getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in BAM does not match the master sequence dictionary.");
        }

        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);

        intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
        intervalMultiset = HashMultiset.create(intervals.size());

        logger.info("Collecting read counts...");
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (currentContig == null || !read.getContig().equals(currentContig)) {
            //if we are on a new contig, create an OverlapDetector covering the contig
            currentContig = read.getContig();
            final List<SimpleInterval> intervalsOnCurrentContig = intervals.stream()
                    .filter(i -> i.getContig().equals(currentContig))
                    .collect(Collectors.toList());
            intervalCachedOverlapDetector = new CachedOverlapDetector<>(intervalsOnCurrentContig);
        }
        final SimpleInterval overlappingInterval = intervalCachedOverlapDetector.getOverlap(
                new SimpleInterval(read.getContig(), read.getStart(), read.getStart()));

        //if read doesn't overlap any of the provided intervals, do nothing
        if (overlappingInterval == null) {
            return;
        }
        intervalMultiset.add(overlappingInterval);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Writing read counts to " + outputCountsFile);
        final SimpleCountCollection readCounts = new SimpleCountCollection(
                metadata,
                ImmutableList.copyOf(intervals.stream()     //making this an ImmutableList avoids a defensive copy in SimpleCountCollection
                        .map(i -> new SimpleCount(i, intervalMultiset.count(i)))
                        .iterator()));

        if (format == Format.HDF5) {
            readCounts.writeHDF5(outputCountsFile);
        } else {
            readCounts.write(outputCountsFile);
        }

        return "SUCCESS";
    }

    /**
     * A simple wrapper around {@link OverlapDetector} to provide naive caching and ensure that overlap sets
     * only contain a single interval.
     */
    private final class CachedOverlapDetector<T extends Locatable> {
        private final OverlapDetector<T> overlapDetector;
        private T cachedResult;

        CachedOverlapDetector(final List<T> intervals) {
            Utils.nonEmpty(intervals);
            this.overlapDetector = OverlapDetector.create(intervals);
            //double check that intervals do not overlap
            Utils.validateArg(intervals.stream().noneMatch(i -> overlapDetector.getOverlaps(i).size() > 1),
                    "Input intervals may not be overlapping.");
            cachedResult = intervals.get(0);
        }

        /**
         * We check the previously cached result first.  Assuming that queries will be made in sorted order,
         * this may slightly save on lookup time.
         * @return {@code null} if no interval overlaps {@code locatable}
         */
        T getOverlap(final Locatable locatable) {
            if (IntervalUtils.overlaps(cachedResult, locatable)) {
                return cachedResult;
            }
            final Set<T> overlaps = overlapDetector.getOverlaps(locatable);
            if (overlaps.size() > 1) {
                //should not reach here since intervals are checked for overlaps;
                //performing a redundant check to protect against future code changes
                throw new GATKException.ShouldNeverReachHereException("Intervals should be non-overlapping, " +
                        "so at most one interval should intersect with the start of a read.");
            }
            final T firstOverlap = overlaps.stream().findFirst().orElse(null);
            if (firstOverlap != null) {
                cachedResult = firstOverlap;
            }
            return firstOverlap;
        }
    }
}