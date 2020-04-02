package org.broadinstitute.hellbender.tools.copynumber;


import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import htsjdk.tribble.readers.LineIterator;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.spark.sql.catalyst.expressions.In;
import org.apache.spark.sql.catalyst.expressions.Sin;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.HDF5SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;

import org.broadinstitute.hellbender.utils.IntervalMergingRule;

import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import scala.Immutable;
import scala.Int;


import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.Arrays;
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
 *         {@link SimpleCountCollection}, and the corresponding entry rows.
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
        summary = "Collects fragment counts at specified intervals",
        oneLineSummary = "Collects fragment counts at specified intervals",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
public final class CollectFragmentCounts extends ReadWalker {
    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 0;


    @Argument(
            doc = "Output file for read counts.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputCountsFile = null;

    @Argument(doc="gff file", shortName = "G")
    private FeatureInput<Gff3Feature> gffFile;

    @Argument(doc="types to group by", shortName = "T")
    private Set<String> grouping_type = new HashSet<>(Collections.singleton("gene"));

    @Argument(doc="overlap type")
    private Set<String> overlap_type = new HashSet<>(Collections.singleton("exon"));

    @Argument(doc="sample name to label counts", shortName = "N")
    private String name;

    @Argument(doc="gene_id key")
    private String gene_id_key = "gene_id";

    @Argument(doc = "which read corresponds to the transcription strand")
    private TrancriptionRead trancriptionRead = TrancriptionRead.R1;

    @Argument(doc = "Whether the rna is spliced.  If spliced, alignments must be from a splice aware aligner (such as star), if unspliced, alignments must be from " +
            "a non-splicing aligner (such as bwa). ")
    private boolean spliced = true;

    final private LinkedHashMap<Gff3BaseData, FeatureCoverage> featureCounts = new LinkedHashMap<>();

    final private OverlapDetector<Pair<FeatureCoverage, Gff3BaseData>> featureOverlapDetector = new OverlapDetector<>(0,0);

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>();
        readFilters.add(ReadFilterLibrary.VALID_ALIGNMENT_START);
        readFilters.add(ReadFilterLibrary.VALID_ALIGNMENT_END);
        readFilters.add(new AlignmentAgreesWithHeaderReadFilter(getHeaderForReads()));
        readFilters.add(ReadFilterLibrary.HAS_READ_GROUP);
        readFilters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        readFilters.add(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH);
        readFilters.add(ReadFilterLibrary.SEQ_IS_STORED);
        readFilters.add(ReadFilterLibrary.MAPPED);
        readFilters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        readFilters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        //readFilters.add(new InwardFragmentsFilter());
        readFilters.add(ReadFilterLibrary.PRIMARY_LINE);
        readFilters.add(new NotMultiMappingReadFilter());
        //readFilters.add(new GoodPairReadFilter());
        return readFilters;
    }

    @Override
    public void onTraversalStart() {
        validateArguments();
        //this check is currently redundant, since the master dictionary is taken from the reads;
        //however, if any other dictionary is added in the future, such a check should be performed
        getHeaderForReads().getSequenceDictionary();
        final SAMSequenceDictionary dict = getBestAvailableSequenceDictionary();
        if (dict == null) {
            throw new GATKException("sequence dictionary must be specified (" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + ").");
        }

        logger.info("collecting list of features");
        for (final SAMSequenceRecord contig : dict.getSequences()) {
            final List<Gff3Feature> contigFeatures = features.getFeatures(gffFile, new SimpleInterval(contig.getSequenceName(), 1, contig.getSequenceLength()));
            logger.info("collecting features in " + contig.getSequenceName());
            for (final Gff3Feature feature : contigFeatures) {
                if(!overlap_type.contains(feature.getType())) {
                    continue;
                }
                final Gff3BaseData overlapBaseData = shrinkBaseData(feature.getBaseData());
                if (grouping_type.contains(feature.getType())) {
                    addGroupingFeature(overlapBaseData, overlapBaseData);
                }

                feature.getAncestors().stream().filter(f -> grouping_type.contains(f.getType())).map(Gff3Feature::getBaseData).map(this::shrinkBaseData).forEach(b -> addGroupingFeature(b, overlapBaseData));
            }
        }

        logger.info("Collecting read counts...");
    }

    private Gff3BaseData shrinkBaseData(final Gff3BaseData baseData) {
        //remove all but gene_id_key attributes
        final Map<String, String> shrunkAttributes = baseData.getAttributes().entrySet().stream().filter(e -> e.getKey().equals(gene_id_key)).collect(Collectors.toMap(Map.Entry::getKey,Map.Entry::getValue));
        final Gff3BaseData shrunkBaseData = new Gff3BaseData(baseData.getContig(), baseData.getSource(), baseData.getType(), baseData.getStart(), baseData.getEnd(), baseData.getStrand(), baseData.getPhase(), shrunkAttributes);
        return shrunkBaseData;
    }

    private void addGroupingFeature(final Gff3BaseData groupingBaseData, final Gff3BaseData overlappingBaseData) {
        final String geneID = groupingBaseData.getAttributes().get(gene_id_key);
        if (geneID == null) {
            throw new UserException("no geneid field " + gene_id_key + " found in feature at " + groupingBaseData.getContig() + ":" + groupingBaseData.getStart() + "-" + groupingBaseData.getEnd());
        }

        final FeatureCoverage featureCoverage = featureCounts.computeIfAbsent(groupingBaseData,
                g -> new FeatureCoverage(g, g.getStrand() != Strand.NONE));
        featureOverlapDetector.addLhs(Pair.of(featureCoverage, overlappingBaseData), overlappingBaseData);


    }

    private void validateArguments() {
        CopyNumberArgumentValidationUtils.validateOutputFiles(outputCountsFile);
    }

    @Override
    public SimpleInterval getReadInterval(final GATKRead read) {
        if (spliced) {
            return super.getReadInterval(read);
        }
        if (read.isUnmapped()) {
            return null;
        }

        if (!inGoodPair(read)) {
            return SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd())? new SimpleInterval(read.getContig(), read.getStart(), read.getEnd()) : null;
            //return SimpleInterval.isValid(read.getContig(), read.getStart(), read.getStart() + 76 - 1)? new SimpleInterval(read.getContig(), read.getStart(), read.getStart() + 76 - 1) : null;
        }

        final int start = Math.min(read.getStart(), read.getMateStart());
        final int end = start + Math.abs(read.getFragmentLength()) - 1;
        //final int end = Math.max(read.getStart() + 76 - 1, read.getMateStart() + 76 - 1);
        return SimpleInterval.isValid(read.getContig(), start, end)? new SimpleInterval(read.getContig(), start, end) : null;
    }

    private static boolean inGoodPair(final GATKRead read) {
        return !(read.mateIsUnmapped() || !read.isProperlyPaired() ||
                !read.getContig().equals(read.getMateContig()) ||
                read.isReverseStrand() == read.mateIsReverseStrand() ||
                (read.isReverseStrand() && read.getEnd() < read.getMateStart()) ||
                (!read.isReverseStrand() && read.getStart() > read.getMateStart() + TextCigarCodec.decode(read.getAttributeAsString("MC")).getReferenceLength()) ||
                read.getAttributeAsInteger("MQ") < DEFAULT_MINIMUM_MAPPING_QUALITY
        );
    }

    private List<Interval> getAlignmentIntervals(final GATKRead read) {

        if (spliced) {
            final SAMFileHeader intervalHeader = new SAMFileHeader(getMasterSequenceDictionary());
            //need to set sortOrder to coordinate in order to avoid inefficient IntervalList.getUniqueIntervals recreating header later
            intervalHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            //final IntervalList alignmentIntervals = new IntervalList(intervalHeader);
            final List<Interval> alignmentIntervals = new ArrayList<>();
            final SAMRecord rec = read.convertToSAMRecord(getHeaderForReads());

            final List<AlignmentBlock> readAlignmentBlocks = rec.getAlignmentBlocks();

            for( final AlignmentBlock block : readAlignmentBlocks) {
                alignmentIntervals.add(new Interval(read.getContig(), block.getReferenceStart(), block.getReferenceStart()+block.getLength() - 1));
            }

            boolean overlapsMate = false;
            if (inGoodPair(read)) {
                if(SAMUtils.getMateCigar(rec) == null) {
                    throw new GATKException("Mate cigar must be present if using spliced reads");
                }
                final List<AlignmentBlock> mateAlignmentBlocks = SAMUtils.getMateAlignmentBlocks(rec);
                for( final AlignmentBlock block : mateAlignmentBlocks) {
                    final Interval alignmentBlockInterval = new Interval(read.getMateContig(), block.getReferenceStart(), block.getReferenceStart()+block.getLength() - 1);
                    alignmentIntervals.add(alignmentBlockInterval);

                    if (!overlapsMate && read.overlaps(alignmentBlockInterval)) {
                        overlapsMate = true;
                    }
                }
            }
            if (overlapsMate) {
                //sort and merge
                Collections.sort(alignmentIntervals);
                final IntervalList.IntervalMergerIterator mergeringIterator = new IntervalList.IntervalMergerIterator(alignmentIntervals.iterator(), true, true, false);
                final List<Interval> merged = new ArrayList<>();
                while (mergeringIterator.hasNext()) {
                    merged.add(mergeringIterator.next());
                }
                return merged;
            } else {
                Collections.sort(alignmentIntervals);
                return alignmentIntervals;
            }
        } else {
            final SimpleInterval fragmentInterval = getReadInterval(read);
            return Arrays.asList(new Interval(fragmentInterval.getContig(), fragmentInterval.getStart(), fragmentInterval.getEnd()));
        }


    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if ((!spliced || !read.isReverseStrand() || !inGoodPair(read))) {
            final List<Interval> alignmentIntervals = getAlignmentIntervals(read);

            //List<Gff3Feature> features = featureContext.getValues(gffFile);
            Set<Pair<FeatureCoverage, Gff3BaseData>> features ;
            if (spliced) {
                features = alignmentIntervals.stream().flatMap(i -> featureOverlapDetector.getOverlaps(i).stream()).collect(Collectors.toSet());
            } else {
                features = featureOverlapDetector.getOverlaps(getReadInterval(read));
            }

            final Strand fragmentStrand = getFragmentStrand(read);



            final int basesOnReference = alignmentIntervals.stream().map(Interval::getLengthOnReference).reduce(0, Integer::sum);
            final Map<FeatureCoverage, Set<Gff3BaseData>> overlapsByGroupingFeatures = new LinkedHashMap<>();
            //count how many bases in alignmentIntervals are not covered by any feature
            final List<Interval> overlappingIntervals = features.stream().map(p -> new Interval (p.getRight().getContig(), p.getRight().getStart(), p.getRight().getEnd())).collect(Collectors.toList());

           // final int uncoveredBases = countUncoveredBases(alignmentIntervals, overlappingIntervals);


            for (final Pair<FeatureCoverage, Gff3BaseData> feature : features) {
                final Set<Gff3BaseData> overlappingFeatures = overlapsByGroupingFeatures.computeIfAbsent(feature.getLeft(), b -> new HashSet<>());
                overlappingFeatures.add(feature.getRight());
            }
//            final Map<FeatureCoverage, Set<Gff3BaseData>> overlapsByGroupingFeatures = features.stream().collect(Collectors.groupingBy(
//                    Pair::getLeft, Collectors.mapping(
//                            Pair::getRight,Collectors.toSet())
//                    )
//            );

            final int nGroupingFeaturesCovered = overlapsByGroupingFeatures.size();
            for (final Map.Entry<FeatureCoverage, Set<Gff3BaseData>> overlapsByGroupingFeature : overlapsByGroupingFeatures.entrySet()) {
                final FeatureCoverage featureCoverage = overlapsByGroupingFeature.getKey();
                if(featureCoverage.baseData.getId().equals("ENSG00000221344")) {
                    System.out.println("hi");
                }
                final boolean isSense = !featureCoverage.isStranded || featureCoverage.baseData.getStrand() == fragmentStrand;
                float maxCounts = (float)1.0/(float)nGroupingFeaturesCovered;
//                float maxCounts = 0;
//                for (final Gff3BaseData overlapBaseData : overlapsByGroupingFeature.getValue()) {
//                    if (!overlapBaseData.getStrand().equals(featureCoverage.baseData.getStrand())) {
//                        throw new GATKException("Grouping feature " + featureCoverage.baseData.getId() + " and subfeature " + overlapBaseData.getId() + " are on different strands (" + overlapBaseData.getStrand() + "," + featureCoverage.baseData.getStrand() + ")");
//                    }
//                    float thisCounts = 0;
//                    for (final SimpleInterval interval : alignmentIntervals) {
//                        if (interval.overlaps(overlapBaseData)) {
//                            thisCounts += (float) interval.intersect(overlapBaseData).size() / (float) basesOnReference;
//                        }
//                    }
//                    if (thisCounts>maxCounts) {
//                        maxCounts = thisCounts;
//                    }
//                }
                if (isSense) {
                    featureCoverage.addSenseCount(maxCounts);
                } else {
                    featureCoverage.addAntiSenseCount(maxCounts);
                }
            }
        }
    }

//    private int countUncoveredBases(final List<Interval> countingIntervals, final List<Interval> coveringIntervals) {
//        if (coveringIntervals.isEmpty()) {
//            return (int)Interval.countBases(countingIntervals);
//        }
//        //make sure sorted
//        Collections.sort(countingIntervals);
//        Collections.sort(coveringIntervals);
//
//        Iterator<Interval> countingIntervalsIterator = countingIntervals.iterator();
//        Iterator<Interval> coveringIntervalsIterator = coveringIntervals.iterator();
//
//        int uncoveredBases = 0;
//        Interval currentCoveringInterval = coveringIntervalsIterator.next();
//        while (countingIntervalsIterator.hasNext()) {
//            final Interval currentCountingInterval = countingIntervalsIterator.next();
//            int currentIntervalCoveredUntil = currentCountingInterval.getStart(); //we will move through the currentCountingInterval
//            //get to correct contig, or last covering interval
//            while (!currentCoveringInterval.getContig().equals(currentCountingInterval.getContig()) && coveringIntervalsIterator.hasNext()) {
//                 currentCoveringInterval = coveringIntervalsIterator.next();
//            }
//
//            //seek until end of currentCoveringInterval is at or past currentCountingIntervalPos
//            while (currentCoveringInterval.getContig().equals(currentCountingInterval.getContig()) && currentCoveringInterval.getEnd() < currentIntervalCoveredUntil && coveringIntervalsIterator.hasNext()) {
//                currentCoveringInterval = coveringIntervalsIterator.next();
//            }
//
//            //if currentCoveringInterval start is past end of currentCountingInterval, continue
//            if (!currentCoveringInterval.getContig().equals(currentCountingInterval.getContig()) || currentCoveringInterval.getStart() > currentCountingInterval.getEnd()) {
//                continue;
//            }
//
//            if (currentCoveringInterval.getStart() > currentIntervalCoveredUntil) {
//                uncoveredBases += currentCoveringInterval.getStart() - currentIntervalCoveredUntil;
//            }
//
//            if (currentIntervalCoveredUntil <= currentCoveringInterval.getEnd()) {
//                currentIntervalCoveredUntil = currentCoveringInterval.getEnd() + 1;
//            }
//
//            if (currentIntervalCoveredUntil > )
//
//
//        }
//
//    }

    @Override
    public Object onTraversalSuccess() {
        logger.info(String.format("Writing read counts to %s...", outputCountsFile.getAbsolutePath()));
        try (final FragmentCountWriter writer = new FragmentCountWriter(outputCountsFile.toPath())) {
            int i=0;
            for (final File input_bam: readArguments.getReadFiles()) {
                writer.writeMetadata("input_bam_"+i, input_bam.toString());
                i++;
            }
            writer.writeMetadata("annotation_file", gffFile.toString());
            for (final FeatureCoverage featureCoverage : featureCounts.values()) {
                for (final SingleStrandFeatureCoverage singleStrandFeatureCoverage : featureCoverage.getSingleStrandCoverages()) {
                    writer.writeRecord(singleStrandFeatureCoverage);
                }
            }
        } catch (IOException ex) {
            throw new UserException(ex.getMessage());
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    public class FragmentCountWriter extends TableWriter<SingleStrandFeatureCoverage> {

        public FragmentCountWriter(final Path file ) throws IOException {
            super(file, new TableColumnCollection("gene_id", "contig", "start", "stop", "strand", "sense_antisense", name+"_counts"));
        }

        protected void composeLine(final SingleStrandFeatureCoverage fragmentCount, final DataLine dataLine) {
            dataLine.set("contig", fragmentCount.baseData.getContig())
                    .set("start", fragmentCount.baseData.getStart())
                    .set("stop", fragmentCount.baseData.getEnd())
                    .set("strand", fragmentCount.baseData.getStrand().encode())
                    .set("sense_antisense", fragmentCount.sense? "sense" : "antisense")
                    .set(name+"_counts", fragmentCount.count, 2)
                    .set("gene_id", fragmentCount.baseData.getAttributes().get(gene_id_key));
        }
    }

    public class FeatureCoverage {
        private final Gff3BaseData baseData;
        private float sense_count;
        private float antisense_count;
        private final boolean isStranded;

        FeatureCoverage(final Gff3BaseData baseData, final boolean isStranded) {
            this.baseData = baseData;
            this.isStranded = isStranded;
        }

        public void addSenseCount(float count) {
            sense_count += count;
        }

        public void addAntiSenseCount(float count) {
            if (!isStranded) {
                throw new GATKException("Should never add an antisense count to a non-stranded feature.");
            }
            antisense_count += count;
        }

        public List<SingleStrandFeatureCoverage> getSingleStrandCoverages() {
            final List<SingleStrandFeatureCoverage> coverages = new ArrayList<>(Collections.singletonList(new SingleStrandFeatureCoverage(baseData, sense_count, true)));
            if (isStranded) {
                coverages.add(new SingleStrandFeatureCoverage(baseData, antisense_count, false));
            }

            return coverages;
        }

        @Override
        public int hashCode() {
            int hash = baseData.hashCode();
            hash = 31 * hash + (int)sense_count;
            hash = 31 * hash + (int)antisense_count;
            hash = 31 * hash + (isStranded? 1 : 0);
            return hash;
        }

        @Override
        public boolean equals(final Object other) {
            if (! (other instanceof FeatureCoverage)) {
                return false;
            }

            return baseData.equals(((FeatureCoverage) other).baseData) && sense_count == ((FeatureCoverage) other).sense_count &&
                    antisense_count == ((FeatureCoverage) other).antisense_count && isStranded == ((FeatureCoverage) other).isStranded;
        }
    }

    public class SingleStrandFeatureCoverage {
        public final Gff3BaseData baseData;
        //public final double count;
        public final float count;
        public final boolean sense;

        SingleStrandFeatureCoverage(final Gff3BaseData baseData, final float count, final boolean sense) {
            this.baseData = baseData;
            this.count = count;
            this.sense = sense;
        }
    }

    public class NotMultiMappingReadFilter extends ReadFilter {
        @Override
        public boolean test(GATKRead read) {
            return read.hasAttribute(SAMTag.NH.toString()) && read.getAttributeAsInteger(SAMTag.NH.toString()) == 1;
        }
    }

    public class GoodPairReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test(GATKRead read) {
            return CollectFragmentCounts.inGoodPair(read);
        }
    }

    public class InwardFragmentsFilter extends ReadFilterLibrary.MateDifferentStrandReadFilter {
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test(GATKRead read) {
            final boolean mateDifferentStrandFilterResult = super.test(read);
            return mateDifferentStrandFilterResult && ((read.isReverseStrand() && read.getEnd() >= read.getMateStart()) ||
                    (!read.isReverseStrand() && read.getStart() <= read.getMateStart() + TextCigarCodec.decode(read.getAttributeAsString("MC")).getReferenceLength()));
        }
    }

    private Strand getFragmentStrand(final GATKRead read) {
        return (read.isFirstOfPair() == (trancriptionRead==TrancriptionRead.R1))? (read.isReverseStrand() ? Strand.NEGATIVE : Strand.POSITIVE) :
                                                                                    (read.isReverseStrand() ? Strand.POSITIVE : Strand.NEGATIVE);
    }

    private enum TrancriptionRead {
        R1,
        R2
    }
}


