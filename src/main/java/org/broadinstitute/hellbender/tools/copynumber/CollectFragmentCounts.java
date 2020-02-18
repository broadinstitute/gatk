package org.broadinstitute.hellbender.tools.copynumber;


import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.spark.sql.catalyst.expressions.Sin;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.HDF5SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;

import org.broadinstitute.hellbender.utils.IntervalMergingRule;

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

    final private LinkedHashSet<FeatureCoverage> featureCounts = new LinkedHashSet<>();

    final private OverlapDetector<Pair<FeatureCoverage, Gff3BaseData>> featureOverlapDetector = new OverlapDetector<>(0,0);

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.MAPPED);
        readFilters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        readFilters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        readFilters.add(new InwardFragmentsFilter());
        readFilters.add(ReadFilterLibrary.PRIMARY_LINE);
        //readFilters.add(new GoodPairReadFilter());
        return readFilters;
    }

    @Override
    public void onTraversalStart() {
        validateArguments();
        //this check is currently redundant, since the master dictionary is taken from the reads;
        //however, if any other dictionary is added in the future, such a check should be performed


        logger.info("Collecting read counts...");

        final SAMSequenceDictionary dict = getMasterSequenceDictionary();
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
                final Gff3BaseData overlapBaseData = feature.getBaseData();
                if (grouping_type.contains(feature.getType())) {
                    addGroupingFeature(overlapBaseData, overlapBaseData);
                }

                feature.getAncestors().stream().filter(f -> grouping_type.contains(f.getType())).map(Gff3Feature::getBaseData).forEach(b -> addGroupingFeature(b, overlapBaseData));
            }
        }
    }

    private void addGroupingFeature(final Gff3BaseData groupingBaseData, final Gff3BaseData overlappingBaseData) {
        final String geneID = groupingBaseData.getAttributes().get(gene_id_key);
        if (geneID == null) {
            throw new UserException("no geneid field " + gene_id_key + " found in feature at " + groupingBaseData.getContig() + ":" + groupingBaseData.getStart() + "-" + groupingBaseData.getEnd());
        }

        final FeatureCoverage featureCoverage = new FeatureCoverage(groupingBaseData, groupingBaseData.getStrand() != Strand.NONE);
        featureOverlapDetector.addLhs(Pair.of(featureCoverage, overlappingBaseData), overlappingBaseData);
        featureCounts.add(featureCoverage);
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

    private List<SimpleInterval> getAlignmentIntervals(final GATKRead read) {

        if (spliced) {
            final IntervalList alignmentIntervals = new IntervalList(getMasterSequenceDictionary());
            final SAMRecord rec = read.convertToSAMRecord(getHeaderForReads());
            if(SAMUtils.getMateCigar(rec) == null) {
                throw new GATKException("Mate cigar must be present if using spliced reads");
            }
            final List<AlignmentBlock> readAlignmentBlocks = rec.getAlignmentBlocks();

            for( final AlignmentBlock block : readAlignmentBlocks) {
                alignmentIntervals.add(new Interval(read.getContig(), block.getReferenceStart(), block.getReferenceStart()+block.getLength()));
            }

            boolean overlapsMate = false;
            if (inGoodPair(read)) {
                final List<AlignmentBlock> mateAlignmentBlocks = SAMUtils.getMateAlignmentBlocks(rec);
                for( final AlignmentBlock block : mateAlignmentBlocks) {
                    final Interval alignmentBlockInterval = new Interval(read.getMateContig(), block.getReferenceStart(), block.getReferenceStart()+block.getLength());
                    alignmentIntervals.add(alignmentBlockInterval);

                    if (!overlapsMate && read.overlaps(alignmentBlockInterval)) {
                        overlapsMate = true;
                    }
                }
            }
            if (overlapsMate) {
                alignmentIntervals.unique();
            }
            return alignmentIntervals.getIntervals().stream().map(i -> new SimpleInterval(i.getContig(), i.getStart(), i.getEnd())).collect(Collectors.toList());
        } else {
            return Arrays.asList(getReadInterval(read));
        }


    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {

        if ((!spliced || !read.isReverseStrand() || !inGoodPair(read))) {
            final List<SimpleInterval> alignmentIntervals = getAlignmentIntervals(read);

            //List<Gff3Feature> features = featureContext.getValues(gffFile);
            Set<Pair<FeatureCoverage, Gff3BaseData>> features = featureOverlapDetector.getOverlaps(getReadInterval(read));

            final Strand fragmentStrand = getFragmentStrand(read);

            final int basesOnReference = alignmentIntervals.stream().map(SimpleInterval::getLengthOnReference).reduce(0, Integer::sum);
            for (final Pair<FeatureCoverage, Gff3BaseData> featureCoverageGff3BaseDataPair : features) {
                final FeatureCoverage featureCoverage = featureCoverageGff3BaseDataPair.getLeft();
                final Gff3BaseData overlapBaseData = featureCoverageGff3BaseDataPair.getRight();
                final boolean isSense = !featureCoverage.isStranded || overlapBaseData.getStrand() == fragmentStrand;

                for (final SimpleInterval interval : alignmentIntervals) {
                    if(interval.overlaps(overlapBaseData)) {
                        if (isSense) {
                            featureCoverage.addSenseCount((float)interval.intersect(overlapBaseData).size() / (float) basesOnReference);
                        } else {
                            featureCoverage.addAntiSenseCount((float) interval.intersect(overlapBaseData).size() / (float) basesOnReference);
                        }
                    }
                }
            }
        }
    }

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
            for (final FeatureCoverage featureCoverage : featureCounts) {
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
                    .set(name+"_counts", fragmentCount.count)
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


