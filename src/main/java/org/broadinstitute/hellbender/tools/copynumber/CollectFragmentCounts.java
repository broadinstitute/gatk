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
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import org.apache.commons.lang3.tuple.MutablePair;
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
    //private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 10;
    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 0;


    @Argument(
            doc = "Output file for read counts.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputCountsFile = null;

    @Argument(doc="gff file", shortName = "G")
    private FeatureInput<Gff3Feature> gffFile;

    @Argument(doc="types to count", shortName = "T")
    private Set<String> type = new HashSet<>(Arrays.asList("CDS"));

    @Argument(doc="sample name to label counts", shortName = "N")
    private String name;

    @Argument(doc="gene_id key")
    private String gene_id_key = "gene_id";

    @Argument(doc = "which read corresponds to the transcription strand")
    private TrancriptionRead trancriptionRead = TrancriptionRead.R1;

    @Argument(doc = "Whether the rna is spliced.  If spliced, alignments must be from a splice aware aligner (such as star), if unspliced, alignments must be from " +
            "a non-splicing aligner (such as bwa). ")
    private boolean spliced = true;

    final private Map<Gff3Feature, MutablePair<Double, Double>> featureCountsMap = new LinkedHashMap<>();

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.MAPPED);
        readFilters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        readFilters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        readFilters.add(new InwardFragmentsFilter());
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

        final Set<String> geneIDSet = new HashSet<>();
        logger.info("collecting list of features");
        for (final SAMSequenceRecord contig : dict.getSequences()) {
            final List<Gff3Feature> contigFeatures = features.getFeatures(gffFile, new SimpleInterval(contig.getSequenceName(), 1, contig.getSequenceLength()));
            logger.info("collecting features in " + contig.getSequenceName());
            for (final Gff3Feature feature : contigFeatures) {
                if(!type.contains(feature.getType())) {
                    continue;
                }
                final String geneID=feature.getAttribute(gene_id_key);
                if (geneID == null) {
                    throw new UserException("no geneid field " + gene_id_key + " found in feature at " + feature.getContig() + ":" + feature.getStart() + "-" + feature.getEnd());
                }
                if (geneIDSet.contains(geneID)) {
                    logger.warn("geneid " + geneID + " is not unique, consists of multiple disjoint regions ");
                }
                featureCountsMap.put(getBareBonesFeature(feature), MutablePair.of(0.0,0.0));
                geneIDSet.add(geneID);
            }
        }
    }

    private void validateArguments() {
        CopyNumberArgumentValidationUtils.validateOutputFiles(outputCountsFile);
    }

    @Override
    public SimpleInterval getReadInterval(final GATKRead read) {
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

        if ((!read.isReverseStrand() || !inGoodPair(read))) {
            if (read.getStart()> 14382 && read.getEnd()<14836) {
                System.out.println();
            }
            final List<SimpleInterval> alignmentIntervals = getAlignmentIntervals(read);

            List<Gff3Feature> features = featureContext.getValues(gffFile);

            final Strand fragmentStrand = getFragmentStrand(read);

            final int basesOnReference = alignmentIntervals.stream().map(SimpleInterval::getLengthOnReference).reduce(0, Integer::sum);
            for (final Gff3Feature feature : features) {
                if(!type.contains(feature.getType())) {
                    continue;
                }
                if (feature.getStart() == 10643) {
                    System.out.println();
                }
                final boolean isSense = feature.getStrand() == Strand.NONE || feature.getStrand() == fragmentStrand;
                final MutablePair<Double, Double> currentCount = featureCountsMap.get(getBareBonesFeature(feature));

                for (final SimpleInterval interval : alignmentIntervals) {
                    if(interval.overlaps(feature)) {
                        if (isSense) {
                            currentCount.left += (float) interval.intersect(feature).size() / (float) basesOnReference;
                        } else {
                            currentCount.right += (float) interval.intersect(feature).size() / (float) basesOnReference;
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
            for (final Map.Entry<Gff3Feature, MutablePair<Double, Double>> entry : featureCountsMap.entrySet()) {
                writer.writeRecord(new FragmentCount(entry.getKey(),(int) Math.round(entry.getValue().left), true));
                if (entry.getKey().getStrand()!= Strand.NONE) {
                    writer.writeRecord(new FragmentCount(entry.getKey(), (int) Math.round(entry.getValue().right), false));
                }
            }
        } catch (IOException ex) {
            throw new UserException(ex.getMessage());
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    public class FragmentCountWriter extends TableWriter<FragmentCount> {

        public FragmentCountWriter(final Path file ) throws IOException {
            super(file, new TableColumnCollection("gene_id", "contig", "start", "stop", "strand", "sense_antisense", name+"_counts"));
        }

        protected void composeLine(final FragmentCount fragmentCount, final DataLine dataLine) {
            dataLine.set("contig", fragmentCount.gff3BaseData.getContig())
                    .set("start", fragmentCount.gff3BaseData.getStart())
                    .set("stop", fragmentCount.gff3BaseData.getEnd())
                    .set("strand", fragmentCount.gff3BaseData.getStrand().encode())
                    .set("sense_antisense", fragmentCount.sense? "sense" : "antisense")
                    .set(name+"_counts", fragmentCount.count)
                    .set("gene_id", fragmentCount.gff3BaseData.getAttribute(gene_id_key));
        }
    }

    private Gff3Feature getBareBonesFeature(final Gff3Feature feature) {
        return new Gff3FeatureImpl(feature.getContig(), "", "", feature.getStart(), feature.getEnd(), feature.getStrand(), 0, ImmutableMap.of(gene_id_key, feature.getAttribute(gene_id_key)));
    }

    public class FragmentCount {
        public final Gff3Feature gff3BaseData;
        //public final double count;
        public final int count;
        public final boolean sense;

//        FragmentCount(final Gff3Feature gtfFeature, final double count, final boolean sense) {
//            this.gtfFeature = gtfFeature;
//            this.count = count;
//            this.sense = sense;
//        }

        FragmentCount(final Gff3Feature gtfFeature, final int count, final boolean sense) {
            this.gff3BaseData = gtfFeature;
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


