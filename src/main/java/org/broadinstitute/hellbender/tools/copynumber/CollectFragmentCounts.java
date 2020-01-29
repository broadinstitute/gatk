package org.broadinstitute.hellbender.tools.copynumber;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Feature;
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


import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
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
        for (final SAMSequenceRecord contig : dict.getSequences()) {
            final List<Gff3Feature> contigFeatures = features.getFeatures(gffFile, new SimpleInterval(contig.getSequenceName(), 1, contig.getSequenceLength()));
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
                featureCountsMap.put(feature, MutablePair.of(0.0,0.0));
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

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {

        if ((!read.isReverseStrand() || !inGoodPair(read))) {
            final SimpleInterval fragment_interval = getReadInterval(read);

            List<Gff3Feature> features = featureContext.getValues(gffFile);
            //final Strand readStrand = read.isFirstOfPair() ? (read.isReverseStrand() ? Strand.NEGATIVE : Strand.POSITIVE) : (read.isReverseStrand() ? Strand.POSITIVE : Strand.NEGATIVE);
            Strand readStrand = read.isFirstOfPair() ? (read.isReverseStrand() ? Strand.NEGATIVE : Strand.POSITIVE) : (read.isReverseStrand() ? Strand.POSITIVE : Strand.NEGATIVE);
            //reverse stand for some reason...
            readStrand = readStrand == Strand.POSITIVE? Strand.NEGATIVE : Strand.POSITIVE;
            /*
            let's figure out the strand backwards (and possibly inconsistenly) yay!!!
             */
//            Strand readStrand = Strand.POSITIVE;
//            if (read.isFirstOfPair()) {
//                if (!read.isReverseStrand() || read.mateIsReverseStrand()) {
//                    readStrand = Strand.NEGATIVE;
//                }
//                if (read.isReverseStrand() || !read.mateIsReverseStrand()) {
//                    readStrand = Strand.POSITIVE;
//                }
//            }
//
//            if (read.isSecondOfPair()) {
//                if (read.isReverseStrand() || !read.mateIsReverseStrand()) {
//                    readStrand = Strand.NEGATIVE;
//                }
//                if (!read.isReverseStrand() || read.mateIsReverseStrand()) {
//                    readStrand = Strand.POSITIVE;
//                }
//            }

            LinkedList<SimpleInterval> codingIntervals = new LinkedList<>();
            int overlappingBases = 0;
            for (final Gff3Feature feature : features) {
                if(!type.contains(feature.getType())) {
                    continue;
                }
                overlappingBases += fragment_interval.intersect(feature).size();

                final List<SimpleInterval> overlappingCodingIntervals = codingIntervals.stream().filter(i -> i.overlaps(fragment_interval.intersect(feature))).collect(Collectors.toList());
                if (overlappingCodingIntervals.size() == 0) {
                    codingIntervals.add(fragment_interval.intersect(feature));
                } else {
                    SimpleInterval overlappingCodingInterval = fragment_interval.intersect(feature);
                    for (final SimpleInterval overlappingInterval : overlappingCodingIntervals) {
                        codingIntervals.remove(overlappingInterval);
                        overlappingCodingInterval = overlappingCodingInterval.mergeWithContiguous(overlappingInterval);
                    }
                    codingIntervals.add(overlappingCodingInterval);
                }
            }
            final int codingBases = codingIntervals.stream().map(SimpleInterval::size).reduce(0, Integer::sum);
            final int utrBases = fragment_interval.size() - codingBases;
            overlappingBases += utrBases;

            for (final Gff3Feature feature : features) {
                if(!type.contains(feature.getType())) {
                    continue;
                }

                // for each assign fraction of fragment feature overlaps
                final SimpleInterval intersection = fragment_interval.intersect(feature);

                final boolean isSense = feature.getStrand() == Strand.NONE || feature.getStrand() == readStrand;
                final MutablePair<Double, Double> currentCount = featureCountsMap.get(feature);

                if (isSense) {

                    //currentCount.left = currentCount.left + (float) intersection.getLengthOnReference() / (float) fragment_interval.getLengthOnReference();
                    currentCount.left += (float) intersection.getLengthOnReference() / (float) overlappingBases;
                } else {
                    //currentCount.right = currentCount.right + (float) intersection.getLengthOnReference() / (float) fragment_interval.getLengthOnReference();
                    currentCount.right += (float) intersection.getLengthOnReference() / (float) overlappingBases;
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
                //writer.writeRecord(new FragmentCount(entry.getKey(),entry.getValue().left, true));
                writer.writeRecord(new FragmentCount(entry.getKey(),(int) Math.round(entry.getValue().left), true));
                //writer.writeRecord(new FragmentCount(entry.getKey(),entry.getValue().right, false));
                writer.writeRecord(new FragmentCount(entry.getKey(),(int) Math.round(entry.getValue().right), false));
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
            dataLine.set("contig", fragmentCount.gtfFeature.getContig())
                    .set("start", fragmentCount.gtfFeature.getStart())
                    .set("stop", fragmentCount.gtfFeature.getEnd())
                    .set("strand", fragmentCount.gtfFeature.getStrand().encode())
                    .set("sense_antisense", fragmentCount.sense? "sense" : "antisense")
                    .set(name+"_counts", fragmentCount.count)
                    .set("gene_id", fragmentCount.gtfFeature.getAttribute(gene_id_key));
        }
    }

    public class FragmentCount {
        public final Gff3Feature gtfFeature;
        //public final double count;
        public final int count;
        public final boolean sense;

//        FragmentCount(final Gff3Feature gtfFeature, final double count, final boolean sense) {
//            this.gtfFeature = gtfFeature;
//            this.count = count;
//            this.sense = sense;
//        }

        FragmentCount(final Gff3Feature gtfFeature, final int count, final boolean sense) {
            this.gtfFeature = gtfFeature;
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
}


