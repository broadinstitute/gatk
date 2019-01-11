package org.broadinstitute.hellbender.tools.spark;

import htsjdk.tribble.Feature;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerContext;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerSpark;
import scala.Tuple3;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Prints read alignments in {@code samtools} pileup format. The tool leverages the Spark framework
 * for faster operation.
 *
 * <p>This tool emulates  the functionality of {@code samtools pileup}. It prints the alignments in a format that is very similar
 * to the {@code samtools} pileup format; see the <a href="http://samtools.sourceforge.net/pileup.shtml">Samtools Pileup format
 * documentation</a> for more details about the original format. The output comprises one line per genomic position, listing the
 * chromosome name, coordinate, reference base, bases from reads, and corresponding base qualities from reads.
 * In addition to these default fields,
 * additional information can be added to the output as extra columns.</p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk PileupSpark \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -O output.txt
 * </pre>
 *
 * <h4>Emulated command:</h4>
 * <pre>
 *  samtools pileup -f reference.fasta input.bam
 * </pre>
 *
 * <h4>Typical output format</h4>
 * <pre>
 *     chr1 257 A CAA '&=
 *     chr1 258 C TCC A:=
 *     chr1 259 C CCC )A=
 *     chr1 260 C ACC (=<
 *     chr1 261 T TCT '44
 *     chr1 262 A AAA '?:
 *     chr1 263 A AGA 1'6
 *     chr1 264 C TCC 987
 *     chr1 265 C CCC (@(
 *     chr1 266 C GCC ''=
 *     chr1 267 T AAT 7%>
 * </pre>
 * <p>
 * This tool can be run without explicitly specifying Spark options.
 * That is to say, the given example command without Spark options will run locally.
 * See <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.
 * </p>
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        summary = "Prints read alignments in samtools pileup format. The tool leverages the Spark framework for " +
                "faster operationThe output comprises one line per genomic position, listing the chromosome name, " +
                "coordinate, reference base, read bases, and read qualities. In addition to these default fields, " +
                "additional information can be added to the output as extra columns.",
        oneLineSummary = "Prints read alignments in samtools pileup format",
        programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class PileupSpark extends LocusWalkerSpark {
    private static final long serialVersionUID = 1L;

    private static final String VERBOSE_DELIMITER = "@"; // it's ugly to use "@" but it's literally the only usable character not allowed in read names

    @Override
    public boolean requiresReads() { return true; }

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc="The output directory to which the scattered output will be written."
    )
    protected String outputFile;

    /**
     * In addition to the standard pileup output, adds 'verbose' output too. The verbose output contains the number of
     * spanning deletions, and for each read in the pileup it has the read name, offset in the base string, read length,
     * and read mapping quality.  These per read items are delimited with an '@' character.
     */
    @Argument(
            fullName = "show-verbose",
            shortName = "verbose",
            doc = "Add extra informative columns to the pileup output. The verbose output contains the number of " +
                    "spanning deletions, and for each read in the pileup it has the read name, offset in the base " +
                    "string, read length, and read mapping quality.  These per read items are delimited with an '@' " +
                    "character.",
            optional = true
    )
    public boolean showVerbose = false;

    /**
     * This enables annotating the pileup to show overlaps with metadata from a Feature file(s). For example, if the
     * user provide a VCF and there is a SNP at a given location covered by the pileup, the pileup output at that
     * position will be annotated with the corresponding source Feature identifier.
     */
    @Argument(
            fullName = "metadata",
            shortName = "metadata",
            doc = "Features file(s) containing metadata. The overlapping sites will be annotated with the corresponding" +
                    " source Feature identifier.",
            optional = true
    )
    public List<FeatureInput<Feature>> metadata = new ArrayList<>();

    /**
     * Adds the length of the insert each base comes from to the output pileup. Here, "insert" refers to the DNA insert
     * produced during library generation before sequencing.
     */
    @Argument(
            fullName = "output-insert-length",
            shortName = "output-insert-length",
            doc = "If enabled, inserts lengths will be added to the output pileup.",
            optional = true
    )
    public boolean outputInsertLength = false;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        List<ReadFilter> filterList = new ArrayList<>(5);
        filterList.add(ReadFilterLibrary.MAPPED);
        filterList.add(ReadFilterLibrary.NOT_DUPLICATE);
        filterList.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filterList.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        filterList.add(new WellformedReadFilter());
        return filterList;
    }

    @Override
    protected void processAlignments(JavaRDD<LocusWalkerContext> rdd, JavaSparkContext ctx) {
        JavaRDD<String> lines = rdd.map(pileupFunction(metadata, outputInsertLength, showVerbose));
        if (numReducers != 0) {
            lines = lines.coalesce(numReducers);
        }
        lines.saveAsTextFile(outputFile);
    }

    private static Function<LocusWalkerContext, String> pileupFunction(List<FeatureInput<Feature>> metadata,
                                                                       boolean outputInsertLength, boolean showVerbose) {
        return (Function<LocusWalkerContext, String>) context -> {
            AlignmentContext alignmentContext = context.getAlignmentContext();
            ReferenceContext referenceContext = context.getReferenceContext();
            FeatureContext featureContext = context.getFeatureContext();
            final String features = getFeaturesString(featureContext, metadata);
            final ReadPileup basePileup = alignmentContext.getBasePileup();
            final StringBuilder s = new StringBuilder();
            s.append(String.format("%s %s",
                    basePileup.getPileupString((referenceContext.hasBackingDataSource()) ? (char) referenceContext.getBase() : 'N'),
                    features));
            if (outputInsertLength) {
                s.append(" ").append(insertLengthOutput(basePileup));
            }
            if (showVerbose) {
                s.append(" ").append(createVerboseOutput(basePileup));
            }
            s.append("\n");
            return s.toString();
        };
    }

    private static String getFeaturesString(final FeatureContext featureContext, List<FeatureInput<Feature>> metadata) {
        String featuresString = featureContext.getValues(metadata).stream()
                .map(Feature::toString).collect(Collectors.joining(", "));
        if (!featuresString.isEmpty()) {
            featuresString = "[Feature(s): " + featuresString + "]";
        }
        return featuresString;
    }

    private static String insertLengthOutput(final ReadPileup pileup) {
        return pileup.getReads().stream()
                .map(r -> String.valueOf(r.getFragmentLength()))
                .collect(Collectors.joining(","));
    }

    private static String createVerboseOutput(final ReadPileup pileup) {
        final StringBuilder sb = new StringBuilder();
        boolean isFirst = true;
        sb.append(pileup.getNumberOfElements(PileupElement::isDeletion));
        sb.append(" ");
        for (final PileupElement p : pileup) {
            if (isFirst) {
                isFirst = false;
            } else {
                sb.append(",");
            }
            sb.append(p.getRead().getName());
            sb.append(VERBOSE_DELIMITER);
            sb.append(p.getOffset());
            sb.append(VERBOSE_DELIMITER);
            sb.append(p.getRead().getLength());
            sb.append(VERBOSE_DELIMITER);
            sb.append(p.getRead().getMappingQuality());
        }
        return sb.toString();
    }
}
