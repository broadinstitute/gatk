package org.broadinstitute.hellbender.tools.walkers.qc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Print read alignments in Pileup-style format
 *
 * <p>This tool emulates the 'samtools pileup' command. It prints the alignment in a format that is very similar to the
 * Samtools pileup format (see the <a href="http://samtools.sourceforge.net/pileup.shtml">Pileup format
 * documentation</a> for more details about the original format). There is one line per genomic position, listing the
 * chromosome name, coordinate, reference base, read bases, and read qualities. In addition to these default fields,
 * additional information can be added to the output as extra columns; see options detailed below.</p>
 *
 * <h4>Emulated command:</h4>
 * <pre>
 *  samtools pileup -f in.ref.fasta -l in.site_list input.bam
 * </pre>
 *
 * <h3>Input</h3> <p> A BAM file and the interval to print. </p>
 *
 * <h3>Output</h3> <p> Alignment of reads formatted in the Pileup style. </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * ./gatk-launch Pileup  \\
 *   -R reference.fasta \
 *   -I my_reads.bam \
 *   -L chr1:257-267
 *   -O output.txt
 * </pre>
 * <h4>Expected output</h4>
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
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
    summary = "This tool emulates the 'samtools pileup' command. It prints the alignment in a format that is very "
        + "similar to the Samtools pileup format (see the documentation in http://samtools.sourceforge.net/pileup.shtml"
        + "for more details about the original format). There is one line per genomic position, listing the chromosome "
        + "name, coordinate, reference base, read bases, and read qualities. In addition to these default fields, "
        + "additional information can be added to the output as extra columns; see options detailed below.",
    oneLineSummary = "Print read alignments in Pileup-style format",
    programGroup = QCProgramGroup.class)
public final class Pileup extends LocusWalker {

    private static final String VERBOSE_DELIMITER = "@"; // it's ugly to use "@" but it's literally the only usable character not allowed in read names

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "An output file created by the walker. Will overwrite contents if file exists")
    public File outFile;

    /**
     * In addition to the standard pileup output, adds 'verbose' output too. The verbose output contains the number of
     * spanning deletions, and for each read in the pileup it has the read name, offset in the base string, read length,
     * and read mapping quality.  These per read items are delimited with an '@' character.
     */
    @Argument(fullName = "showVerbose", shortName = "verbose", doc = "Add an extra verbose section to the pileup output", optional = true)
    public boolean showVerbose = false;

    /**
     * This enables annotating the pileup to show overlaps with metadata from a Feature file(s). For example, if you provide a
     * VCF and there is a SNP at a given location covered by the pileup, the pileup output at that position will be
     * annotated with the corresponding source Feature identifier.
     */
    @Argument(fullName = "metadata", shortName = "metadata", doc = "Features file(s) containing metadata", optional = true)
    public List<FeatureInput<Feature>> metadata = new ArrayList<>();

    /**
     * Adds the length of the insert each base comes from to the output pileup. Here, "insert" refers to the DNA insert
     * produced during library generation before sequencing.
     */
    @Hidden
    @Argument(fullName = "outputInsertLength", shortName = "outputInsertLength", doc = "Output insert length", optional = true)
    public boolean outputInsertLength = false;

    private PrintStream out;

    @Override
    public boolean includeDeletions() {
        return false;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = super.getDefaultReadFilters();
        defaultFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        defaultFilters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        defaultFilters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        return defaultFilters;
    }

    @Override
    public void onTraversalStart() {
        try {
            out = new PrintStream(outFile);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final String features = getFeaturesString(featureContext);
        final ReadPileup basePileup = alignmentContext.getBasePileup();
        final StringBuilder s = new StringBuilder();
        s.append(String.format("%s %s",
                basePileup.getPileupString((hasReference()) ? (char) referenceContext.getBase() : 'N'),
                features));
        if (outputInsertLength) {
            s.append(" ").append(insertLengthOutput(basePileup));
        }
        if (showVerbose) {
            s.append(" ").append(createVerboseOutput(basePileup));
        }
        s.append("\n");
        out.print(s.toString());
    }

    /**
     * Get a string representation for the metadata
     *
     * @param featureContext Context for the features.
     *
     * @return String representation of the metadata
     */
    private String getFeaturesString(final FeatureContext featureContext) {
        String featuresString = featureContext.getValues(metadata).stream()
                .map(Feature::toString).collect(Collectors.joining(", "));
        if (!featuresString.isEmpty()) {
            featuresString = "[Feature(s): " + featuresString + "]";
        }
        return featuresString;
    }

    /**
     * Format the insert length for a pileup
     * @param pileup the pileup to format
     * @return a comma-separated string with insert lengths
     */
    @VisibleForTesting
    static String insertLengthOutput(final ReadPileup pileup) {
        return pileup.getReads().stream()
                .map(r -> String.valueOf(r.getFragmentLength()))
                .collect(Collectors.joining(","));
    }

    /**
     * Collect information for the overlapping reads, delimited by {@link #VERBOSE_DELIMITER}
     * @param pileup the pileup to format
     * @return formatted string with read information for the pileup
     */
    @VisibleForTesting
    static String createVerboseOutput(final ReadPileup pileup) {
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

    @Override
    public void closeTool() {
        if (out!=null) {
            out.close();
        }
    }
}
