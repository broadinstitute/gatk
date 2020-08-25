package org.broadinstitute.hellbender.tools.walkers.qc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
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
 * Prints read alignments in {@code samtools} pileup format.
 *
 * <p>This tool emulates the functionality of {@code samtools pileup}. It prints the alignments in a format that is very similar
 * to the {@code samtools} pileup format. For more details about the original format,
 * see the <a href="http://samtools.sourceforge.net/pileup.shtml">Samtools Pileup format
 * documentation</a>. The output comprises one line per genomic position, listing the
 * chromosome name, coordinate, reference base, bases from reads, and and corresponding  base qualities from reads.
 * In addition to these default fields,
 * additional information can be added to the output as extra columns.</p>
 *
 * Note that if the reference is omitted from the command, reference bases in the output will be <i>N</i>s.
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk Pileup \
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
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
    summary = "Prints read alignments in samtools pileup format. The output comprises one line per genomic position, " +
            "listing the chromosome name, coordinate, reference base, read bases, and read qualities. In addition to " +
            "these default fields, additional information can be added to the output as extra columns.",
    oneLineSummary = "Prints read alignments in samtools pileup format",
    programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
public final class Pileup extends LocusWalker {

    private static final String VERBOSE_DELIMITER = "@"; // it's ugly to use "@" but it's literally the only usable character not allowed in read names

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "An output file created to be created by the walker. Will overwrite the contents if the file exists."
    )
    public File outFile;

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

    private PrintStream out;

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
        final ReadPileup basePileup = alignmentContext.getBasePileup().makeFilteredPileup(pe -> !pe.isDeletion());
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
