package org.broadinstitute.hellbender.tools.walkers.pacbio;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.nio.file.Path;

/**
 * Quickly count errors
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A collection of BAM files where care has been keep reads from the same ZMW in the same file</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Quickly count errors</h4>
 * <pre>
 *   gatk QuicklyCountErrors \
 *     -I input.bam \
 *     -O output.txt
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Quickly count errors",
        oneLineSummary = "Quickly count errors",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class QuicklyCountErrors extends ReadWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Path to which variants should be written")
    public String outPathName;
    //final Path outPath = IOUtils.getPath(outPathName);

    private long numMismatches = 0;
    private long numInsertions = 0;
    private long numDeletions = 0;
    private long numBases = 0;

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        boolean xseen = false;
        for (CigarElement ce : read.getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.X)) {
                xseen = true;
                numMismatches++;
            } else if (ce.getOperator().equals(CigarOperator.D)) {
                numDeletions += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.I)) {
                numInsertions += ce.getLength();
            }
        }

        if (!xseen && read.hasAttribute("NM")) {
            numMismatches += read.getAttributeAsInteger("NM");
        }

        numBases += read.getLength();
    }

    @Override
    public void closeTool() {
        try {
            //PrintStream out = new BufferedOutputStream(new PrintStream(outPathName));
            PrintStream out = new PrintStream(outPathName);

            out.println(Joiner.on(" ").join("numMismatches", "numDeletions", "numInsertions", "numBases"));
            out.println(Joiner.on(" ").join(numMismatches, numDeletions, numInsertions, numBases));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        logger.info("numMismatches: " + numMismatches + " numDeletions: " + numDeletions + " numInsertions: " + numInsertions + " numBases:" + numBases);
    }
}
