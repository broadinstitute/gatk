package org.broadinstitute.hellbender.tools;


import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.*;

import java.io.File;

/**
 * Left-aligns indels in read data
 *
 * <p>
 * This tool left-aligns any indels in the read data contained in a BAM or CRAM file. The same indel can often be
 * placed at multiple positions and still represent the same haplotype.  While it is a commonly used convention to place
 * an indel at the left-most position, this doesn't always happen (either because upstream tools broke ties between
 * equivalent representations randomly or used different aligning conventions), so this tool can be used to left-align
 * them according to convention. </p>
 *
 * <h3>Input</h3>
 * <p>
 * A BAM or CRAM file to left-align.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A left-aligned BAM or CRAM file.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk LeftAlignIndels \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -O output.bam
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Left-aligns indels from reads in a SAM/BAM/CRAM file.",
        oneLineSummary = "Left-aligns indels from reads in a SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class LeftAlignIndels extends ReadWalker {

    @Argument(doc="Output BAM")
    private String OUTPUT;

    private SAMFileGATKReadWriter outputWriter = null;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(OUTPUT), true);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext ref, FeatureContext featureContext ) {
        // we can not deal with screwy records
        if ( read.isUnmapped() || read.numCigarElements() == 0 ) {
            outputWriter.addRead(read);
            return;
        }

        // move existing indels (for 1 indel reads only) to leftmost position within identical sequence
        int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
        if ( numBlocks == 2 ) {
            // We checked in onTraversalStart() that a reference is present, so ref.get() is safe
            Cigar newCigar = AlignmentUtils.leftAlignIndel(CigarUtils.trimReadToUnclippedBases(read.getCigar()), ref.getBases(), read.getBases(), 0, 0, true);
            newCigar = CigarUtils.reclipCigar(newCigar, read);
            read.setCigar(newCigar);
        }

        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}