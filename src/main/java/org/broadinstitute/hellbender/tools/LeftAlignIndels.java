package org.broadinstitute.hellbender.tools;


import htsjdk.samtools.*;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;

/**
 * Left-aligns indels from reads in a bam file.
 *
 * <p>
 * LeftAlignIndels is a tool that takes a bam file and left-aligns any indels inside it.  The same indel can often be
 * placed at multiple positions and still represent the same haplotype.  While a standard convention is to place an
 * indel at the left-most position this doesn't always happen, so this tool can be used to left-align them.
 *
 * <h3>Input</h3>
 * <p>
 * A bam file to left-align.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A left-aligned bam.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx3g -jar GenomeAnalysisTK.jar \
 *   -T LeftAlignIndels \
 *   -I input.bam \
 *   -o output.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        usage = "LeftAlignIndels is a tool that takes a bam file and left-aligns any indels inside it.  The same indel can often be\n" +
                "placed at multiple positions and still represent the same haplotype.  While a standard convention is to place an\n" +
                "indel at the left-most position this doesn't always happen, so this tool can be used to left-align them.",
        usageShort = "Left-aligns indels from reads in a bam file.",
        programGroup = ReadProgramGroup.class
)
public class LeftAlignIndels extends ReadWalker {

    @Argument(doc="Output BAM")
    private File OUTPUT;

    private SAMFileWriter outputWriter = null;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final SAMFileHeader outputHeader = ReadUtils.clone(getHeaderForReads());
        outputWriter = new SAMFileWriterFactory().makeWriter(outputHeader, true, OUTPUT, referenceArguments.getReferenceFile());
    }

    @Override
    public void apply( SAMRecord read, ReferenceContext ref, FeatureContext featureContext ) {
        // we can not deal with screwy records
        if ( read.getReadUnmappedFlag() || read.getCigar().numCigarElements() == 0 ) {
            outputWriter.addAlignment(read);
            return;
        }

        // move existing indels (for 1 indel reads only) to leftmost position within identical sequence
        int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
        if ( numBlocks == 2 ) {
            // We checked in onTraversalStart() that a reference is present, so ref.get() is safe
            Cigar newCigar = AlignmentUtils.leftAlignIndel(CigarUtils.unclipCigar(read.getCigar()), ref.getBases(), read.getReadBases(), 0, 0, true);
            newCigar = CigarUtils.reclipCigar(newCigar, read);
            read.setCigar(newCigar);
        }

        outputWriter.addAlignment(read);
    }
}