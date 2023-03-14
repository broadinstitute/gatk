package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler.Path;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler.PathPart;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler.Traversal;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.tools.walkers.ErrorCorrectHiFi.CallIterator;
import org.broadinstitute.hellbender.utils.read.UnalignedRead;
import org.broadinstitute.hellbender.utils.read.UnalignedRead.ByteSequence;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Assemble read fragments near SV breakpoint.",
        oneLineSummary = "Assemble read fragments near SV breakpoint.",
        programGroup = VariantEvaluationProgramGroup.class
)
public final class AssembleSVBreakpoint extends VariantWalker {
    public static final int PADDING = 500;
    public static final String SCRIPT_TEXT =
        "#!/bin/sh\n" +
        "    minimap2 -axmap-hifi $1 correctedReads.fq | samtools sort -OBAM - > correctedReads.bam &&\n" +
        "    samtools index correctedReads.bam &&\n" +
        "    minimap2 -axmap-hifi $1 assembledSequences.fa | samtools sort -OBAM - > assembledSequences.bam &&\n" +
        "    samtools index assembledSequences.bam &&\n" +
        "    gfaviz --no-groups --seg-as-arrow --seg-outline-width 0.1 assembly.gfa &&\n" +
        "    igv -g hg38 -l $2 assembledSequences.bam correctedReads.bam $3\n";

    public final File tmpDir = IOUtils.createTempDir("asvb");
    public final File readsFile = new File(tmpDir, "correctedReads.fq");
    public final File assembledSequencesFile = new File(tmpDir, "assembledSequences.fa");
    public final File scriptFile = new File(tmpDir, "run.sh");
    public final ProcessBuilder scriptRunner = new ProcessBuilder();

    @Argument(fullName = "mmi-index", shortName = "M", doc = "A minimap2 reference index")
    public File mmiIndex;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        IOUtils.writeByteArrayToFile(SCRIPT_TEXT.getBytes(), scriptFile);
        if ( !scriptFile.setExecutable(true) ) {
            throw new UserException("Can't make minimap2/igv script executable.");
        }
        scriptRunner.directory(tmpDir).inheritIO();
    }

    @Override
    public void apply( final VariantContext variant,
                       final ReadsContext readsContext,
                       final ReferenceContext refContext,
                       final FeatureContext featureContext ) {
        if ( variant.hasAttribute(GATKSVVCFConstants.SVTYPE) ) {
            final int start = variant.getStart();
            final int paddedStart = Math.max(1, start - PADDING);
            final int paddedEnd = variant.getEnd() + PADDING - 1;
            final SimpleInterval window =
                    new SimpleInterval(variant.getContig(), paddedStart, paddedEnd);
            final List<UnalignedRead> correctedReads = errorCorrect(window, readsContext);
            writeFASTQ(readsFile.getAbsolutePath(), correctedReads);
            final List<UnalignedRead> assembledSequences = assemble(correctedReads);
            writeFASTA(assembledSequencesFile.getAbsolutePath(), assembledSequences);
            alignAndDisplay(window);
        }
    }

    private List<UnalignedRead> errorCorrect( final SimpleInterval interval,
                                              final ReadsContext readsContext ) {
        final Iterator<GATKRead> readsIterator = readsContext.iterator(interval);
        final List<GATKRead> reads = new ArrayList<>();
        while ( readsIterator.hasNext() ) {
            reads.add(readsIterator.next());
        }
        final List<CallIterator> callIteratorList =
                ErrorCorrectHiFi.errorCorrectWindow(reads, interval.getStart(), interval.getEnd());
        final List<UnalignedRead> correctedReads = new ArrayList<>(callIteratorList.size());
        for ( final CallIterator callIterator : callIteratorList ) {
            correctedReads.add(callIterator.getUnalignedRead());
        }
        return correctedReads;
    }

    private List<UnalignedRead> assemble( final List<UnalignedRead> reads ) {
        final LocalAssembler assembler = new LocalAssembler(2*PADDING, reads);
        final List<Traversal> traversals = assembler.getPathedTraversals();
        assembler.writeGFA(new GATKPath(new File(tmpDir, "assembly.gfa").getAbsolutePath()), traversals);
        final List<UnalignedRead> assembledSequences = new ArrayList<>();
        int idx = 0;
        for ( final Traversal trav : traversals ) {
            final ByteSequence calls = new ByteSequence(trav.getSequence().getBytes());
            assembledSequences.add(new UnalignedRead("t" + ++idx, calls, null));
        }
        return assembledSequences;
    }

    private void writeFASTQ( final String filename, final List<UnalignedRead> reads ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(filename)) ) {
            for ( final UnalignedRead read : reads ) {
                read.writeFASTQ(writer);
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write to " + filename, ioe);
        }
    }

    private void writeFASTA( final String filename, final List<UnalignedRead> reads ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(filename)) ) {
            for ( final UnalignedRead read : reads ) {
                read.writeFASTA(writer);
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write to " + filename, ioe);
        }
    }

    private void alignAndDisplay( final SimpleInterval window ) {
        try {
            scriptRunner.command(scriptFile.getAbsolutePath(),
                                    mmiIndex.getAbsolutePath(),
                                    window.toString(),
                                    readArguments.getReadPaths().get(0).toAbsolutePath().toString());
            final Process process = scriptRunner.start();
            final int exitValue = process.waitFor();
            if ( exitValue != 0 ) {
                throw new UserException("Script that runs minimap2 and igv returned " + exitValue);
            }
        } catch ( final IOException | InterruptedException ex ) {
            throw new UserException("Script that runs minimap2 and igv failed.", ex);
        }
    }
}
