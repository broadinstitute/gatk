package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** This LocalAssemblyHandler aligns assembly contigs with BWA, along with some optional writing of intermediate results. */
public final class FermiLiteAssemblyHandler implements FindBreakpointEvidenceSpark.LocalAssemblyHandler {
    private static final long serialVersionUID = 1L;
    private final String alignerIndexFile;
    private final int maxFastqSize;
    private final String fastqDir;
    private final boolean writeGFAs;

    public FermiLiteAssemblyHandler( final String alignerIndexFile, final int maxFastqSize,
                                     final String fastqDir, final boolean writeGFAs ) {
        this.alignerIndexFile = alignerIndexFile;
        this.maxFastqSize = maxFastqSize;
        this.fastqDir = fastqDir;
        this.writeGFAs = writeGFAs;
    }

    @Override
    public AlignedAssemblyOrExcuse apply( final Tuple2<Integer, List<SVFastqUtils.FastqRead>> intervalAndReads ) {
        final int intervalID = intervalAndReads._1();
        final List<SVFastqUtils.FastqRead> readsList = intervalAndReads._2();

        final int fastqSize = readsList.stream().mapToInt(FastqRead -> FastqRead.getBases().length).sum();
        if (fastqSize > maxFastqSize) {
            return new AlignedAssemblyOrExcuse(intervalID, "no assembly -- too big (" + fastqSize + " bytes).");
        }

        if ( fastqDir != null ) {
            final String fastqName = String.format("%s/%s.fastq", fastqDir, AlignedAssemblyOrExcuse.formatAssemblyID(intervalID));
            final ArrayList<SVFastqUtils.FastqRead> sortedReads = new ArrayList<>(readsList);
            sortedReads.sort(Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
            SVFastqUtils.writeFastqFile(fastqName, sortedReads.iterator());
        }
        final long timeStart = System.currentTimeMillis();
        final FermiLiteAssembly assembly = new FermiLiteAssembler().createAssembly(readsList);
        final int secondsInAssembly = (int)((System.currentTimeMillis() - timeStart + 500)/1000);
        if ( fastqDir != null && writeGFAs ) {
            final String gfaName =  String.format("%s/%s.gfa",fastqDir, AlignedAssemblyOrExcuse.formatAssemblyID(intervalID));
            try ( final OutputStream os = BucketUtils.createFile(gfaName) ) {
                assembly.writeGFA(os);
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write "+gfaName, ioe);
            }
        }
        final List<byte[]> tigSeqs =
                assembly.getContigs().stream()
                        .map(FermiLiteAssembly.Contig::getSequence)
                        .collect(SVUtils.arrayListCollector(assembly.getNContigs()));
        try ( final BwaMemAligner aligner = new BwaMemAligner(BwaMemIndexCache.getInstance(alignerIndexFile)) ) {
            aligner.setIntraCtgOptions();
            final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(tigSeqs);
            return new AlignedAssemblyOrExcuse(intervalID, assembly, secondsInAssembly, alignments);
        }
    }
}
