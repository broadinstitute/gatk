package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import scala.Tuple2;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collector;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.RunSGAViaProcessBuilderOnSpark.ContigsCollection;
import static org.broadinstitute.hellbender.tools.spark.sv.RunSGAViaProcessBuilderOnSpark.ContigsCollection.ContigID;
import static org.broadinstitute.hellbender.tools.spark.sv.RunSGAViaProcessBuilderOnSpark.ContigsCollection.ContigSequence;

public class ContigAligner implements Closeable {

    static String referencePath;

    final BwaIndex index;
    final BwaMem bwaMem;

    private static final Logger log = LogManager.getLogger(ContigAligner.class);

    static {
        BWANativeLibrary.load();
    }

    public ContigAligner(final String referenceFilePath) throws IOException {
        referencePath = referenceFilePath;

        try {
            index = new BwaIndex(LocalizedReference.INSTANCE);
            log.info("Created BWA index");
        } catch (final IOException e) {
            throw new GATKException("Unable to load reference", e);
        }
        bwaMem = new BwaMem(index);
        log.info("Created BWA MEM");
    }

    public List<AssembledBreakpoint> alignContigs (final ContigsCollection contigsCollection) {
        final List<AssembledBreakpoint> assembledBreakpoints = new ArrayList<>();
        try {
            for(final Tuple2<ContigID, ContigSequence> contigInfo : contigsCollection.getContents()) {
                final String contigId = contigInfo._1.toString();
                final byte[] sequence = contigInfo._2.toString().getBytes();
                final AlnRgn[] alnRgns = bwaAlignSequence(bwaMem, contigId, sequence);

                log.info("alnRgns : " + (alnRgns == null ? "null" : alnRgns.length));
                // todo: parse cigar for internal indels
                if (alnRgns.length > 1) {

                    // filter out secondary alignments, convert to AlignmentRegion objects and sort by alignment start pos
                    final List<AlignmentRegion> alignmentRegionList = Arrays.stream(alnRgns)
                            .filter(a -> a.getSecondary() < 0)
                            .map(AlignmentRegion::new)
                            .sorted(Comparator.comparing(a -> a.startInAssembledContig))
                            .collect(arrayListCollector(alnRgns.length));


                    final Iterator<AlignmentRegion> iterator = alignmentRegionList.iterator();
                    if ( iterator.hasNext() ) {
                        AlignmentRegion current = iterator.next();
                        while ( iterator.hasNext() ) {
                            final AlignmentRegion previous = current;
                            current = iterator.next();
                            String homology = "NA";
                            if (previous.endInAssembledContig >= current.startInAssembledContig) {
                                homology = new String(Arrays.copyOfRange(sequence, current.startInAssembledContig - 1, previous.endInAssembledContig));
                            }

                            String insertedSequence = "NA";
                            if (previous.endInAssembledContig < current.startInAssembledContig - 1) {
                                insertedSequence = new String(Arrays.copyOfRange(sequence, previous.endInAssembledContig, current.startInAssembledContig));
                            }
                            final AssembledBreakpoint assembledBreakpoint = new AssembledBreakpoint(contigId, previous, current, homology, insertedSequence);
                            assembledBreakpoints.add(assembledBreakpoint);
                        }
                    }
                }
            }
        } catch (final IOException e) {
            throw new GATKException("could not execute BWA");
        }

        return assembledBreakpoints;
    }

    private Collector<AlignmentRegion, ?, ArrayList<AlignmentRegion>> arrayListCollector(final int size) {
        return Collectors.toCollection( () -> new ArrayList<>(size));
    }

    /**
     * Wrap a contig sequence in a ShortRead object and pass it to BWA to align
     */
    private AlnRgn[] bwaAlignSequence(final BwaMem bwaMem, final String contigId, final byte[] sequence) throws IOException {
        final ShortRead contigShortRead = new ShortRead(contigId, sequence, qualSequence(sequence.length));
        log.debug("Calling bwaMem.align");
        return bwaMem.align(contigShortRead);
    }

    /**
     * Generate a bogus base quality sequence to pass in for the aligned contig (since the jBWA API requires that reads have qualities)
     */
    private byte[] qualSequence(final int length) {
        final byte[] quals = new byte[length];
        Arrays.fill(quals, (byte)'A');
        return quals;
    }

    public void close() {
        log.info("closing BWA mem and index");
        bwaMem.dispose();
        index.close();
    }

    /**
     * Holds information about a split alignment of a contig, which may represent an SV breakpoint. Each AssembledBreakpoint
     * represents the junction on the contig of two aligned regions. For example, if a contig aligns to three different regions
     * of the genome (with one primary and two supplementary alignment records), there will be two AssembledBreakpoint
     * objects created, one to represent each junction between alignment regions:
     *
     * Example Contig:
     * ACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTG
     * Alignment regions:
     * |---------1:100-200------------|
     *                                 |----------2:100-200------------------|
     *                                                                       |----------3:100-200-----------------|
     * Assmbled breakpoints:
     * 1) links 1:100-200 to 2:100-200
     * 2) links 2:100-200 to 3:100-200
     *
     * Inserted sequence contains portions of the contig that are aligned to neither region, and therefore may be inserted in
     * the sample. For example, a translocation breakpoint with a microinsertion:
     *
     * Contig:
     * ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
     * Alignment regions:
     * |-----1:100-200-------|
     *                          |----2:100-200-----|
     * Inserted sequence:
     *  GA
     *
     * Homology represents ambiguity about the exact location of the breakpoint. For example, in this case one alignment
     * region ends with "AC" and the next begins with AC, so we don't know if the AC truly belongs with the first or
     * second alignment region.
     *
     * Contig:
     * ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
     * Alignment regions:
     * |-----1:100-200-------|
     *                    |-----2:100-200----------|
     * Homology:
     *  AC
     */
    static class AssembledBreakpoint {
        String contigId;
        AlignmentRegion region1;
        AlignmentRegion region2;
        String insertedSequence;
        String homology;

        public AssembledBreakpoint(final String contigId, final AlignmentRegion region1, final AlignmentRegion region2, final String homology, final String insertedSequence) {
            this.contigId = contigId;
            this.region1 = region1;
            this.region2 = region2;
            this.insertedSequence = insertedSequence;
            this.homology = homology;
        }

        @Override
        public String toString() {
            return contigId +
                    "\t" +
                    region1.toString() +
                    "\t" +
                    region2.toString() +
                    "\t" +
                    insertedSequence +
                    "\t" +
                    homology;
        }
    }

    static class AlignmentRegion {

        final Cigar cigar;
        final boolean forwardStrand;
        final SimpleInterval referenceInterval;
        final int mqual;
        final int startInAssembledContig;
        final int endInAssembledContig;
        final int assembledContigLength;


        public AlignmentRegion(final AlnRgn alnRgn) {
            this.forwardStrand = alnRgn.getStrand() == '+';
            final Cigar alignmentCigar = TextCigarCodec.decode(alnRgn.getCigar());
            this.cigar = forwardStrand ? alignmentCigar : CigarUtils.invertCigar(alignmentCigar);
            this.referenceInterval = new SimpleInterval(alnRgn.getChrom(), (int) alnRgn.getPos() + 1, (int) (alnRgn.getPos() + 1 + cigar.getReferenceLength()));
            this.mqual = alnRgn.getMQual();
            this.assembledContigLength = cigar.getReadLength();
            this.startInAssembledContig = startOfAlignmentInContig(cigar);
            this.endInAssembledContig = endOfAlignmentInContig(assembledContigLength, cigar);
        }

        private static int startOfAlignmentInContig(final Cigar cigar) {
            return getClippedBases(true, cigar) + 1;
        }

        private static int endOfAlignmentInContig(final int assembledContigLength, final Cigar cigar) {
            return assembledContigLength - getClippedBases(false, cigar);
        }

        private static int getClippedBases(final boolean fromStart, final Cigar cigar) {
            int posInContig = 0;
            int j = fromStart ? 0 : cigar.getCigarElements().size() - 1;
            final int offset = fromStart ? 1 : -1;
            CigarElement ce = cigar.getCigarElement(j);
            while(ce.getOperator().isClipping()){
                posInContig += ce.getLength();
                j += offset;
                ce = cigar.getCigarElement(j);
            }
            return posInContig;
        }

        @Override
        public String toString() {
            return referenceInterval.getContig() +
                    "\t" +
                    referenceInterval.getStart() +
                    "\t" +
                    referenceInterval.getEnd() +
                    "\t" +
                    (forwardStrand ? "+" : "-") +
                    "\t" +
                    cigar.toString() +
                    "\t" +
                    mqual +
                    "\t" +
                    startInAssembledContig +
                    "\t" +
                    endInAssembledContig;
        }
    }

    private static class LocalizedReference {
        static File INSTANCE;

        static {
            try {
                INSTANCE = BucketUtils.isHadoopUrl(ContigAligner.referencePath) ? BwaSparkEngine.localizeReferenceAndBwaIndexFiles(ContigAligner.referencePath) : new File(ContigAligner.referencePath);
            } catch (final IOException e) {
                throw new GATKException("unable to localize reference", e);
            }
        }
    }

}
