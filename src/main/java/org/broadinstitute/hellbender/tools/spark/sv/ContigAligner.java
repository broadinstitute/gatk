package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
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

    public List<AlignmentRegion> alignContigs(String breakpointId, final ContigsCollection contigsCollection) {
        final List<AlignmentRegion> alignedContigs = new ArrayList<>();
        try {
            for(final Tuple2<ContigID, ContigSequence> contigInfo : contigsCollection.getContents()) {
                final String contigId = contigInfo._1.toString();
                final byte[] sequence = contigInfo._2.toString().getBytes();
                final AlnRgn[] alnRgns = bwaAlignSequence(bwaMem, contigId, sequence);

                log.info("alnRgns : " + (alnRgns == null ? "null" : alnRgns.length));
                // filter out secondary alignments, convert to AlignmentRegion objects and sort by alignment start pos
                final List<AlignmentRegion> alignmentRegionList = Arrays.stream(alnRgns)
                        .filter(a -> a.getSecondary() < 0)
                        .map(a -> new AlignmentRegion(breakpointId, contigId, a))
                        .sorted(Comparator.comparing(a -> a.startInAssembledContig))
                        .collect(arrayListCollector(alnRgns.length));
                 alignedContigs.addAll(alignmentRegionList);
            }
        } catch (final IOException e) {
            throw new GATKException("could not execute BWA");
        }

        return alignedContigs;
    }

    @VisibleForTesting
    public static List<AssembledBreakpoint> getAssembledBreakpointsFromAlignmentRegions(final byte[] sequence, final List<AlignmentRegion> alignmentRegionList, final Integer minAlignLength) {
        final List<AssembledBreakpoint> results = new ArrayList<>(alignmentRegionList.size() - 1);
        final Iterator<AlignmentRegion> iterator = alignmentRegionList.iterator();
        final List<String> insertionAlignmentRegions = new ArrayList<>();
        if ( iterator.hasNext() ) {
            AlignmentRegion current = iterator.next();
            while (treatAlignmentRegionAsInsertion(current) && iterator.hasNext()) {
                current = iterator.next();
            }
            while ( iterator.hasNext() ) {
                final AlignmentRegion next = iterator.next();
                if (currentAlignmentRegionIsTooSmall(current, next, minAlignLength)) {
                    continue;
                }

                if (treatNextAlignmentRegionInPairAsInsertion(current, next, minAlignLength)) {
                    if (iterator.hasNext()) {
                        insertionAlignmentRegions.add(next.toPackedString());
                        // todo: track alignments of skipped regions for classification as duplications, mei's etc.
                        continue;
                    } else {
                        break;
                    }
                }

                final AlignmentRegion previous = current;
                current = next;

                final byte[] sequenceCopy = Arrays.copyOf(sequence, sequence.length);

                String homology = "";
                if (previous.endInAssembledContig >= current.startInAssembledContig) {
                    final byte[] homologyBytes = Arrays.copyOfRange(sequenceCopy, current.startInAssembledContig - 1, previous.endInAssembledContig);
                    if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                        SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
                    }
                    homology = new String(homologyBytes);
                }

                String insertedSequence = "";
                if (previous.endInAssembledContig < current.startInAssembledContig - 1) {

                    final int insertionStart;
                    final int insertionEnd;

                    insertionStart = previous.endInAssembledContig + 1;
                    insertionEnd = current.startInAssembledContig - 1;

                    final byte[] insertedSequenceBytes = Arrays.copyOfRange(sequenceCopy, insertionStart - 1, insertionEnd);
                    if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                        SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
                    }
                    insertedSequence = new String(insertedSequenceBytes);
                }
                final AssembledBreakpoint assembledBreakpoint = new AssembledBreakpoint(current.contigId, previous, current, insertedSequence, homology, insertionAlignmentRegions);

                results.add(assembledBreakpoint);
            }
        }
        return results;
    }

    private static boolean currentAlignmentRegionIsTooSmall(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return current.referenceInterval.size() - current.overlapOnContig(next) < minAlignLength;
    }

    protected static boolean treatNextAlignmentRegionInPairAsInsertion(AlignmentRegion current, AlignmentRegion next, final Integer minAlignLength) {
        return treatAlignmentRegionAsInsertion(next) ||
                (next.referenceInterval.size() - current.overlapOnContig(next) < minAlignLength) ||
                current.referenceInterval.contains(next.referenceInterval) ||
                next.referenceInterval.contains(current.referenceInterval);
    }

    private static boolean treatAlignmentRegionAsInsertion(final AlignmentRegion next) {
        return next.mqual < 60;
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

    static class BreakpointAllele {
        SimpleInterval leftAlignedLeftBreakpoint;
        SimpleInterval leftAlignedRightBreakpoint;
        String insertedSequence;
        String homology;
        boolean fiveToThree;
        boolean threeToFive;

        public BreakpointAllele(final SimpleInterval leftAlignedLeftBreakpoint, final SimpleInterval leftAlignedRightBreakpoint, final String insertedSequence, final String homology, final boolean fiveToThree, final boolean threeToFive) {
            this.leftAlignedLeftBreakpoint = leftAlignedLeftBreakpoint;
            this.leftAlignedRightBreakpoint = leftAlignedRightBreakpoint;
            this.insertedSequence = insertedSequence;
            this.homology = homology;
            this.fiveToThree = fiveToThree;
            this.threeToFive = threeToFive;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            final BreakpointAllele that = (BreakpointAllele) o;
            return fiveToThree == that.fiveToThree &&
                    threeToFive == that.threeToFive &&
                    Objects.equals(leftAlignedLeftBreakpoint, that.leftAlignedLeftBreakpoint) &&
                    Objects.equals(leftAlignedRightBreakpoint, that.leftAlignedRightBreakpoint) &&
                    Objects.equals(insertedSequence, that.insertedSequence) &&
                    Objects.equals(homology, that.homology);
        }

        @Override
        public int hashCode() {
            return Objects.hash(leftAlignedLeftBreakpoint, leftAlignedRightBreakpoint, insertedSequence, homology, fiveToThree, threeToFive);
        }
    }

    static class AlignmentRegion {

        final String contigId;
        final String breakpointId;
        final Cigar cigar;
        final boolean forwardStrand;
        final SimpleInterval referenceInterval;
        final int mqual;
        final int startInAssembledContig;
        final int endInAssembledContig;
        final int assembledContigLength;
        final int mismatches;

        public AlignmentRegion(final String breakpointId, final String contigId, final AlnRgn alnRgn) {
            this.contigId = contigId;
            this.breakpointId = breakpointId;
            this.forwardStrand = alnRgn.getStrand() == '+';
            final Cigar alignmentCigar = TextCigarCodec.decode(alnRgn.getCigar());
            this.cigar = forwardStrand ? alignmentCigar : CigarUtils.invertCigar(alignmentCigar);
            this.referenceInterval = new SimpleInterval(alnRgn.getChrom(), (int) alnRgn.getPos() + 1, (int) (alnRgn.getPos() + 1 + cigar.getReferenceLength()));
            this.mqual = alnRgn.getMQual();
            this.assembledContigLength = cigar.getReadLength();
            this.startInAssembledContig = startOfAlignmentInContig(cigar);
            this.endInAssembledContig = endOfAlignmentInContig(assembledContigLength, cigar);
            this.mismatches = alnRgn.getNm();
        }

        public AlignmentRegion(final String breakpointId, final String contigId, final Cigar cigar, final boolean forwardStrand, final SimpleInterval referenceInterval, final int mqual, final int startInAssembledContig, final int endInAssembledContig, final int mismatches) {
            this.contigId = contigId;
            this.breakpointId = breakpointId;
            this.cigar = cigar;
            this.forwardStrand = forwardStrand;
            this.referenceInterval = referenceInterval;
            this.mqual = mqual;
            this.startInAssembledContig = startInAssembledContig;
            this.endInAssembledContig = endInAssembledContig;
            this.assembledContigLength = cigar.getReadLength();
            this.mismatches = mismatches;
        }

        public AlignmentRegion(final GATKRead read) {
            this.breakpointId = null;
            this.contigId = read.getName();
            this.forwardStrand = ! read.isReverseStrand();
            this.cigar = forwardStrand ? read.getCigar() : CigarUtils.invertCigar(read.getCigar());
            this.referenceInterval = new SimpleInterval(read);
            this.assembledContigLength = cigar.getReadLength() + getHardClipping(cigar);
            this.startInAssembledContig = startOfAlignmentInContig(cigar);
            this.endInAssembledContig = endOfAlignmentInContig(assembledContigLength, cigar);
            this.mqual = read.getMappingQuality();
            if (read.hasAttribute("NM")) {
                this.mismatches = read.getAttributeAsInteger("NM");
            } else {
                this.mismatches = 0;
            }
        }

        public int overlapOnContig(final AlignmentRegion other) {
            return Math.max(0, Math.min(endInAssembledContig, other.endInAssembledContig) - Math.max(startInAssembledContig, other.startInAssembledContig));
        }

        private static int getHardClipping(final Cigar cigar) {
            final List<CigarElement> cigarElements = cigar.getCigarElements();
            return (cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0) +
                    (cigarElements.get(cigarElements.size() - 1).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(cigarElements.size() - 1).getLength() : 0);
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
            while (ce.getOperator().isClipping()) {
                posInContig += ce.getLength();
                j += offset;
                ce = cigar.getCigarElement(j);
            }
            return posInContig;
        }

        @Override
        public String toString() {
            return breakpointId +
                    "\t" +
                    contigId +
                    "\t" +
                    referenceInterval.getContig() +
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
                    endInAssembledContig +
                    "\t" +
                    mismatches;
        }

        public static AlignmentRegion fromString(final String[] fields) {
            final String breakpointId = fields[0];
            final String contigId = fields[1];
            final String refContig = fields[2];
            final Integer refStart = Integer.valueOf(fields[3]);
            final Integer refEnd = Integer.valueOf(fields[4]);
            final SimpleInterval refInterval = new SimpleInterval(refContig, refStart, refEnd);
            final boolean refStrand = ("+".equals(fields[5]));
            final Cigar cigar = TextCigarCodec.decode(fields[6]);
            final int mqual = Integer.valueOf(fields[7]);
            final int contigStart = Integer.valueOf(fields[8]);
            final int contigEnd = Integer.valueOf(fields[9]);
            final int mismatches = Integer.valueOf(fields[10]);
            return new AlignmentRegion(breakpointId, contigId, cigar, refStrand, refInterval, mqual, contigStart, contigEnd, mismatches);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            AlignmentRegion that = (AlignmentRegion) o;
            return forwardStrand == that.forwardStrand &&
                    mqual == that.mqual &&
                    startInAssembledContig == that.startInAssembledContig &&
                    endInAssembledContig == that.endInAssembledContig &&
                    assembledContigLength == that.assembledContigLength &&
                    mismatches == that.mismatches &&
                    Objects.equals(contigId, that.contigId) &&
                    Objects.equals(breakpointId, that.breakpointId) &&
                    Objects.equals(cigar, that.cigar) &&
                    Objects.equals(referenceInterval, that.referenceInterval);
        }

        @Override
        public int hashCode() {
            return Objects.hash(contigId, breakpointId, cigar, forwardStrand, referenceInterval, mqual, startInAssembledContig, endInAssembledContig, assembledContigLength, mismatches);
        }

        public String toPackedString() {
            return "" + startInAssembledContig + "-" + endInAssembledContig + ":" + referenceInterval.getContig() + ',' + referenceInterval.getStart() + ',' + (forwardStrand ? '+' : '-') + ',' + TextCigarCodec.encode(cigar) + ',' + mqual + ',' + mismatches;
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
