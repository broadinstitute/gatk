package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.google.common.annotations.VisibleForTesting;
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

                    assembledBreakpoints.addAll(getAssembledBreakpointsFromAlignmentRegions(contigId, sequence, alignmentRegionList));
                }
            }
        } catch (final IOException e) {
            throw new GATKException("could not execute BWA");
        }

        return assembledBreakpoints;
    }

    @VisibleForTesting
    public static List<AssembledBreakpoint> getAssembledBreakpointsFromAlignmentRegions(final String contigId, final byte[] sequence, final List<AlignmentRegion> alignmentRegionList) {
        final List<AssembledBreakpoint> results = new ArrayList<>(alignmentRegionList.size() / 2);
        final Iterator<AlignmentRegion> iterator = alignmentRegionList.iterator();
        if ( iterator.hasNext() ) {
            AlignmentRegion current = iterator.next();
            while ( iterator.hasNext() ) {
                final AlignmentRegion next = iterator.next();
                if (next.mqual < 60 && iterator.hasNext()) {
                    continue;
                }

                final AlignmentRegion previous = current;
                current = next;

                String homology = "NA";
                if (previous.endInAssembledContig >= current.startInAssembledContig) {
                    homology = new String(Arrays.copyOfRange(sequence, current.startInAssembledContig - 1, previous.endInAssembledContig));
                }

                String insertedSequence = "NA";
                if (previous.endInAssembledContig < current.startInAssembledContig - 1) {
                    insertedSequence = new String(Arrays.copyOfRange(sequence, previous.endInAssembledContig, current.startInAssembledContig));
                }
                final AssembledBreakpoint assembledBreakpoint = new AssembledBreakpoint(contigId, previous, current, insertedSequence, homology);
                results.add(assembledBreakpoint);
            }
        }
        return results;
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

        public AssembledBreakpoint(final String contigId, final AlignmentRegion region1, final AlignmentRegion region2, final String insertedSequence, final String homology) {
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

        /**
         *  Parses a tab-delimited assembled breakpoint line into an AssembledBreakpoint object
         */
        public static AssembledBreakpoint fromString(String assembledBreakpointLine) {
            final String[] fields = assembledBreakpointLine.split("\t");
            return fromFields(fields);
        }

        public static AssembledBreakpoint fromFields(final String[] fields) {
            try {
                final String contigId = fields[0].replaceFirst("^>","");
                final String[] alignmentRegion1Fields = Arrays.copyOfRange(fields, 1, 9);
                final AlignmentRegion alignmentRegion1 = AlignmentRegion.fromString(alignmentRegion1Fields);
                final String[] alignmentRegion2Fields = Arrays.copyOfRange(fields, 9, 17);
                final AlignmentRegion alignmentRegion2 = AlignmentRegion.fromString(alignmentRegion2Fields);
                final String insertedSequence = fields[17].equals("NA") ? "" : fields[17];
                final String homology = fields[18].equals("NA") ? "" : fields[18];
                return new AssembledBreakpoint(contigId, alignmentRegion1, alignmentRegion2, insertedSequence, homology);
            } catch (final NumberFormatException nfe) {
                throw new GATKException(Arrays.toString(fields), nfe);
            }

        }

        public SimpleInterval getLeftAlignedLeftBreakpointOnAssembledContig() {
            final int alignmentEnd = region1.forwardStrand ? region1.referenceInterval.getEnd() : region1.referenceInterval.getStart();
            final int position = region1.forwardStrand ? alignmentEnd - homology.length() : alignmentEnd + homology.length();
            return new SimpleInterval(region1.referenceInterval.getContig(), position, position);
        }

        public SimpleInterval getLeftAlignedRightBreakpointOnAssembledContig() {
            final int position = region2.forwardStrand ? region2.referenceInterval.getStart() : region2.referenceInterval.getEnd();
            return new SimpleInterval(region2.referenceInterval.getContig(), position, position);
        }

        public BreakpointAllele getBreakpointAllele() {
            final SimpleInterval leftAlignedLeftBreakpointOnAssembledContig = getLeftAlignedLeftBreakpointOnAssembledContig();
            final SimpleInterval leftAlignedRightBreakpointOnAssembledContig = getLeftAlignedRightBreakpointOnAssembledContig();
            if (! leftAlignedLeftBreakpointOnAssembledContig.getContig().equals(leftAlignedRightBreakpointOnAssembledContig.getContig())) {
                return new BreakpointAllele(leftAlignedLeftBreakpointOnAssembledContig, leftAlignedRightBreakpointOnAssembledContig, insertedSequence, region1.forwardStrand, region2.forwardStrand);
            } else if ( leftAlignedLeftBreakpointOnAssembledContig.getStart() < leftAlignedRightBreakpointOnAssembledContig.getStart()) {
                return new BreakpointAllele(leftAlignedLeftBreakpointOnAssembledContig, leftAlignedRightBreakpointOnAssembledContig, insertedSequence, region1.forwardStrand, region2.forwardStrand);
            } else {
                return new BreakpointAllele(leftAlignedRightBreakpointOnAssembledContig, leftAlignedLeftBreakpointOnAssembledContig, insertedSequence, region1.forwardStrand, region2.forwardStrand);
            }
        }
    }

    static class BreakpointAllele {
        SimpleInterval leftAlignedLeftBreakpoint;
        SimpleInterval leftAlignedRightBreakpoint;
        String insertedSequence;
        boolean left5Prime;
        boolean right5Prime;

        public BreakpointAllele(final SimpleInterval leftAlignedLeftBreakpoint, final SimpleInterval leftAlignedRightBreakpoint, final String insertedSequence, final boolean left5Prime, final boolean right5Prime) {
            this.leftAlignedLeftBreakpoint = leftAlignedLeftBreakpoint;
            this.leftAlignedRightBreakpoint = leftAlignedRightBreakpoint;
            this.insertedSequence = insertedSequence;
            this.left5Prime = left5Prime;
            this.right5Prime = right5Prime;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            final BreakpointAllele that = (BreakpointAllele) o;
            return left5Prime == that.left5Prime &&
                    right5Prime == that.right5Prime &&
                    Objects.equals(leftAlignedLeftBreakpoint, that.leftAlignedLeftBreakpoint) &&
                    Objects.equals(leftAlignedRightBreakpoint, that.leftAlignedRightBreakpoint) &&
                    Objects.equals(insertedSequence, that.insertedSequence);
        }

        @Override
        public int hashCode() {
            return Objects.hash(leftAlignedLeftBreakpoint, leftAlignedRightBreakpoint, insertedSequence, left5Prime, right5Prime);
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

        public AlignmentRegion(final Cigar cigar, final boolean forwardStrand, final SimpleInterval referenceInterval, final int mqual, final int startInAssembledContig, final int endInAssembledContig) {
            this.cigar = cigar;
            this.forwardStrand = forwardStrand;
            this.referenceInterval = referenceInterval;
            this.mqual = mqual;
            this.startInAssembledContig = startInAssembledContig;
            this.endInAssembledContig = endInAssembledContig;
            this.assembledContigLength = cigar.getReadLength();
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

        public static AlignmentRegion fromString(final String[] fields) {
            final String alignmentRegion1RefContig = fields[0];
            final Integer alignmentRegion1RefStart = Integer.valueOf(fields[1]);
            final Integer alignmentRegion1RefEnd = Integer.valueOf(fields[2]);
            final SimpleInterval alignmentRegion1Interval = new SimpleInterval(alignmentRegion1RefContig, alignmentRegion1RefStart, alignmentRegion1RefEnd);
            final boolean alignmentRegion1ForwardStrand = ("+".equals(fields[3]));
            final Cigar alignmentRegion1Cigar = TextCigarCodec.decode(fields[4]);
            final int alignmentRegion1MQual = Integer.valueOf(fields[5]);
            final int alignmentRegion1ContigStart = Integer.valueOf(fields[6]);
            final int alignmentRegion1ContigEnd = Integer.valueOf(fields[7]);
            return new AlignmentRegion(alignmentRegion1Cigar, alignmentRegion1ForwardStrand, alignmentRegion1Interval, alignmentRegion1MQual, alignmentRegion1ContigStart, alignmentRegion1ContigEnd);
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
