package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Gathers paired-end and split read evidence files",
        oneLineSummary = "Gathers paired-end and split read evidence files",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class PairedEndAndSplitReadEvidenceCollection extends ReadWalker {

    @Argument(shortName = "p", fullName = "pe-file", doc = "Output file for paired end evidence", optional=false)
    public String peFile;

    @Argument(shortName = "s", fullName = "sr-file", doc = "Output file for split read evidence", optional=false)
    public String srFile;

    @Argument(fullName = "sample-name", doc = "Sample name")
    String sampleName = null;

    @Argument(fullName = "max-split-dist", doc = "Maxiumum split distance", optional = true)
    Integer maxSplitDist = 300;


    Set<String> observedDiscordantNames = new HashSet<>();
    int currentDiscordantPosition = -1;
    int prevSplitPos = -1;
    String currentChrom = null;
    private OutputStreamWriter peWriter;
    private OutputStreamWriter srWriter;

    SortedMap<SplitPos, Integer> splitCounts;
    SortedSet<DiscordantRead> discordantPairs;

    @Override
    public boolean requiresReads() {
        return true;
    }


    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        peWriter = new OutputStreamWriter(new BlockCompressedOutputStream(peFile));
        srWriter = new OutputStreamWriter(new BlockCompressedOutputStream(srFile));

        final Comparator<SplitPos> comparator = (o1, o2) -> {
            if (o1.pos != o2.pos) {
                return new Integer(o1.pos).compareTo(o2.pos);
            } else {
                return o1.direction.compareTo(o2.direction);
            }
        };

        splitCounts = new TreeMap<>(comparator);

        final Comparator<DiscordantRead> discReadComparator =
                Comparator.comparing((DiscordantRead r) -> getBestAvailableSequenceDictionary().getSequenceIndex(r.getContig()))
                        .thenComparing(DiscordantRead::getStart)
                        .thenComparing((DiscordantRead r) -> getBestAvailableSequenceDictionary().getSequenceIndex(r.getMateContig()))
                        .thenComparing(DiscordantRead::getMateStart)
                        .thenComparing(DiscordantRead::getName);

        discordantPairs = new TreeSet<>(discReadComparator);

    }

    public boolean isExcluded(final GATKRead read) {
        return read.isUnmapped() || read.mateIsUnmapped() || read.isSecondaryAlignment() || read.isDuplicate() || read.isSupplementaryAlignment();
    }

    static class DiscordantRead {
        boolean readReverseStrand;
        boolean mateReverseStrand;
        String contig;
        int start;
        String mateContig;
        int mateStart;
        String name;

        public DiscordantRead(final GATKRead read) {
            this.readReverseStrand = read.isReverseStrand();
            this.mateReverseStrand = read.mateIsReverseStrand();
            this.contig = read.getContig();
            this.start = read.getStart();
            this.mateContig = read.getMateContig();
            this.mateStart = read.getMateStart();
            this.name = read.getName();
        }

        public boolean isReadReverseStrand() {
            return readReverseStrand;
        }

        public void setReadReverseStrand(final boolean readReverseStrand) {
            this.readReverseStrand = readReverseStrand;
        }

        public boolean isMateReverseStrand() {
            return mateReverseStrand;
        }

        public void setMateReverseStrand(final boolean mateReverseStrand) {
            this.mateReverseStrand = mateReverseStrand;
        }

        public String getContig() {
            return contig;
        }

        public void setContig(final String contig) {
            this.contig = contig;
        }

        public int getStart() {
            return start;
        }

        public void setStart(final int start) {
            this.start = start;
        }

        public String getMateContig() {
            return mateContig;
        }

        public void setMateContig(final String mateContig) {
            this.mateContig = mateContig;
        }

        public int getMateStart() {
            return mateStart;
        }

        public void setMateStart(final int mateStart) {
            this.mateStart = mateStart;
        }

        public String getName() {
            return name;
        }

        public void setName(final String name) {
            this.name = name;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final DiscordantRead that = (DiscordantRead) o;

            if (readReverseStrand != that.readReverseStrand) return false;
            if (mateReverseStrand != that.mateReverseStrand) return false;
            if (start != that.start) return false;
            if (mateStart != that.mateStart) return false;
            if (contig != null ? !contig.equals(that.contig) : that.contig != null) return false;
            if (mateContig != null ? !mateContig.equals(that.mateContig) : that.mateContig != null) return false;
            return name != null ? name.equals(that.name) : that.name == null;
        }

        @Override
        public int hashCode() {
            int result = (readReverseStrand ? 1 : 0);
            result = 31 * result + (mateReverseStrand ? 1 : 0);
            result = 31 * result + (contig != null ? contig.hashCode() : 0);
            result = 31 * result + start;
            result = 31 * result + (mateContig != null ? mateContig.hashCode() : 0);
            result = 31 * result + mateStart;
            result = 31 * result + (name != null ? name.hashCode() : 0);
            return result;
        }
    }

    static class SplitPos {
        POSITION direction;
        int pos;

        public SplitPos(final int start, final POSITION direction) {
            this.pos = start;
            this.direction = direction;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final SplitPos splitPos = (SplitPos) o;

            if (pos != splitPos.pos) return false;
            return direction.ordinal() == splitPos.direction.ordinal();
        }

        @Override
        public int hashCode() {
            int result = direction != null ? direction.ordinal() : 0;
            result = 31 * result + pos;
            return result;
        }
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (isExcluded(read)) {
            return;
        }

        if (isSoftClipped(read)) {
            countSplitRead(read);
        }

        if (! read.isProperlyPaired()) {
            reportDiscordantReadPair(read);
        }
    }

    private void reportDiscordantReadPair(final GATKRead read) {
        if (read.getStart() != currentDiscordantPosition) {
            flushDiscordantReadPairs();
            currentDiscordantPosition = read.getStart();
        }

        final int readSeqId = getBestAvailableSequenceDictionary().getSequenceIndex(read.getContig());
        final int mateSeqId = getBestAvailableSequenceDictionary().getSequenceIndex(read.getMateContig());
        if (readSeqId < mateSeqId) {
            discordantPairs.add(new DiscordantRead(read));
        } else if (readSeqId == mateSeqId) {
            if (read.getStart() < read.getMateStart()) {
                discordantPairs.add(new DiscordantRead(read));
            } else if (read.getStart() == read.getMateStart()) {

                final boolean seenBefore = observedDiscordantNames.remove(read.getName());
                if (! seenBefore) {
                    discordantPairs.add(new DiscordantRead(read));
                    observedDiscordantNames.add(read.getName());
                }
            }
        }
    }

    private void flushDiscordantReadPairs() {
        discordantPairs.forEach(this::writeDiscordantPair);
        discordantPairs.clear();
    }

    private void writeDiscordantPair(final DiscordantRead r) {
        final String strandA = r.isReadReverseStrand() ? "-" : "+";
        final String strandB = r.isMateReverseStrand() ? "-" : "+";

        try {
            // subtract 1 from positions to match pysam output
            peWriter.write(r.getContig() + "\t" + (r.getStart() - 1) + "\t" + strandA + "\t" + r.getMateContig() + "\t" + (r.getMateStart() - 1) + "\t" + strandB + "\t" + sampleName + "\n");
        } catch (IOException e) {
            throw new GATKException("Could not write to PE file", e);
        }
    }

    private void countSplitRead(final GATKRead read) {
        final SplitPos splitPosition = getSplitPosition(read);
        if (splitPosition.direction == POSITION.MIDDLE) {
            return;
        }
        final int dist;
        if (prevSplitPos == -1) {
            dist = 0;
        } else {
            dist = Math.abs(splitPosition.pos - prevSplitPos);
        }
        prevSplitPos = splitPosition.pos;
        if (currentChrom == null) {
            currentChrom = read.getContig();
        }
        if (dist > maxSplitDist) {
            flushSplitCounts();
            currentChrom = read.getContig();
        }

        if (! splitCounts.containsKey(splitPosition)) {
            splitCounts.put(splitPosition, 1);
        } else {
            splitCounts.put(splitPosition, splitCounts.get(splitPosition) + 1);
        }
    }

    private void flushSplitCounts() {
        splitCounts.forEach((position, count) -> {
            try {
                // subtract one from pos to match pysam results
                srWriter.write(currentChrom + "\t" + (position.pos - 1) + "\t" + position.direction.getDescription() + "\t" + count + "\t" + sampleName + "\n");
            } catch (IOException e) {
                throw new GATKException("Could not write to sr file", e);
            }
        });
        splitCounts.clear();
    }

    private SplitPos getSplitPosition(GATKRead read) {
        if (read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.M) {
            final int matchLength = read.getCigar().getCigarElements().stream().filter(e -> e.getOperator() == CigarOperator.M).mapToInt(CigarElement::getLength).sum();
            return new SplitPos(read.getStart() + matchLength, POSITION.RIGHT);
        } else if (read.getCigar().getLastCigarElement().getOperator() == CigarOperator.M) {
            return new SplitPos(read.getStart(), POSITION.LEFT);
        }

        return new SplitPos(-1, POSITION.MIDDLE);
    }

    private boolean isSoftClipped(final GATKRead read) {
        final CigarOperator firstOperator = read.getCigar().getFirstCigarElement().getOperator();
        final CigarOperator lastOperator = read.getCigar().getLastCigarElement().getOperator();
        return (firstOperator == CigarOperator.SOFT_CLIP && lastOperator != CigarOperator.SOFT_CLIP) ||
                (firstOperator != CigarOperator.SOFT_CLIP && lastOperator == CigarOperator.SOFT_CLIP);
    }

    @Override
    public Object onTraversalSuccess() {
        flushSplitCounts();
        return null;
    }

    @Override
    public void closeTool() {
        super.closeTool();
        try {
            peWriter.close();
            srWriter.close();
        } catch (IOException e) {
            throw new GATKException("error closing output file", e);
        }
    }

    enum POSITION {
        LEFT ("left"),
        MIDDLE ("middle"),
        RIGHT ("right");

        private String description;

        POSITION(final String description) {
            this.description = description;
        }

        public String getDescription() {
            return description;
        }
    }
}
