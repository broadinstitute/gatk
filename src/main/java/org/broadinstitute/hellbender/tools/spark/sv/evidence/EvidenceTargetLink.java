package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervals;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Set;
import java.util.function.Function;

/**
 * This class holds information about pairs of intervals on the reference that are connected by one or more
 * BreakpointEvidence objects that have distal targets. Source and target are the stranded intervals pointed at by
 * the breakpoint evidence, with strand defined as whether the putative breakpoint exists upstream of the interval
 * start or downstream of the interval end, for + and -. The class also tracks the number of split reads and read pairs
 * that contributed to the link, as well as the template names of the reads.
 */
@DefaultSerializer(EvidenceTargetLink.Serializer.class)
public final class EvidenceTargetLink {
    private static final StrandedInterval.Serializer intervalSerializer = new StrandedInterval.Serializer();

    final StrandedInterval source;
    final StrandedInterval target;
    final int splitReads;
    final int readPairs;
    private final Set<String> readPairTemplateNames;
    private final Set<String> splitReadTemplateNames;
    private static final int NUM_BEDPE_TOKENS = 12;

    public EvidenceTargetLink(final StrandedInterval source, final StrandedInterval target,
                              final int splitReads, final int readPairs,
                              final Set<String> readPairTemplateNames, final Set<String> splitReadTemplateNames) {
        this.splitReadTemplateNames = splitReadTemplateNames;
        Utils.validateArg(source != null, "Can't construct EvidenceTargetLink with null source interval");
        if (source.getInterval().isUpstreamOf(target.getInterval())) {
            this.source = source;
            this.target = target;
        } else {
            this.source = target;
            this.target = source;
        }
        this.splitReads = splitReads;
        this.readPairs = readPairs;
        this.readPairTemplateNames = readPairTemplateNames;
    }

    public Set<String> getReadPairTemplateNames() {
        return readPairTemplateNames;
    }

    public Set<String> getSplitReadTemplateNames() {
        return splitReadTemplateNames;
    }

    @SuppressWarnings("unchecked")
    public EvidenceTargetLink(final Kryo kryo, final Input input) {
        this.source = intervalSerializer.read(kryo, input, StrandedInterval.class);
        this.target = intervalSerializer.read(kryo, input, StrandedInterval.class);

        this.splitReads = input.readInt();
        this.readPairs = input.readInt();

        this.readPairTemplateNames = (Set<String>) kryo.readClassAndObject(input);
        this.splitReadTemplateNames = (Set<String>) kryo.readClassAndObject(input);
    }

    protected void serialize(final Kryo kryo, final Output output ) {
        intervalSerializer.write(kryo, output, source);
        intervalSerializer.write(kryo, output, target);

        output.writeInt(splitReads);
        output.writeInt(readPairs);

        kryo.writeClassAndObject(output, readPairTemplateNames);
        kryo.writeClassAndObject(output, splitReadTemplateNames);
    }

    public String toBedpeString(ReadMetadata readMetadata) {
        final SVInterval sourceInterval = source.getInterval();
        final SVInterval targetInterval = target.getInterval();
        return readMetadata.getContigName(sourceInterval.getContig()) + "\t" + (sourceInterval.getStart() - 1) + "\t" + sourceInterval.getEnd() +
                "\t" + readMetadata.getContigName(targetInterval.getContig()) + "\t" + (targetInterval.getStart() - 1) + "\t" + targetInterval.getEnd() +
                "\t"  + getId(readMetadata) + "\t" +
                (readPairs + splitReads) + "\t" + (source.getStrand() ? "+" : "-") + "\t" + (target.getStrand() ? "+" : "-")
                + "\t" + "SR:" + Utils.join(",", splitReadTemplateNames) + "\t" + "RP:" + Utils.join(",", readPairTemplateNames);
    }

    public static EvidenceTargetLink fromBedpeString(final String str, final SAMSequenceDictionary dictionary) {
        return fromBedpeString(str, dictionary::getSequenceIndex);
    }
    public static EvidenceTargetLink fromBedpeString(final String str, final ReadMetadata readMetadata) {
        return fromBedpeString(str, readMetadata::getContigID);
    }

    public static EvidenceTargetLink fromBedpeString(final String str, Function<String,Integer> sequenceNameToIndexFunction) {
        final String[] tokens = str.split("\t");
        if (tokens.length != NUM_BEDPE_TOKENS) {
            throw new IllegalArgumentException("Could not create " + EvidenceTargetLink.class.getSimpleName() + " because " + NUM_BEDPE_TOKENS + " tab-delimited tokens were expected but found " + tokens.length + " in the bedpe string: " + str);
        }
        final int sourceContig = sequenceNameToIndexFunction.apply(tokens[0]);
        final int sourceStart = parseInteger(tokens[1]);
        final int sourceEnd = parseInteger(tokens[2]);
        final int targetContig = sequenceNameToIndexFunction.apply(tokens[3]);
        final int targetStart = parseInteger(tokens[4]);
        final int targetEnd = parseInteger(tokens[5]);
        final int numReadPairsAndSplitReads = parseInteger(tokens[7]);
        final boolean sourceStrand = parseStrand(tokens[8]);
        final boolean targetStrand = parseStrand(tokens[9]);
        final Set<String> splitReadNames = parseTemplateNames(tokens[10]);
        final Set<String> readPairNames = parseTemplateNames(tokens[11]);

        if (numReadPairsAndSplitReads != splitReadNames.size() + readPairNames.size()) {
            throw new IllegalArgumentException("Sum of split read (" + splitReadNames.size() + ") and read pair (" + readPairNames.size() + ") template names does not equal the listed sum " + numReadPairsAndSplitReads + " in record: " + str);
        }

        final StrandedInterval source = new StrandedInterval(new SVInterval(sourceContig, sourceStart, sourceEnd), sourceStrand);
        final StrandedInterval target = new StrandedInterval(new SVInterval(targetContig, targetStart, targetEnd), targetStrand);
        return new EvidenceTargetLink(source, target, splitReadNames.size(), readPairNames.size(), readPairNames, splitReadNames);

    }

    private static int parseInteger(final String str) {
        try {
            return Integer.parseInt(str);
        } catch (final NumberFormatException e) {
            throw new IllegalArgumentException("Could not parse string as integer: " + str);
        }
    }

    private static boolean parseStrand(final String str) {
        if (str.equals("+")) {
            return true;
        } else if (str.equals("-")) {
            return false;
        } else {
            throw new IllegalArgumentException("Unrecognized strand token:" + str);
        }
    }

    private static Set<String> parseTemplateNames(final String str) {
        final String trimmedStr = str.substring(3); //remove "RP:" or "SR:" prefix
        if (trimmedStr.isEmpty()) return Collections.emptySet();
        final String[] tokens = trimmedStr.split(",");
        return Sets.newHashSet(tokens);
    }

    private String getId(final ReadMetadata readMetadata) {
        final SVInterval sourceInterval = source.getInterval();
        final SVInterval targetInterval = target.getInterval();

        return readMetadata.getContigName(sourceInterval.getContig()) + "_" + (sourceInterval.getStart() - 1) + "_" + sourceInterval.getEnd() +
                "_" + readMetadata.getContigName(targetInterval.getContig()) + "_" + (targetInterval.getStart() - 1) + "_" + targetInterval.getEnd() +
                "_" + (source.getStrand() ? "P" : "M")  + (target.getStrand() ? "P" : "M") + "_" + splitReads + "_" + readPairs;
    }

    public int getSplitReads() {
        return splitReads;
    }

    public int getReadPairs() {
        return readPairs;
    }

    public PairedStrandedIntervals getPairedStrandedIntervals() {
        return new PairedStrandedIntervals(source, target);
    }

    public boolean isImpreciseDeletion() {
        return getPairedStrandedIntervals().getLeft().getInterval().getContig() == getPairedStrandedIntervals().getRight().getInterval().getContig()
                && getPairedStrandedIntervals().getLeft().getStrand() && !getPairedStrandedIntervals().getRight().getStrand();
    }

    /**
     * Distance between the two intervals. Defined to be -1 if the intervals are on different contigs. Otherwise, returns the
     * inner distance: the distance between the right boundary of the left interval and the left boundary of the right interval
     */
    public int getDistance() {
        if (getPairedStrandedIntervals().getLeft().getInterval().getContig() != getPairedStrandedIntervals().getRight().getInterval().getContig()) {
            return -1;
        } else {
            return getPairedStrandedIntervals().getRight().getInterval().getStart() - getPairedStrandedIntervals().getLeft().getInterval().getEnd();
        }

    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<EvidenceTargetLink> {
        @Override
        public void write( final Kryo kryo, final Output output, final EvidenceTargetLink evidence ) {
            evidence.serialize(kryo, output);
        }

        @Override
        public EvidenceTargetLink read(final Kryo kryo, final Input input, final Class<EvidenceTargetLink> klass ) {
            return new EvidenceTargetLink(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final EvidenceTargetLink link = (EvidenceTargetLink) o;

        if (splitReads != link.splitReads) return false;
        if (readPairs != link.readPairs) return false;
        if (!source.equals(link.source)) return false;
        return target.equals(link.target);
    }

    @Override
    public int hashCode() {
        int result = source.hashCode();
        result = 31 * result + target.hashCode();
        result = 31 * result + splitReads;
        result = 31 * result + readPairs;
        return result;
    }

    @Override
    public String toString() {
        return "EvidenceTargetLink{" +
                "source=" + source +
                ", target=" + target +
                ", splitReads=" + splitReads +
                ", readPairs=" + readPairs +
                '}';
    }
}
