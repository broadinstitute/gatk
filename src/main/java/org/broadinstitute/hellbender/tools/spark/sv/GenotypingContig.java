package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by valentin on 8/25/17.
 */
class GenotypingContig {

    class HaplotypeAlignment {
        public final int score;
        public final List<AlignmentInterval> intervals;
        public final int numberOfMatches;
        public final int numberOfMismatches;
        public final int numberOfIndels;
        public final int numberOfReversals;
        public final int numberOfIndelBases;

        public HaplotypeAlignment(final String scoreString, final String alingmentString) {
            final String[] parts = scoreString.split(",");
            int nextIndex = 0;
            score = Integer.parseInt(parts[nextIndex++]);
            numberOfMatches = Integer.parseInt(parts[nextIndex++]);
            numberOfMismatches = Integer.parseInt(parts[nextIndex++]);
            numberOfIndels = Integer.parseInt(parts[nextIndex++]);
            numberOfIndelBases = Integer.parseInt(parts[nextIndex++]);
            numberOfReversals = Integer.parseInt(parts[nextIndex]);
            intervals = parseAlignmentString(alingmentString);
        }

        private List<AlignmentInterval> parseAlignmentString(final String alingmentString) {
            final String[] alignments = alingmentString.replaceAll(";$","").split(";");
            // TODO AI(null) -> AI(s);
            return Stream.of(alignments).map(AlignmentInterval::new).collect(Collectors.toList());
        }
    }

    private final byte[] bases;

    private final HaplotypeAlignment[] alignments;

    public GenotypingContig(final GATKRead read) {
        assertRecordIsValid(read);
        this.bases = read.getBases();
        final String refScore = getMantarodyAttributeAsString(read, "RS");
        final String altScore = getMantarodyAttributeAsString(read, "XS");
        final String refAlign = getMantarodyAttributeAsString(read, "RA");
        final String altAlign = getMantarodyAttributeAsString(read, "AA");
        alignments = new HaplotypeAlignment[] { new HaplotypeAlignment(refScore, refAlign), new HaplotypeAlignment(altScore, altAlign)};
    }

    public HaplotypeAlignment getHaplotypeAlignment(final int index) {
        ParamUtils.inRange(index, 0, 1, "input index is out of range");
        return alignments[index];
    }

    public byte[] getBases() {
        return bases.clone();
    }

    private static GATKRead assertRecordIsValid(final GATKRead read) {
        Utils.nonNull(read);
        if (!read.isUnmapped() || read.isReverseStrand() || !read.getCigar().containsOperator(CigarOperator.H)) {
            throw new IllegalArgumentException("the input read does is not a genotyping-contig record");
        }
        return read;
    }

    private static String getMantarodyAttributeAsString(final GATKRead read, final String name) {
        final String result = read.getAttributeAsString(name);
        if (result == null) {
            throw new GATKException.MissingReadField(name, read);
        }
        return result;
    }
}
