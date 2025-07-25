package org.broadinstitute.hellbender.tools.sv;


import com.google.common.collect.ImmutableSet;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;


/**
 * TODO docs
 */

public class SelectSVPairs {

    // SV pair table column names
    public static final String VID_A_COLUMN = "VID_A";
    public static final String VID_B_COLUMN = "VID_B";
    public static final String SCORE_COLUMN = "SCORE";
    protected static final Set<String> COLUMN_NAMES = ImmutableSet.of(VID_A_COLUMN, VID_B_COLUMN, SCORE_COLUMN);

    public HashMap<String, String> vidAtoB = new HashMap<>();
    public HashMap<String, String> vidBtoA = new HashMap<>();

    public SelectSVPairs(final GATKPath svPairFilePath) {
        final List<SVPair> sortedPairs = loadSVPairs(svPairFilePath, this);
        final HashSet<String> visitedA = new HashSet<>();
        final HashSet<String> visitedB = new HashSet<>();

        for (SVPair pair : sortedPairs) {
            if (!visitedA.contains(pair.getVidA()) && !visitedB.contains(pair.getVidB())) {
                // keep pair if neither SV is already in another pair with a higher score
                vidAtoB.put(pair.getVidA(), pair.getVidB());
                vidBtoA.put(pair.getVidB(), pair.getVidA());
                visitedA.add(pair.getVidA());
                visitedB.add(pair.getVidB());
            }
        }
    }

    public HashMap<String, String> getVidAToBMap() { return vidAtoB; }
    public HashMap<String, String> getVidBToAMap() { return vidBtoA; }

    protected List<SVPair> loadSVPairs(final GATKPath svPairFilePath, final SelectSVPairs selector) {
        Utils.nonNull(svPairFilePath);
        final HashSet<SVPair> pairSet = new HashSet<>();  // use HashSet to remove duplicates
        try (final TableReader<SVPair> tableReader = TableUtils.reader(svPairFilePath.toPath(), selector::tableParser)) {
            for (final SVPair svPair : tableReader) {
                pairSet.add(svPair);
            }
        } catch (final IOException e) {
            throw new GATKException("IO error while reading SV pair table", e);
        }
        // convert to list and sort in descending order by score
        // don't implement as compareTo() because inconsistent with equals()
        return pairSet.stream()
                .sorted((a, b) -> b.getScore().compareTo(a.getScore()))
                .collect(Collectors.toList());
    }

    protected Function<DataLine,SVPair> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
        // Check for expected columns
        for (final String column : COLUMN_NAMES) {
            if (!columns.contains(column)) {
                throw exceptionFactory.apply("Missing column " + column);
            }
        }
        // Check there are no extra columns
        if (columns.columnCount() != COLUMN_NAMES.size()) {
            throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
        }
        return this::parseTableLine;
    }

    protected SVPair parseTableLine(final DataLine dataLine) {
        final Double score = Double.parseDouble(dataLine.get(SCORE_COLUMN));
        return new SVPair(dataLine.get(VID_A_COLUMN), dataLine.get(VID_B_COLUMN), score);
    }

    public static class SVPair {
        final String vidA;
        final String vidB;
        final Double score;

        public SVPair(final String vidA, final String vidB, final Double score) {
            this.vidA = Utils.nonNull(vidA);
            this.vidB = Utils.nonNull(vidB);
            this.score = Utils.nonNull(score);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            final SVPair that = (SVPair) o;
            return getVidA().equals(that.getVidA()) && getVidB().equals(that.getVidB()) &&
                    getScore().equals(that.getScore());
        }

        @Override
        public int hashCode() {
            return Objects.hash(vidA, vidB, score);
        }

        public String getVidA() { return vidA; }
        public String getVidB() { return vidB; }
        public Double getScore() {return score; }
    }

}
