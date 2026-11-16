package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * The tables of histograms of the F1R2 counts of a sample, which {@link F1R2CountsCollector} writes: of sites without alt
 * reads, <sample>.ref_histogram, and of sites with one alt read, <sample>.alt_histogram.  Each row has the sites of one
 * block of the genome (see {@link F1R2FilterConstants#BLOCK_SIZE}) and reference context, and for one alt read its base
 * and orientation, and there is a column for the number of these sites in each depth bin:
 *
 * <pre>
 * #&lt;METADATA&gt;SAMPLE=tumor
 * contig  block_start  context  1  2  ...  24   30    40     ...  200
 * chr1    1            AAA      0  1  ...  134  7050  20369  ...  0
 * </pre>
 *
 * with columns alt and orientation after context in an alt histogram table.  The rows are in order of block, and only
 * rows with sites are written.
 */
public final class F1R2HistogramTable {
    public static final String CONTIG_COLUMN = "contig";
    public static final String BLOCK_START_COLUMN = "block_start";
    public static final String CONTEXT_COLUMN = "context";
    public static final String ALT_COLUMN = "alt";
    public static final String ORIENTATION_COLUMN = "orientation";

    private static final List<String> REF_KEY_COLUMNS = List.of(CONTIG_COLUMN, BLOCK_START_COLUMN, CONTEXT_COLUMN);
    private static final List<String> ALT_KEY_COLUMNS = List.of(CONTIG_COLUMN, BLOCK_START_COLUMN, CONTEXT_COLUMN, ALT_COLUMN, ORIENTATION_COLUMN);

    private F1R2HistogramTable() { }

    /**
     * The sites of one row: the block of the genome, by its contig and 1-based start; the reference context; for sites
     * with one alt read, its base and orientation, both null for sites without alt reads; and the number of sites of
     * each depth
     */
    public record Row(String contig, int blockStart, String context, Nucleotide altBase, ReadOrientation orientation,
                      int[] depths, int[] counts) {
        // the label of the histogram of these sites in LearnReadOrientationModel, eg ACG, or ACG_T_F1R2 for sites with one alt read
        public String histogramLabel() {
            return altBase == null ? context : F1R2FilterUtils.tripletToLabel(context, altBase, orientation);
        }
    }

    /**
     * A writer of a table of a sample's sites without alt reads, or with one alt read if alt is true, with a column for
     * each of the given depths
     */
    public static TableWriter<Row> writer(final Path path, final String sample, final boolean alt, final int[] depths) throws IOException {
        final List<String> columns = new ArrayList<>(alt ? ALT_KEY_COLUMNS : REF_KEY_COLUMNS);
        Arrays.stream(depths).forEach(depth -> columns.add(Integer.toString(depth)));
        final TableWriter<Row> writer = TableUtils.writer(path, new TableColumnCollection(columns), (row, dataLine) -> {
            dataLine.append(row.contig()).append(row.blockStart()).append(row.context());
            if (alt) {
                dataLine.append(row.altBase().toString()).append(row.orientation().toString());
            }
            dataLine.append(row.counts());
        });
        writer.writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
        return writer;
    }

    /**
     * A reader of a ref or alt histogram table, whose metadata name its sample
     */
    public static TableReader<Row> reader(final Path path) throws IOException {
        return TableUtils.reader(path, (columns, formatException) -> {
            final boolean alt = columns.contains(ALT_COLUMN);
            final List<String> keyColumns = alt ? ALT_KEY_COLUMNS : REF_KEY_COLUMNS;
            if (!columns.matchesAll(0, keyColumns.toArray(new String[0]))) {
                throw formatException.apply("a histogram table of F1R2 counts must begin with the columns " + keyColumns);
            }
            final int[] depths;
            try {
                depths = columns.names().subList(keyColumns.size(), columns.columnCount()).stream().mapToInt(Integer::parseInt).toArray();
            } catch (final NumberFormatException e) {
                throw formatException.apply("the columns of a histogram table of F1R2 counts after " + keyColumns + " must be depths");
            }
            if (depths.length == 0 || depths[0] < 1 || IntStream.range(1, depths.length).anyMatch(n -> depths[n] <= depths[n - 1])) {
                throw formatException.apply("the depths of a histogram table of F1R2 counts must be positive and increasing");
            }
            return dataLine -> new Row(dataLine.get(CONTIG_COLUMN), dataLine.getInt(BLOCK_START_COLUMN), dataLine.get(CONTEXT_COLUMN),
                    alt ? Nucleotide.valueOf(dataLine.get(ALT_COLUMN)) : null,
                    alt ? ReadOrientation.valueOf(dataLine.get(ORIENTATION_COLUMN)) : null,
                    depths, IntStream.range(keyColumns.size(), columns.columnCount()).map(dataLine::getInt).toArray());
        });
    }
}
