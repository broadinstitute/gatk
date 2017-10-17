package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.apache.hadoop.yarn.webapp.hamlet.HamletSpec;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.sql.Ref;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by tsato on 10/11/17.
 */
public class RefSiteHistogram {
    // depth greater than or equal to this value will be grouped into the last column
    public final static int MAX_DEPTH = 200;

    // counts starts at 1 at index 0. e.g. for MAX_DEPTH = 200:
    // index: [0 1 2 ... 199]
    // depth: [1 2 3 ... 200+]
    private final int[] counts;
    private final String referenceContext;

    public RefSiteHistogram(final String referenceContext){
        this.referenceContext = referenceContext;
        this.counts = new int[MAX_DEPTH];
    }

    public RefSiteHistogram(final String referenceContext, final int[] counts) {
        this.referenceContext = referenceContext;
        this.counts = Arrays.copyOf(counts, counts.length);
    }

    public void increment(final int depth){
        // we don't count depth 0 sites so the 0th index will have depth 1
        Utils.validateArg(depth > 0, String.format("depth must be greater than 0 but got %d", depth));
        if (depth > MAX_DEPTH) {
            counts[MAX_DEPTH - 1] += 1;
            return;
        }

        counts[depth - 1] += 1;
    }

    public String getReferenceContext(){ return referenceContext; }

    public int[] getCounts() { return counts; }

    // writer
    public static class RefSiteHistogramWriter extends TableWriter<RefSiteHistogram> {
        public RefSiteHistogramWriter(final File output) throws IOException {
            super(output, RefSiteHistogram.RefSiteHistogramTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final RefSiteHistogram record, final DataLine dataLine) {
            dataLine.set(RefSiteHistogramTableColumn.CONTEXT.toString(), record.getReferenceContext())
                    .set(RefSiteHistogramTableColumn.COUNTS.toString(), Utils.intArrayToString(record.getCounts()));
        }
    }

    public static void writeRefSiteHistograms(final List<RefSiteHistogram> histograms, final File output){
        try (RefSiteHistogramWriter writer = new RefSiteHistogramWriter(output)) {
            writer.writeAllRecords(histograms);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing ref site histograms to %s.", output), e);
        }
    }

    // reader
    public static List<RefSiteHistogram> readRefSiteHistograms(final File table) {
        try (RefSiteHistogramReader reader = new RefSiteHistogramReader(table)) {
            return reader.toList();
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", table), e);
        }
    }

    private static class RefSiteHistogramReader extends TableReader<RefSiteHistogram> {
        private RefSiteHistogramReader(final File table) throws IOException {
            super(table);
        }

        @Override
        protected RefSiteHistogram createRecord(final DataLine dataLine) {
            final String referenceContext = dataLine.get(RefSiteHistogramTableColumn.CONTEXT);
            final int[] counts = Arrays.stream(dataLine.get(RefSiteHistogramTableColumn.COUNTS).split(","))
                    .mapToInt(Integer::parseInt).toArray();
            return new RefSiteHistogram(referenceContext, counts);
        }
    }

    enum RefSiteHistogramTableColumn {
        CONTEXT("context"),
        COUNTS("counts"); // {@code MAX_DEPTH dimensional array of counts}

        private String columnName;

        RefSiteHistogramTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
}
