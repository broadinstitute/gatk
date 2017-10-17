package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by tsato on 10/11/17.
 */
public class AltSiteRecord {
    private String referenceContext;
    private int[] baseCounts;
    private int[] f1r2Counts;
    private int depth;
    private Nucleotide altAllele;


    public AltSiteRecord(final String referenceContext, final int[] baseCounts, final int[] f1r2Counts, final int depth, final Nucleotide altAllele){
        this.referenceContext = referenceContext;
        this.baseCounts = baseCounts;
        this.f1r2Counts = f1r2Counts;
        this.depth = depth;
        this.altAllele = altAllele;
    }

    public String getReferenceContext(){ return referenceContext; }

    public int[] getBaseCounts(){ return baseCounts; }

    public int[] getF1R2Counts(){ return f1r2Counts; }

    public int getDepth(){ return depth; }

    public Nucleotide getAltAllele(){ return altAllele; }


    /*** Table IO code ***/
    private enum AltSiteRecordTableColumn {
        CONTEXT("context"),
        BASE_COUNTS("base-counts"),
        F1R2_COUNTS("f1r2-counts"),
        DEPTH("depth"),
        ALT_ALLELE("alt");

        private String columnName;

        AltSiteRecordTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }


    /*** Writer ***/
    public static class AltSiteRecordTableWriter extends TableWriter<AltSiteRecord> {
        public AltSiteRecordTableWriter(final File output) throws IOException {
            super(output, AltSiteRecord.AltSiteRecordTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final AltSiteRecord record, final DataLine dataLine) {
            // it'd be nice to set() less manually...
            // Note that allele fraction f is not allele-specific, thus the same f array will be printed
            // four times for each context
            dataLine.set(AltSiteRecordTableColumn.CONTEXT.toString(), record.getReferenceContext())
                    .set(AltSiteRecordTableColumn.BASE_COUNTS.toString(), Utils.intArrayToString(record.getBaseCounts()))
                    .set(AltSiteRecordTableColumn.F1R2_COUNTS.toString(), Utils.intArrayToString(record.getF1R2Counts()))
                    .set(AltSiteRecordTableColumn.DEPTH.toString(), record.getDepth())
                    .set(AltSiteRecordTableColumn.ALT_ALLELE.toString(), record.getAltAllele().toString());
        }
    }

    /*** Reader ***/
    private static class AltSiteRecordTableReader extends TableReader<AltSiteRecord> {
        private AltSiteRecordTableReader(final File table) throws IOException {
            super(table);
        }

        @Override
        protected AltSiteRecord createRecord(final DataLine dataLine) {
            final String referenceContext = dataLine.get(AltSiteRecord.AltSiteRecordTableColumn.CONTEXT);
            final int[] baseCounts = Arrays.stream(dataLine.get(AltSiteRecordTableColumn.BASE_COUNTS).split(","))
                    .mapToInt(Integer::parseInt).toArray();
            final int[] f1r2Counts = Arrays.stream(dataLine.get(AltSiteRecordTableColumn.F1R2_COUNTS).split(","))
                    .mapToInt(Integer::parseInt).toArray();
            final int depth = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.DEPTH));
            final Nucleotide altAllele = Nucleotide.valueOf(dataLine.get(AltSiteRecordTableColumn.ALT_ALLELE));
            return new AltSiteRecord(referenceContext, baseCounts, f1r2Counts, depth, altAllele);
        }
    }

    /** Code for reading hyperparameters from a table **/
    public static List<AltSiteRecord> readAltSiteRecords(final File table, final int initialListSize) {
        List<AltSiteRecord> designMatrix = new ArrayList<>(initialListSize);
        try (AltSiteRecordTableReader reader = new AltSiteRecordTableReader(table)) {
            for (AltSiteRecord record : reader){
                designMatrix.add(record);
            }
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", table), e);
        }
        return designMatrix;
    }

    public static List<AltSiteRecord> readAltSiteRecords(final File table) {
        // arbitrarily initialize the list to size 100
        return readAltSiteRecords(table, 100);
    }
}
