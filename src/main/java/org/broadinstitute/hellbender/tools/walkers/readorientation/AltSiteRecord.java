package org.broadinstitute.hellbender.tools.walkers.readorientation;

import avro.shaded.com.google.common.primitives.Ints;
import htsjdk.samtools.util.SequenceUtil;
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
import java.util.Map;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.CollectDataForReadOrientationFilter.REGULAR_BASES;

/**
 * Created by tsato on 10/11/17.
 */
public class AltSiteRecord {
    private String contig;
    private int position;
    private String referenceContext;
    private int[] baseCounts;
    private int[] f1r2Counts;
    private int depth;
    private Nucleotide altAllele;


    public AltSiteRecord(final String contig, final int position, final String referenceContext, final int[] baseCounts,
                         final int[] f1r2Counts, final int depth, final Nucleotide altAllele){
        this.contig = contig;
        this.position = position;
        this.referenceContext = referenceContext;
        this.baseCounts = baseCounts;
        this.f1r2Counts = f1r2Counts;
        this.depth = depth;
        this.altAllele = altAllele;
    }

    public String getContig(){ return contig; }

    public int getPosition(){ return position; }

    public String getReferenceContext(){ return referenceContext; }

    public int[] getBaseCounts(){ return baseCounts; }

    public int getAltCount() {
        return baseCounts[altAllele.ordinal()];
    }

    public int[] getF1R2Counts(){ return f1r2Counts; }

    public int getDepth(){ return depth; }

    public Nucleotide getAltAllele(){ return altAllele; }

    public AltSiteRecord getReverseComplementOfRecord(){
        final int indexA = Nucleotide.A.ordinal();
        final int indexC = Nucleotide.C.ordinal();
        final int indexG = Nucleotide.G.ordinal();
        final int indexT = Nucleotide.T.ordinal();

        final int[] revCompBaseCounts = new int[]{ baseCounts[indexT], baseCounts[indexG], baseCounts[indexC], baseCounts[indexA] };
        final int[] revCompF1R2Counts = new int[]{ f1r2Counts[indexT], f1r2Counts[indexG], f1r2Counts[indexC], f1r2Counts[indexA] };

        final Nucleotide[] revCompMap = new Nucleotide[REGULAR_BASES.size()];
        revCompMap[indexA] = Nucleotide.T;
        revCompMap[indexC] = Nucleotide.G;
        revCompMap[indexG] = Nucleotide.C;
        revCompMap[indexT] = Nucleotide.A;

        return new AltSiteRecord(contig, position, SequenceUtil.reverseComplement(referenceContext), revCompBaseCounts,
                revCompF1R2Counts, depth, revCompMap[altAllele.ordinal()]);
    }



    /*** Table IO code ***/
    private enum AltSiteRecordTableColumn {
        CONTIG("contig"),
        POSITION("position"),
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
            dataLine.set(AltSiteRecordTableColumn.CONTIG.toString(), record.getContig())
                    .set(AltSiteRecordTableColumn.POSITION.toString(), record.getPosition())
                    .set(AltSiteRecordTableColumn.CONTEXT.toString(), record.getReferenceContext())
                    .set(AltSiteRecordTableColumn.BASE_COUNTS.toString(), Ints.join(",", record.getBaseCounts()))
                    .set(AltSiteRecordTableColumn.F1R2_COUNTS.toString(), Ints.join(",", record.getF1R2Counts()))
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
            final String contig = dataLine.get(AltSiteRecordTableColumn.CONTIG);
            final int position = Integer.valueOf(dataLine.get(AltSiteRecordTableColumn.POSITION));
            final String referenceContext = dataLine.get(AltSiteRecordTableColumn.CONTEXT);
            // TODO: eventually split basecounts into ref and alt counts
            final int[] baseCounts = Arrays.stream(dataLine.get(AltSiteRecordTableColumn.BASE_COUNTS).split(","))
                    .mapToInt(Integer::parseInt).toArray();
            final int[] f1r2Counts = Arrays.stream(dataLine.get(AltSiteRecordTableColumn.F1R2_COUNTS).split(","))
                    .mapToInt(Integer::parseInt).toArray();
            final int depth = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.DEPTH));
            final Nucleotide altAllele = Nucleotide.valueOf(dataLine.get(AltSiteRecordTableColumn.ALT_ALLELE));
            return new AltSiteRecord(contig, position, referenceContext, baseCounts, f1r2Counts, depth, altAllele);
        }
    }

    /** Code for reading alt site records from a table **/
    public static List<AltSiteRecord> readAltSiteRecords(final File table, final int initialListSize) {
        List<AltSiteRecord> designMatrix = new ArrayList<>(initialListSize);
        try (AltSiteRecordTableReader reader = new AltSiteRecordTableReader(table)) {
            reader.forEach(designMatrix::add);
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
