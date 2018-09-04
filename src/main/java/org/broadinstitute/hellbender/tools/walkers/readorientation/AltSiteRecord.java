package org.broadinstitute.hellbender.tools.walkers.readorientation;

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
import java.util.List;

/**
 * Created by tsato on 10/11/17.
 */
public class AltSiteRecord {
    private String referenceContext;
    private int refCount;
    private int altCount;
    private int refF1R2;
    private int altF1R2;
    private int depth;
    private Nucleotide altAllele;


    public AltSiteRecord(final String referenceContext, final int refCount, final int altCount,
                         final int refF1R2, final int altF1R2, final Nucleotide altAllele){
        Utils.validateArg(altAllele != null && altAllele.isStandard(), "altAllele must be one of {A,C,G,T} but got " + altAllele);

        this.referenceContext = referenceContext;
        this.refCount = refCount;
        this.altCount = altCount;
        this.refF1R2 = refF1R2;
        this.altF1R2 = altF1R2;
        this.depth = refCount + altCount;
        this.altAllele = altAllele;
    }

    public String getReferenceContext(){ return referenceContext; }

    public int getRefCount() { return refCount; }

    public int getAltCount() { return altCount; }

    public int getRefF1R2() { return refF1R2; }

    public int getAltF1R2(){ return altF1R2; }

    public int getDepth(){ return depth; }

    public Nucleotide getAltAllele(){ return altAllele; }


    /**
     * Contract: only call this method on an {@link AltSiteRecord} whose reference context is *not* in the
     * @{link F1R2FilterConstants.CANONICAL_KMERS}. The idea is that this method should be called only to make standard
     * the non-standard representation of an alt site record.
     */
    public AltSiteRecord getReverseComplementOfRecord(){
        Utils.validate(!F1R2FilterConstants.CANONICAL_KMERS.contains(referenceContext), "for consistency, don't make the " +
                        "revcomp record of a canonical reference context");
        final Nucleotide revCompOfAlt = altAllele.complement();
        final int newRefF1R2 = refCount - refF1R2;
        final int newAltF1R2 = altCount - altF1R2;
        return new AltSiteRecord(SequenceUtil.reverseComplement(referenceContext), refCount, altCount, newRefF1R2,
                newAltF1R2, revCompOfAlt);
    }

    /*** Table IO code ***/
    private enum AltSiteRecordTableColumn {
        CONTEXT("context"),
        REF_COUNT("ref_count"),
        ALT_COUNT("alt_count"),
        REF_F1R2("ref_f1r2"),
        ALT_F1R2("alt_f1r2"),
        DEPTH("depth"),
        ALT_BASE("alt");

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
                    .set(AltSiteRecordTableColumn.REF_COUNT.toString(), record.getRefCount())
                    .set(AltSiteRecordTableColumn.ALT_COUNT.toString(), record.getAltCount())
                    .set(AltSiteRecordTableColumn.REF_F1R2.toString(), record.getRefF1R2())
                    .set(AltSiteRecordTableColumn.ALT_F1R2.toString(), record.getAltF1R2())
                    .set(AltSiteRecordTableColumn.DEPTH.toString(), record.getDepth())
                    .set(AltSiteRecordTableColumn.ALT_BASE.toString(), record.getAltAllele().toString());
        }
    }

    /*** Reader ***/
    private static class AltSiteRecordTableReader extends TableReader<AltSiteRecord> {
        private AltSiteRecordTableReader(final File table) throws IOException {
            super(table);
        }

        @Override
        protected AltSiteRecord createRecord(final DataLine dataLine) {
            final String referenceContext = dataLine.get(AltSiteRecordTableColumn.CONTEXT);
            final int refCount = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.REF_COUNT));
            final int altCount = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.ALT_COUNT));
            final int refF1R2 = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.REF_F1R2));
            final int altF1R2 = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.ALT_F1R2));
            final int depth = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.DEPTH));
            final Nucleotide altAllele = Nucleotide.valueOf(dataLine.get(AltSiteRecordTableColumn.ALT_BASE));
            return new AltSiteRecord(referenceContext, refCount, altCount, refF1R2, altF1R2, altAllele);
        }
    }

    /** Code for reading alt site records from a table **/
    public static List<AltSiteRecord> readAltSiteRecords(final File table, final int initialListSize) {
        List<AltSiteRecord> records = new ArrayList<>(initialListSize);
        try (AltSiteRecordTableReader reader = new AltSiteRecordTableReader(table)) {
            reader.forEach(records::add);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", table), e);
        }
        return records;
    }

    public static List<AltSiteRecord> readAltSiteRecords(final File table) {
        // arbitrarily initialize the list to size 100
        return readAltSiteRecords(table, 100);
    }
}
