package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import java.nio.file.Path;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.*;

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
    private String contig;
    private int position;


    public AltSiteRecord(final String referenceContext, final int refCount, final int altCount,
                         final int refF1R2, final int altF1R2, final Nucleotide altAllele,
                         final String contig, final int position){
        Utils.validateArg(altAllele != null && altAllele.isStandard(), "altAllele must be one of {A,C,G,T} but got " + altAllele);

        // records with the same context share one String to save memory
        this.referenceContext = referenceContext.intern();
        this.refCount = refCount;
        this.altCount = altCount;
        this.refF1R2 = refF1R2;
        this.altF1R2 = altF1R2;
        this.depth = refCount + altCount;
        this.altAllele = altAllele;
        // likewise for the contig
        this.contig = Utils.nonNull(contig).intern();
        this.position = position;
    }

    public String getReferenceContext(){ return referenceContext; }

    public int getRefCount() { return refCount; }

    public int getAltCount() { return altCount; }

    public int getRefF1R2() { return refF1R2; }

    public int getAltF1R2(){ return altF1R2; }

    public int getDepth(){ return depth; }

    public Nucleotide getAltAllele(){ return altAllele; }

    public String getContig() { return contig; }

    public int getPosition() { return position; }


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
                newAltF1R2, revCompOfAlt, contig, position);
    }

    /*** Table IO code ***/
    private enum AltSiteRecordTableColumn {
        CONTEXT("context"),
        REF_COUNT("ref_count"),
        ALT_COUNT("alt_count"),
        REF_F1R2("ref_f1r2"),
        ALT_F1R2("alt_f1r2"),
        ALT_BASE("alt"),
        CONTIG("contig"),
        POSITION("position");

        private String columnName;

        AltSiteRecordTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() {
            return columnName;
        }

        // The columns of a table of one reference context, which is in the table's name.  The depth is the ref count
        // plus the alt count.
        public static final TableColumnCollection COLUMNS = new TableColumnCollection(
                REF_COUNT, ALT_COUNT, REF_F1R2, ALT_F1R2, ALT_BASE, CONTIG, POSITION);

        // the columns of a table of any reference contexts
        public static final TableColumnCollection COLUMNS_WITH_CONTEXT = new TableColumnCollection((Object[]) values());
    }


    /*** Writer ***/
    public static class AltSiteRecordTableWriter extends TableWriter<AltSiteRecord> {
        private final boolean writeContexts;

        // a writer of a table of the records of one reference context, which is in the table's name
        public AltSiteRecordTableWriter(final Path output, final String sample) throws IOException {
            this(output, sample, false);
        }

        // a writer of a table with a column of reference contexts if writeContexts is true
        AltSiteRecordTableWriter(final Path output, final String sample, final boolean writeContexts) throws IOException {
            super(output, writeContexts ? AltSiteRecordTableColumn.COLUMNS_WITH_CONTEXT : AltSiteRecordTableColumn.COLUMNS);
            this.writeContexts = writeContexts;
            writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
        }

        @Override
        protected void composeLine(final AltSiteRecord record, final DataLine dataLine) {
            // it'd be nice to set() less manually...
            if (writeContexts) {
                dataLine.set(AltSiteRecordTableColumn.CONTEXT.toString(), record.getReferenceContext());
            }
            dataLine.set(AltSiteRecordTableColumn.REF_COUNT.toString(), record.getRefCount())
                    .set(AltSiteRecordTableColumn.ALT_COUNT.toString(), record.getAltCount())
                    .set(AltSiteRecordTableColumn.REF_F1R2.toString(), record.getRefF1R2())
                    .set(AltSiteRecordTableColumn.ALT_F1R2.toString(), record.getAltF1R2())
                    .set(AltSiteRecordTableColumn.ALT_BASE.toString(), record.getAltAllele().toString())
                    .set(AltSiteRecordTableColumn.CONTIG.toString(), record.getContig())
                    .set(AltSiteRecordTableColumn.POSITION.toString(), record.getPosition());
        }
    }

    /*** Reader ***/
    static class AltSiteRecordTableReader extends TableReader<AltSiteRecord> {
        // the reference context of every record of a table of one context, or null for a table with a column of contexts
        private final String referenceContext;

        // a reader of a table with a column of reference contexts
        AltSiteRecordTableReader(final Path table) throws IOException {
            this(table, null);
        }

        // a reader of a table of one reference context, or, if it is null, of a table with a column of contexts
        AltSiteRecordTableReader(final Path table, final String referenceContext) throws IOException {
            super(table);
            this.referenceContext = referenceContext;
        }

        @Override
        protected AltSiteRecord createRecord(final DataLine dataLine) {
            final String context = referenceContext != null ? referenceContext : dataLine.get(AltSiteRecordTableColumn.CONTEXT);
            final int refCount = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.REF_COUNT));
            final int altCount = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.ALT_COUNT));
            final int refF1R2 = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.REF_F1R2));
            final int altF1R2 = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.ALT_F1R2));
            final Nucleotide altAllele = Nucleotide.valueOf(dataLine.get(AltSiteRecordTableColumn.ALT_BASE));
            final String contig = dataLine.get(AltSiteRecordTableColumn.CONTIG);
            final int position = Integer.parseInt(dataLine.get(AltSiteRecordTableColumn.POSITION));
            return new AltSiteRecord(context, refCount, altCount, refF1R2, altF1R2, altAllele, contig, position);
        }
    }

    /** Code for reading sample name and alt site records from a table of one reference context **/
    public static Pair<String, List<AltSiteRecord>> readAltSiteRecords(final Path table, final String referenceContext, final int initialListSize) {
        List<AltSiteRecord> records = new ArrayList<>(initialListSize);
        final String sample;
        try (AltSiteRecordTableReader reader = new AltSiteRecordTableReader(table, Utils.nonNull(referenceContext))) {
            sample = reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG);
            reader.forEach(records::add);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", table), e);
        }
        return ImmutablePair.of(sample, records);
    }

    public static Pair<String, List<AltSiteRecord>> readAltSiteRecords(final Path table, final String referenceContext) {
        // arbitrarily initialize the list to size 100
        return readAltSiteRecords(table, referenceContext, 100);
    }
}
