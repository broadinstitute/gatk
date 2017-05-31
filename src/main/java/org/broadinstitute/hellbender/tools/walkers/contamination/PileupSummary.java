package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.OptionalDouble;
import java.util.stream.StreamSupport;

/**
 * Created by David Benjamin on 2/14/17.
 */
public class PileupSummary implements Locatable {
    private final String contig;
    private final int position;
    private final int refCount;
    private final int altCount;
    private final int otherAltsCount;
    private final int totalCount;

    private final double alleleFrequency;

    public PileupSummary(String contig, int position, int refCount, int altCount, int otherAltsCount, double alleleFrequency) {
        this.contig = contig;
        this.position = position;
        this.altCount = altCount;
        this.refCount = refCount;
        this.otherAltsCount = otherAltsCount;
        this.totalCount = refCount + altCount + otherAltsCount;
        this.alleleFrequency = alleleFrequency;
    }

    public PileupSummary(final VariantContext vc, final ReadPileup pileup) {
        contig = vc.getContig();
        position = vc.getStart();
        alleleFrequency = vc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0);
        final byte altBase = vc.getAlternateAllele(0).getBases()[0];
        final byte refBase = vc.getReference().getBases()[0];
        final int[] baseCounts = pileup.getBaseCounts();
        altCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(altBase)];
        refCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(refBase)];
        totalCount = (int) MathUtils.sum(baseCounts);
        otherAltsCount = totalCount - altCount - refCount;
    }

    @Override
    public String getContig() { return contig; }

    @Override
    public int getStart() { return position; }

    @Override
    public int getEnd() { return position; }

    public int getAltCount() {
        return altCount;
    }
    public int getRefCount() {
        return refCount;
    }
    public int getOtherAltCount() {
        return otherAltsCount;
    }
    public int getTotalCount() {
        return totalCount;
    }
    public double getAlleleFrequency() {
        return alleleFrequency;
    }
    public double getRefFrequency() {
        return 1 - alleleFrequency;
    }
    public double getAltFraction() {
        return totalCount == 0 ? 0 : (double) altCount / totalCount;
    }


    //----- The following two public static methods read and write pileup summary files
    public static void writePileupSummaries(final List<PileupSummary> records, final File outputTable) {
        try ( PileupSummaryTableWriter writer = new PileupSummaryTableWriter(outputTable) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<PileupSummary> readPileupSummaries(final File tableFile) {
        try( PileupSummaryTableReader reader = new PileupSummaryTableReader(tableFile) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing pileup summary tables
    private static class PileupSummaryTableWriter extends TableWriter<PileupSummary> {
        private PileupSummaryTableWriter(final File output) throws IOException {
            super(output, PileupSummaryTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final PileupSummary record, final DataLine dataLine) {
            dataLine.set(PileupSummaryTableColumn.CONTIG.toString(), record.getContig())
                    .set(PileupSummaryTableColumn.POSITION.toString(), record.getStart())
                    .set(PileupSummaryTableColumn.REF_COUNT.toString(), record.getRefCount())
                    .set(PileupSummaryTableColumn.ALT_COUNT.toString(), record.getAltCount())
                    .set(PileupSummaryTableColumn.OTHER_ALT_COUNT.toString(), record.getOtherAltCount())
                    .set(PileupSummaryTableColumn.ALT_ALLELE_FREQUENCY.toString(), record.getAlleleFrequency());
        }
    }

    private static class PileupSummaryTableReader extends TableReader<PileupSummary> {
        public PileupSummaryTableReader(final File file) throws IOException { super(file); }

        @Override
        protected PileupSummary createRecord(final DataLine dataLine) {
            final String contig = dataLine.get(PileupSummaryTableColumn.CONTIG);
            final int position = dataLine.getInt(PileupSummaryTableColumn.POSITION);
            final int refCount = dataLine.getInt(PileupSummaryTableColumn.REF_COUNT);
            final int altCount = dataLine.getInt(PileupSummaryTableColumn.ALT_COUNT);
            final int otherAltCount = dataLine.getInt(PileupSummaryTableColumn.OTHER_ALT_COUNT);
            final double alleleFrequency = dataLine.getDouble(PileupSummaryTableColumn.ALT_ALLELE_FREQUENCY);

            return new PileupSummary(contig, position, refCount, altCount, otherAltCount, alleleFrequency);
        }
    }

    private enum PileupSummaryTableColumn {
        CONTIG("contig"),
        POSITION("position"),
        REF_COUNT("ref_count"),
        ALT_COUNT("alt_count"),
        OTHER_ALT_COUNT("other_alt_count"),
        ALT_ALLELE_FREQUENCY("allele_frequency");

        private final String columnName;

        PileupSummaryTableColumn(final String columnName) { this.columnName = Utils.nonNull(columnName); }

        @Override
        public String toString() { return columnName; }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

}
