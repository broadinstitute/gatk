package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Reference context and ounts of substitution and indel errors at a ReadArtifact
 * Created by David Benjamin on 6/13/17.
 */
public class ArtifactSummary implements Locatable {
    private final String contig;
    private final int position;
    private final byte[] refBases;

    private final int[] baseCounts;
    private int insertionStartCount = 0;
    private int deletionStartCount = 0;

    public ArtifactSummary(final String contig, final int position, final byte[] refBases, final int[] baseCounts, final int insertionStartCount, final int deletionStartCount) {
        this.contig = contig;
        this.position = position;
        this.refBases = refBases;
        this.baseCounts = baseCounts;
        this.insertionStartCount = insertionStartCount;
        this.deletionStartCount = deletionStartCount;
    }

    public ArtifactSummary(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final int numBasesOnEitherSide, final byte minBaseQuality) {
        contig = alignmentContext.getContig();
        position = (int) alignmentContext.getPosition();

        final int basesToDiscardInFront = Math.max(position - referenceContext.getWindow().getStart() - numBasesOnEitherSide, 0);
        final byte[] allBases = referenceContext.getBases();
        final int lastBaseIndex = Math.min(basesToDiscardInFront + 2 * numBasesOnEitherSide + 1, allBases.length);
        refBases = Arrays.copyOfRange(allBases, basesToDiscardInFront, lastBaseIndex);

        baseCounts = new int[4];
        for (final PileupElement pe : alignmentContext.getBasePileup()) {
            if (pe.isDeletion()) {
                continue;
            } else if (pe.isBeforeDeletionStart()) {   // if context is ACGTTTT and we have a deletion GTTTT -> G, then the G is the "before deletion start" base
                deletionStartCount++;
            } else if (pe.isBeforeInsertion()) {
                insertionStartCount++;
            } else if ( pe.getQual() >= minBaseQuality) {
                final int index = BaseUtils.simpleBaseToBaseIndex(pe.getBase());
                if (index != -1) {
                    baseCounts[index]++;
                }
            }
        }
    }

    @Override
    public String getContig() { return contig; }

    @Override
    public int getStart() { return position; }

    @Override
    public int getEnd() { return position; }

    public byte[] getRefBases() {
        return refBases;
    }

    public int[] getBaseCounts() {
        return baseCounts;
    }

    public int getDeletionStartCount() {
        return deletionStartCount;
    }

    public int getInsertionStartCount() {
        return insertionStartCount;
    }

    //----- The following two public static methods read and write tables
    public static void writeArtifactSummaries(final List<ArtifactSummary> records, final File outputTable) {
        try ( ArtifactSummaryTableWriter writer = new ArtifactSummaryTableWriter(outputTable) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<ArtifactSummary> readArtifactSummaries(final File tableFile) {
        try( ArtifactSummaryTableReader reader = new ArtifactSummaryTableReader(tableFile) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing pileup summary tables
    public static class ArtifactSummaryTableWriter extends TableWriter<ArtifactSummary> {
        public ArtifactSummaryTableWriter(final File output) throws IOException {
            super(output, ArtifactSummaryTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final ArtifactSummary record, final DataLine dataLine) {
            final int[] baseCounts = record.getBaseCounts();
            dataLine.set(ArtifactSummaryTableColumn.CONTIG.toString(), record.getContig())
                    .set(ArtifactSummaryTableColumn.POSITION.toString(), record.getStart())
                    .set(ArtifactSummaryTableColumn.CONTEXT.toString(), new String(record.refBases))
                    .set(ArtifactSummaryTableColumn.A.toString(), baseCounts[0])
                    .set(ArtifactSummaryTableColumn.C.toString(), baseCounts[1])
                    .set(ArtifactSummaryTableColumn.G.toString(), baseCounts[2])
                    .set(ArtifactSummaryTableColumn.T.toString(), baseCounts[3])
                    .set(ArtifactSummaryTableColumn.INSERTION.toString(), record.getInsertionStartCount())
                    .set(ArtifactSummaryTableColumn.DELETION.toString(), record.getDeletionStartCount());
        }
    }

    private static class ArtifactSummaryTableReader extends TableReader<ArtifactSummary> {
        public ArtifactSummaryTableReader(final File file) throws IOException { super(file); }

        @Override
        protected ArtifactSummary createRecord(final DataLine dataLine) {
            final String contig = dataLine.get(ArtifactSummaryTableColumn.CONTIG);
            final int position = dataLine.getInt(ArtifactSummaryTableColumn.POSITION);
            final byte[] refBases = dataLine.get(ArtifactSummaryTableColumn.CONTEXT).getBytes();

            final int[] baseCounts = new int[] {dataLine.getInt(ArtifactSummaryTableColumn.A),
                    dataLine.getInt(ArtifactSummaryTableColumn.C),
                    dataLine.getInt(ArtifactSummaryTableColumn.G),
                    dataLine.getInt(ArtifactSummaryTableColumn.T)};
            final int insertionCount = dataLine.getInt(ArtifactSummaryTableColumn.INSERTION);
            final int deletionCount = dataLine.getInt(ArtifactSummaryTableColumn.DELETION);

            return new ArtifactSummary(contig, position, refBases, baseCounts, insertionCount, deletionCount);
        }
    }

    private enum ArtifactSummaryTableColumn {
        CONTIG("contig"),
        POSITION("position"),
        CONTEXT("context"),
        A("A"),
        C("C"),
        G("G"),
        T("T"),
        INSERTION("deletion"),
        DELETION("insertion");

        private final String columnName;

        ArtifactSummaryTableColumn(final String columnName) { this.columnName = Utils.nonNull(columnName); }

        @Override
        public String toString() { return columnName; }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

}
