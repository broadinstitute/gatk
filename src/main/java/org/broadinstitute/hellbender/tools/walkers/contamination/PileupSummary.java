package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import java.nio.file.Path;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

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
    public double getMinorAlleleFraction() {
        final double altFraction = getAltFraction();
        return FastMath.min(altFraction, 1 - altFraction);
    }


    //----- The following two public static methods read and write pileup summary files
    public static void writeToFile(final String sample, final List<PileupSummary> records, final File outputTable) {
        try ( PileupSummaryTableWriter writer = new PileupSummaryTableWriter(IOUtils.fileToPath(outputTable)) ) {
            writer.writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    // Takes a list of PileupSummaryTable files and write them all in one output file in order
    public static void writeToFile(final List<File> inputFiles, final File output) {
        boolean headerWritten = false;
        String sample = "";
        try ( PileupSummaryTableWriter writer = new PileupSummaryTableWriter(output.toPath()) ) {
            for (final File inputFile : inputFiles){
                try ( PileupSummaryTableReader reader = new PileupSummaryTableReader(inputFile.toPath())){
                    if (! headerWritten){
                        sample = reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG);
                        writer.writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
                        headerWritten = true;
                    }

                    final String thisSample = reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG);
                    if (! thisSample.equals(sample)){
                        throw new UserException.BadInput(String.format("Combining PileupSummaryTables from different samples is not supported. Got samples %s and %s",
                                sample, thisSample));
                    }

                    final List<PileupSummary> pileupSummaries = reader.toList();
                    writer.writeAllRecords(pileupSummaries);
                } catch (IOException e){
                    throw new UserException(String.format("Encountered an IO exception while reading from %s.", inputFile));
                }
            }

        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", output));
        }
    }

    public static ImmutablePair<String, List<PileupSummary>> readFromFile(final File tableFile) {
        try( PileupSummaryTableReader reader = new PileupSummaryTableReader(IOUtils.fileToPath(tableFile)) ) {
            final List<PileupSummary> pileupSummaries = reader.toList();
            return ImmutablePair.of(reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG), pileupSummaries);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    public static class PileupSummaryComparator implements Comparator<PileupSummary> {
        final SAMSequenceDictionary sequenceDictionary;
        final List<String> contigsInOrder;

        public PileupSummaryComparator(final SAMSequenceDictionary sequenceDictionary){
            this.sequenceDictionary = sequenceDictionary;
            contigsInOrder = sequenceDictionary.getSequences().stream()
                    .map(ssr -> ssr.getSequenceName())
                    .collect(Collectors.toList());
        }

        @Override
        public int compare(PileupSummary ps1, PileupSummary ps2) {
            // Use Contig Index in case the contig name is e.g. chr2
            final int contigIndex1 = contigsInOrder.indexOf(ps1.getContig());
            final int contigIndex2 = contigsInOrder.indexOf(ps2.getContig());

            if (contigIndex1 != contigIndex2){
                return Integer.compare(contigIndex1, contigIndex2);
            } else {
                return Integer.compare(ps1.getStart(), ps2.getStart());
            }
        }
    }

    //-------- The following methods are boilerplate for reading and writing pileup summary tables
    public static class PileupSummaryTableWriter extends TableWriter<PileupSummary> {
        public PileupSummaryTableWriter(final Path output) throws IOException {
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
        public PileupSummaryTableReader(final Path path) throws IOException { super(path); }

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
