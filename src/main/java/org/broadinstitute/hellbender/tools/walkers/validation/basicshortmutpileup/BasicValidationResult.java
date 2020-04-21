package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

public class BasicValidationResult implements Locatable {
    private final int minValidationReadCount;
    private final boolean isEnoughValidationReads;
    private final boolean isOutOfNoiseFloor;
    private final double power;
    private final int validationAltCount;
    private final int validationRefCount;
    private final int discoveryAltCount;
    private final int discoveryRefCount;
    private final Locatable interval;
    private final Allele reference;
    private final Allele alternate;
    private final String filters;
    private final long numAltSupportingReadsInNormal;

    public BasicValidationResult(final Locatable interval, final int minValidationReadCount, final boolean isEnoughValidationReads,
                                 final boolean isOutOfNoiseFloor, final double power, final int validationAltCount,
                                 final int validationRefCount, final int discoveryAltCount, final int discoveryRefCount, final Allele ref,
                                 final Allele alt, final String filters, final long numAltSupportingReadsInNormal) {
        this.minValidationReadCount = minValidationReadCount;
        this.isEnoughValidationReads = isEnoughValidationReads;
        this.isOutOfNoiseFloor = isOutOfNoiseFloor;
        this.power = power;
        this.validationAltCount = validationAltCount;
        this.validationRefCount = validationRefCount;
        this.discoveryAltCount = discoveryAltCount;
        this.discoveryRefCount = discoveryRefCount;
        this.interval = interval;
        this.alternate = alt;
        this.reference = ref;
        this.filters = filters;
        this.numAltSupportingReadsInNormal = numAltSupportingReadsInNormal;
    }

    public int getValidationAltCount() {
        return validationAltCount;
    }

    public int getValidationRefCount() {
        return validationRefCount;
    }

    public int getDiscoveryAltCount() {
        return discoveryAltCount;
    }

    public int getDiscoveryRefCount() {
        return discoveryRefCount;
    }

    public int getMinValidationReadCount() {

        return minValidationReadCount;
    }

    public boolean isEnoughValidationReads() {
        return isEnoughValidationReads;
    }

    public boolean isOutOfNoiseFloor() {
        return isOutOfNoiseFloor;
    }

    public double getPower() {
        return power;
    }

    public Locatable getInterval() {
        return interval;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public Allele getReference() {
        return reference;
    }

    public Allele getAlternate() {
        return alternate;
    }

    public String getFilters() {
        return filters;
    }

    public long getNumAltSupportingReadsInNormal() {
        return numAltSupportingReadsInNormal;
    }

    //----- The following two public static methods read and write files
    public static void write(final List<BasicValidationResult> records, final File file) {
        try (final BasicValidationResult.BasicValidationResultTableWriter writer = new BasicValidationResult.BasicValidationResultTableWriter(IOUtils.toPath(file))) {
            writer.writeHeaderIfApplies();
            writer.writeAllRecords(records);
        } catch (final IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", file));
        }
    }

    public static List<BasicValidationResult> read(final File file) {
        try(final BasicValidationResult.BasicValidationResultTableReader reader = new BasicValidationResult.BasicValidationResultTableReader(IOUtils.toPath(file))) {
            return reader.toList();
        } catch (final IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", file));
        }
    }

    //-------- The following methods are boilerplate for reading and writing tables
    private static class BasicValidationResultTableWriter extends TableWriter<BasicValidationResult> {
        public BasicValidationResultTableWriter(final Path output) throws IOException {
            super(output, BasicValidationResult.BasicValidationResultTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final BasicValidationResult record, final DataLine dataLine) {
            dataLine.set(BasicValidationResultTableColumn.CONTIG.toString(), record.getContig());
            dataLine.set(BasicValidationResultTableColumn.START.toString(), record.getStart());
            dataLine.set(BasicValidationResultTableColumn.END.toString(), record.getEnd());
            dataLine.set(BasicValidationResultTableColumn.REF.toString(), record.getReference().getBaseString());
            dataLine.set(BasicValidationResultTableColumn.ALT.toString(), record.getAlternate().getBaseString());
            dataLine.set(BasicValidationResultTableColumn.DISCOVERY_ALT_COVERAGE.toString(), record.getDiscoveryAltCount());
            dataLine.set(BasicValidationResultTableColumn.DISCOVERY_REF_COVERAGE.toString(), record.getDiscoveryRefCount());
            dataLine.set(BasicValidationResultTableColumn.VALIDATION_ALT_COVERAGE.toString(), record.getValidationAltCount());
            dataLine.set(BasicValidationResultTableColumn.VALIDATION_REF_COVERAGE.toString(), record.getValidationRefCount());
            dataLine.set(BasicValidationResultTableColumn.MIN_VAL_COUNT.toString(), record.getMinValidationReadCount());
            dataLine.set(BasicValidationResultTableColumn.POWER.toString(), record.getPower());
            dataLine.set(BasicValidationResultTableColumn.IS_NOT_NOISE.toString(), record.isOutOfNoiseFloor());
            dataLine.set(BasicValidationResultTableColumn.IS_ENOUGH_VALIDATION_COVERAGE.toString(), record.isEnoughValidationReads());
            dataLine.set(BasicValidationResultTableColumn.DISCOVERY_VCF_FILTER.toString(), record.getFilters() == null ? "" : record.getFilters());
            dataLine.set(BasicValidationResultTableColumn.NUM_ALT_READS_IN_VALIDATION_NORMAL.toString(), record.getNumAltSupportingReadsInNormal());
        }
    }

    private static class BasicValidationResultTableReader extends TableReader<BasicValidationResult> {
        public BasicValidationResultTableReader(final Path path) throws IOException {
            super(path);
        }

        @Override
        protected BasicValidationResult createRecord(final DataLine dataLine) {
            final String contig = dataLine.get(BasicValidationResultTableColumn.CONTIG);
            final int start = dataLine.getInt(BasicValidationResultTableColumn.START);
            final int end = dataLine.getInt(BasicValidationResultTableColumn.END);
            final Allele ref = Allele.create(dataLine.get(BasicValidationResultTableColumn.REF).getBytes(), true);
            final Allele alt = Allele.create(dataLine.get(BasicValidationResultTableColumn.ALT).getBytes(), false);
            final int discoveryAltCount = dataLine.getInt(BasicValidationResultTableColumn.DISCOVERY_ALT_COVERAGE);
            final int discoveryRefCount = dataLine.getInt(BasicValidationResultTableColumn.DISCOVERY_REF_COVERAGE);
            final int validationAltCount = dataLine.getInt(BasicValidationResultTableColumn.VALIDATION_ALT_COVERAGE);
            final int validationRefCount = dataLine.getInt(BasicValidationResultTableColumn.VALIDATION_REF_COVERAGE);
            final int minValidationReadCount = dataLine.getInt(BasicValidationResultTableColumn.MIN_VAL_COUNT);
            final double power = dataLine.getDouble(BasicValidationResultTableColumn.POWER);
            final boolean isOutOfNoiseFloor = dataLine.getBoolean(BasicValidationResultTableColumn.IS_NOT_NOISE.toString());
            final boolean isEnoughValidationReads = dataLine.getBoolean(BasicValidationResultTableColumn.IS_ENOUGH_VALIDATION_COVERAGE.toString());
            final String filters = dataLine.get(BasicValidationResultTableColumn.DISCOVERY_VCF_FILTER);
            final long numAltSupportingReadsInNormal = dataLine.getLong(BasicValidationResultTableColumn.NUM_ALT_READS_IN_VALIDATION_NORMAL);

            final Locatable interval = new SimpleInterval(contig, start, end);

            return new BasicValidationResult(interval, minValidationReadCount, isEnoughValidationReads, isOutOfNoiseFloor,
                    power, validationAltCount, validationRefCount, discoveryAltCount, discoveryRefCount, ref, alt,
                    filters, numAltSupportingReadsInNormal);
        }
    }

    private enum BasicValidationResultTableColumn {
        CONTIG("CONTIG"),
        START("START"),
        END("END"),
        REF("ref_allele"),
        ALT("alt_allele"),
        DISCOVERY_ALT_COVERAGE("t_alt_count"),
        DISCOVERY_REF_COVERAGE("t_ref_count"),
        VALIDATION_ALT_COVERAGE("tv_alt_count"),
        VALIDATION_REF_COVERAGE("tv_ref_count"),
        MIN_VAL_COUNT("min_val_count"),
        POWER("power"),
        IS_NOT_NOISE("validated"),
        IS_ENOUGH_VALIDATION_COVERAGE("sufficient_tv_alt_coverage"),
        DISCOVERY_VCF_FILTER("discovery_vcf_filter"),
        NUM_ALT_READS_IN_VALIDATION_NORMAL("num_alt_reads_in_validation_normal");

        private final String columnName;

        BasicValidationResultTableColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(CONTIG, START, END, REF, ALT, DISCOVERY_ALT_COVERAGE, DISCOVERY_REF_COVERAGE,
                VALIDATION_ALT_COVERAGE, VALIDATION_REF_COVERAGE, MIN_VAL_COUNT, POWER, IS_NOT_NOISE, IS_ENOUGH_VALIDATION_COVERAGE,
                DISCOVERY_VCF_FILTER, NUM_ALT_READS_IN_VALIDATION_NORMAL);
    }
}
