package org.broadinstitute.hellbender.tools.spark.longread;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.StringUtils;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.Arrays;
import java.util.List;

/**
 * Collect some metrics on alignments produced from long read alignments.
 *
 * <ul>
 * <li>Total clipped bases</li>
 * <li>Total count of unmapped queries</li>
 * <li>Total count of supplementary (SAM flag 2048) alignment records</li>
 * <li>Total count of secondary (SAM flag 256) alignment records</li>
 * <li>average number CIGAR operations per 100 base pair (mapped queries)</li>
 * </ul>
 * <p>
 * In addition to summary metrics, we also collect a metric that aims to measure contiguity of the alignments.
 * That is, for one alignment, there typically are several operations in the CIGAR.
 * We construct a vector of lengths of the alignment blocks, then normalize it against the total non-clipped bases.
 * This is output as a vector of 2-digit integers (that is, the fraction with precision up to 2 decimal places then * 100)
 * for each SAM record.
 * The output format is
 * "read_name\tclipped_start\tclipped_end\tSAM_flag\tvector_separated_by_,"
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Collect some metrics on alignments produced from long read alignments",
        summary = "Collect some metrics on alignments produced from long read alignments",
        programGroup = LongReadAnalysisProgramGroup.class)
public class CollectLongReadAlignmentMetricsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String output;

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * This is overriden because the alignment might be not "WellFormed".
     */
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(new ReadFilterLibrary.AllowAllReadsReadFilter());
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> allAlignmentRecords = getReads().cache();

        reportHighLevelCounts(allAlignmentRecords);

        final SAMFileHeader headerForReads = getHeaderForReads();
        final JavaRDD<ReadInformation> cache = allAlignmentRecords
                .filter(read -> ! read.isUnmapped())
                .filter(read -> ! read.isSecondaryAlignment())
                .map(read -> new ReadInformation(read, headerForReads));
        allAlignmentRecords.unpersist(false);

        reportSummaryMetricsOnCIGARs(cache);

        writePerAlignmentReport(cache);
    }

    private void reportHighLevelCounts(final JavaRDD<GATKRead> allAlignmentRecords) {
        final long totalCount = allAlignmentRecords.count();
        final long unmappedCount = allAlignmentRecords.filter(GATKRead::isUnmapped).count();
        final long secondaryMappingCounts = allAlignmentRecords.filter(GATKRead::isSecondaryAlignment).count();
        final long supplementaryMappingCounts = allAlignmentRecords.filter(GATKRead::isSupplementaryAlignment).count();
        formattedLogger("Total number of alignment records: ", totalCount);
        formattedLogger("Total number of unmapped queries: ", unmappedCount);
        formattedLogger("Total number of secondary (256) mappings: ", secondaryMappingCounts);
        formattedLogger("Total number of supplementary (2048) mappings: ", supplementaryMappingCounts);
    }

    private static final class ReadInformation implements Serializable {
        private static final long serialVersionUID = 1L;

        final String name;
        final int start;
        final int end;
        final int samFlag;
        final double clippedBases;
        final double clippedLength;
        final double numOperations;
        final int[] alnBlockLength;

        ReadInformation(final GATKRead read, final SAMFileHeader headerForReads) {
            name = read.getName();
            start = read.getStart();
            end = read.getEnd();
            samFlag = read.convertToSAMRecord(headerForReads).getFlags();

            int nOperations = 0;
            int softClippedBases = 0;
            int hardClippedBases = 0;
            for (final CigarElement cigarElement : read.getCigarElements()) {
                if (cigarElement.getOperator() == CigarOperator.H) {
                    hardClippedBases += cigarElement.getLength();
                } else if (cigarElement.getOperator() == CigarOperator.S) {
                    softClippedBases += cigarElement.getLength();
                } else {
                    ++nOperations;
                }
            }
            numOperations = nOperations;
            clippedBases = softClippedBases + hardClippedBases;
            clippedLength = read.getLength() - softClippedBases;

            alnBlockLength = read.getCigarElements().stream()
                    .filter(ele -> ele.getOperator().isAlignment())
                    .mapToInt(     ele -> (int) Math.round( (100.0 * ele.getLength()) / clippedLength )     )
                    .toArray();
        }

        @Override
        public String toString() {
            final String arr = Arrays.toString(alnBlockLength)
                    .replaceFirst("\\[", "").replaceFirst("]", "");
            return String.format("%s\t%d\t%d\t%d\t%s", name, start, end, samFlag, arr);
        }
    }

    private void reportSummaryMetricsOnCIGARs( final JavaRDD<ReadInformation> cache) {
        final Double totalClippedBases = cache.map(info -> info.clippedBases).reduce(Double::sum);
        formattedLogger("Total number of clipped bases (soft + hard): ", totalClippedBases);

        final Double numerator = cache.map(info -> info.numOperations).reduce(Double::sum);
        final Double denominator = cache.map(info -> info.clippedLength).reduce(Double::sum);
        double frac = numerator / denominator;
        final double meanOperationPerHundredBases = Math.round(100. * (100. * frac) ) / 100.; // keep 2 decimals
        final double meanOperationPerTwoHundredBases = Math.round(100. * (200 * frac) ) / 100.;
        final double meanOperationPerFiveHundredBases = Math.round(100. * (500. * frac) ) / 100.;
        final double meanOperationPerThousandBases = Math.round(100. * (1000. * frac) ) / 100.;
        formattedLogger("Mean number of CIGAR operations per  100 base pair: ", meanOperationPerHundredBases);
        formattedLogger("Mean number of CIGAR operations per  200 base pair: ", meanOperationPerTwoHundredBases);
        formattedLogger("Mean number of CIGAR operations per  500 base pair: ", meanOperationPerFiveHundredBases);
        formattedLogger("Mean number of CIGAR operations per 1000 base pair: ", meanOperationPerThousandBases);
    }

    private void writePerAlignmentReport(final JavaRDD<ReadInformation> cache) {
        final List<String> perAlnReport = cache.map(ReadInformation::toString).collect();
        try ( final Writer writer = new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(output))) ) {
            writer.write("READ_NAME\tSTART\tEND\tSAM_FLAG\tALN_LENGTHS\n");
            for ( final String line : perAlnReport ) {
                writer.write(line);
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Unable to write reports to " + output, ioe);
        }
    }

    private void formattedLogger(final String prefix, final Number numerical) {
        final int rightPaddingLength = 53;
        logger.info(StringUtils.rightPad(prefix, rightPaddingLength) + numerical);
    }
}
