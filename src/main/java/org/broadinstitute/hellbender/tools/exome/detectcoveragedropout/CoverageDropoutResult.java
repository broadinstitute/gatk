package org.broadinstitute.hellbender.tools.exome.detectcoveragedropout;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class for capturing the results of a coverage dropout.  A very simple class.
 *
 * If there is a coverage dropout (i.e. isCoverageDropout() is true), then the case sample and PoN combination is likely
 *  not a good match.
 *
 *  goodSegmentThreshold -- proportion of segments marked as "good" in order for isCoverageDropout to be true
 *  minWeight -- minimum
 *
 * @author lichtens &lt;lichtens@broadinstitute.org&gt;
 */
public final class CoverageDropoutResult {

    private boolean isCoverageDropout;
    private double goodSegmentThreshold;
    private double minWeight;
    private double minTargetProportion;
    private long numGoodSegments;
    private long numSegments;
    private long numSegmentsAfterFiltering;
    private double thresholdDistancePerSegment;
    private String comment;

    private static final String COVERAGE_DROPOUT_COLUMN = "isCoverageDropout";
    private static final String GOOD_SEGMENT_THRESHOLD_COLUMN = "goodSegmentThreshold";
    private static final String MININUM_WEIGHT_COLUMN = "minWeight";
    private static final String MININUM_TARGET_PROPORTION_COLUMN = "minTargetProportion";
    private static final String NUM_GOOD_SEGMENTS_COLUMN = "numGoodSegments";
    private static final String NUM_SEGMENTS_AFTER_FILTERING_COLUMN = "numSegmentsAfterFiltering";
    private static final String NUM_SEGMENTS_COLUMN = "numSegments";
    private static final String THRESHOLD_DISTANCE_PER_SEGMENTS_COLUMN = "thresholdDistancePerSegment";
    private static final String COMMENT_COLUMN = "comment";

    /** All columns that are required */
    private static final TableColumnCollection outputColumnNames =
            new TableColumnCollection(COVERAGE_DROPOUT_COLUMN, GOOD_SEGMENT_THRESHOLD_COLUMN, MININUM_WEIGHT_COLUMN,
            MININUM_TARGET_PROPORTION_COLUMN, NUM_GOOD_SEGMENTS_COLUMN, NUM_SEGMENTS_AFTER_FILTERING_COLUMN, NUM_SEGMENTS_COLUMN,
            THRESHOLD_DISTANCE_PER_SEGMENTS_COLUMN, COMMENT_COLUMN);

    public CoverageDropoutResult(final boolean isCoverageDropout, final double goodSegmentThreshold,
                                 final double minWeight, final double minTargetProportion, final long numGoodSegments,
                                 final long numSegments, final long numSegmentsAfterFiltering,
                                 final double thresholdDistancePerSegment, final String comment) {
        this.isCoverageDropout = isCoverageDropout;
        this.goodSegmentThreshold = goodSegmentThreshold;
        this.minWeight = minWeight;
        this.minTargetProportion = minTargetProportion;
        this.numGoodSegments = numGoodSegments;
        this.numSegments = numSegments;
        this.numSegmentsAfterFiltering = numSegmentsAfterFiltering;
        this.thresholdDistancePerSegment = thresholdDistancePerSegment;
        this.comment = comment;
    }

    public boolean isCoverageDropout() {
        return isCoverageDropout;
    }
    
    public double getGoodSegmentThreshold() {
        return goodSegmentThreshold;
    }

    public double getMinWeight() {
        return minWeight;
    }

    public double getMinTargetProportion() {
        return minTargetProportion;
    }

    public long getNumGoodSegments() {
        return numGoodSegments;
    }

    public long getNumSegments() {
        return numSegments;
    }

    public double getThresholdDistancePerSegment() {
        return thresholdDistancePerSegment;
    }

    public long getNumSegmentsAfterFiltering() {
        return numSegmentsAfterFiltering;
    }

    public String getComment() {
        return comment;
    }

    public static void writeCoverageDropoutResults(final List<CoverageDropoutResult> coverageDropoutResults, final File outFile) {
        try (final TableWriter<CoverageDropoutResult> writer = TableUtils.writer(outFile, outputColumnNames, (cdr, dataLine) -> {
                dataLine.set(COVERAGE_DROPOUT_COLUMN, cdr.isCoverageDropout())
                        .set(GOOD_SEGMENT_THRESHOLD_COLUMN, cdr.getGoodSegmentThreshold())
                        .set(MININUM_WEIGHT_COLUMN, cdr.getMinWeight())
                        .set(MININUM_TARGET_PROPORTION_COLUMN, cdr.getMinTargetProportion())
                        .set(NUM_GOOD_SEGMENTS_COLUMN, cdr.getNumGoodSegments())
                        .set(NUM_SEGMENTS_AFTER_FILTERING_COLUMN, cdr.getNumSegmentsAfterFiltering())
                        .set(NUM_SEGMENTS_COLUMN, cdr.getNumSegments())
                        .set(THRESHOLD_DISTANCE_PER_SEGMENTS_COLUMN, cdr.getThresholdDistancePerSegment())
                        .set(COMMENT_COLUMN, cdr.getComment());
            }))
        {
            writer.writeHeaderIfApplies();
            for (CoverageDropoutResult coverageDropoutResult: coverageDropoutResults) {
                writer.writeRecord(coverageDropoutResult);
            }
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile(outFile, ioe.getMessage());
        }
    }

    public static List<CoverageDropoutResult> readCoverageDropoutResultsFromTsv(final File inFile) throws IOException {
        try (final TableReader<CoverageDropoutResult> reader = TableUtils.reader(inFile,
                (columns, formatExceptionFactory) -> {
                    TableUtils.checkMandatoryColumns(columns, outputColumnNames, formatExceptionFactory);

                    // return the lambda to translate dataLines into coverage dropout result.
                    return (dataLine) -> new CoverageDropoutResult(dataLine.getBoolean(COVERAGE_DROPOUT_COLUMN),
                            dataLine.getDouble(GOOD_SEGMENT_THRESHOLD_COLUMN),
                            dataLine.getDouble(MININUM_WEIGHT_COLUMN),
                            dataLine.getDouble(MININUM_TARGET_PROPORTION_COLUMN),
                            dataLine.getLong(NUM_GOOD_SEGMENTS_COLUMN),
                            dataLine.getLong(NUM_SEGMENTS_COLUMN),
                            dataLine.getLong(NUM_SEGMENTS_AFTER_FILTERING_COLUMN),
                            dataLine.getDouble(THRESHOLD_DISTANCE_PER_SEGMENTS_COLUMN),
                            dataLine.get(COMMENT_COLUMN));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final UncheckedIOException e) {
            throw e.getCause();
        }
    }
}
