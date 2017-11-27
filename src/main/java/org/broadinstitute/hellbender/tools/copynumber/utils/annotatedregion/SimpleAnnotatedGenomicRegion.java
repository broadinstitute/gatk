package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Sets;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Simple class that just has an interval and name value pairs.
 *
 * When reading a TSV file of simple annotated genomic regions, the genomic region columns for the header when reading
 *  are specified in {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
 * {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and {@link SimpleAnnotatedGenomicRegion::END_HEADER}
 */
final public class SimpleAnnotatedGenomicRegion implements Locatable {
    public final static String CONTIG_HEADER = "CONTIG";
    public final static String START_HEADER = "START";
    public final static String END_HEADER = "END";

    private SimpleInterval interval;
    private final SortedMap<String, String> annotations;

    public SimpleAnnotatedGenomicRegion(final SimpleInterval interval, final SortedMap<String, String> annotations) {
        this.interval = interval;
        this.annotations = annotations;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public SortedMap<String, String> getAnnotations() {
        return annotations;
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

    public void setInterval(final Locatable interval) {
        this.interval = new SimpleInterval(interval);
    }

    public void setEnd(final int end) {
        this.interval = new SimpleInterval(this.interval.getContig(), this.interval.getStart(), end);
    }

    public String getAnnotationValue(final String annotationName) {
        return annotations.get(annotationName);
    }

    public boolean hasAnnotation(final String annotationName) {
        return annotations.containsKey(annotationName);
    }

    /**
     * Creates the annotation if it does not exist.
     *
     * @param annotationName the name for the annotation
     * @param annotationValue the value
     * @return the previous value or null
     */
    public String setAnnotation(final String annotationName, final String annotationValue) {
        return annotations.put(annotationName, annotationValue);
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final SimpleAnnotatedGenomicRegion that = (SimpleAnnotatedGenomicRegion) o;
        return this.interval.equals(that.getInterval()) && this.getAnnotations().equals(that.getAnnotations());
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + annotations.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return interval.toString() + " :: " + annotations.entrySet().stream()
                .map(e -> e.getKey() + "->" + e.getValue()).collect(Collectors.joining(","));
    }

    /**
     *
     * This method will likely be removed in future releases.
     *
     *  See {@link SimpleAnnotatedGenomicRegion#readAnnotatedRegions(File, Set)}, except that all columns (except those to define the interval)
     *   are considered headers of interest.
     *
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @return annotated regions with one line in the input file for each entry of the list.  Never {<@UNVERIFIED|@code> null}
     */
    public static List<SimpleAnnotatedGenomicRegion> readAnnotatedRegions(final File tsvRegionFile) {
        try (final TableReader<SimpleAnnotatedGenomicRegion> reader = new TableReader<SimpleAnnotatedGenomicRegion>(tsvRegionFile) {

            @Override
            protected boolean isCommentLine(final String[] line) {
                return super.isCommentLine(line) || line[0].startsWith("@");
            }

            @Override
            protected SimpleAnnotatedGenomicRegion createRecord(final DataLine dataLine) {
                final Set<String> headersOfInterestPresent = Sets.difference(new HashSet<>(this.columns().names()),
                        Sets.newHashSet(CONTIG_HEADER, START_HEADER, END_HEADER));
                final Map<String, String> annotationMap = headersOfInterestPresent.stream()
                        .collect(Collectors.toMap(Function.identity(), s -> dataLine.get(s)));
                return new SimpleAnnotatedGenomicRegion( new SimpleInterval(dataLine.get(CONTIG_HEADER), dataLine.getInt(START_HEADER), dataLine.getInt(END_HEADER)),
                        new TreeMap<>(annotationMap));
            }
        }) {
            return reader.toList();
        } catch (final IOException ioe) {
            throw new UserException.CouldNotReadInputFile("Cannot read input file: " + tsvRegionFile.getAbsolutePath(), ioe);
        }
    }

    /**
     *  This method will likely be removed in future releases.
     *
     *  Reads entire TSV file in one command and stores in RAM.  Please see {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} for defining the genomic region.
     *
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @param headersOfInterest -- The columns that should be in the resulting annotated region.  This list should
     *                          not include any headers that are used to define the region (e.g. contig, start, end)
     * @return annotated regions with one line in the input file for each entry of the list.  Never {<@UNVERIFIED|@code> null}
     */
    public static List<SimpleAnnotatedGenomicRegion> readAnnotatedRegions(final File tsvRegionFile, final Set<String> headersOfInterest) {
        try (final TableReader<SimpleAnnotatedGenomicRegion> reader = new TableReader<SimpleAnnotatedGenomicRegion>(tsvRegionFile) {

            @Override
            protected boolean isCommentLine(final String[] line) {
                return super.isCommentLine(line) || line[0].startsWith("@");
            }

            @Override
            protected SimpleAnnotatedGenomicRegion createRecord(final DataLine dataLine) {
                final Set<String> headersOfInterestPresent = Sets.intersection(headersOfInterest, new HashSet<>(this.columns().names()));
                final Map<String, String> annotationMap = headersOfInterestPresent.stream()
                        .collect(Collectors.toMap(Function.identity(), dataLine::get));
                return new SimpleAnnotatedGenomicRegion( new SimpleInterval(dataLine.get(CONTIG_HEADER), dataLine.getInt(START_HEADER), dataLine.getInt(END_HEADER)),
                        new TreeMap<>(annotationMap));
            }
        }) {
            return reader.toList();
        } catch (final IOException ioe) {
            throw new UserException.CouldNotReadInputFile("Cannot read input file: " + tsvRegionFile.getAbsolutePath(), ioe);
        }
    }

    /**
     * @return list of all of the annotations (not the values) in this annotated region.
     */
    public List<String> getAnnotationNames() {
        return getAnnotations().keySet().stream().sorted().collect(Collectors.toList());
    }
}
