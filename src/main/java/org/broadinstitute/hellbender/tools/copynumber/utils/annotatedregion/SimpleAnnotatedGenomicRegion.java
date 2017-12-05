package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

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

    private final SimpleInterval interval;
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
     *  Reads entire TSV file in one command and stores in RAM.  Please see {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} for defining the genomic region.
     *
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @param headersOfInterest -- should not include any headers that are used to define the region (e.g. contig, start, end)
     * @return annotated regions with one line in the input file for each entry of the list.  Never {@code null}
     */
    public static List<SimpleAnnotatedGenomicRegion> readAnnotatedRegions(final File tsvRegionFile, final Set<String> headersOfInterest) {
        try (final TableReader<SimpleAnnotatedGenomicRegion> reader = new TableReader<SimpleAnnotatedGenomicRegion>(tsvRegionFile) {

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
     *  Write a tsv file of annotated regions with the specified strings for the position headers.
     *
     * @param regions list of locatables and associated annotations for writing.  One per line.
     * @param outputFile location of the output file.
     */
    public static void writeAnnotatedRegionsAsTsv(final List<SimpleAnnotatedGenomicRegion> regions, final File outputFile) {

        // Can't do Arrays.asList, since we will be adding to this list.
        final List<String> columnHeaders = Lists.newArrayList(CONTIG_HEADER, START_HEADER, END_HEADER);
        final List<String> otherHeaders = new ArrayList<>(regions.get(0).getAnnotations().keySet());
        otherHeaders.sort(String::compareTo);

        columnHeaders.addAll(otherHeaders);
        TableColumnCollection outputColumns = new TableColumnCollection(columnHeaders);

        try (final TableWriter<SimpleAnnotatedGenomicRegion> writer = new TableWriter<SimpleAnnotatedGenomicRegion>(outputFile, outputColumns) {
            @Override
            protected void composeLine(SimpleAnnotatedGenomicRegion record, DataLine dataLine) {
                dataLine.set(CONTIG_HEADER, record.getContig());
                dataLine.set(START_HEADER, record.getStart());
                dataLine.set(END_HEADER, record.getEnd());
                otherHeaders.forEach(h -> dataLine.set(h, record.getAnnotations().getOrDefault(h, "")));
            }
        }) {

            writer.writeHeaderIfApplies();
            writer.writeAllRecords(regions);
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile("Cannot write file: " + outputFile.getAbsolutePath(), ioe);
        }
    }
}
