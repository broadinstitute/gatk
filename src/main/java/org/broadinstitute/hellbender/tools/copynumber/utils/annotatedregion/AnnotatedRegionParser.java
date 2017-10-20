package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.collections.Lists;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;


public class AnnotatedRegionParser {
    protected final static String CONTIG_HEADER = "CONTIG";
    protected final static String START_HEADER = "START";
    protected final static String END_HEADER = "END";

    /**
     *  Reads entire file in one command and stores in RAM.
     *
     *  Comments can be obtained after reading the file with a call to {@link AnnotatedRegionParser ::getComments}
     *
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @param headersOfInterest -- should not include any headers that are used to define the region (e.g. contig, start, end)
     * @return annotated regions with one line in the input file for each entry of the list.  Never {@code null}
     * @throws IOException
     */
    public static List<SimpleAnnotatedGenomicRegion> readAnnotatedRegions(final File tsvRegionFile, final Set<String> headersOfInterest) {
        try (final TableReader<SimpleAnnotatedGenomicRegion> reader = new TableReader<SimpleAnnotatedGenomicRegion>(tsvRegionFile) {

            @Override
            protected SimpleAnnotatedGenomicRegion createRecord(DataLine dataLine) {
                final Set<String> headersOfInterestPresent = Sets.intersection(headersOfInterest, new HashSet<>(this.columns().names()));
                final Map<String, String> annotationMap = headersOfInterestPresent.stream()
                        .collect(Collectors.toMap(Function.identity(), s -> dataLine.get(s)));
                return new SimpleAnnotatedGenomicRegion( new SimpleInterval(dataLine.get(CONTIG_HEADER), dataLine.getInt(START_HEADER), dataLine.getInt(END_HEADER)),
                        annotationMap);
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
     * @param comments list of comments to prepend to the tsv.  Each line will be prepended with "#"
     */
    public static void writeAnnotatedRegionsAsTsv(final List<SimpleAnnotatedGenomicRegion> regions, final File outputFile,
                                           final List<String> comments) {

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

            for (final String comment : comments) {
                writer.writeComment(comment);
            }
            writer.writeHeaderIfApplies();
            writer.writeAllRecords(regions);
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile("Cannot write file: " + outputFile.getAbsolutePath(), ioe);
        }
    }
}
