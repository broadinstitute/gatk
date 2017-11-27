package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion.*;

/**
 * Represents a collection of annotated regions.  The annotations do not need to be known ahead of time, if reading from a file.
 * Though {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
 *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
 *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} are always expected for defining the genomic region.
 */
public class SimpleAnnotatedGenomicRegionCollection extends AbstractLocatableCollection<LocatableMetadata, SimpleAnnotatedGenomicRegion> {
    SimpleAnnotatedGenomicRegionCollection(LocatableMetadata metadata, List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions, TableColumnCollection mandatoryColumns, Function<DataLine, SimpleAnnotatedGenomicRegion> recordFromDataLineDecoder, BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> recordToDataLineEncoder) {
        super(metadata, simpleAnnotatedGenomicRegions, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    SimpleAnnotatedGenomicRegionCollection(File inputFile, TableColumnCollection mandatoryColumns, Function<DataLine, SimpleAnnotatedGenomicRegion> recordFromDataLineDecoder, BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
    }

    /**
     *  Reads entire TSV file in one command and stores in RAM.  Please see {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} for defining the genomic region.
     *
     *  A sequence dictionary must be included above the column headers.
     *
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @param headersOfInterest -- should not include any headers that are used to define the region (e.g. contig, start, end)
     * @return annotated regions with one line in the input file for each entry of the list.  Never {@code null}
     */
    public static SimpleAnnotatedGenomicRegionCollection readAnnotatedRegions(final File tsvRegionFile, final Set<String> headersOfInterest) {
        IOUtils.canReadFile(tsvRegionFile);
        Utils.nonNull(headersOfInterest);

        headersOfInterest.remove(CONTIG_HEADER);
        headersOfInterest.remove(START_HEADER);
        headersOfInterest.remove(END_HEADER);

        final Function<DataLine, SimpleAnnotatedGenomicRegion> datalineToRecord = getDataLineToRecordFunction(headersOfInterest);

        final List<String> otherHeaders = new ArrayList<>(headersOfInterest);
        otherHeaders.sort(String::compareTo);

        final BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> recordToDataLine = getRecordToDataLineBiConsumer(otherHeaders);

        return new SimpleAnnotatedGenomicRegionCollection(tsvRegionFile, new TableColumnCollection(Lists.newArrayList(CONTIG_HEADER, START_HEADER, END_HEADER)), datalineToRecord, recordToDataLine);
    }

    private static Function<DataLine, SimpleAnnotatedGenomicRegion> getDataLineToRecordFunction(final Set<String> headersOfInterest) {
        return dataLine -> {

                final Map<String, String> annotationMap = headersOfInterest.stream()
                        .filter(h -> dataLine.columns().contains(h))
                        .collect(Collectors.toMap(Function.identity(), dataLine::get));

                return new SimpleAnnotatedGenomicRegion( new SimpleInterval(dataLine.get(CONTIG_HEADER), dataLine.getInt(START_HEADER), dataLine.getInt(END_HEADER)),
                        new TreeMap<>(annotationMap));
            };
    }

    private static BiConsumer<SimpleAnnotatedGenomicRegion, DataLine> getRecordToDataLineBiConsumer(final List<String> otherHeaders) {
        final List<String> finalHeaders = new ArrayList<>(otherHeaders);
        finalHeaders.remove(CONTIG_HEADER);
        finalHeaders.remove(START_HEADER);
        finalHeaders.remove(END_HEADER);

        return (record, dataLine) -> {

                dataLine.set(CONTIG_HEADER, record.getContig());
                dataLine.set(START_HEADER, record.getStart());
                dataLine.set(END_HEADER, record.getEnd());

                finalHeaders.stream()
                        .filter(h -> dataLine.columns().contains(h))
                        .forEach(h -> dataLine.set(h, record.getAnnotations().getOrDefault(h, "")));
            };
    }

    /**
     * Read in a collection of simple annotated genomic regions from a tsv.
     *
     * Please see {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} for defining the genomic region.  These headers must be present.
     *
     *  Additionally, a sequence dictionary at the top of the file must be present.
     *
     *  *.interval_list is the common extension.
     *
     * @param tsvRegionFile input file. Never {@code null} and must be readable.
     * @return a collection of annotated regions
     */
    public static SimpleAnnotatedGenomicRegionCollection readAnnotatedRegions(final File tsvRegionFile) {
        IOUtils.canReadFile(tsvRegionFile);
        try (final TableReader<SimpleAnnotatedGenomicRegion> reader = new TableReader<SimpleAnnotatedGenomicRegion>(tsvRegionFile) {

                @Override
                protected boolean isCommentLine(final String[] line) {
                    // TODO: Fix magic constant
                    return line.length > 0 && (line[0].startsWith(TableUtils.COMMENT_PREFIX) || line[0].startsWith("@"));
                }
                @Override
                protected SimpleAnnotatedGenomicRegion createRecord(final DataLine dataLine) {
                    // no op
                    return null;
                }
            }) {
            return readAnnotatedRegions(tsvRegionFile, new HashSet<>(reader.columns().names()));
        } catch (final IOException ioe) {
            throw new UserException.CouldNotReadInputFile("Cannot read input file: " + tsvRegionFile.getAbsolutePath(), ioe);
        }
    }

    /** Creates a collection with the same metadata as the given collection, but with the regions specified
     * @param simpleAnnotatedGenomicRegions new regions to use
     * @param collection existing collection to draw the metadata
     * @return a new collection.  Note that it is created with references to the
     */
    public static SimpleAnnotatedGenomicRegionCollection createCollectionFromExistingCollection(final List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions,
                                                                                                final SimpleAnnotatedGenomicRegionCollection collection, final List<String> annotations) {
        return create(simpleAnnotatedGenomicRegions, collection.getMetadata().getSequenceDictionary(), annotations);
    }

    /** Create a collection from a list of annotated regions, a sequence dictionary, and a list of annotations to be included in the collection.
     * The locatable columns {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
     *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
     *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} should not be included.
     *
     * @param simpleAnnotatedGenomicRegions annotated genomic regions
     * @param dictionary a sequence dictionary that includes the contigs in the annotated genomic regions
     * @param annotations annotations of interest that must be present in the annotated regions
     * @return collection of annotated regions
     */
    public static SimpleAnnotatedGenomicRegionCollection create(final List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions, final SAMSequenceDictionary dictionary,
                                                                final List<String> annotations) {
        final List<String> finalColumnList = Lists.newArrayList(SimpleAnnotatedGenomicRegion.CONTIG_HEADER,
                SimpleAnnotatedGenomicRegion.START_HEADER,
                SimpleAnnotatedGenomicRegion.END_HEADER);
        finalColumnList.addAll(annotations);
        final TableColumnCollection annotationColumns = new TableColumnCollection(finalColumnList);
        return new SimpleAnnotatedGenomicRegionCollection(new SimpleLocatableMetadata(dictionary), simpleAnnotatedGenomicRegions, annotationColumns,
                getDataLineToRecordFunction(new HashSet<>(finalColumnList)), getRecordToDataLineBiConsumer(finalColumnList));
    }
}
