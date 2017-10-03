package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Sets;
import org.apache.commons.collections4.list.UnmodifiableList;
import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.supercsv.io.CsvMapReader;
import org.supercsv.io.ICsvMapReader;
import org.supercsv.prefs.CsvPreference;
import org.testng.collections.Lists;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Generally a good idea to create a new instance for each file being read or written.
 */
public class VersatileAnnotatedRegionParser {
    protected final static Set<String> CONTIG_HEADERS = Sets.newHashSet("Chromosome", "contig", "chr", "CONTIG");
    protected final static Set<String> START_HEADERS = Sets.newHashSet("Start_position", "Start_Position", "start", "START");
    protected final static Set<String> END_HEADERS = Sets.newHashSet("End_position", "End_Position", "end", "stop", "END");

    protected final HashCommentMatcher hashCommentMatcher = new HashCommentMatcher();

    /**
     *  Reads entire file in one command and stores in RAM.
     *
     *  Comments can be obtained after reading the file with a call to {@link VersatileAnnotatedRegionParser::getComments}
     *
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @param headersOfInterest -- should not include any headers that are used to define the region (e.g. contig, start, end)
     * @return annotated regions with one line in the input file for each entry of the list.  Never {@code null}
     * @throws IOException
     */
    public List<SimpleAnnotatedGenomicRegion> readAnnotatedRegions(final File tsvRegionFile, final Set<String> headersOfInterest) throws IOException {

        Utils.nonNull(tsvRegionFile);
        Utils.nonNull(headersOfInterest);

        ICsvMapReader mapReader = null;
        final List<SimpleAnnotatedGenomicRegion> result = new ArrayList<>();
        try {
            mapReader = new CsvMapReader(new FileReader(tsvRegionFile), getTsvPreference());

            // the header columns are used as the keys to the Map
            final String[] headers = mapReader.getHeader(true);

            // Create mapping to get the contig, start, and end.
            final Set<String> headerSet = new HashSet<>(Arrays.asList(headers));
            final Set<String> contigHeaders = Sets.intersection(headerSet, CONTIG_HEADERS);
            final Set<String> startHeaders = Sets.intersection(headerSet, START_HEADERS);
            final Set<String> endHeaders = Sets.intersection(headerSet, END_HEADERS);

            validateInputHeader(contigHeaders, startHeaders, endHeaders);

            final Set<String> otherHeaders = Sets.intersection(headerSet, headersOfInterest);
            final String contigHeader = contigHeaders.iterator().next();
            final String startHeader = startHeaders.iterator().next();
            final String endHeader = endHeaders.iterator().next();

            Map<String, String> segmentInfoMap;
            while ((segmentInfoMap = mapReader.read(headers)) != null) {
                final Map<String, String> annotations = segmentInfoMap.entrySet().stream()
                        .filter(e -> otherHeaders.contains(e.getKey()))
                        .collect(Collectors.toMap(e -> e.getKey(), e -> (e.getValue() == null? "" : e.getValue()) ));

                result.add(new SimpleAnnotatedGenomicRegion(
                        new SimpleInterval(segmentInfoMap.get(contigHeader),
                                Integer.parseInt(segmentInfoMap.get(startHeader)),
                                Integer.parseInt(segmentInfoMap.get(endHeader))),
                        annotations));
            }

        } finally {
            if (mapReader != null) {
                mapReader.close();
            }
        }

        return result;
    }

    private CsvPreference getTsvPreference() {
        // Tab-delimited text file with comments on lines starting with "#"
        return new CsvPreference.Builder('"', '\t', "\n")
                .skipComments(hashCommentMatcher).build();
    }

    private void validateInputHeader(Set<String> contigHeaders, Set<String> startHeaders, Set<String> endHeaders) {
        if (contigHeaders.size() != 1) {
            throw new UserException.BadInput("Wrong number of headers that could be used for contig.  Must be one.  Found: " + contigHeaders.toString());
        }
        if (startHeaders.size() != 1) {
            throw new UserException.BadInput("Wrong number of headers that could be used for start.  Must be one.  Found: " + startHeaders.toString());
        }
        if (endHeaders.size() != 1) {
            throw new UserException.BadInput("Wrong number of headers that could be used for end.  Must be one.  Found: " + endHeaders.toString());
        }
    }

    public UnmodifiableList<String> getComments() {
        return hashCommentMatcher.getComments();
    }

    /**
     *  Write a tsv file of annotated regions with the specified strings for the position headers.
     *
     * @param regions list of locatables and associated annotations for writing.  One per line.
     * @param outputFile location of the output file.
     * @param comments list of comments to prepend to the tsv.  Each line will be prepended with "#"
     * @param contigHeader header name to use for the contig
     * @param startHeader header name to use for the start position (inclusive)
     * @param endHeader header name to use for the end position (inclusive)
     */
    public void writeAnnotatedRegionsAsTsv(final List<SimpleAnnotatedGenomicRegion> regions, final File outputFile,
                                           final List<String> comments, final String contigHeader, final String startHeader, final String endHeader) {
        final List<String> headersAsList = Lists.newArrayList(contigHeader, startHeader, endHeader);
        final List<String> additionalAnnotations = new ArrayList<>(regions.get(0).getAnnotations().keySet());
        additionalAnnotations.sort(String::compareTo);
        headersAsList.addAll(additionalAnnotations);

        final List<Map<String, String>> lines = regions.stream()
                .map(r -> VersatileAnnotatedRegionParser.convertSimpleAnnotatedGenomicRegionToStringMap(r, contigHeader, startHeader, endHeader))
                .collect(Collectors.toList());
        BasicTsvWriter.writeLines(headersAsList.toArray(new String[headersAsList.size()]), lines, outputFile, comments, getTsvPreference());
    }

    private static Map<String, String> convertSimpleAnnotatedGenomicRegionToStringMap(final SimpleAnnotatedGenomicRegion region, final String contigHeader, final String startHeader, final String endHeader) {
        final Map<String, String> result = new HashMap<>(region.getAnnotations());
        result.put(contigHeader, region.getContig());
        result.put(startHeader, String.valueOf(region.getStart()));
        result.put(endHeader, String.valueOf(region.getEnd()));
        return result;
    }
}
