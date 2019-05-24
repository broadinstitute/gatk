package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalHeader;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Properties;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Read AnnotatedIntervals from a xsv file (see {@link XsvLocatableTableCodec}.
 *
 * This class also depends on the config file detailed in {@link XsvLocatableTableCodec}.
 *
 * For writing, see {@link AnnotatedIntervalCodec#createHeaderForWriter(Path, List, SAMFileHeader)}.  In the case of
 *  config files that have multiple columns specified (comma delimited), the first value will be selected as the one to
 *  use for output.
 *
 */
public class AnnotatedIntervalCodec extends AsciiFeatureCodec<AnnotatedInterval> {
    public static final String ANNOTATED_INTERVAL_DEFAULT_CONFIG_RESOURCE = "org/broadinstitute/hellbender/tools/copynumber/utils/annotatedinterval/annotated_region_default.config";

    public static final String CONTIG_COL_COMMENT = "_ContigHeader=";
    public static final String START_COL_COMMENT = "_StartHeader=";
    public static final String END_COL_COMMENT = "_EndHeader=";

    private Path configFilePath;
    private XsvLocatableTableCodec xsvLocatableTableCodec;
    private AnnotatedIntervalHeader header;

    /** Use the default config for tribble. */
    public AnnotatedIntervalCodec() {
        super(AnnotatedInterval.class);
        final String resourcePath = ANNOTATED_INTERVAL_DEFAULT_CONFIG_RESOURCE;
        try {
            configFilePath = Resource.getResourceContentsAsFile(resourcePath).toPath();
            xsvLocatableTableCodec = new XsvLocatableTableCodec(configFilePath);
        } catch (final IOException ioe) {
            throw new GATKException.ShouldNeverReachHereException("Could not read config file: " + resourcePath,
                    ioe);
        }
    }

    public AnnotatedIntervalCodec(final Path configFilePath) {
        super(AnnotatedInterval.class);
        this.configFilePath = configFilePath;
        xsvLocatableTableCodec = new XsvLocatableTableCodec(configFilePath);
    }

    @Override
    public AnnotatedInterval decode(final String s) {
        final XsvTableFeature feature = xsvLocatableTableCodec.decode(s);
        if (feature == null) {
            return null;
        }

        final List<String> annotationFields = header.getAnnotations();
        final SortedMap<String, String> annotations = new TreeMap<>();
        IntStream.range(0, annotationFields.size()).boxed()
                .forEach(i -> annotations.put(annotationFields.get(i), feature.getValuesWithoutLocationColumns().get(i)));

        return new AnnotatedInterval(
                new SimpleInterval(feature.getContig(), feature.getStart(), feature.getEnd()),
                annotations);
    }

    @Override
    public AnnotatedIntervalHeader readActualHeader(final LineIterator reader) {
        xsvLocatableTableCodec.readActualHeader(reader);
        header = new AnnotatedIntervalHeader(xsvLocatableTableCodec.getContigColumn(), xsvLocatableTableCodec.getStartColumn(),
                xsvLocatableTableCodec.getEndColumn(), xsvLocatableTableCodec.getHeaderWithoutLocationColumns(),
                xsvLocatableTableCodec.renderSamFileHeader());
        return header;
    }

    @Override
    public boolean canDecode(final String pathString) {
        return (pathString.endsWith(".seg") || pathString.endsWith(".maf") || pathString.endsWith(".maf.annotated")) && xsvLocatableTableCodec.canDecodeFileChecks(configFilePath.toUri().toString(), pathString);
    }

    /**
     * Create an annotated interval header based on a config file (for locatable field names only) and a list of annotations (the rest of the fields).
     *
     * @param outputConfigFile config path for determining the locatable column headers in the output file.  If comma separated lists are present as values, the first entry will be chosen.
     *                         If the config file contains numeric indexes as output column name, then an exception is thrown.  Never {@code null}.
     * @param annotations  Names of the annotations to render.  If any of the locatable columns are in the annotation, those columns will be removed from the annotations list in the header.
     *                     Never {@code null}.
     * @param samFileHeader SAM FileHeader to prepend to the data.  {@code null} is allowed.
     * @return a header that can be used in an AnnotatedFileWriter.  Structured comments will be updated here.  Never {@code null}.
     */
    public static AnnotatedIntervalHeader createHeaderForWriter(final Path outputConfigFile, final List<String> annotations, final SAMFileHeader samFileHeader) {

        Utils.nonNull(annotations);
        Utils.nonNull(outputConfigFile);

        final Pair<Boolean, Properties> validityAndPropertiesPair = XsvLocatableTableCodec.getAndValidateConfigFileContentsOnPath(outputConfigFile, true);
        final boolean                   isValid                   = validityAndPropertiesPair.getLeft();
        final Properties                headerNameProperties      = validityAndPropertiesPair.getRight();
        if ( !isValid ) {
            throw new UserException.BadInput("Error: invalid configuration file given: " + outputConfigFile.toUri().toString());
        }
        final String                    contigColumnName          = determineOutputColumnFromList(headerNameProperties.getProperty(DataSourceUtils.CONFIG_FILE_FIELD_NAME_CONTIG_COLUMN));
        final String                    startColumnName           = determineOutputColumnFromList(headerNameProperties.getProperty(DataSourceUtils.CONFIG_FILE_FIELD_NAME_START_COLUMN));
        final String                    endColumnName             = determineOutputColumnFromList(headerNameProperties.getProperty(DataSourceUtils.CONFIG_FILE_FIELD_NAME_END_COLUMN));

        XsvLocatableTableCodec.validateLocatableColumnName(contigColumnName);
        XsvLocatableTableCodec.validateLocatableColumnName(startColumnName);
        XsvLocatableTableCodec.validateLocatableColumnName(endColumnName);

        final List<String> finalAnnotations = annotations.stream()
                .filter(a -> !a.equals(contigColumnName) && !a.equals(startColumnName) && !a.equals(endColumnName))
                .collect(Collectors.toList());

        //TODO: Test NIO?  See https://github.com/broadinstitute/gatk/issues/4579
        final List<String> commentsToWrite = AnnotatedIntervalCodec.updateStructuredComments(contigColumnName, startColumnName, endColumnName, samFileHeader.getComments());

        // A bit more manual to write the SAM Header
        final SAMFileHeader finalSamHeader = samFileHeader.clone();
        finalSamHeader.setComments(commentsToWrite);

        return new AnnotatedIntervalHeader(contigColumnName, startColumnName, endColumnName, finalAnnotations, finalSamHeader);
    }

    private static String determineOutputColumnFromList(final String rawColumnNamesAsString) {
        final String result = StringUtils.split(rawColumnNamesAsString, ",")[0];
        if (StringUtils.isNumeric(result)) {
            throw new UserException.BadInput("Index numbers were used in output configuration.  This is not allowed: " + rawColumnNamesAsString);
        }
        return result;
    }

    /**
     *  See {@link #createHeaderForWriter(Path, List, SAMFileHeader)}
     *
     *  This will use the default headers for annotated regions.  Call this method if no config file is available.
     *
     * @param annotations Annotations that should be used in the header.  Never {@code null}.
     * @param samFileHeader SAM File header to use for this header.  {@code null} is allowed.
     * @return A header to be used for output.  Never {@code null}
     */
    public static AnnotatedIntervalHeader createHeaderForWriter(final List<String> annotations, final SAMFileHeader samFileHeader) {
        Utils.nonNull(annotations);

        try {
            final Path resourceFile = Resource.getResourceContentsAsFile(ANNOTATED_INTERVAL_DEFAULT_CONFIG_RESOURCE).toPath();
            return createHeaderForWriter(resourceFile, annotations, samFileHeader);
        } catch (final IOException ioe) {
            throw new GATKException.ShouldNeverReachHereException("Could not load the default config file for annotated intervals.", ioe);
        }
    }

    /**
     * Update the comments to represent the new field names for the locatable information.
     *
     * @param finalContigColumnName Never {@code null} and never blank.
     * @param finalStartColumnName Never {@code null} and never blank.
     * @param finalEndColumnName Never {@code null} and never blank.
     * @param currentComments Never {@code null}.  Can be empty.
     * @return list of the comments.  Never {@code null} and never empty.
     */
    private static List<String> updateStructuredComments(final String finalContigColumnName, final String finalStartColumnName, final String finalEndColumnName, final List<String> currentComments) {

        Utils.nonNull(currentComments);
        Utils.validateArg(!StringUtils.isEmpty(finalContigColumnName), "Contig column name was null or blank, which is invalid.");
        Utils.validateArg(!StringUtils.isEmpty(finalStartColumnName), "Start column name was null or blank, which is invalid.");
        Utils.validateArg(!StringUtils.isEmpty(finalEndColumnName), "End column name was null or blank, which is invalid.");

        // Remove old structured comments, if present.
        final List<String> commentsToWrite = currentComments.stream()
                .filter(c -> !c.startsWith(SAMTextHeaderCodec.COMMENT_PREFIX + CONTIG_COL_COMMENT))
                .filter(c -> !c.startsWith(SAMTextHeaderCodec.COMMENT_PREFIX + START_COL_COMMENT))
                .filter(c -> !c.startsWith(SAMTextHeaderCodec.COMMENT_PREFIX + END_COL_COMMENT)).collect(Collectors.toList());

        // Write out the column headers as a comment
        commentsToWrite.add(CONTIG_COL_COMMENT + finalContigColumnName);
        commentsToWrite.add(START_COL_COMMENT + finalStartColumnName);
        commentsToWrite.add(END_COL_COMMENT + finalEndColumnName);
        return commentsToWrite;
    }
}
