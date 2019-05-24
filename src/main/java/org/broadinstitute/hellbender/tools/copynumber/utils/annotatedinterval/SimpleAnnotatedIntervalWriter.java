package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.AnnotatedIntervalCodec;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Path;
import java.util.List;


/**
 * Callers must call {@link #writeHeader} before {@link #add}.
 *
 * This class is not thread-safe.
 */
public class SimpleAnnotatedIntervalWriter implements AnnotatedIntervalWriter {

    private SimpleTableWriter writer;
    private FileWriter fileWriter;
    private String contigColumnHeader;
    private String startColumnHeader;
    private String endColumnHeader;
    private File outputFile;
    private boolean hasHeaderBeenWritten = false;
    private static final Logger logger = LogManager.getLogger(SimpleAnnotatedIntervalWriter.class);

    private class SimpleTableWriter extends TableWriter<AnnotatedInterval> {

        SimpleTableWriter(final Writer writer, TableColumnCollection tableColumns) throws IOException {
            super(writer, tableColumns);
        }

        @Override
        protected void composeLine(final AnnotatedInterval record, final DataLine dataLine) {
            // First the Locatable info
            dataLine.set(contigColumnHeader, record.getContig());
            dataLine.set(startColumnHeader, record.getStart());
            dataLine.set(endColumnHeader, record.getEnd());

            // Now everything else.
            record.getAnnotations().keySet().forEach(k -> dataLine.set(k, record.getAnnotationValue(k)));
        }
    }

    /**
     * Initialize this writer to the given output file.
     *
     * @param outputFile destination file.  Must be writeable.
     */
    public SimpleAnnotatedIntervalWriter(final File outputFile) {
        IOUtil.assertFileIsWritable(outputFile);
        this.outputFile = outputFile;
    }

    private void initializeForWriting(final String contigColumnName, final String startColumnName, final String endColumnName, final List<String> annotations) {
        final List<String> finalColumnList = Lists.newArrayList(contigColumnName, startColumnName, endColumnName);
        finalColumnList.addAll(annotations);
        try {
            fileWriter = new FileWriter(outputFile);

            // By initializing writer to be based on fileWriter, writer.close will close the fileWriter as well.
            writer = new SimpleTableWriter(fileWriter, new TableColumnCollection(finalColumnList));
        } catch (final IOException ioe) {
            throw new GATKException("Could not create: " + outputFile.getAbsolutePath(), ioe);
        }

        this.contigColumnHeader = contigColumnName;
        this.startColumnHeader = startColumnName;
        this.endColumnHeader = endColumnName;
    }

    /**
     * @param annotatedIntervalHeader Header to write.  Recommended to create the header using
     *  {@link AnnotatedIntervalCodec#createHeaderForWriter(Path, List, SAMFileHeader)} or one of its overloaded
     *  versions.  Never {@code null}
     */
    @Override
    public void writeHeader(final AnnotatedIntervalHeader annotatedIntervalHeader) {
        Utils.nonNull(annotatedIntervalHeader);
        if (!hasHeaderBeenWritten) {

            final String contigColumnName = annotatedIntervalHeader.getContigColumnName();
            final String startColumnName = annotatedIntervalHeader.getStartColumnName();
            final String endColumnName = annotatedIntervalHeader.getEndColumnName();

            initializeForWriting(contigColumnName, startColumnName, endColumnName, annotatedIntervalHeader.getAnnotations());
            try {
                final SAMFileHeader samFileHeader = annotatedIntervalHeader.getSamFileHeader();
                if (samFileHeader != null) {
                    fileWriter.write(samFileHeader.getSAMString());
                }
                writer.writeHeaderIfApplies();
                hasHeaderBeenWritten = true;
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile("Could not write to file.", e);
            }
        } else {
            logger.warn("Attempted to write header more than once.  Ignoring this request.");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() {
        try {
            // Only close the SimpleTableWriter, since it will close fileWriter
            writer.close();
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not close file writing.", e);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void add(final AnnotatedInterval annotatedInterval) {
        if (!hasHeaderBeenWritten) {
            throw new GATKException.ShouldNeverReachHereException("This is an invalid code path, since the header of the output file should already be written.  Please post this error to the GATK forum (https://gatkforums.broadinstitute.org/gatk)");
        }
        try {
            writer.writeRecord(annotatedInterval);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not write to file.", e);
        }
    }
}
