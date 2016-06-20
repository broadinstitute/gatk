package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Table writer for target tables.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetWriter extends TableWriter<Target> {

    private List<TargetAnnotation> annotations;

    /**
     * Creates a new writer given the file name with no annotations.
     * @param file the destination file.
     * @throws IOException if such an exception is throw when creating or writing in the file provided.
     */
    public TargetWriter(final File file) throws IOException {
        super(file, composeTableColumns(Collections.emptySet()));
        this.annotations = Collections.emptyList();
    }

    /**
     * Creates a new writer given the file name and the set of annotations to be output.
     * @param file the destination file.
     * @param annotations the output annotations.
     * @throws IOException if such an exception is throw when creating or writing in the file provided.
     */
    public TargetWriter(final File file, final Set<TargetAnnotation> annotations) throws IOException {
        super(file, composeTableColumns(annotations));
        this.annotations = new ArrayList<>(annotations);
    }

    private static TableColumnCollection composeTableColumns(final Set<TargetAnnotation> annotations) {
        Utils.nonNull(annotations);

        final List<TargetTableColumn> columns = new ArrayList<>(4 + annotations.size());
        columns.add(TargetTableColumn.CONTIG);
        columns.add(TargetTableColumn.START);
        columns.add(TargetTableColumn.END);
        columns.add(TargetTableColumn.NAME);
        annotations.stream()
                .map(annotation -> annotation.column)
                .forEach(columns::add);

        return new TableColumnCollection(columns.toArray());
    }

    /**
     * {@inheritDoc}
     * @throws IllegalArgumentException if the record does not have a value
     *  for any of the annotations output by this writer.
     * dataLine the destination data-line object.
     */
    @Override
    protected void composeLine(final Target record, final DataLine dataLine) {
        dataLine.append(record.getContig())
                .append(record.getStart())
                .append(record.getEnd())
                .append(record.getName());
        final TargetAnnotationCollection annotationCollection = record.getAnnotations();
        for (final TargetAnnotation annotation : annotations) {
            dataLine.append(annotationCollection.get(annotation));
        }
    }

    /**
     * Utility method that wraps the try. . . catch block of writing targets to file.
     *
     * @param outputFile    file for output
     * @param targets       targets to write to file
     */
    public static void writeTargetsToFile(final File outputFile, final Collection<Target> targets) {
        try (final TargetWriter writer = new TargetWriter(outputFile)) {
            writer.writeAllRecords(targets);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("Could not write output targets file", ex);
        }
    }
}
