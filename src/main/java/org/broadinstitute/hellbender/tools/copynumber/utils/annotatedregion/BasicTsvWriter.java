package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.supercsv.io.CsvMapWriter;
import org.supercsv.io.ICsvMapWriter;
import org.supercsv.prefs.CsvPreference;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

public class BasicTsvWriter {

    private BasicTsvWriter() {
    }

    public static  void writeLines(final String[] headers, final List<Map<String, String>> lines, final File outputFile, final List<String> comments, final CsvPreference pref)  {

        try ( ICsvMapWriter mapWriter = new CsvMapWriter(new FileWriter(outputFile), pref)) {

            for (final String comment: comments) {
                mapWriter.writeComment(comment);
            }

            mapWriter.writeHeader(headers);

            for (final Map<String, String> line: lines) {
                mapWriter.write(line, headers);
            }

        } catch (final IOException ioe) {
            throw new UserException.BadInput("Cannot write the file: " + ioe.getMessage());
        }
    }
}
