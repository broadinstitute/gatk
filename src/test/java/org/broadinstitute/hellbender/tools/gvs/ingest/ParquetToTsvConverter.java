package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.example.data.Group;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.hadoop.example.GroupReadSupport;
import org.apache.parquet.schema.GroupType;
import org.apache.parquet.schema.PrimitiveType;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

/**
 * Test utility for converting Parquet files to TSV format so that Parquet output can be
 * compared against TSV golden files.
 *
 * <p>Conversion rules:
 * <ul>
 *   <li>Columns are written in the caller-specified order, which may differ from the
 *       field order in the Parquet schema (e.g. ref_ranges golden files use
 *       {@code location, sample_id, length, state} while the schema orders them
 *       {@code sample_id, location, length, state}).</li>
 *   <li>INT64 Parquet values are converted to their decimal string representation.</li>
 *   <li>BINARY/UTF8 Parquet values are passed through as-is.</li>
 *   <li>Optional Parquet fields that are absent from a record are written as empty strings.</li>
 *   <li>Lines are terminated with {@code \n} (not {@code System.lineSeparator()}) to keep
 *       golden-file comparisons platform-independent.</li>
 * </ul>
 */
public class ParquetToTsvConverter {

    /**
     * Reads a Parquet file and writes its contents as a tab-separated TSV to {@code outputTsvFile}.
     * A header row containing the given {@code columns} is written first, followed by one data row
     * per Parquet record, preserving the original record order.
     *
     * @param parquetFile   the Parquet file to read
     * @param columns       column names to include in the TSV output, in the desired column order;
     *                      every name must exist as a field in the Parquet file's schema
     * @param outputTsvFile destination TSV file (will be created or overwritten)
     * @throws IOException if an error occurs reading the Parquet file or writing the TSV
     */
    public static void convert(File parquetFile, List<String> columns, File outputTsvFile) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(outputTsvFile.toPath(), StandardCharsets.UTF_8);
             ParquetReader<Group> reader = ParquetReader
                     .builder(new GroupReadSupport(), new Path(parquetFile.toURI()))
                     .withConf(new Configuration())
                     .build()) {

            // Write header — use \n to match CSVWriter's DEFAULT_LINE_END
            writer.write(String.join("\t", columns));
            writer.write("\n");

            Group record;
            while ((record = reader.read()) != null) {
                final GroupType schema = record.getType();
                final List<String> row = new ArrayList<>(columns.size());

                for (final String col : columns) {
                    final int fieldIndex = schema.getFieldIndex(col);

                    if (record.getFieldRepetitionCount(fieldIndex) == 0) {
                        // Absent optional field → empty string
                        row.add("");
                    } else {
                        final PrimitiveType.PrimitiveTypeName typeName =
                                schema.getType(fieldIndex).asPrimitiveType().getPrimitiveTypeName();
                        if (typeName == PrimitiveType.PrimitiveTypeName.INT64) {
                            row.add(String.valueOf(record.getLong(fieldIndex, 0)));
                        } else {
                            // BINARY (UTF8) and any other type: use string representation
                            row.add(record.getString(fieldIndex, 0));
                        }
                    }
                }

                writer.write(String.join("\t", row));
                writer.write("\n");
            }
        }
    }
}

