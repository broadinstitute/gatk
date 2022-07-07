package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;

/**
 *
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReferenceProgramGroup.class
)
public class CompareReferences extends GATKTool {

    public static final String MISSING_ENTRY = "---";
    public static final String MD5_COLUMN_NAME = "MD5";
    public static final String LENGTH_COLUMN_NAME = "Length";

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "", optional = true)
    private GATKPath output;

    @Argument(fullName = "md5-calculation-mode", shortName = "md5-calculation-mode", doc = "", optional = true)
    private MD5CalculationMode md5CalculationMode = MD5CalculationMode.USE_DICT;

    public enum MD5CalculationMode {
          USE_DICT,
          RECALCULATE_IF_MISSING,
          ALWAYS_RECALCULATE;
    }

    private Map<GATKPath, ReferenceDataSource> referenceSources;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        // add data source for -R reference
        referenceSources = new LinkedHashMap<>();
        referenceSources.put(getReferencePath(), directlyAccessEngineReferenceDataSource());

        // add data sources for remaining references
        for(GATKPath path : references){
            referenceSources.put(path, ReferenceDataSource.of(path.toPath()));
        }
    }

    @Override
    public void traverse(){
        ReferenceSequenceTable table = new ReferenceSequenceTable(referenceSources, md5CalculationMode);
        table.build();

        if(output == null){
            writeTableToStdOutput(table);
        }
        else{
            writeTableToFileOutput(table);
        }
    }

    private void writeTableToStdOutput(ReferenceSequenceTable table){
        // print header --> change to use indices to get rid of tab at the end
        for(String str : table.getColumnNames()){
            System.out.print(str + "\t");
        }
        System.out.println();

        // use string format to output as a table
        for(ReferenceSequenceTable.TableRow row : table){
            String currMd5 = row.getMd5();
            System.out.printf("%s\t%d", currMd5, row.getLength());
            for(ReferenceSequenceTable.TableEntry currEntry : row.getEntries()){
                if(currEntry == null){
                    System.out.print("\t" + MISSING_ENTRY);
                }
                else{
                    System.out.printf("\t%s", currEntry.getSequenceName());
                }
            }
            System.out.println();
        }
    }

    private void writeTableToFileOutput(ReferenceSequenceTable table) {
        TableColumnCollection columns = new TableColumnCollection(table.getColumnNames());
        try(CompareReferences.CompareReferencesOutputTableWriter writer = new CompareReferences.CompareReferencesOutputTableWriter(output.toPath(), columns)){
            writer.writeHeaderIfApplies();
            for(ReferenceSequenceTable.TableRow row : table){
                writer.writeRecord(row);
            }
        }
        catch(IOException exception){
            throw new UserException.CouldNotCreateOutputFile(output, "Failed to write output table.", exception);
        }
    }

    protected static String getReferenceDisplayName(GATKPath reference){
        return reference.toPath().getFileName().toString();
    }

    @Override
    public Object onTraversalSuccess() {
        return null;
    }

    @Override
    public void closeTool() {
        for(Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()){
            if(!entry.getKey().equals(getReferencePath())){
                entry.getValue().close();
            }
        }
    }


    public static class CompareReferencesOutputTableWriter extends TableWriter<ReferenceSequenceTable.TableRow> {
        private TableColumnCollection columnCollection;

        public CompareReferencesOutputTableWriter(final Path table, TableColumnCollection columns) throws IOException {
            super(table, columns);
            columnCollection = columns;
        }

        @Override
        protected void composeLine(final ReferenceSequenceTable.TableRow record, final DataLine dataLine) {
            // change this !!! (TableRow getColumns() method and iterate over that directly)
            /*List<String> columnNames = record.getColumnNames();
            List<TableEntry> entries = record.getEntries();
            for(int i = 0; i < columnNames.size(); i++){
                dataLine.set(columnNames.get(i), entries.get(i) == null ? MISSING_ENTRY : entries.get(i).getSequenceName());
            }*/

           dataLine.set(MD5_COLUMN_NAME, record.getMd5())
                    .set(LENGTH_COLUMN_NAME, record.getLength());

            List<String> columnNames = columnCollection.names();
            List<ReferenceSequenceTable.TableEntry> entries = record.getEntries();
            for (int i = 2; i < columnNames.size(); i++) {
                dataLine.set(columnNames.get(i), entries.get(i-2) == null ? MISSING_ENTRY : entries.get(i-2).getSequenceName());
            }

        }
    }
}
