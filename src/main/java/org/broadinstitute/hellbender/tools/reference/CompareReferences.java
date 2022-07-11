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

        analyzeTable(table);
    }

    private void writeTableToStdOutput(ReferenceSequenceTable table){
        // print header
        List<String> columnNames = table.getColumnNames();
        for(int i = 0 ; i < columnNames.size(); i++){
            if(i == 0){
                System.out.print(columnNames.get(i));
            }
            else{
                System.out.print("\t" + columnNames.get(i));
            }
        }
        System.out.println();

        // use string format to output as a table
        for(ReferenceSequenceTable.TableRow row : table){
            for(ReferenceSequenceTable.TableEntry currEntry : row.getEntries()){
                System.out.printf("%s\t", currEntry.getColumnValue());
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

    public void analyzeTable(ReferenceSequenceTable table) {
        boolean noMissingEntries = true;
        String output = "";
        int currRow = 0;

        for(ReferenceSequenceTable.TableRow row : table){
            currRow++;
            ReferenceSequenceTable.TableEntry[] entries = row.getEntries();

            for(int i = 0; i < entries.length; i++){
                if(entries[i].getColumnValue().equals(ReferenceSequenceTable.MISSING_ENTRY)){
                    output += String.format("Row %d: Missing entry in %s column.\n", currRow, entries[i].getColumnName());
                    noMissingEntries = false;
                }
            }
        }

        if(noMissingEntries){
            output += "No missing entries. References are interchangeable.\n";
        }
        else{
            output += "References are not an exact match. See table output for details.\n";
        }
        System.out.println(output);
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
            List<String> columnNames = record.getColumnNames();
            ReferenceSequenceTable.TableEntry[] entries = record.getEntries();
            for(int i = 0; i < columnNames.size(); i++){
                dataLine.set(entries[i].getColumnName(), entries[i].getColumnValue());
            }
        }
    }
}
