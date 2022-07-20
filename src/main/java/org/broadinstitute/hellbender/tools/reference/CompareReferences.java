package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
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
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Path;
import java.util.*;

/**
 * Display reference comparison as a tab-delimited table and summarize reference differences.
 *
 * <p>This tool generates a MD5-keyed table comparing specified references and does an analysis to summarize the differences
 * between the references provided. Comparisons are made against a "special" reference, specified with the -R
 * argument. Subsequent references to be compared may be specified using the -references-to-compare argument.
 * The table can be directed to a file or standard output using provided command-line arguments.
 * A supplementary table keyed by sequence name can be displayed using the -display-sequences-by-name argument;
 * to display only sequence names for which the references are not consistent, run with the -display-only-differing-sequences
 * argument as well.</p>
 *
 * <h3>Input</h3>
 * <p>
 * The references and MD5 calculation mode.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A TSV file or output stream displaying the MD5-keyed table comparison, and a table analysis displayed to standard output.
 * <pre>
 * MD5	Length	hg19mini.fasta	hg19mini_1renamed.fasta
 * 8c0c38e352d8f3309eabe4845456f274	16000	1	chr1
 * 5f8388fe3fb34aa38375ae6cf5e45b89	16000	2	2
 * 94de808a3a2203dbb02434a47bd8184f	16000	3	3
 * 7d397ee919e379328d8f52c57a54c778	16000	4	4
 *
 * REFERENCE PAIR: hg19mini.fasta, hg19mini_1renamed.fasta
 * Status:
 * 	DIFFER_IN_SEQUENCE_NAMES
 * </pre>
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CompareReferences \
 *   -R reference1.fasta \
 *   -ref-comp reference2.fasta
 *   -O output.fasta \
 *   -md5-calculation-mode USE_DICT
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = "Compare multiple references and output a tab-delimited table detailing what the differences are and a summarized analysis of each pair of references.",
        oneLineSummary = "Display reference comparison as a tab-delimited table and summarize reference differences.",
        programGroup = ReferenceProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public class CompareReferences extends GATKTool {

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    /**
     * Output file will be written here.
     *
     * Note: If no output file provided, table will print to standard output.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "If provided, file to output table.", optional = true)
    private GATKPath output;

    @Argument(fullName = "md5-calculation-mode", shortName = "md5-calculation-mode", doc = "MD5CalculationMode indicating method of MD5 calculation.", optional = true)
    private MD5CalculationMode md5CalculationMode = MD5CalculationMode.RECALCULATE_IF_MISSING;

    @Argument(fullName = "display-sequences-by-name", doc = "If provided, the table by sequence name will be printed.", optional = true)
    private boolean displaySequencesByName = false;

    @Argument(fullName = "display-only-differing-sequences", doc = "If provided, only display sequence names that differ in their actual sequence.", optional = true)
    private boolean onlyDisplayDifferingSequences = false;

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
        logger.warn("Building reference sequence table.");
        table.build();
        logger.warn("Finished building table.");

        writeTable(table);

        if(displaySequencesByName){
            writeTableBySequenceName(table);
        }

        logger.warn("Analyzing table.");
        List<ReferencePair> referencePairs = table.analyzeTable();
        logger.warn("Finished analyzing table.");
        System.out.println("*********************************************************");
        for(ReferencePair pair : referencePairs){
            System.out.println(pair);
        }
    }


    /**
     * Given a table, write table to output.
     *
     * Note: if no output file specified, displays to standard output.
     *
     * @param table
     */
    private void writeTable(ReferenceSequenceTable table){
        TableColumnCollection columns = new TableColumnCollection(table.getColumnNames());
        try(CompareReferences.CompareReferencesOutputTableWriter writer = output == null
            ? new CompareReferences.CompareReferencesOutputTableWriter(new OutputStreamWriter(System.out), columns)
            : new CompareReferences.CompareReferencesOutputTableWriter(output.toPath(), columns)
        ){
            writer.writeHeaderIfApplies();
            for(ReferenceSequenceTable.TableRow row : table){
                writer.writeRecord(row);
            }
        }
        catch(IOException exception){
             throw (output == null)
             ? new UserException.CouldNotCreateOutputFile("System.out", "Failed to write output table.", exception)
             : new UserException.CouldNotCreateOutputFile(output, "Failed to write output table.", exception);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        return null;
    }

    /**
     * Given a table, write table by sequence name to standard output
     *
     * @param table
     */
    public void writeTableBySequenceName(ReferenceSequenceTable table){
        List<String> output = new ArrayList<>();
        output.add("Sequence \tMD5 \tReference\n");

        for(String sequenceName : table.getAllSequenceNames()){
            Set<ReferenceSequenceTable.TableRow> rows = table.queryBySequenceName(sequenceName);
            if(onlyDisplayDifferingSequences) {
                if(rows.size() > 1) {
                    output.addAll(displayBySequenceName(rows, sequenceName));
                }
            } else {
                output.addAll(displayBySequenceName(rows, sequenceName));
            }
        }

        for(String str : output){
            System.out.print(str);
        }
        System.out.println();
    }

    private List<String> displayBySequenceName(Set<ReferenceSequenceTable.TableRow> rows, String sequenceName){
        List<String> output = new ArrayList<>();
        output.add(sequenceName);
        for(ReferenceSequenceTable.TableRow row : rows) {
            ReferenceSequenceTable.TableEntry[] entries = row.getEntries();
            output.add("\n\t" + row.getMd5() + "\t");

            for(int i = 2; i < entries.length; i++) {
                if(entries[i].getColumnValue().equals(sequenceName)) {
                    output.add(entries[i].getColumnName() + "\t");
                }
            }
        }
        output.add("\n");
        return output;
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

        public CompareReferencesOutputTableWriter(final Path table, TableColumnCollection columns) throws IOException {
            super(table, columns);
        }

        public CompareReferencesOutputTableWriter(final Writer writer, TableColumnCollection columns) throws IOException {
            super(writer, columns);
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
