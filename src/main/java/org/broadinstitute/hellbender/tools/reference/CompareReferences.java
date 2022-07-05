package org.broadinstitute.hellbender.tools.reference;

import com.google.flatbuffers.Table;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.readorientation.AltSiteRecord;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
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

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "", optional = true)
    private GATKPath output;

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
    public void traverse() {
        Map<String, TableRow> tableOutput = buildTable();

        if(output == null){
            writeTableToStdOutput(tableOutput);
        }
        else{

            writeTableToFileOutput(tableOutput);
        }
    }

    private Map<String, TableRow> buildTable(){
        Map<String, TableRow> table = new LinkedHashMap<>();

        for(Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()){
            SAMSequenceDictionary dictionary = entry.getValue().getSequenceDictionary();
            for(SAMSequenceRecord record : dictionary.getSequences()){
                String name = record.getSequenceName();
                int length = record.getSequenceLength();
                String md5 = record.getMd5();
                GATKPath reference = entry.getKey();
                TableEntry newEntry = new TableEntry(name, length, md5, reference);
                TableRow newRow = null;

                // map each MD5 to List of TableEntry objects containing length, md5, and name
                if (!table.containsKey(md5)) {
                    newRow = new TableRow(md5, length, new ArrayList<>(referenceSources.keySet()));
                    table.put(md5, newRow);
                }
                table.get(md5).add(newEntry);
            }
        }
        return table;
    }

    private void writeTableToStdOutput(Map<String, TableRow> table){
        // MD5  Length  Ref1  Ref2 ...
        System.out.printf("%s\t%s\t", "MD5", "Length");
        for(GATKPath file : referenceSources.keySet()){
            System.out.printf("%s\t", getReferenceDisplayName(file));
        }
        System.out.println();

        // use string format to output as a table
        for(TableRow row : table.values()){
            String currMd5 = row.getMd5();
            System.out.printf("%s\t%d", currMd5, row.getLength());
            for(TableEntry currEntry : row.getEntries()){
                if(currEntry == null){
                    System.out.print("\t" + MISSING_ENTRY);
                }
                else{
                    System.out.printf("\t%s", currEntry.sequenceName);
                }
            }
            System.out.println();
        }
    }

    private void writeTableToFileOutput(){


    }

    private String getReferenceDisplayName(GATKPath reference){
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

    private static class TableEntry {
        private final String sequenceName;
        private final int sequenceLength;
        private final String md5;
        private final GATKPath reference;

        public TableEntry(String name, int length, String md5, GATKPath reference){
            this.md5 = md5;
            sequenceLength = length;
            this.reference = reference;
            sequenceName = name;
        }

        public String getSequenceName() {
            return sequenceName;
        }

        public int getSequenceLength() {
            return sequenceLength;
        }

        public String getMd5() {
            return md5;
        }

        public GATKPath getReference() {
            return reference;
        }
    }

    private static class TableRow {

        private final String md5;
        private final List<TableEntry> entries;
        private final int length;
        private final Map<GATKPath, Integer> columnIndices;

        public TableRow(String md5, int length, List<GATKPath> references){
            this.md5 = md5;
            this.length = length;
            entries = new ArrayList<>(references.size());
            for(int i = 0; i < references.size(); i++){
                entries.add(null);
            }
            columnIndices = new HashMap<>();
            for(int i = 0; i < references.size(); i++){
                columnIndices.put(references.get(i), i);
            }
        }

        public void add(TableEntry entry){
            int idx = columnIndices.get(entry.getReference());
            entries.set(idx, entry);
        }

        public String getMd5() {
            return md5;
        }

        public List<TableEntry> getEntries() {
            return entries;
        }

        public int getLength() {
            return length;
        }
    }

    public static class CompareReferencesOutputTableWriter extends TableWriter<TableRow> {
        public CompareReferencesOutputTableWriter(final Path output, TableRow row) throws IOException {

            super(output, row.getEntries());
            //writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
        }

        @Override
        protected void composeLine(final AltSiteRecord record, final DataLine dataLine) {
            // it'd be nice to set() less manually...
            // Note that allele fraction f is not allele-specific, thus the same f array will be printed
            // four times for each context
            dataLine.set(AltSiteRecord.AltSiteRecordTableColumn.CONTEXT.toString(), record.getReferenceContext())
                    .set(AltSiteRecord.AltSiteRecordTableColumn.REF_COUNT.toString(), record.getRefCount())
                    .set(AltSiteRecord.AltSiteRecordTableColumn.ALT_COUNT.toString(), record.getAltCount())
                    .set(AltSiteRecord.AltSiteRecordTableColumn.REF_F1R2.toString(), record.getRefF1R2())
                    .set(AltSiteRecord.AltSiteRecordTableColumn.ALT_F1R2.toString(), record.getAltF1R2())
                    .set(AltSiteRecord.AltSiteRecordTableColumn.DEPTH.toString(), record.getDepth())
                    .set(AltSiteRecord.AltSiteRecordTableColumn.ALT_BASE.toString(), record.getAltAllele().toString());
        }
    }


}
