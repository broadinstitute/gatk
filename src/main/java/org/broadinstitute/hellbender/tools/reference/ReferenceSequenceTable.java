package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.sql.Ref;
import java.util.*;

public class ReferenceSequenceTable implements Iterable<ReferenceSequenceTable.TableRow> {

    public static final String MISSING_ENTRY = "---";
    public static final String MD5_COLUMN_NAME = "MD5";
    public static final String LENGTH_COLUMN_NAME = "Length";
    private static final int MD5_COLUMN_INDEX = 0;
    private static final int LENGTH_COLUMN_INDEX = 1;
    private Map<String, TableRow> tableByMD5;
    private Map<String, Set<TableRow>> tableBySequenceName;
    private Map<GATKPath, ReferenceDataSource> referenceSources;
    private List<String> columnNames;
    private CompareReferences.MD5CalculationMode md5Mode;
    private final Map<String, Integer> columnIndices;
    private List<GATKPath> references;

    public ReferenceSequenceTable(Map<GATKPath, ReferenceDataSource> referenceSources, CompareReferences.MD5CalculationMode mode) {
        this.referenceSources = referenceSources;
        columnNames = generateColumnNames();
        md5Mode = mode;
        references = new ArrayList<>(referenceSources.keySet());

        columnIndices = new HashMap<>();
        for (int i = 0; i < columnNames.size(); i++) {
            if(i == MD5_COLUMN_INDEX){
                columnIndices.put(MD5_COLUMN_NAME, i);
            }
            else if(i == LENGTH_COLUMN_INDEX){
                columnIndices.put(LENGTH_COLUMN_NAME, i);
            }
            else{
                columnIndices.put(getReferenceDisplayName(references.get(i-2)), i);
            }
        }
    }

    private List<String> generateColumnNames() {
        List<String> columns = new ArrayList<>();
        columns.add(MD5_COLUMN_NAME);
        columns.add(LENGTH_COLUMN_NAME);
        for (GATKPath path : referenceSources.keySet()) {
            columns.add(getReferenceDisplayName(path));
        }
        return columns;
    }

    public List<String> getColumnNames() {
        return columnNames;
    }

    public static String getReferenceDisplayName(GATKPath reference){
        return reference.toPath().getFileName().toString();
    }

    public Map<String, Integer> getColumnIndices() {
        return columnIndices;
    }

    public void build() {
        tableByMD5 = new LinkedHashMap<>();
        tableBySequenceName = new LinkedHashMap<>();

        for (Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()) {
            SAMSequenceDictionary dictionary = entry.getValue().getSequenceDictionary();
            for (SAMSequenceRecord record : dictionary.getSequences()) {
                String name = record.getSequenceName();
                int length = record.getSequenceLength();
                String md5 = calculateMD5(record, entry.getValue());
                GATKPath reference = entry.getKey();
                TableEntry newEntry = new TableEntry(getReferenceDisplayName(reference), name);
                TableRow newRow = null;

                // map each MD5 to List of TableEntry objects containing length, md5, and name
                if (!tableByMD5.containsKey(md5)) {
                    newRow = new TableRow(md5, length, new ArrayList<>(referenceSources.keySet()));
                    tableByMD5.put(md5, newRow);
                }
                tableByMD5.get(md5).add(newEntry);

                if(!tableBySequenceName.containsKey(name)){
                    Set<TableRow> rowSet = new HashSet<>();
                    tableBySequenceName.put(name, rowSet);
                }
                tableBySequenceName.get(name).add(tableByMD5.get(md5));
            }
        }
    }

    public Set<String> getAllSequenceNames(){
        return tableBySequenceName.keySet();
    }

    // number of rows in md5 keyed table
    public int size(){
        int size = 0;
        for(TableRow row : this){
            size++;
        }
        return size;
    }

    public TableRow queryByMD5(String md5){
        return tableByMD5.get(md5);
    }

    public Set<TableRow> queryBySequenceName(String sequenceName){
        return tableBySequenceName.get(sequenceName) == null ? Collections.emptySet() : tableBySequenceName.get(sequenceName);
    }

    private String calculateMD5(SAMSequenceRecord record, ReferenceDataSource source) {
        final String md5FromDict = record.getMd5();
        final SimpleInterval referenceInterval = new SimpleInterval(record.getSequenceName(), 1, record.getSequenceLength());
        String md5 = null;

        switch (md5Mode) {
            case USE_DICT:
                if (md5FromDict == null || md5FromDict.isEmpty()) {
                    throw new UserException.BadInput("Running in USE_DICT mode, but MD5 missing for sequence " + record.getSequenceName() + ". Run --md5-calculation-mode with a different mode to recalculate MD5.");
                }
                md5 = md5FromDict;
                break;
            case RECALCULATE_IF_MISSING:
                if (md5FromDict == null || md5FromDict.isEmpty()) {
                    md5 = ReferenceUtils.calculateMD5(source, referenceInterval);
                } else {
                    md5 = md5FromDict;
                }
                break;
            case ALWAYS_RECALCULATE:
                md5 = ReferenceUtils.calculateMD5(source, referenceInterval);
                break;
        }
        return md5;
    }

    public List<ReferencePair> generateReferencePairs(List<GATKPath> references){
        List<ReferencePair> referencePairs = new ArrayList<>();
        for(int i = 0; i < references.size(); i++){
            for(int j = i + 1; j < references.size(); j++){
                referencePairs.add(new ReferencePair(this, references.get(i), references.get(j)));
            }
        }
        return referencePairs;
    }

    public void analyzeTable(List<ReferencePair> refPairs){

        for(TableRow row : tableByMD5.values()) {
            for (ReferencePair pair : refPairs) {
                int ref1Index = pair.getRef1ColumnIndex();
                int ref2Index = pair.getRef2ColumnIndex();
                TableEntry[] entries = row.getEntries();

                TableEntry ref1Value = entries[ref1Index];
                TableEntry ref2Value = entries[ref2Index];
                if (!ref1Value.equals(ref2Value)) {
                    pair.removeStatus(ReferencePair.Status.EXACT_MATCH);
                } if(!ref1Value.getColumnValue().equals(ref2Value.getColumnValue()) &&
                        !(ref1Value.getColumnValue().equals(MISSING_ENTRY) ^ ref2Value.getColumnValue().equals(MISSING_ENTRY))){
                    pair.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE_NAMES);
                }
            }
        }

        for(String sequenceName : tableBySequenceName.keySet()){
            for(ReferencePair pair : refPairs){
                int ref1Index = pair.getRef1ColumnIndex();
                int ref2Index = pair.getRef2ColumnIndex();
                Set<ReferenceSequenceTable.TableRow> rows = queryBySequenceName(sequenceName);
                int sequenceNameFound = 0;

                for(TableRow row : rows) {
                    TableEntry[] entries = row.getEntries();
                    TableEntry ref1Value = entries[ref1Index];
                    TableEntry ref2Value = entries[ref2Index];

                    if((ref1Value.getColumnValue().equals(sequenceName) ^ ref2Value.getColumnValue().equals(sequenceName)) && (ref1Value.getColumnValue().equals(MISSING_ENTRY) ^ ref2Value.getColumnValue().equals(MISSING_ENTRY))){
                        sequenceNameFound++;
                        continue;
                    }
                }

                if(sequenceNameFound == 2){
                    pair.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE);
                }
                if(sequenceNameFound == 1){
                    pair.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCES_PRESENT);
                }
            }
        }
        for(ReferencePair pair : refPairs){
            System.out.println(pair.toString());
        }
    }

    @Override
    public Iterator<TableRow> iterator() {
        return tableByMD5.values().iterator();
    }

    public class TableEntry {

        private final String columnName;
        private final String columnValue;

        public TableEntry(String columnName, String columnValue){
            this.columnName = columnName;
            this.columnValue = columnValue;
        }

        public String getColumnName() { return columnName; }
        public String getColumnValue() { return columnValue; }

        public boolean isEmpty(){
            return columnValue.equals(MISSING_ENTRY);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            TableEntry that = (TableEntry) o;
            return columnName.equals(that.columnName) && columnValue.equals(that.columnValue);
        }

        @Override
        public int hashCode() {
            return Objects.hash(columnName, columnValue);
        }
    }

    public class TableRow {

//        private static final int MD5_COLUMN_INDEX = 0;
//        private static final int LENGTH_COLUMN_INDEX = 1;
        private final String md5;
        private final TableEntry[] entries;
        private final int length;

        public TableRow(String md5, int length, List<GATKPath> references) {
            this.md5 = md5;
            this.length = length;
            entries = new TableEntry[references.size()+2];

            add(new TableEntry(MD5_COLUMN_NAME, md5));
            add(new TableEntry(LENGTH_COLUMN_NAME, Integer.toString(length)));

            for(int i = LENGTH_COLUMN_INDEX + 1; i < entries.length; i++){
                entries[i] = new TableEntry(getColumnNames().get(i), MISSING_ENTRY);
            }
        }

        public void add(TableEntry entry) {
            int idx = columnIndices.get(entry.getColumnName());
            entries[idx] = entry;
        }

        public String getMd5() {
            return md5;
        }

        public TableEntry[] getEntries() {
            return entries;
        }

        public int getLength() {
            return length;
        }

        public List<String> getColumnNames() {
            return columnNames;
        }

        public String toString(){
            String row = String.format("%s\t", this.md5);
            for(int i = 0; i < entries.length; i++){
                row += entries[i].columnValue;
            }
            return row;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            TableRow tableRow = (TableRow) o;
            return md5.equals(tableRow.md5);
        }

        @Override
        public int hashCode() {
            return Objects.hash(md5);
        }
    }
}
