package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.*;

/**
 * Table utilized by CompareReferences tool to compare and analyze sequences found in specified references.
 * Underlying table is keyed by MD5, with the ability to query by sequence name across references as well.
 */
public class ReferenceSequenceTable implements Iterable<ReferenceSequenceTable.TableRow> {

    public static final Logger logger = LogManager.getLogger(ReferenceSequenceTable.class);
    public static final String MISSING_ENTRY_DISPLAY_STRING = "---";
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
    private boolean tableBuilt;

    /**
     * Returns a ReferenceSequenceTable object.
     *
     * Note: table is not constructed until a call to build() is made. All methods which attempt to use
     * the table before it has been constructed will crash.
     *
     * @param referenceSources map containing references and their data sources
     * @param mode MD5 calculation mode
     */
    public ReferenceSequenceTable(Map<GATKPath, ReferenceDataSource> referenceSources, CompareReferences.MD5CalculationMode mode) {
        this.referenceSources = referenceSources;
        columnNames = generateColumnNames();
        md5Mode = mode;
        references = new ArrayList<>(referenceSources.keySet());
        tableBuilt = false;

        columnIndices = new HashMap<>();
        for (int i = 0; i < columnNames.size(); i++) {
            if(i == MD5_COLUMN_INDEX){
                columnIndices.put(MD5_COLUMN_NAME, i);
            }
            else if(i == LENGTH_COLUMN_INDEX){
                columnIndices.put(LENGTH_COLUMN_NAME, i);
            }
            else{
                columnIndices.put(getReferenceColumnName(references.get(i-2)), i);
            }
        }
    }

    private List<String> generateColumnNames() {
        List<String> columns = new ArrayList<>();
        columns.add(MD5_COLUMN_NAME);
        columns.add(LENGTH_COLUMN_NAME);
        for (GATKPath path : referenceSources.keySet()) {
            columns.add(getReferenceColumnName(path));
        }
        return columns;
    }

    public List<String> getColumnNames() {
        return columnNames;
    }

    /**
     * Given a GATKPath, return the name of the file as a String.
     *
     * TODO: may need to add functionality to uniquify reference names
     *
     * @param reference The path to a reference.
     * @return the name of the reference as a String.
     */
    public static String getReferenceColumnName(GATKPath reference){
        return reference.toPath().getFileName().toString();
    }

    public Map<String, Integer> getColumnIndices() {
        return columnIndices;
    }

    /**
     * Construct 2 tables of references: one keyed by MD5, one keyed by sequence name.
     */
    public void build() {
        tableByMD5 = new LinkedHashMap<>();
        tableBySequenceName = new LinkedHashMap<>();

        for (Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()) {
            SAMSequenceDictionary dictionary = entry.getValue().getSequenceDictionary();
            for (SAMSequenceRecord record : dictionary.getSequences()) {
                String name = record.getSequenceName();
                int length = record.getSequenceLength();
                String md5 = getMD5ForRecord(record, entry.getKey());
                GATKPath reference = entry.getKey();
                TableEntry newEntry = new TableEntry(getReferenceColumnName(reference), name);
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
        tableBuilt = true;
    }

    public Set<String> getAllSequenceNames(){
        if(!tableBuilt){
            throw new IllegalStateException("Table has not been built yet.");
        }

        return tableBySequenceName.keySet();
    }

    /**
     * Given an MD5, returns its corresponding row
     *
     * @param md5 The MD5 as a String
     * @return the corresponding TableRow from the tableByMD5
     */
    public TableRow queryByMD5(String md5){
        if(!tableBuilt){
            throw new IllegalStateException("Table has not been built yet.");
        }

        return tableByMD5.get(md5);
    }

    /**
     * Given a sequence name, returns the set of its corresponding rows or an empty set if no such rows found.
     *
     * @param sequenceName The sequence name as a String
     * @return the set of TableRows that contain the sequence name
     */
    public Set<TableRow> queryBySequenceName(String sequenceName){
        if(!tableBuilt){
            throw new IllegalStateException("Table has not been built yet.");
        }

        Set<TableRow> rows = tableBySequenceName.get(sequenceName);
        return rows == null ? Collections.emptySet() : rows;
    }

    private String getMD5ForRecord(SAMSequenceRecord record, GATKPath referencePath) {
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
                    logger.warn(String.format("%s in sequence dictionary for %s missing MD5. MD5 recalculated.", record.getSequenceName(), getReferenceColumnName(referencePath)));
                    md5 = ReferenceUtils.calculateMD5(referencePath, referenceInterval);
                }
                else {
                    md5 = md5FromDict;
                }
                break;
            case ALWAYS_RECALCULATE:
                md5 = ReferenceUtils.calculateMD5(referencePath, referenceInterval);
                if(md5FromDict != null && !md5FromDict.equals(md5)){
                    logger.warn(String.format("MD5 Mismatch for sequence %s. Found '%s', but calculated '%s'", record.getSequenceName(), md5FromDict, md5));
                }
                break;
        }
        return md5;
    }

    /**
     * Generate ReferencePairs for pairwise comparison of all references present in the table
     *
     * @return the list of ReferencePairs for every pair of references
     */
    public List<ReferencePair> generateReferencePairs(){
        List<ReferencePair> referencePairs = new ArrayList<>();
        for(int i = 0; i < references.size(); i++){
            for(int j = i + 1; j < references.size(); j++){
                referencePairs.add(new ReferencePair(this, references.get(i), references.get(j)));
            }
        }
        return referencePairs;
    }

    /**
     * Analyze the table by doing a pairwise comparison for all table references. Generates all ReferencePairs, then
     * analyzes each pair. First, traverses the table by MD5 and checks for any discrepancies across the column values
     * that violate an exact match or a difference in sequence name. Then traverses by sequence name and checks
     * each sequence name across the references to determine if the sequence name is present in both but corresponding to
     * different sequences or missing in one, and if so, if it constitutes a superset/subset case.
     *
     * Each ReferencePair is analyzed as above and assigned an analysis as a set of the following statuses:
     *      EXACT_MATCH,
     *      DIFFER_IN_SEQUENCE_NAMES,
     *      DIFFER_IN_SEQUENCE,
     *      DIFFER_IN_SEQUENCES_PRESENT,
     *      SUPERSET,
     *      SUBSET
     *
     * @return list of ReferencePairs with updated status sets
     */
    public List<ReferencePair> analyzeTable(){
        if(!tableBuilt){
            throw new IllegalStateException("Table has not been built yet.");
        }

        List<ReferencePair> refPairs = generateReferencePairs();

        for(TableRow row : tableByMD5.values()) {
            for(ReferencePair pair : refPairs) {
                int ref1Index = pair.getRef1ColumnIndex();
                int ref2Index = pair.getRef2ColumnIndex();
                TableEntry[] entries = row.getEntries();

                TableEntry ref1Value = entries[ref1Index];
                TableEntry ref2Value = entries[ref2Index];
                if (!ref1Value.getColumnValue().equals(ref2Value.getColumnValue())) {
                    pair.removeStatus(ReferencePair.Status.EXACT_MATCH);
                }
                if(!ref1Value.getColumnValue().equals(ref2Value.getColumnValue()) &&
                        !(ref1Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING) || ref2Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING))){
                    pair.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE_NAMES);
                }
            }
        }

        for(ReferencePair pair : refPairs){
            boolean superset = false;
            boolean subset = false;
            for(String sequenceName : tableBySequenceName.keySet()){
                int ref1Index = pair.getRef1ColumnIndex();
                int ref2Index = pair.getRef2ColumnIndex();
                Set<ReferenceSequenceTable.TableRow> rows = queryBySequenceName(sequenceName);
                int sequenceNameFoundInOneRef = 0;

                for(TableRow row : rows) {
                    TableEntry[] entries = row.getEntries();
                    TableEntry ref1Value = entries[ref1Index];
                    TableEntry ref2Value = entries[ref2Index];

                    if((ref1Value.getColumnValue().equals(sequenceName) ^ ref2Value.getColumnValue().equals(sequenceName))
                            && (ref1Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING) ^ ref2Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING))){
                        sequenceNameFoundInOneRef++;
                    }
                    if(ref1Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING) ^ ref2Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING)){
                        if(ref1Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING)){
                            subset = true;
                        }
                        else if(ref2Value.getColumnValue().equals(MISSING_ENTRY_DISPLAY_STRING)){
                            superset = true;
                        }
                    }
                }

                // if 2 difference sequences with same name found, mismatching sequences for that name
                if(sequenceNameFoundInOneRef == 2){
                    pair.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE);
                // if a sequence name is found only once, it is present in one reference but not the other
                } else if(sequenceNameFoundInOneRef == 1){
                    pair.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCES_PRESENT);
                } else if (sequenceNameFoundInOneRef > 2){
                    throw new UserException.BadInput(String.format("Duplicate of sequence '%s' found in %s or %s.", sequenceName, pair.getRef1(), pair.getRef2()));
                }
            }

            // if a pair has exclusively superset or subset marked as true and no naming discrepancies found,
            // mark it as super/sub set and remove the 'DIFFER_IN_SEQUENCES_PRESENT' tag
            if(superset ^ subset && !pair.getAnalysis().contains(ReferencePair.Status.DIFFER_IN_SEQUENCE_NAMES)){
                pair.removeStatus(ReferencePair.Status.DIFFER_IN_SEQUENCES_PRESENT);

                if(superset && !subset){
                    pair.addStatus(ReferencePair.Status.SUPERSET);
                } else if(subset && !superset) {
                    pair.addStatus(ReferencePair.Status.SUBSET);
                }
            }
        }
        return refPairs;
    }

    @Override
    public Iterator<TableRow> iterator() {
        if(!tableBuilt){
            throw new IllegalStateException("Table has not been built yet.");
        }

        return tableByMD5.values().iterator();
    }

    /**
     * Minimal class representing a table entry within a row of the ReferenceSequenceTable.
     * Stores the entry's value and the associated column.
     */
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
            return columnValue.equals(MISSING_ENTRY_DISPLAY_STRING);
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

    /**
     * Class representing a row of the ReferenceSequenceTable. Rows are keyed by MD5 and the remaining columns
     * contain the sequence length, followed by the sequence name across each reference.
     */
    public class TableRow {
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
                entries[i] = new TableEntry(getColumnNames().get(i), MISSING_ENTRY_DISPLAY_STRING);
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

        // TableRow keyed by MD5 - equals() and hashCode() check only MD5 to allow for TableRow equality based on MD5 only
        // 2 TableRows can be considered equal when they have matching MD5s
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
