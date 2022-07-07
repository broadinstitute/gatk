package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.*;

public class ReferenceSequenceTable implements Iterable<ReferenceSequenceTable.TableRow> {
    private Map<String, TableRow> table;
    private Map<GATKPath, ReferenceDataSource> referenceSources;
    private List<String> columnNames;
    private CompareReferences.MD5CalculationMode md5Mode;

    public ReferenceSequenceTable(Map<GATKPath, ReferenceDataSource> referenceSources, CompareReferences.MD5CalculationMode mode) {
        this.referenceSources = referenceSources;
        columnNames = generateColumnNames();
        md5Mode = mode;
    }

    private List<String> generateColumnNames() {
        List<String> columns = new ArrayList<>();
        columns.add(CompareReferences.MD5_COLUMN_NAME);
        columns.add(CompareReferences.LENGTH_COLUMN_NAME);
        for (GATKPath path : referenceSources.keySet()) {
            columns.add(CompareReferences.getReferenceDisplayName(path));
        }
        return columns;
    }

    public List<String> getColumnNames() {
        return columnNames;
    }

    public void build() {
        table = new LinkedHashMap<>();

        for (Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()) {
            SAMSequenceDictionary dictionary = entry.getValue().getSequenceDictionary();
            for (SAMSequenceRecord record : dictionary.getSequences()) {
                String name = record.getSequenceName();
                int length = record.getSequenceLength();
                String md5 = calculateMD5(record, entry.getValue());
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

    @Override
    public Iterator<TableRow> iterator() {
        return table.values().iterator();
    }

    public class TableEntry {
        private final String sequenceName;
        private final int sequenceLength;
        private final String md5;
        private final GATKPath reference;

        public TableEntry(String name, int length, String md5, GATKPath reference) {
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

    public class TableRow {

        private final String md5;
        private final List<TableEntry> entries;
        private final int length;
        //private List<String> columnNames;
        private final Map<GATKPath, Integer> columnIndices;

        public TableRow(String md5, int length, List<GATKPath> references) {
            this.md5 = md5;
            this.length = length;

            entries = new ArrayList<>(references.size());
            for (int i = 0; i < references.size(); i++) {
                entries.add(null);
            }

            columnIndices = new HashMap<>();
            for (int i = 0; i < references.size(); i++) {
                columnIndices.put(references.get(i), i);
            }
        }

        public void add(TableEntry entry) {
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

        public Map<GATKPath, Integer> getColumnIndices() {
            return columnIndices;
        }

        public List<String> getColumnNames() {
            return columnNames;
        }
    }
}
