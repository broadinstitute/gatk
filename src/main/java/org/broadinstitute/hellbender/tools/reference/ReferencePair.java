package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.engine.GATKPath;

import java.util.*;

/**
 * Class representing a pair of references and their differences.
 *
 */
public class ReferencePair {

    private final GATKPath ref1;

    private final GATKPath ref2;

    private final int ref1ColumnIndex;

    private final int ref2ColumnIndex;

    public enum Status{
        // references match exactly (table entries identical across references)
        EXACT_MATCH,
        // found sequences with same MD5s but inconsistent names between the 2 references
        DIFFER_IN_SEQUENCE_NAMES,
        // found sequences with same names but inconsistent MD5s between the 2 references
        DIFFER_IN_SEQUENCE,
        // at least one sequence is present in one reference but not the other
        DIFFER_IN_SEQUENCES_PRESENT,
        // ref1 is a superset of ref2
        SUPERSET,
        // ref1 is a subset of ref2
        SUBSET;
    }

    private EnumSet<Status> analysis;

    public ReferencePair(ReferenceSequenceTable table, GATKPath reference1, GATKPath reference2){
        ref1 = reference1;
        ref2 = reference2;
        ref1ColumnIndex = table.getColumnIndices().get(ReferenceSequenceTable.getReferenceColumnName(ref1));
        ref2ColumnIndex = table.getColumnIndices().get(ReferenceSequenceTable.getReferenceColumnName(ref2));
        // assume EXACT_MATCH until proven otherwise - if any discrepancies found, EXACT_MATCH status removed
        analysis = EnumSet.of(Status.EXACT_MATCH);
    }

    /**
     * Given a Status, add it to the ReferencePair's analysis
     *
     * @param status
     */
    public void addStatus(Status status) {
        analysis.add(status);
    }

    /**
     * Given a Status, remove it from the ReferencePair's analysis, if present
     *
     * @param status
     */
    public void removeStatus(Status status){
        analysis.remove(status);
    }

    public int getRef1ColumnIndex() {
        return ref1ColumnIndex;
    }

    public int getRef2ColumnIndex() {
        return ref2ColumnIndex;
    }

    public String getRef1(){
        return ReferenceSequenceTable.getReferenceColumnName(ref1);
    }
    public String getRef2(){
        return ReferenceSequenceTable.getReferenceColumnName(ref2);
    }

    /**
     * Displays a ReferencePair's status set
     *
     * @return the status set as a formatted String
     */
    public String statusAsString(){
        String output = "";
        for(Status status : analysis){
            output += String.format("\t%s\n", status.name());
        }
        return output;
    }

    public EnumSet<Status> getAnalysis(){
        return analysis;
    }

    public String toString(){
        return String.format("REFERENCE PAIR: %s, %s\nStatus:\n%s",
                ReferenceSequenceTable.getReferenceColumnName(ref1),
                ReferenceSequenceTable.getReferenceColumnName(ref2),
                statusAsString());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ReferencePair that = (ReferencePair) o;
        return Objects.equals(ref1, that.ref1) && Objects.equals(ref2, that.ref2);
    }

    @Override
    public int hashCode() {
        return Objects.hash(ref1, ref2);
    }
}
