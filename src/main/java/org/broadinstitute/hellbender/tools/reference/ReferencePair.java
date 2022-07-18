package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.engine.GATKPath;

import java.util.*;

public class ReferencePair {

    private final GATKPath ref1;
    private final GATKPath ref2;
    private final int ref1ColumnIndex;
    private final int ref2ColumnIndex;
    public enum Status{
        EXACT_MATCH,
        DIFFER_IN_SEQUENCE_NAMES,
        DIFFER_IN_SEQUENCE,
        DIFFER_IN_SEQUENCES_PRESENT,
        SUPERSET,
        SUBSET;
    }
    private EnumSet<Status> analysis;

    public ReferencePair(ReferenceSequenceTable table, GATKPath reference1, GATKPath reference2){
        ref1 = reference1;
        ref2 = reference2;
        ref1ColumnIndex = table.getColumnIndices().get(ref1.toPath().getFileName().toString());
        ref2ColumnIndex = table.getColumnIndices().get(ref2.toPath().getFileName().toString());
        analysis = EnumSet.of(Status.EXACT_MATCH);
    }

    public void addStatus(Status status) {
        analysis.add(status);
    }

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
        return ref1.toPath().getFileName().toString();
    }
    public String getRef2(){
        return ref2.toPath().getFileName().toString();
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

    public EnumSet<Status> getStatus(){
        return analysis;
    }

    public String toString(){
        return String.format("REFERENCE PAIR: %s, %s\nStatus:\n%s",
                ReferenceSequenceTable.getReferenceDisplayName(ref1),
                ReferenceSequenceTable.getReferenceDisplayName(ref2),
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
