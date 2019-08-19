package org.broadinstitute.hellbender.utils.read.markduplicates.sparksql;

public class NonDuplicatesSparkSqlRecord {
    private Integer partitionIndex;
    private String name;
    private Integer numOpticalDuplicates;

    public NonDuplicatesSparkSqlRecord() {
    }

    public NonDuplicatesSparkSqlRecord(Integer partitionIndex, String name, Integer numOpticalDuplicates) {
        this.partitionIndex = partitionIndex;
        this.name = name;
        this.numOpticalDuplicates = numOpticalDuplicates;
    }

    public Integer getPartitionIndex() {
        return partitionIndex;
    }

    public void setPartitionIndex(Integer partitionIndex) {
        this.partitionIndex = partitionIndex;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Integer getNumOpticalDuplicates() {
        return numOpticalDuplicates;
    }

    public void setNumOpticalDuplicates(Integer numOpticalDuplicates) {
        this.numOpticalDuplicates = numOpticalDuplicates;
    }
}
