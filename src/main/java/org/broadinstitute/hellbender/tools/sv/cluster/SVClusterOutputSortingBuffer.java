package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;

import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

public final class SVClusterOutputSortingBuffer {
    private final TreeSet<SVCallRecord> buffer;
    private final SVClusterEngine<SVCallRecord> engine;
    private final Comparator<SVCallRecord> recordComparator;

    public SVClusterOutputSortingBuffer(final SVClusterEngine<SVCallRecord> engine, final SAMSequenceDictionary dictionary) {
        this.buffer = new TreeSet<>(SVCallRecordUtils.getCallComparator(dictionary));
        this.recordComparator = SVCallRecordUtils.getCallComparator(dictionary);
        this.engine = engine;
    }

    public List<SVCallRecord> flush(final String currentContig) {
        buffer.addAll(engine.getOutput());
        final Integer minActiveStart = engine.getMinActiveStartingPosition();
        final int minPos = minActiveStart == null ? Integer.MAX_VALUE : minActiveStart;
        final List<SVCallRecord> result = buffer.stream()
                .filter(record -> !record.getContigA().equals(currentContig) || record.getPositionA() < minPos)
                .sorted(recordComparator)
                .collect(Collectors.toList());
        buffer.removeAll(result);
        return result;
    }

    public List<SVCallRecord> forceFlush() {
        buffer.addAll(engine.forceFlushAndGetOutput());
        final List<SVCallRecord> result = buffer.stream().sorted(recordComparator).collect(Collectors.toList());
        buffer.clear();
        return result;
    }
}
