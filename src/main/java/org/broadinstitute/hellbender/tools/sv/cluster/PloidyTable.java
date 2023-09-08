package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

public class PloidyTable {

    private final Map<String, PloidyRecord> samplePloidyMap;

    public PloidyTable(final Path path) {
        try {
            final TableReader<PloidyRecord> reader = TableUtils.reader(path, (columns, exceptionFactory) ->
                (dataLine) -> {
                    final String sample = dataLine.get(0);
                    final Map<String, Integer> dataMap = new HashMap<>();
                    for (int i = 1; i < columns.columnCount(); i++) {
                        final String contig = columns.names().get(i);
                        final Integer ploidy = dataLine.getInt(i);
                        dataMap.put(contig, ploidy);
                    }
                    return new PloidyRecord(sample, dataMap);
                }
            );
            samplePloidyMap = reader.stream().collect(Collectors.toMap(PloidyRecord::getSample, r -> r));
            reader.close();
        } catch (final IOException e) {
            throw new GATKException("IO error while reading ploidy table", e);
        }
    }

    public PloidyTable(final Map<String, Map<String, Integer>> samplePloidyMap) {
        this.samplePloidyMap = samplePloidyMap.entrySet().stream()
                .collect(Collectors.toMap(p -> p.getKey(), p -> new PloidyRecord(p.getKey(), p.getValue())));
    }

    public Integer get(final String sample, final String contig) {
        Utils.validateArg(samplePloidyMap.containsKey(sample), "Sample " + sample + " not found in ploidy records");
        return samplePloidyMap.get(sample).getPloidy(contig);
    }

    private static final class PloidyRecord {
        private final String sample;
        private final Map<String, Integer> ploidyMap;

        public PloidyRecord(final String sample, final Map<String, Integer> ploidyMap) {
            this.sample = sample;
            this.ploidyMap = ploidyMap;
        }

        public String getSample() {
            return sample;
        }

        public Integer getPloidy(final String contig) {
            Utils.validateArg(ploidyMap.containsKey(contig), "No ploidy entry for sample " + sample +
                    " at contig " + contig);
            return ploidyMap.get(contig);
        }
    }
}
