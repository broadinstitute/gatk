package org.broadinstitute.hellbender.tools.exome.eval;

import org.apache.parquet.it.unimi.dsi.fastutil.objects.Object2IntLinkedOpenHashMap;
import org.apache.parquet.it.unimi.dsi.fastutil.objects.Object2IntMap;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.stream.Collectors;

/**
 * Represents a sample evaluation class count record.
 */
public final class EvaluationSampleSummaryRecord {

    private final String sample;
    private final Object2IntMap<EvaluationClass> countsByClass;
    private int total;

    public EvaluationSampleSummaryRecord(final String sampleName) {
        sample = Utils.nonNull(sampleName);
        countsByClass = new Object2IntLinkedOpenHashMap<>(EvaluationClass.values().length);
        for (final EvaluationClass ec : EvaluationClass.values()) {
            countsByClass.put(ec, 0);
        }
        total = 0;
    }

    @Override
    public String toString() {
        final String countsString = countsByClass.entrySet().stream()
                .map(entry -> String.format("%s = %d", entry.getKey(), entry.getValue()))
                .collect(Collectors.joining(", "));
        return String.format("Overall stats: ALL = %d, %s", getTotal(), countsString);
    }

    public void increase(final EvaluationClass evalClass) {
        Utils.nonNull(evalClass);
        countsByClass.put(evalClass, countsByClass.getInt(evalClass) + 1);
        total = getTotal() + 1;
    }

    public void set(final EvaluationClass evalClass, final int newValue) {
        Utils.nonNull(evalClass);
        ParamUtils.isPositiveOrZero(newValue, "the new value cannot be negative");
        final int previousCount = countsByClass.put(evalClass, newValue);
        total += newValue - previousCount;
    }

    public void aggregate(final EvaluationSampleSummaryRecord record) {
        Utils.nonNull(record);
        for (final EvaluationClass evalClass : countsByClass.keySet()) {
            countsByClass.put(evalClass,
                    countsByClass.get(evalClass) + record.countsByClass.get(evalClass));
        }
        total = getTotal() + record.getTotal();
    }

    public int getTotal() {
        return total;
    }

    public int get(final EvaluationClass clazz) {
        Utils.nonNull(clazz);
        return countsByClass.getInt(clazz);
    }

    public String getSample() {
        return sample;
    }
}
