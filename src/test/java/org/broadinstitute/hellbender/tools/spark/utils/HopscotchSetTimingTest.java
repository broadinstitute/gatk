package org.broadinstitute.hellbender.tools.spark.utils;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;

import java.util.*;

/**
 * Timing tests for HopscotchSet.
 */
public final class HopscotchSetTimingTest {

    @FunctionalInterface
    public interface Action {
        void execute();
    }
    private static double time( final Action action ) {
        final long nanosecs = System.nanoTime();
        action.execute();
        return (System.nanoTime() - nanosecs)/1.E9;
    }

    private static final int N_TRIALS = 5;
    private static final int N_VALUES = 25000000;

    public static void main( final String[] args ) {

        final Random rng = new Random(0xdeadbeef);
        final List<Integer[]> trials = new ArrayList<>(N_TRIALS);
        for ( int trialId = 0; trialId != N_TRIALS; ++trialId ) {
            final Integer[] values = new Integer[N_VALUES];
            for ( int valueId = 0; valueId != N_VALUES; ++valueId ) {
                values[valueId] = rng.nextInt();
            }
            trials.add(values);
        }

        final List<Set<Integer>> hashSets = new ArrayList<>(N_TRIALS);
        System.out.println("HashSet construction: "+time( () -> {
            for ( final Integer[] values : trials ) {
                final Set<Integer> hashSet = new HashSet<>(SVUtils.hashMapCapacity(N_VALUES));
                for (final Integer value : values) {
                    hashSet.add(value);
                }
                hashSets.add(hashSet);
            }
        }));

        final List<Set<Integer>> hopscotchHashSets = new ArrayList<>(N_TRIALS);
        System.out.println("HopscotchSet construction: "+time( () -> {
            for (final Integer[] values : trials) {
                final Set<Integer> hashSet = new HopscotchSet<>(N_VALUES);
                for (final Integer value : values) {
                    hashSet.add(value);
                }
                hopscotchHashSets.add(hashSet);
            }
        }));

        System.out.println("HashSet +retrieval: "+time( () -> {
            for (int trialId = 0; trialId != N_TRIALS; ++trialId) {
                final Set<Integer> hashSet = hashSets.get(trialId);
                for (final Integer value : trials.get(trialId))
                    hashSet.contains(value);
            }
        }));

        System.out.println("HopscotchSet +retrieval: "+time( () -> {
            for (int trialId = 0; trialId != N_TRIALS; ++trialId) {
                final Set<Integer> hashSet = hopscotchHashSets.get(trialId);
                for (final Integer value : trials.get(trialId))
                    hashSet.contains(value);
            }
        }));

        final Integer[] missingValues = new Integer[N_VALUES];
        for ( int valueId = 0; valueId != N_VALUES; ++valueId ) {
            missingValues[valueId] = rng.nextInt();
        }

        System.out.println("HashSet -retrieval: "+time( () -> {
            for (int trialId = 0; trialId != N_TRIALS; ++trialId) {
                final Set<Integer> hashSet = hashSets.get(trialId);
                for (final Integer value : missingValues)
                    hashSet.contains(value);
            }
        }));

        System.out.println("HopscotchSet -retrieval: "+time( () -> {
            for (int trialId = 0; trialId != N_TRIALS; ++trialId) {
                final Set<Integer> hashSet = hopscotchHashSets.get(trialId);
                for (final Integer value : missingValues)
                    hashSet.contains(value);
            }
        }));

        System.out.println("HashSet removal: "+time( () -> {
            for (int trialId = 0; trialId != N_TRIALS; ++trialId) {
                final Set<Integer> hashSet = hashSets.get(trialId);
                for (final Integer value : trials.get(trialId))
                    hashSet.remove(value);
            }
        }));

        System.out.println("HopscotchSet removal: "+time( () -> {
            for (int trialId = 0; trialId != N_TRIALS; ++trialId) {
                final Set<Integer> hashSet = hopscotchHashSets.get(trialId);
                for (final Integer value : trials.get(trialId))
                    hashSet.remove(value);
            }
        }));
    }
}
