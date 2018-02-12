package org.broadinstitute.hellbender.utils.read.markduplicates;

import org.testng.annotations.Test;

import java.util.Random;

public class ReadsKeyUnitTest {

    @Test (enabled = false)
    // This test makes the argument that using a StringBuilder in ReadsKey is faster than an equivalent call to String.format()
    /* an example of the results:
        With Format:	13909867289
        With Builder:	1130077051
        With String +:	1155343164
        With Joiner:	2013227121
     */
    public void testKeyForFragmentPerformanceEvidence() throws Exception {
        Random rand = new Random();
        int numTrials = 10000000;
        for (int i = 0; i < numTrials; i++) {
            rand.nextInt();
        }

        long startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = String.format(
                    "%s|%d|%d|%s",
                    rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-",
                    rand.nextInt(30),
                    rand.nextInt(100),
                    rand.nextBoolean() ? "r" : "f");
        }
        long endTime = System.nanoTime();
        System.out.println("With Format:\t" + (endTime - startTime));

        startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = new StringBuilder().append(rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-")
                    .append("|")
                    .append(rand.nextInt(30))
                    .append("|")
                    .append(rand.nextInt(100))
                    .append("|")
                    .append( rand.nextBoolean() ? "r" : "f")
                    .toString();

        }
        endTime = System.nanoTime();
        System.out.println("With Builder:\t" + (endTime - startTime));

        startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = (rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-") +
                    "|" +
                    rand.nextInt(30) +
                    "|" +
                    rand.nextInt(100) +
                    "|" +
                    (rand.nextBoolean() ? "r" : "f");

        }
        endTime = System.nanoTime();
        System.out.println("With String +:\t" + (endTime - startTime));

        startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = String.join("|", rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-",
                    Integer.toString(rand.nextInt(30)),
                    Integer.toString(rand.nextInt(100)),
                    rand.nextBoolean() ? "r" : "f");

        }
        endTime = System.nanoTime();
        System.out.println("With Joiner:\t" + (endTime - startTime));
    }

}