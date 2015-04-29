package org.broadinstitute.hellbender.utils.illumina;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Just some simple tests for the IlluminaUtil.getTileFromReadName() method for now.
 *
 * @author Tim Fennell
 */
public class IlluminaUtilTest {
    @Test(dataProvider="readNames") public void testFindTileInReadName(final String readName, final Integer tile) {
        final Integer otherTile = IlluminaUtil.getTileFromReadName(readName);
        Assert.assertEquals(otherTile, tile, "Tile numbers do not match for read name: " + readName);
    }

    @Test
    public void performanceTestGetTileFromReadName() {
        final int ITERATIONS = 5000000;

        final long startTime = System.currentTimeMillis();
        for (int i=0; i<ITERATIONS; ++i) {
            final Integer tile = IlluminaUtil.getTileFromReadName("300WFAAXX090909:1:1:1024:978#0/1");
            if (tile == null || tile != 1) throw new RuntimeException("WTF?");
        }
        final long endTime = System.currentTimeMillis();

        System.out.println("Time taken: " + (endTime-startTime) + "ms.");
    }

    @DataProvider(name="readNames")
    public Object[][] readNames() {
        return new Object[][] {
                new Object[] {"300WFAAXX:1:119:1024:978#0/1", 119},
                new Object[] {"300WFAAXX090909:1:1:1024:978#0/1", 1},
                new Object[] {"FOO", null},
                new Object[] {"FOO:BAR_splat", null}
        };
    }

    public static byte [] iToB(int [] intVals) {
        byte [] byteVals = new byte[intVals.length];
        for(int i = 0; i < byteVals.length; i++) {
            byteVals[i] = (byte) intVals[i];
        }
        return byteVals;
    }

    @DataProvider(name = "solexaQualStrToPhreds")
    public Object[][] solexaQualStrToPhredsToPhreds() {
        return new Object[][] {
            new Object[] {
                "x@Axy" + ((char)156) + ((char) 157) + ((char) 0) + ((char) 1) + "?",
                new byte[][] {
                    new byte[]{}, iToB(new int[]{56}),  iToB(new int[]{0, 1, 56, 57, 92, 93}),
                        iToB(new int[]{-64, -63, -1}), new byte[]{}
                }
            },
            new Object[] {
                "nNpZo" + ((char)250) + ((char) 255) + ((char) 1) + ((char) 1) + "CCaB",
                new byte[][] {
                    iToB(new int[]{46, 14}),  iToB(new int[]{48, 26, 47, -70, -65, -63, -63, 3, 3, 33}), iToB(new int[]{2})
                }
            }
        };
    }
}
