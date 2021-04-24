package org.broadinstitute.hellbender.tools.dragstr;

import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.Utils;   
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class InterleavingListSpliteratorUnitTest {

    enum Operation {
        // trySplit
        SPLIT('^'),
        // tryAdvance
        ADVANCE('.'),
        // forEachRemaining.
        READTHRU('*'),
        // popEmpty.
        POPEMPTY('$');

        final char code;

        Operation(final char code) {
            this.code  = code;
        }

        static Operation of(char ch) {
            if (ch == SPLIT.code) {
                return SPLIT;
            } else if (ch == ADVANCE.code) {
                return ADVANCE;
            } else if (ch == READTHRU.code) {
                return READTHRU;
            } else if (ch == POPEMPTY.code) {
                return POPEMPTY;
            } else if (Character.isDigit(ch)) {
                return ADVANCE;
            } else {
                throw new IllegalArgumentException();
            }
        }
    }

    @Test(dataProvider = "lengthAndNumberOfThreadsData")
    public void testParallel(final int length, final int numberOfThreads) {
        final List<Integer> values = firstNaturalIntegers(length);
        final List<Integer> result = Utils.runInParallel(numberOfThreads,
                () -> StreamSupport.stream(new InterleavingListSpliterator<>(values), true)
                    .collect(Collectors.toList()));

        Collections.sort(result);
        Assert.assertEquals(result, values);
    }

    @Test(dataProvider = "lengthAndSplitLevelsData")
    public void testInterleaveness(final int length, final int splitLevels) {
        final Spliterator<Integer> subject = new InterleavingListSpliterator<>(firstNaturalIntegers(length));
        final Deque<Pair<Spliterator<Integer>, List<Integer>>> spliteratorsAndExpected = new ArrayDeque<>(1 << splitLevels);
        spliteratorsAndExpected.add(new Pair<>(subject, firstNaturalIntegers(length)));
        // we split all spliterators splitLevels times recursively using a dequeue to implement it with a single iteration.
        for (int i = 0; i < splitLevels; i++) {
            final int count = spliteratorsAndExpected.size();
            for (int j = 0; j < count; j++) {
                final Pair<Spliterator<Integer>, List<Integer>> pair = spliteratorsAndExpected.removeFirst();
                final Spliterator<Integer> secondSplit = pair.getFirst().trySplit();
                if (secondSplit == null) {
                    Assert.assertTrue(pair.getFirst().estimateSize() <= 1);
                    spliteratorsAndExpected.addLast(pair);
                } else {
                    final List<Integer> oldExpected = pair.getSecond();
                    final List<Integer> firstExpected = new ArrayList<>(oldExpected.size());
                    final List<Integer> secondExpected = new ArrayList<>(oldExpected.size());
                    for (int k = 0; k < oldExpected.size(); k++) {
                        ((k & 1) == 0 ? firstExpected : secondExpected).add(oldExpected.get(k));
                    }
                    spliteratorsAndExpected.addLast(new Pair<>(pair.getFirst(), firstExpected));
                    spliteratorsAndExpected.addLast(new Pair<>(secondSplit, secondExpected));
                }
            }
        }
        // Now we check actual and expected sequences of all spliterators.
        for (final Pair<Spliterator<Integer>, List<Integer>> pair : spliteratorsAndExpected) {
            final List<Integer> actual = StreamSupport.stream(pair.getFirst(), false).collect(Collectors.toList());
            Assert.assertEquals(actual, pair.getSecond());
        }
    }

    /**
     * Tests a particular length list given the list of operations to perform on it.
     * <p>
     *     The operation sequence describe the set of operations to perfom on the initial
     *     full list splitterator.
     * </p>
     * <p>
     *     Each operator indicate what to do with the spliterator that seat at the top of
     *     a stack.
     *     <table border="1" align="center">
     *      <tr><td width="10%" align="center">code</td><td align="center">action</td><td align="center" width="35%" >pseudo-code</td></tr>
     *          <td style="text-align: center">^<br/>(hat)</td><td>we try-split the stack top (popped); the new spliterator is pushed followed by the original one; if try-split fails (cannot be splitted)
     *                        the spliterator is push back as it was.<p/> Fails if there is no spliterator left on the stack.</td>
     *          <td style="padding-left: 10px"><code>
     *              A = pop; <br/>
     *              B = A.trySplit();<br/>
     *              push B if B != null; <br/>
     *              push A;
     *          </pre></td></tr>
     *      <tr><td style="text-align: center">.<br/>(dot)<br/>or integer<br/>&ge;1</td>
     *
     *          <td>we advance once ('.' or 1) or several times (2, 3, etc) the stack top (left on the stack or "peeked").
     *              The new spliterator is pushed followed by the original one.
     *              <p/> Fails if there is no splitterator or if it has no enough remaining elements.
     *                              the spliterator is push back as it was.</td>
     *      <td style="padding-left: 10px"><code>A = peek;<br/> A.tryAdvance;</td></tr>
     *      <tr><td style="text-align: center">*<br/>(star)</td>
     *          <td>we advance all the way to the end of the top (not popped) spliterator. <p/>Fails if there are no more spliterators on the stack.</td>
     *          <td style="padding-left: 10px"><code>
     *              A = peek;<br/>
     *              A.forEachRemaining();
     *          </code></td>
     *      </tr>
     *      <tr><td style="text-align: center">$<br/>(dollar)</td>
     *                <td>we pop an empty spliterator from the stack
     *                <p/>Fails if there are no more spliterators on the stack or the top one is not empty</td>
     *                <td style="padding-left: 10px"><code>
     *                    A = pop;<br/>
     *                    A.isEmpty or fail;
     *                </code></td>
     *            </tr>
     *     </table>
     * </p>
     *
     * <p>
     *     The test succeds if on operation failed and it the end we have an empty stack.
     * </p>
     *
     * <h3>Examples</h3>
     * <ul>
     *     <li>with a list of arbitrary length we simply go thru all elements with a for-each-remaining (don't forget to pop the empty spliterator at the end: <code>*$</code></li>
     *     <li>with a list with 5 elements we advance thru all of the one at a time: <code>.....$</code> or <code>5$</code> or <code>.2.1$</code></li>
     *     <li>with a list of arbitrary length we split once at the beginning and we read all elements of both splits: <code>^*$*$</code></li>
     *     <li>Assuming the list has 10 elements we split once and we advace exactly 5 explicitly for each split: <code>^5$5$</code></li>
     *     <li>We have a 10 element list we split twice resulting in ((3, 2), 5) splitertors, the we advance only the first element of each and we do a foreach with the reminding: <code>^^.*$.*$.*$</code></li>
     *     <li>Same as above but we advance thru explicitly using . or numbers: <code>^^..$.2$.3.$</code>
     *     <li>With a list with 5 or more elements we avance 3 then we split and  read the rest: <code>...^*$*$</code> or <code>3^*$*$</code></li>
     * </ul>
     * @param length the lenght of the list to traverse.
     * @param operations the operation to perform on it.
     */
    @Test(dataProvider="lengthAndOperationsSequenceData")
    public void testOperations(final int length, final String operations) {
        final List<Integer> subject = firstNaturalIntegers(length);
        final Set<Integer> visited = new HashSet<>();
        final Stack<Spliterator<Integer>> spliteratorStack = new Stack<>();
        spliteratorStack.push(new InterleavingListSpliterator<>(subject));
        for (int pos = 0; pos < operations.length();) {
            char code = operations.charAt(pos++);
            if (Character.isDigit(code)) {
                int times = code - '0';
                while (Character.isDigit(code = operations.charAt(pos))) {
                    times = times * 10 + (code - '0');
                    pos++;
                }
                performAdvance(visited, spliteratorStack, times);
            } else {
                final Spliterator<Integer> top = spliteratorStack.peek();
                switch (Operation.of(code)) {
                    case SPLIT:
                        Assert.assertFalse(spliteratorStack.isEmpty());
                        final long topPrevSize = top.estimateSize();
                        final Spliterator<Integer> other = top.trySplit();
                        Assert.assertEquals(topPrevSize > 1, other != null);
                        Assert.assertEquals((other == null ? 0 : other.estimateSize()) + top.estimateSize(), topPrevSize);
                        if (other != null) {
                            Assert.assertEquals(((top.estimateSize() - other.estimateSize()) & ~0x01), 0, "unbalanced split"); // checking this subtraction is 0 or 1 by checing every other bit except the first one.
                            spliteratorStack.pop();
                            spliteratorStack.push(other);
                            spliteratorStack.push(top);
                        }
                        break;
                    case ADVANCE:
                        performAdvance(visited, spliteratorStack, 1);
                        break;
                    case READTHRU:
                        top.forEachRemaining(i -> {
                            Assert.assertTrue(visited.add(i), "element visited twice: " + i);
                        });
                        Assert.assertEquals(top.estimateSize(), 0);
                        break;
                    case POPEMPTY:
                        Assert.assertEquals(spliteratorStack.peek().estimateSize(), 0);
                        spliteratorStack.pop();
                        break;
                }
            }
        }
        Assert.assertTrue(spliteratorStack.isEmpty());
        Assert.assertEquals(visited.size(), length);
        Assert.assertTrue(visited.containsAll(subject));
    }

    private void performAdvance(Set<Integer> visited, Stack<Spliterator<Integer>> spliteratorStack, int times) {
        long prevSize = spliteratorStack.peek().estimateSize();
        if (prevSize < times) {
            Assert.fail("advance beyond end");
        }
        for (int i = 0; i < times; i++) {
            spliteratorStack.peek().tryAdvance(val -> {
                Assert.assertNotNull(val, "bad argument, null not expected");
                Assert.assertTrue(visited.add(val), "same element added twice: " + val);
            });
            long newSize = spliteratorStack.peek().estimateSize();
            Assert.assertEquals(newSize, prevSize - 1);
            prevSize = newSize;
        }
    }

    /**
     * returns a list from 0 to N-1
     * @param N the lenght of the list.
     * @return never {@code null}.
     */
    private List<Integer> firstNaturalIntegers(final int N) {
        return new AbstractList<Integer>() {
            @Override
            public Integer get(final int index) {
                if (index < 0 || index >= N) {
                    throw new IndexOutOfBoundsException("" + index);
                }
                return index;
            }

            @Override
            public int size() {
                return N;
            }
        };
    }

    @DataProvider(name = "lengthAndSplitLevelsData")
    private Object[][] lengthAndSplitLevels() {
        final List<Object[]> result = new ArrayList<>();
        final int[] sizes = {0, 1, 11, 13, 10, 100, 2313};
        final int[] splitLevels = {0, 1, 2, 3, 4};
        for (final int size : sizes) {
            for (final int splitLevel : splitLevels) {
                result.add(new Object[] {size, splitLevel});
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "lengthAndNumberOfThreadsData")
    private Object[][] lengthAndNumberOfThreadsData() {
        final List<Object[]> result = new ArrayList<>();
        final int[] sizes = {0, 500, 517, 31, 10001};
        final int[] splitLevels = {0, 1, 2, 3, 4, 5};
        for (final int size : sizes) {
            for (final int splitLevel : splitLevels) {
                result.add(new Object[] {size, splitLevel});
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider
    private Object[][] lengthAndOperationsSequenceData() {
        return new Object[][]{
                {0, "$"}, {0, "*$"}, {0, "***$"},
                {0, "^$"}, {0, "^^^^^^^^^^$"}, // splits on an unesplitable spliterator don't do anything.
                {10, "10$"}, {10, "*$"}, {10, "...*$"}, {10, "...5..$"},
                {10, "^5$5$"}, {10, "^5$*$"}, {10, "^*$*$"},
                {11, "^6$5$"}, {11, "^^3$3$5$"},
                // 20 elm split into ^(10,10) advance 2 to 2^(8,10)
                // then split the 8 into 2x4 to 2^(4, 4, 10)
                // read the first 4 and split the other 4 in 2x2 (
                // 6^(2, 2, 10) -> 8^(2, 10) -> 10^10 -> 20^
                {20, "^..^4$^2$2$10$"},
                {32, "^^^^^.$.$..$....$8$16$"},
                {32, "^16$^8$^4$^..$^.$.$"}
        };
    }
}
