package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.readers.LineIterator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit test for {@link LineIteratorReader}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class LineIteratorReaderUnitTest {

    private static final Random RANDOM = new Random();

    private static final int[] TEXT_SIZE = new int[] { 1, 100, 200, 400, 500 };

    private static final int[] LINE_SIZE = new int[] { 10, 11, 34, 60 };

    private static final int[] BUFFER_SIZE = new int[] { 10, 100 };


    @Test(dataProvider="testParameter")
    public void testRead(@SuppressWarnings("unused") final int textSize, @SuppressWarnings("unused") final int lineSize, final int bufferSize, final String text) {
        final LineIterator lineIterator = new TestLineIterator(text);
        final LineIteratorReader reader = new LineIteratorReader(lineIterator);
        final StringBuilder sb = new StringBuilder();
        final char[] buffer = new char[bufferSize];
        int readResult;
        while ((readResult = reader.read(buffer, 0, bufferSize)) >= 0) {
            sb.append(buffer, 0, readResult);
        }
        final String actualText = sb.toString();
        Assert.assertEquals(actualText, text);
    }

    @DataProvider(name="testParameter")
    private Object[][] testParameters() {
        final List<Object[]> result = new ArrayList<>();
        for (final int textSize : TEXT_SIZE) {
            for (final int lineSize : LINE_SIZE) {
                for (final int bufferSize : BUFFER_SIZE) {
                    result.add(new Object[] { textSize, lineSize, bufferSize, createText(textSize, lineSize, true) });
                }
            }
        }
        return result.toArray(new Object[result.size()][]);
    }
    private static String createText(final int textSize, final int lineSize, final boolean finishWithNewLine) {
        final int fullLineCount = textSize % lineSize;
        final int shortLineSize = textSize - lineSize * fullLineCount;
        final StringBuilder sb = new StringBuilder(textSize);
        for (int i = 0; i < fullLineCount; i++) {
            sb.append(createRandomLine(lineSize)).append('\n');
        }
        if (shortLineSize > 0) {
            sb.append(createRandomLine(shortLineSize)).append('\n');
        }
        if (!finishWithNewLine && sb.length() > 0) {
            sb.setLength(sb.length() - 1);
        } else if (finishWithNewLine && sb.length() == 0) {
            sb.append('\n');
        }
        return sb.toString();
    }

    private static char[] createRandomLine(final int lineSize) {
        final char[] result = new char[lineSize];
        for (int i = 0; i < result.length; i++) {
            result[i] = (char) ((RANDOM.nextInt(127 - 32)) + 32);
        }
        return result;
    }

    private class TestLineIterator implements LineIterator {
        private final ListIterator<String> lines;
        public TestLineIterator(final String text) {
            lines = Arrays.asList(text.split("\\n")).listIterator();
        }

        @Override
        public String peek() {
            final String result = lines.next();
            lines.previous();
            return result;
        }

        @Override
        public boolean hasNext() {
            return lines.hasNext();
        }

        @Override
        public String next() {
            return lines.next();
        }
    }
}
