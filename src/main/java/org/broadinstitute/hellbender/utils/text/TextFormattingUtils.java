/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.text;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Common utilities for dealing with text formatting.
 */
public final class TextFormattingUtils {
    private TextFormattingUtils(){}

    /**
     * The default line width, for GATK output written to the screen.
     */
    public static final int DEFAULT_LINE_WIDTH = 120;

    /**
     * Simple implementation of word-wrap for a line of text.  Idea and
     * regexp shamelessly stolen from http://joust.kano.net/weblog/archives/000060.html.
     * Regexp can probably be simplified for our application.
     * @param text Text to wrap.
     * @param width Maximum line width.
     * @return A list of word-wrapped lines.
     */
    public static List<String> wordWrap( String text, int width ) {
        Pattern wrapper = Pattern.compile(String.format(".{0,%d}(?:\\S(?:[\\s|]|$)|$)", width - 1));
        Matcher matcher = wrapper.matcher( text );

        List<String> wrapped = new ArrayList<>();
        while( matcher.find() ) {
            // Regular expression is supersensitive to whitespace.
            // Assert that content is present before adding the line.
            String line = matcher.group().trim();
            if( line.length() > 0 )
                wrapped.add( matcher.group() );
        }
        return wrapped;
    }

    /**
     * Returns the word starting positions within line, excluding the first position 0.
     * The returned list is compatible with splitFixedWidth.
     * @param line Text to parse.
     * @return the word starting positions within line, excluding the first position 0.
     */
    public static List<Integer> getWordStarts(String line) {
        if (line == null)
            throw new GATKException("line is null");
        List<Integer> starts = new ArrayList<>();
        int stop = line.length();
        for (int i = 1; i < stop; i++)
            if (Character.isWhitespace(line.charAt(i - 1)))
                if(!Character.isWhitespace(line.charAt(i)))
                    starts.add(i);
        return starts;
    }

    /**
     * Parses a fixed width line of text.
     * @param line Text to parse.
     * @param columnStarts the column starting positions within line, excluding the first position 0.
     * @return The parsed string array with each entry trimmed.
     */
    public static String[] splitFixedWidth(String line, List<Integer> columnStarts) {
        if (line == null)
            throw new GATKException("line is null");
        if (columnStarts == null)
            throw new GATKException("columnStarts is null");
        int startCount = columnStarts.size();
        String[] row = new String[startCount + 1];
        if (startCount == 0) {
            row[0] = line.trim();
        } else {
            row[0] = line.substring(0, columnStarts.get(0)).trim();
            for (int i = 1; i < startCount; i++)
                row[i] = line.substring(columnStarts.get(i - 1), columnStarts.get(i)).trim();
            row[startCount] = line.substring(columnStarts.get(startCount - 1)).trim();
        }
        return row;
    }

    /**
     * Parses a line of text by whitespace.
     * @param line Text to parse.
     * @return The parsed string array.
     */
    public static String[] splitWhiteSpace(String line) {
        if (line == null)
            throw new GATKException("line is null");
        return line.trim().split("\\s+");
    }
}
