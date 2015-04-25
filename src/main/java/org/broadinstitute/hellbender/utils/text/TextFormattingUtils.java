package org.broadinstitute.hellbender.utils.text;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.List;

/**
 * Common utilities for dealing with text formatting.
 */
public final class TextFormattingUtils {
    private TextFormattingUtils(){}

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
