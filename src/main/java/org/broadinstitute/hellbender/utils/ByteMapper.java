package org.broadinstitute.hellbender.utils;

/**
 * Util to map between byte values using a file based lookup table
 *
 *  The lookup table file should contain a single mapping on each line, using the following format
 *  from1 to1
 *  from2 to2
 *  ...
 *  * defautlTo (optional)
 *
 *  for example
 *  0 0
 *  1 0
 *  2 0
 *  3 3
 *  4 3
 *  * 10
 *
 *  This will map 0->0, 1->0, 2->0, 3->3, 4->3 and all other values to 10
 */

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedHashMap;
import java.util.Map;

public class ByteMapper {

    // locals
    private final Map<Byte, Byte> mapping = new LinkedHashMap<>();
    private Byte defaultMapping;

    public ByteMapper(final GATKPath path) {

        // get mapping from file
        try (final BufferedReader reader = new BufferedReader(new InputStreamReader(path.getInputStream()))) {

            String line;
            while ( (line = reader.readLine()) != null ) {
                line = line.trim();
                if ( line.length() == 0 || line.startsWith("#") ) {
                    continue;
                }
                final String toks[] = line.split(" ");
                Utils.validate(toks.length == 2, "each line should contain exactly two tokens separated by a single space");
                if ( toks[0].equals("*") ) {
                    defaultMapping = Byte.parseByte(toks[1]);
                } else {
                    mapping.put(Byte.parseByte(toks[0]), Byte.parseByte(toks[1]));
                }

            }
        } catch (IOException e) {
            throw new GATKException("", e);
        }
    }

    public Byte map(Byte from) {

        if ( mapping.containsKey(from) ) {
            return mapping.get(from);
        } else if ( defaultMapping != null ) {
            return defaultMapping;
        } else {
            return from;
        }
    }
}
