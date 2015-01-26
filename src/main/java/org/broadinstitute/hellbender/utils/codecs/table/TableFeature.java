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

package org.broadinstitute.hellbender.utils.codecs.table;


import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * A feature representing a single row out of a text table
 */
public class TableFeature implements Feature {
    // stores the values for the columns seperated out
    private final List<String> values;

    // if we have column names, we store them here
    private final List<String> keys;

    // our location
    private final GenomeLoc position;

    public TableFeature(GenomeLoc position, List<String> values, List<String> keys) {
        this.values = values;
        this.keys = keys;
        this.position = position;
    }

    @Override
    public String getChr() {
        return position.getContig();
    }

    @Override
    public int getStart() {
        return position.getStart();
    }

    @Override
    public int getEnd() {
        return position.getStop();
    }

    public String getValue(int columnPosition) {
        if (columnPosition >= values.size()) throw new IllegalArgumentException("We only have " + values.size() + "columns, the requested column = " + columnPosition);
        return values.get(columnPosition);
    }

    public String toString() {
        return String.format("%s\t%s",position.toString(), Utils.join("\t", values));
    }

    public String get(String columnName) {
        int position = keys.indexOf(columnName);
        if (position < 0) throw new IllegalArgumentException("We don't have a column named " + columnName);
        return values.get(position);
    }

    public GenomeLoc getLocation() {
        return this.position;
    }

    public List<String> getAllValues() {
        return getValuesTo(values.size());
    }

    public List<String> getValuesTo(int columnPosition) {
        return values.subList(0,columnPosition);
    }

    public List<String> getHeader() {
        return keys;
    }
}
