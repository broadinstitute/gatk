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

package org.broadinstitute.hellbender.utils.diffengine;

public class Difference implements Comparable<Difference> {
    final String path; // X.Y.Z
    final String[] parts;
    int count = 1;
    DiffElement master = null , test = null;

    public Difference(String path) {
        this.path = path;
        this.parts = DiffEngine.diffNameToPath(path);
    }

    public Difference(DiffElement master, DiffElement test) {
        this(createPath(master, test), master, test);
    }

    public Difference(String path, DiffElement master, DiffElement test) {
        this(path);
        this.master = master;
        this.test = test;
    }

    public String[] getParts() {
        return parts;
    }

    public void incCount() { count++; }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    /**
     * The fully qualified path object A.B.C etc
     * @return
     */
    public String getPath() {
        return path;
    }

    /**
     * @return the length of the parts of this summary
     */
    public int length() {
        return this.parts.length;
    }

    /**
     * Returns true if the string parts matches this summary.  Matches are
     * must be equal() everywhere where this summary isn't *.
     * @param otherParts
     * @return
     */
    public boolean matches(String[] otherParts) {
        if ( otherParts.length != length() )
            return false;

        // TODO optimization: can start at right most non-star element
        for ( int i = 0; i < length(); i++ ) {
            String part = parts[i];
            if ( ! part.equals("*") && ! part.equals(otherParts[i]) )
                return false;
        }

        return true;
    }

    @Override
    public String toString() {
        return String.format("%s:%d:%s", getPath(), getCount(), valueDiffString());
    }

    @Override
    public int compareTo(Difference other) {
        // sort first highest to lowest count, then by lowest to highest path
        int countCmp = Integer.valueOf(count).compareTo(other.count);
        return countCmp != 0 ? -1 * countCmp : path.compareTo(other.path);
    }

    public String valueDiffString() {
        if ( hasSpecificDifference() ) {
            return String.format("%s!=%s", getOneLineString(master), getOneLineString(test));
        } else {
            return "N/A";
        }
    }

    private static String createPath(DiffElement master, DiffElement test) {
        return (master == null ? test : master).fullyQualifiedName();
    }

    private static String getOneLineString(DiffElement elt) {
        return elt == null ? "MISSING" : elt.getValue().toOneLineString();
    }

    public boolean hasSpecificDifference() {
        return master != null || test != null;
    }

    public DiffElement getMaster() {
        return master;
    }

    public DiffElement getTest() {
        return test;
    }
}
