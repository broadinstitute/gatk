package org.broadinstitute.hellbender.utils.diffengine;

final class Difference implements Comparable<Difference> {
    private final String path; // X.Y.Z
    private final String[] parts;
    private int count = 1;
    private DiffElement master = null;
    private DiffElement test = null;

    public Difference(final String path) {
        this.path = path;
        this.parts = DiffEngine.diffNameToPath(path);
    }

    public Difference(final DiffElement master, final DiffElement test) {
        this(createPath(master, test), master, test);
    }

    public Difference(final String path, final DiffElement master, final DiffElement test) {
        this(path);
        this.master = master;
        this.test = test;
    }

    /**
     * Returns this difference's parts.
     */
    public String[] getParts() {
        return parts;
    }

    /**
     * Increases this difference's count.
     */
    public void incCount() { count++; }

    /**
     * Returns this difference's count.
     */
    public int getCount() {
        return count;
    }

    /**
     * Reset the count to 0.
     */
    public void resetCount() {
        this.count = 0;
    }

    /**
     * The fully qualified path object A.B.C etc
     */
    public String getPath() {
        return path;
    }

    /**
     * Returns the length of the parts of this summary.
     */
    public int length() {
        return this.parts.length;
    }

    /**
     * Returns true if the string parts matches this summary.  Matches are
     * must be equal() everywhere where this summary isn't *.
     */
    public boolean matches(final String[] otherParts) {
        if ( otherParts.length != length() ) {
            return false;
        }

        // TODO optimization: can start at right most non-star element
        for ( int i = 0; i < length(); i++ ) {
            final String part = parts[i];
            if ( ! part.equals("*") && ! part.equals(otherParts[i]) ) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String toString() {
        return String.format("%s:%d:%s", getPath(), getCount(), valueDiffString());
    }

    @Override
    public int compareTo(final Difference other) {
        // sort first highest to lowest count, then by lowest to highest path
        final int countCmp = Integer.valueOf(count).compareTo(other.count);
        return countCmp != 0 ? -1 * countCmp : path.compareTo(other.path);
    }

    public String valueDiffString() {
        if (master != null || test != null) {
            return String.format("%s!=%s", getOneLineString(master), getOneLineString(test));
        } else {
            return "N/A";
        }
    }

    private static String createPath(final DiffElement master, final DiffElement test) {
        return (master == null ? test : master).fullyQualifiedName();
    }

    private static String getOneLineString(final DiffElement elt) {
        return elt == null ? "MISSING" : elt.getValue().toOneLineString();
    }

    public DiffElement getMaster() {
        return master;
    }

    public DiffElement getTest() {
        return test;
    }
}
