package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * Represents one set of cycles in an ReadStructure (e.g. if the ReadStructure is 36TB836T then
 * 36T, 8B, and 36T are invidually represented internally as a ReadDescriptor).
 */
public final class ReadDescriptor {
    public final int length;
    public final ReadType type;

    public ReadDescriptor(final int length, final ReadType type) {
        this.length = length;
        this.type = type;
    }

    @Override
    public String toString() {
        return this.length + this.type.name();
    }

    public boolean equals(final Object other) {
        if (this == other) return true;
        if (other.getClass() != this.getClass()) return false;

        final ReadDescriptor that = (ReadDescriptor) other;
        return this.length == that.length && this.type == that.type;
    }

    @Override
    public int hashCode() {
        return 31 * this.type.ordinal() + 379 * length;
    }
}
