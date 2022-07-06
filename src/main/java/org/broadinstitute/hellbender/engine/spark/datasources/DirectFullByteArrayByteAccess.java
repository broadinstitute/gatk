package org.broadinstitute.hellbender.engine.spark.datasources;

import org.bdgenomics.utils.io.ByteArrayByteAccess;

/**
 * A version of org.bdgenomics.utils.io.ByteArrayByteAccess that makes no copies of byte array used for initialization.
 * However DirectFullByteArrayByteAccess.readFully can only return a reference to the full underlying byte array.
 * Therefore, the user should exercise caution that the underlying data does not get mutated.
 */
class DirectFullByteArrayByteAccess extends ByteArrayByteAccess {
    private static final long serialVersionUID = 1L;

    DirectFullByteArrayByteAccess(byte[] bytes) {
        super(bytes);
    }

    @Override
    public byte[] readFully(long offset, int length) {
        if ((offset != 0) || (length != this.length())) {
            throw new IllegalArgumentException("readFully can only return a reference to the full underlying byte array");
        }
        return this.bytes();
    }
}
