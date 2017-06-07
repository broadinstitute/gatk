package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A Duplicable wrapper around INDArray
 *
 * NOTE: INDArray already has the dup() method. We simply refer to that.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class DuplicableNDArray implements Duplicable {

    final private INDArray value;

    public DuplicableNDArray() {
        this.value = null;
    }

    public DuplicableNDArray(final INDArray value) {
        this.value = value;
    }

    @Override
    public DuplicableNDArray duplicate() {
        return value == null
                ? new DuplicableNDArray()
                : new DuplicableNDArray(value.dup());
    }

    @Override
    public boolean hasValue() {
        return value != null;
    }

    @Override
    public INDArray value() {
        return value;
    }

    @Override
    public String toString() {
        return value == null ? "null" : value.toString();
    }
}