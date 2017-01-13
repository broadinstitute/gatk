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
    public DuplicableNDArray deepCopy() {
        return new DuplicableNDArray(value.dup());
    }

    @Override
    public boolean isNull() {
        return value == null;
    }
    public INDArray value() {
        return value;
    }

    public static INDArray of(final Duplicable obj) {
        if (obj == null) {
            throw new NullPointerException("The input duplicable object is null.");
        }
        if (obj instanceof DuplicableNDArray) {
            return ((DuplicableNDArray)obj).value();
        } else {
            throw new ClassCastException("Can not cast " + obj + " to an INDArray.");
        }
    }

    @Override
    public String toString() {
        if (value == null) {
            return "null";
        } else {
            return value.toString();
        }
    }

}