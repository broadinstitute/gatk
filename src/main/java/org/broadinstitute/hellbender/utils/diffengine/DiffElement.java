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

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 12:55 PM
 *
 * An interface that must be implemented to allow us to calculate differences
 * between structured objects
 */
@Invariant({
        "name != null",
        "value != null",
        "parent != null || name.equals(\"ROOT\")",
        "value == null || value.getBinding() == this"})
public class DiffElement {
    public final static DiffElement ROOT = new DiffElement();

    final private String name;
    final private DiffElement parent;
    final private DiffValue value;

    /**
     * For ROOT only
     */
    private DiffElement() {
        this.name = "ROOT";
        this.parent = null;
        this.value = new DiffValue(this, "ROOT");
    }

    @Requires({"name != null", "parent != null", "value != null"})
    public DiffElement(String name, DiffElement parent, DiffValue value) {
        if ( name.equals("ROOT") ) throw new IllegalArgumentException("Cannot use reserved name ROOT");
        this.name = name;
        this.parent = parent;
        this.value = value;
        this.value.setBinding(this);
    }

    @Ensures({"result != null"})
    public String getName() {
        return name;
    }

    public DiffElement getParent() {
        return parent;
    }

    @Ensures({"result != null"})
    public DiffValue getValue() {
        return value;
    }

    public boolean isRoot() { return this == ROOT; }

    @Ensures({"result != null"})
    @Override
    public String toString() {
        return getName() + "=" + getValue().toString();
    }

    public String toString(int offset) {
        return (offset > 0 ? Utils.dupString(' ', offset) : 0) + getName() + "=" + getValue().toString(offset);
    }

    @Ensures({"result != null"})
    public final String fullyQualifiedName() {
        if ( isRoot() )
            return "";
        else if ( parent.isRoot() )
            return name;
        else
            return parent.fullyQualifiedName() + "." + name;
    }

    @Ensures({"result != null"})
    public String toOneLineString() {
        return getName() + "=" + getValue().toOneLineString();
    }

    @Ensures({"result != null"})
    public DiffNode getValueAsNode() {
        if ( getValue().isCompound() )
            return (DiffNode)getValue();
        else
            throw new ReviewedGATKException("Illegal request conversion of a DiffValue into a DiffNode: " + this);
    }

    public int size() {
        return 1 + getValue().size();
    }
}
