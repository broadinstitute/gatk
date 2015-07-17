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

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 12:55 PM
 *
 * An interface that must be implemented to allow us to calculate differences
 * between structured objects
 */
public class DiffValue {
    private DiffElement binding = null;
    final private Object value;

    public DiffValue(Object value) {
        this.value = value;
    }

    public DiffValue(DiffElement binding, Object value) {
        this.binding = binding;
        this.value = value;
    }

    public DiffValue(DiffValue parent, Object value) {
        this(parent.getBinding(), value);
    }

    public DiffValue(String name, DiffElement parent, Object value) {
        this.binding = new DiffElement(name, parent, this);
        this.value = value;
    }

    public DiffValue(String name, DiffValue parent, Object value) {
        this(name, parent.getBinding(), value);
    }

    public DiffElement getBinding() {
        return binding;
    }

    protected void setBinding(DiffElement binding) {
        this.binding = binding;
    }

    public Object getValue() {
        return value;
    }

    public String toString() {
        return getValue().toString();
    }

    public String toString(int offset) {
        return toString();
    }

    public String toOneLineString() {
        return getValue().toString();
    }

    public boolean isAtomic() { return true; }
    public boolean isCompound() { return ! isAtomic(); }
    public int size() { return 1; }
}
