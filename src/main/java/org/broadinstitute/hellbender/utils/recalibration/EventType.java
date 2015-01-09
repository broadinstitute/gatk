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

package org.broadinstitute.hellbender.utils.recalibration;

public enum EventType {
    BASE_SUBSTITUTION("M", "Base Substitution"),
    BASE_INSERTION("I", "Base Insertion"),
    BASE_DELETION("D", "Base Deletion");

    private final String representation;
    private final String longRepresentation;

    private EventType(String representation, String longRepresentation) {
        this.representation = representation;
        this.longRepresentation = longRepresentation;
    }

    /**
     * Get the EventType corresponding to its ordinal index
     * @param index an ordinal index
     * @return the event type corresponding to ordinal index
     */
    public static EventType eventFrom(int index) {
        return EventType.values()[index];
    }

    /**
     * Get the EventType with short string representation
     * @throws IllegalArgumentException if representation doesn't correspond to one of EventType
     * @param representation short string representation of the event
     * @return an EventType
     */
    public static EventType eventFrom(String representation) {
        for (EventType eventType : EventType.values())
            if (eventType.representation.equals(representation))
                return eventType;

        throw new IllegalArgumentException(String.format("Event %s does not exist.", representation));
    }

    @Override
    public String toString() {
        return representation;
    }

    public String prettyPrint() {
        return longRepresentation;
    }
}