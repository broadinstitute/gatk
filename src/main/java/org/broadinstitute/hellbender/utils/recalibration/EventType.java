package org.broadinstitute.hellbender.utils.recalibration;

public enum EventType {
    BASE_SUBSTITUTION("M", "Base Substitution"),
    BASE_INSERTION("I", "Base Insertion"),
    BASE_DELETION("D", "Base Deletion");

    private final String representation;
    private final String longRepresentation;

    /**
     * Returns a cached value of the EventType.values() call (more precisely - an unmodifiable list view of that array).
     *
     * Every call to EventType.values() (or any enum type) creates a new array instance but they are all equal (ie contain identical elements).
     * This is very expensive and wasteful when this array is created billions of times as in the case of BQSR.
     *
     * The solution is to create it once and reuse.
     * However, we can't expose this array in an API because we can't make an array immutable.
     * Exposing this array as list also does not work because performance of Collections.UnmodifiableCollection.iterator() is very bad and ruins our performance.
     * The solution is to expose this array via read only calls and have clients iterate explicitly.
     */
    private static final EventType[] cachedValues = EventType.values();

    private EventType(final String representation, final String longRepresentation) {
        this.representation = representation;
        this.longRepresentation = longRepresentation;
    }

    /**
     * Get the EventType corresponding to its ordinal index
     * @param index an ordinal index
     * @return the event type corresponding to ordinal index
     */
    public static EventType eventFrom(final int index) {
        return cachedValues[index];
    }

    /**
     * Get the EventType with short string representation
     * @throws IllegalArgumentException if representation doesn't correspond to one of EventType
     * @param representation short string representation of the event
     * @return an EventType
     */
    public static EventType eventFrom(final String representation) {
        for (EventType eventType : cachedValues) {
            if (eventType.representation.equals(representation)) {
                return eventType;
            }
        }

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