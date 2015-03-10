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