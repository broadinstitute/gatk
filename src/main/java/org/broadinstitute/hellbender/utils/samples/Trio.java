package org.broadinstitute.hellbender.utils.samples;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * A class for imposing a trio structure on three samples; a common paradigm
 */
public final class Trio {
    private final Sample mother;
    private final Sample father;
    private final Sample child;

    public Trio(final Sample mom, final Sample dad, final Sample child) {
        Utils.nonNull(mom, "mom is null");
        Utils.nonNull(dad, "dad is null");
        Utils.nonNull(child, "child is null");
        if (! mom.getID().equals(child.getMaternalID())){
            throw new IllegalArgumentException("incorrect mother");
        }
        if (! dad.getID().equals(child.getPaternalID())){
            throw new IllegalArgumentException("incorrect father");
        }
        this.mother = mom;
        this.father = dad;
        this.child = child;
    }

    public Sample getMother() {
        return mother;
    }

    public String getMaternalID() {
        return mother.getID();
    }

    public Sample getFather() {
        return father;
    }

    public String getPaternalID() {
        return father.getID();
    }

    public Sample getChild() {
        return child;
    }

    public String getChildID() {
        return child.getID();
    }
}
