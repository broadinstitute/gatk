package org.broadinstitute.hellbender.utils.samples;

import org.broadinstitute.hellbender.utils.Utils;

/**
 *  Represents an individual under study.
 */
public final class Sample {
    private final String familyID;
    private final String paternalID;
    private final String maternalID;
    private final Sex gender;
    private final String ID;

    public Sample(final String ID,
                  final String familyID,
                  final String paternalID,
                  final String maternalID,
                  final Sex gender) {
        Utils.nonNull(ID, "ID is null");
        Utils.nonNull(gender, "sex is null");
        this.ID = ID;
        this.familyID = familyID;
        this.paternalID = paternalID;
        this.maternalID = maternalID;
        this.gender = gender;
    }

    public String getID() {
        return ID;
    }

    public String getFamilyID() {
        return familyID;
    }

    public String getPaternalID() {
        return paternalID;
    }

    public String getMaternalID() {
        return maternalID;
    }

    public Sex getSex() {
        return gender;
    }

    @Override
    public String toString() {
        return String.format("Sample %s fam=%s dad=%s mom=%s gender=%s", ID, familyID, paternalID, maternalID, gender);
    }

    @Override
    public int hashCode() {
        return ID.hashCode();
    }

    @Override
    public boolean equals(final Object o) {
        if(o == null) {
            return false;
        }
        return (o instanceof Sample && ID.equals(((Sample)o).ID));
    }
}
