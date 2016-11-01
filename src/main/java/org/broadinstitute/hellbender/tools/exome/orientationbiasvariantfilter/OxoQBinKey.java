package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;


import scala.Tuple2;

import java.io.Serializable;

public class OxoQBinKey implements Serializable {
    static final long serialVersionUID = 33733733742L;

    private ArtifactMode artifactMode;
    private Tuple2<Byte, Byte> prePostBases;
    private boolean isOriginalMolecule;

    public OxoQBinKey(ArtifactMode artifactMode, Tuple2<Byte, Byte> prePostBases, boolean isOriginalMolecule) {
        this.artifactMode = artifactMode;
        this.prePostBases = prePostBases;
        this.isOriginalMolecule = isOriginalMolecule;
    }

    public ArtifactMode getArtifactMode() {
        return artifactMode;
    }

    public Tuple2<Byte, Byte> getPrePostBases() {
        return prePostBases;
    }

    public boolean isOriginalMolecule() {
        return isOriginalMolecule;
    }

    @Override
    public String toString() {
        return (char) artifactMode.getRef() + ">" + (char) artifactMode.getAlt() + "  context: " +
                (char) prePostBases._1().byteValue() + "_" + (char) prePostBases._2().byteValue() + "   is original molecule: " + isOriginalMolecule;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        final OxoQBinKey that = (OxoQBinKey) obj;
        return (this.getArtifactMode().equals(that.getArtifactMode())) && (this.getPrePostBases().equals((that.getPrePostBases())) &&
                (this.isOriginalMolecule() == that.isOriginalMolecule()));
    }
}
