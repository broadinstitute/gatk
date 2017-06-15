package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.io.Serializable;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by valentin on 6/9/17.
 */
public class TemplateMappingInformation implements Serializable {

    private static final long serialVersionUID = 1L;

    public String name;
    public final ReadPairOrientation pairOrientation;
    public final OptionalDouble firstAlignmentScore;
    public final OptionalDouble secondAlignmentScore;
    public final OptionalInt insertSize;

    public TemplateMappingInformation(final double firstAlignment, final double secondAlignment, final ReadPairOrientation orientation) {
        Utils.nonNull(firstAlignment);
        Utils.nonNull(secondAlignment);
        Utils.nonNull(orientation);
        if (orientation.isProper()) {
            throw new IllegalArgumentException("you cannot create a mapping information object with proper orientation without indicating the insert size");
        }
        firstAlignmentScore = OptionalDouble.of(firstAlignment);
        secondAlignmentScore = OptionalDouble.of(secondAlignment);
        pairOrientation = orientation;
        insertSize = OptionalInt.empty();
    }

    public TemplateMappingInformation(final double firstAlignment, final double secondAlignment, final int insertSize) {
        Utils.nonNull(firstAlignment);
        Utils.nonNull(secondAlignment);
        if (insertSize < 1) {
            throw new IllegalArgumentException("the input insert size cannot be less than 1");
        }
        firstAlignmentScore = OptionalDouble.of(firstAlignment);
        secondAlignmentScore = OptionalDouble.of(secondAlignment);
        this.insertSize = OptionalInt.of(insertSize);
        pairOrientation = ReadPairOrientation.PROPER;
    }

    public TemplateMappingInformation(final double firstAlignment) {
        Utils.nonNull(firstAlignment);
        firstAlignmentScore = OptionalDouble.of(firstAlignment);
        secondAlignmentScore = OptionalDouble.empty();
        insertSize = OptionalInt.empty();
        pairOrientation = ReadPairOrientation.XX;
    }

    public TemplateMappingInformation() {
        firstAlignmentScore = OptionalDouble.empty();
        secondAlignmentScore = OptionalDouble.empty();
        insertSize = OptionalInt.empty();
        pairOrientation = ReadPairOrientation.XX;
    }


}
