package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.util.Map;
import java.util.Objects;

public final class VCFAnnotationMeta {

    private final VCFAnnotationTarget target;
    private final String id;
    private final VCFHeaderLineType type;
    private final VCFHeaderLineCount count;
    private final int number;
    private final String description;
    private final String format;


    private  VCFAnnotationMeta(final VCFAnnotationTarget target,
                                         final String id,
                                         final VCFHeaderLineType type,
                                         final VCFHeaderLineCount count,
                                         final int number,
                                         final String description,
                                         final String format) {
        this.target = target;
        this.id = id;
        this.type = type;
        this.count = count;
        this.number = number;
        this.description = description;
        this.format = format;
    }

    public static VCFAnnotationMeta flag(final VCFAnnotationTarget target, final String id, final String description) {
        Objects.requireNonNull(target);
        Objects.requireNonNull(id);
        Objects.requireNonNull(description);
        return new VCFAnnotationMeta(target, id, VCFHeaderLineType.Flag, VCFHeaderLineCount.INTEGER, 0, description, "");
    }

    public static VCFAnnotationMeta make(final VCFAnnotationTarget target, final String id, final VCFHeaderLineType type, final int integerNumber, final String description, final String format) {
        Objects.requireNonNull(target);
        Objects.requireNonNull(id);
        Objects.requireNonNull(description);
        Objects.requireNonNull(type);
        Objects.requireNonNull(format);
        if (integerNumber <= 0) {
            throw new IllegalArgumentException("the integer number/count must be 1 or greater");
        } else if (type == VCFHeaderLineType.Flag) {
            throw new IllegalArgumentException("when a integer number is provided the type cannot be Flag");
        }
        return new VCFAnnotationMeta(target, id, type, VCFHeaderLineCount.INTEGER, integerNumber, description, format);
    }

    public static VCFAnnotationMeta make(final VCFAnnotationTarget target, final String id, final VCFHeaderLineType type, final VCFHeaderLineCount number, final String description, final String format) {
        Objects.requireNonNull(target);
        Objects.requireNonNull(id);
        Objects.requireNonNull(description);
        Objects.requireNonNull(type);
        Objects.requireNonNull(number);
        Objects.requireNonNull(format);
        if (number != VCFHeaderLineCount.INTEGER) {
            throw new IllegalArgumentException("the integer number/count cannot be INTEGER");
        } else if (type == VCFHeaderLineType.Flag) {
            throw new IllegalArgumentException("when a number is provided the type cannot be Flag");
        }
        return new VCFAnnotationMeta(target, id, type, number, -1, description, format);
    }

    public VCFAnnotationTarget target() {
        return target;
    }

    public String id() {
        return id;
    }

    public int integerNumber() {
        return number;
    }

    public VCFHeaderLineCount number() {
        return count;
    }

    public String description() {
        return description;
    }

    public String format()  {
        return format;
    }

    public void requiresType(final VCFHeaderLineType type) {
        if (this.type != type) {
            throw new IllegalArgumentException("invalid type " + this.type + " when requires " + type);
        }
    }
}
