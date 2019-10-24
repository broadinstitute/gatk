package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.util.Objects;
import java.util.OptionalInt;
import java.util.function.Function;

public class VCFAnnotations {

    public static final class Builder {

        private VCFAnnotationTarget target;
        private String id;
        private VCFHeaderLineType type;
        private VCFHeaderLineCount count;
        private int numericNumber = -1;
        private String description = "";
        private String format = "";

        private Builder() {}

        public Builder target(final VCFAnnotationTarget target) {
            this.target = target;
            return this;
        }

        public Builder description(final String description) {
            this.description = Objects.requireNonNull(description);
            return this;
        }

        public Builder type(final VCFHeaderLineType type) {
            this.type = type;
            return this;
        }

        public Builder number(final int value) {
            if ((numericNumber = value) < 0) {
                throw new IllegalArgumentException("bad numeric value " + value);
            }
            count = VCFHeaderLineCount.INTEGER;
            return this;
        }

        public Builder number(final VCFHeaderLineCount number) {
            Objects.requireNonNull(number);
            if ((count = number) == VCFHeaderLineCount.INTEGER) {
                if (type == VCFHeaderLineType.Flag) {
                    throw new IllegalStateException("cannot set number to anything but zero when type is Flag");
                }
                if (numericNumber < 0) {
                    numericNumber = 1;
                }
            }
            return this;
        }

        public VCFAnnotationMeta meta() {
            if (target == null) {
                throw new IllegalStateException("the target must have been specified");
            } else if (id == null) {
                throw new IllegalStateException("the id must have been specified");
            } else if (type == null) {
                throw new IllegalStateException("the type must have been specified");
            } else if (type == VCFHeaderLineType.Flag) {
                return VCFAnnotationMeta.flag(target, id, description);
            } else if (count == null) {
                throw new IllegalStateException("the count must have been specified");
            } else if (count == VCFHeaderLineCount.INTEGER) {
                return VCFAnnotationMeta.make(target, id, type, numericNumber, description, format);
            } else {
                return VCFAnnotationMeta.make(target, id, type, count, description, format);
            }
        }

        public <V extends VCFAnnotation<?>> V make(final Function<VCFAnnotationMeta, V> maker) {
            return maker.apply(meta());
        }
    }

    public static Builder builder() {
        return new Builder();
    }
}
