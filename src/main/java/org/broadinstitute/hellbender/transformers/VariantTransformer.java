package org.broadinstitute.hellbender.transformers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.SerializableFunction;

import java.util.Objects;
import java.util.function.UnaryOperator;

/**
 * Classes which perform transformations from {@link htsjdk.variant.variantcontext.VariantContext} -> {@link htsjdk.variant.variantcontext.VariantContext}
 * should implement this interface by overriding {@link SerializableFunction < VariantContext ,VariantContext>#apply(VariantContext)}
 * Created by jonn on 6/26/18.
 */
@FunctionalInterface
public interface VariantTransformer extends UnaryOperator<VariantContext>, SerializableFunction<VariantContext, VariantContext> {
    public static final long serialVersionUID = 1L;

    //HACK: These methods are a hack to get to get the type system to accept compositions of VariantTransformers.
    @SuppressWarnings("overloads")
    default VariantTransformer andThen(final VariantTransformer after) {
        Objects.requireNonNull(after);
        return (VariantContext vc) -> after.apply(apply(vc));
    }

    @SuppressWarnings("overloads")
    default VariantTransformer compose(final VariantTransformer before) {
        Objects.requireNonNull(before);
        return (VariantContext vc) -> apply(before.apply(vc));
    }

    static VariantTransformer identity(){
        return variantContext -> variantContext;
    }
}
