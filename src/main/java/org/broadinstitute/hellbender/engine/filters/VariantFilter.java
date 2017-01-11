package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;

import java.io.Serializable;
import java.util.function.Predicate;

@FunctionalInterface
public interface VariantFilter extends Predicate<VariantContext>, Serializable {
    static final long serialVersionUID = 1L;

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadFilters.
    default VariantFilter and(VariantFilter filter ) { return Predicate.super.and(filter)::test; }

    default VariantFilter or(VariantFilter filter ) {
        return Predicate.super.or(filter)::test;
    }

    @Override
    default VariantFilter negate(){
        return Predicate.super.negate()::test;
    }}
