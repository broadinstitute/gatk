package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.function.Predicate;

@FunctionalInterface
public interface VariantFilter extends Predicate<VariantContext>{

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadFilters.
    default VariantFilter and(VariantFilter filter ) { return Predicate.super.and(filter)::test; }

    default VariantFilter or(VariantFilter filter ) {
        return Predicate.super.or(filter)::test;
    }

    default VariantFilter negate(){
        return Predicate.super.negate()::test;
    }}
