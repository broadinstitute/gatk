package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.util.Lazy;
import htsjdk.variant.variantcontext.VariantContextUtils;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.function.BiFunction;

/**
 * Keep only reads that the attributes of meet a given set of jexl expressions
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_READFILTERS, groupSummary = HelpConstants.DOC_CAT_READFILTERS_SUMMARY,
        summary = "Keep only reads that meet all given jexl expressions (on their attributes)")
public final class JexlExpressionReadTagValueFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;
    private static final Logger logger = LogManager.getLogger(JexlExpressionReadTagValueFilter.class);

    @Argument(fullName=ReadFilterArgumentDefinitions.READ_FILTER_EXPRESSION_LONG_NAME, shortName="filter", doc="One or more JEXL expressions used to filter", optional=false)
    public List<String> filterExpressions = new ArrayList<>();

    private Lazy<List<Expression>> jexlExprs = new Lazy<>(() -> {
        List<Expression>        l = new LinkedList<>();
        for ( String expr : filterExpressions ) {
            final Expression    jexl =  VariantContextUtils.engine.get().createExpression(expr);
            logger.info("created jexl: " + jexl);
            l.add(jexl);
        }
        return l;
    });

    private static class GATKReadJexlContext implements JexlContext {

        final private GATKRead        read;

        GATKReadJexlContext(final GATKRead read) {
            this.read = read;
        }

        @Override
        public Object get(final String name) {
            return read.getAttributeAsString(name);
        }

        @Override
        public void set(final String name, final Object value) {
            throw new IllegalArgumentException("setting attributes is not allowed");
        }

        @Override
        public boolean has(final String name) {
            return read.hasAttribute(name);
        }
    }

    public JexlExpressionReadTagValueFilter() {
    }

    // convenience constructor for using a single jexl expression
    public JexlExpressionReadTagValueFilter(final String jexlExpr) {
        this.filterExpressions = Collections.singletonList(jexlExpr);
    }

    // convenience constructor for using a multiple jexl expressions
    public JexlExpressionReadTagValueFilter(final List<String> jexlExprs) {
        this.filterExpressions = jexlExprs;
    }

    @Override
    public boolean test(final GATKRead read) {

        // loop over expressions. At this point expressions are ANDed
        for ( Expression expr : jexlExprs.get() ) {
            Object v = expr.evaluate(new GATKReadJexlContext(read));
            if (!v.equals(Boolean.TRUE)) {
                return false;
            }
        }

        // if here, all expressions matched
        return true;
    }
}
