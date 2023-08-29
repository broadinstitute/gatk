package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.SnpEffUtil;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies variants as genes or coding regions, according to the effect modifier, as indicated by snpEff.
 * The 'gene' category includes category 'coding region', and additionally includes introns. 'Coding regions'
 * includes transcripts and, implicitly, UTRs.
 */
public class SnpEffPositionModifier extends VariantStratifier {

	public static final String EFFECT_KEY = "SNPEFF_EFFECT";

	public enum PositionModifier {
		GENE,          // EXON
		CODING_REGION, // CDS
		SPLICE_SITE,   // not a straight translation -- see getRelevantStates
		STOP_GAINED,   // STOP_GAINED
		STOP_LOST      // STOP_LOST
	}

	public SnpEffPositionModifier(VariantEvalEngine engine) {
		super(engine);

		for (final PositionModifier type : PositionModifier.values()) states.add(type.name());
	}

	@Override
	public List<Object> getRelevantStates(
			final VariantEvalContext context,
			final VariantContext comp,
			final String compName,
			final VariantContext eval,
			final String evalName,
			final String sampleName,
			final String FamilyName)
	{
		final List<Object> relevantStates = new ArrayList<>();
		if (eval != null && eval.isVariant() && eval.hasAttribute(EFFECT_KEY)) {
			final SnpEffUtil.EffectType effectType = SnpEffUtil.EffectType.valueOf(
					eval.getAttribute(EFFECT_KEY).toString());

			if (SnpEffUtil.isSubTypeOf(effectType, SnpEffUtil.EffectType.EXON))        relevantStates.add(PositionModifier.GENE.name());
			if (SnpEffUtil.isSubTypeOf(effectType, SnpEffUtil.EffectType.CDS))         relevantStates.add(PositionModifier.CODING_REGION.name());
			if (SnpEffUtil.isSubTypeOf(effectType, SnpEffUtil.EffectType.STOP_GAINED)) relevantStates.add(PositionModifier.STOP_GAINED.name());
			if (SnpEffUtil.isSubTypeOf(effectType, SnpEffUtil.EffectType.STOP_LOST))   relevantStates.add(PositionModifier.STOP_LOST.name());

			if (SnpEffUtil.isSubTypeOf(effectType, SnpEffUtil.EffectType.SPLICE_SITE_ACCEPTOR) ||
				SnpEffUtil.isSubTypeOf(effectType, SnpEffUtil.EffectType.SPLICE_SITE_DONOR))
					relevantStates.add(PositionModifier.SPLICE_SITE.name());
		}

		return relevantStates;
	}
}
