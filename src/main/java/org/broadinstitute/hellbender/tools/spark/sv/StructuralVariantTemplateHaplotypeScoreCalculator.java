package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.List;

/**
 * Created by valentin on 5/18/17.
 */
public interface StructuralVariantTemplateHaplotypeScoreCalculator {

    default TemplateHaplotypeScoreTable calculate(final List<Template> templates, final List<Haplotype> haplotypes, final InsertSizeDistribution dist) {
        final TemplateHaplotypeScoreTable table = new TemplateHaplotypeScoreTable(templates, haplotypes);
        calculate(table);
        return table;
    }

    void calculate(final TemplateHaplotypeScoreTable table);
}
