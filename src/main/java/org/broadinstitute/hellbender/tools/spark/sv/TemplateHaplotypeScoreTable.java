package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.base.Functions;
import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.io.Serializable;
import java.util.*;

/**
 * Created by valentin on 5/18/17.
 */
public class TemplateHaplotypeScoreTable implements Serializable {

    private static long serialVersionUID = 1L;

    private final double[][] values;

    private final List<Template> templates;

    private final Map<String, Integer> templateIndex;

    private final List<Haplotype> haplotypes;

    public TemplateHaplotypeScoreTable(final Iterable<Template> templates, final Iterable<Haplotype> haplotypes)

    {
        this.templates = Collections.unmodifiableList(CollectionUtils.collect(templates, t -> t, new ArrayList<>(1000)));
        this.haplotypes = Collections.unmodifiableList(CollectionUtils.collect(haplotypes, t -> t, new ArrayList<>()));
        values = new double[this.haplotypes.size()][this.templates.size()];
        this.templateIndex = composeTemplateIndex(this.templates);
    }

    private Map<String,Integer> composeTemplateIndex(final List<Template> templates) {
        final Map<String, Integer> result = new HashMap<>(templates.size());
        for (int i = 0; i < templates.size(); i++) {
            result.put(templates.get(i).name(), i);
        }
        return result;
    }

    public double get(final Haplotype allele, final Template template) {
        return get(indexOf(allele), indexOf(template));
    }

    public int indexOf(final Template template) {
        Utils.nonNull(template);
        return templateIndex.getOrDefault(template.name(), -1);
    }

    public int indexOf(final Haplotype allele) {
        Utils.nonNull(allele);
        int i = 0;
        while (i < haplotypes.size()) {
            if (Arrays.equals(haplotypes.get(i).getBases(), allele.getBases())) {
                return i;
            }
            i++;
        }
        return -1;
    }

    public double get(final int alleleIndex, final int templateIndex) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        return values[alleleIndex][templateIndex];
    }

    public void set(final int alleleIndex, final int templateIndex, final double value) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        values[alleleIndex][templateIndex] = value;
    }

    public int numberOfTemplates() {
        return templates.size();
    }

    public int numberOfHaplotypes() {
        return haplotypes.size();
    }

    public List<Template> templates() {
        return templates;
    }

    public List<Haplotype> haplotypes() {
        return haplotypes;
    }

    public String toString() {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numberOfHaplotypes(); i++) {
            sb.append(Arrays.toString(values[i]));
            sb.append(", ");
        }
        if (sb.length() > 0) {
            sb.setLength(sb.length() - 1);
        }
        return sb.toString();
    }
}
