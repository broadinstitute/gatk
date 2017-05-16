package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.variant.variantcontext.VariantContext;
import org.bdgenomics.formats.avro.StructuralVariantType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.util.Collections;
import java.util.List;

/**
 * Created by valentin on 4/26/17.
 */
public class StructuralVariantContext extends VariantContext {

    private static final long serialVersionUID = 1L;

    private transient List<Integer> assemblyIDs;

    protected StructuralVariantContext(final VariantContext other) {
        super(other);
    }

    public boolean isStructural() {
        return getNAlleles() > 1 && getAlternateAlleles().stream()
                .anyMatch(StructuralVariantAllele::isStructural);
    }

    /**
     * Returns the assembly ids for this context's structural variant.
     * <p>
     *     The list returned is an immutable list.
     * </p>
     * @throws IllegalStateException if the {@link GATKSVVCFHeaderLines#ASSEMBLY_IDS} annotation contains invalid
     *
     * @return never {@code null}, an empty list if no structural variant is specified.
     */
    public List<Integer> assemblyIDs() {
        if (assemblyIDs == null) {
            if (hasAttribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS)) {
                assemblyIDs = Collections.emptyList();
            } else {
                try {
                    final List<Integer> result = getAttributeAsIntList(GATKSVVCFHeaderLines.ASSEMBLY_IDS, -1);
                    if (result.stream().anyMatch(i -> i < 0)) {
                        throw new IllegalStateException("this variant context contains unspecified or negative assembly ids");
                    }
                    assemblyIDs = Collections.unmodifiableList(result);
                } catch (final NumberFormatException ex) {
                    throw new IllegalArgumentException(String.format("the variant context '%s' annotation contains non-integer entries: %s", GATKSVVCFHeaderLines.ASSEMBLY_IDS,
                            String.valueOf(getAttribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS))));
                }

            }
        }
        return assemblyIDs;
    }
}
