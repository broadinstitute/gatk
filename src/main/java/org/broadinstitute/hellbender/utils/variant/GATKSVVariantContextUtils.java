package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class GATKSVVariantContextUtils {

    public static final Allele BND_ALLELE = Allele.create("<BND>");
    public static final Allele CPX_ALLELE = Allele.create("<CPX>");
    public static final Allele CTX_ALLELE = Allele.create("<CTX>");

    /**
     * Build the list of called alleles based on reference and called copy numbers
     * For CNVs only, i.e. will only assign reference, DEL, and DUP alleles
     * @param copyNumberCall
     * @param refCopyNumber
     * @param refAllele reference allele representation for the position of interest
     * @return a list of alleles appropriate to pass to a GenotypeBuilder
     */
    public static List<Allele> makeGenotypeAllelesFromCopyNumber(final int copyNumberCall, final int refCopyNumber, final Allele refAllele) {
        final List<Allele> returnAlleles = new ArrayList<>();
        final Allele genotypeAllele = getAlleleForCopyNumber(copyNumberCall, refCopyNumber, refAllele);
        //some allosomes like Y can have ref copy number zero, in which case we just no-call
        if (refCopyNumber == 0) {
            return GATKVariantContextUtils.noCallAlleles(1);
        }
        //for only one haplotype we know which allele it has
        if (refCopyNumber == 1) {
           return Collections.singletonList(genotypeAllele);
        //can't determine counts per haplotypes if there is a duplication
        } if (genotypeAllele.equals(GATKSVVCFConstants.DUP_ALLELE)) {
            return GATKVariantContextUtils.noCallAlleles(refCopyNumber);
        //for homDels, hetDels or homRefs
        } if (refCopyNumber == 2) {
            if (copyNumberCall == 0) {
                returnAlleles.add(genotypeAllele);
            } else {
                returnAlleles.add(refAllele);
            }
            returnAlleles.add(genotypeAllele);
            return returnAlleles;
        }
        //multiploid dels
        for (int i = 0; i < copyNumberCall; i++) {
            returnAlleles.add(refAllele);
        }
        for (int i = copyNumberCall; i < refCopyNumber; i++) {
            returnAlleles.add(GATKSVVCFConstants.DEL_ALLELE);
        }
        return returnAlleles;

    }

    /**
     *
     * @param copyNumberCall
     * @param refCopyNumber
     * @param refAllele
     * @return variant allele if copyNumberCall != refCopyNumber, else refAllele
     */
    public static Allele getAlleleForCopyNumber(final int copyNumberCall, final int refCopyNumber, final Allele refAllele) {
        if (copyNumberCall > refCopyNumber) {
            return GATKSVVCFConstants.DUP_ALLELE;
        } else if (copyNumberCall < refCopyNumber) {
            return GATKSVVCFConstants.DEL_ALLELE;
        } else {
            return refAllele;
        }
    }

    /**
     * Parses allele for the base type string, e.g. "INS" for "&lt;INS:MEI&gt;"
     */
    public static final String[] getSymbolicAlleleSymbols(final Allele allele) {
        return allele.getDisplayString()
                .replace("<", "")
                .replace(">", "")
                .split(":");
    }

    /**
     * Determines whether a given SV type represents a CNV.
     */
    public static boolean isCnvType(final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        return type == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL
                || type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP
                || type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
    }

    /**
     * Determines whether a given SV type represents a CNV.
     */
    public static boolean isCnvType(final StructuralVariantType type) {
        return type == StructuralVariantType.DEL
                || type == StructuralVariantType.DUP
                || type == StructuralVariantType.CNV;
    }
}
