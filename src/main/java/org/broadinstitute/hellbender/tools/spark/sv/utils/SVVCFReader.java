package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;

public class SVVCFReader {
    public static SVIntervalTree<String> readBreakpointsFromTruthVCF(final String truthVCF,
                                                                     final SAMSequenceDictionary dictionary,
                                                                     final int padding ) {
        SVIntervalTree<String> breakpoints = new SVIntervalTree<>();
        try ( final FeatureDataSource<VariantContext> dataSource =
                      new FeatureDataSource<>(truthVCF, null, 0, VariantContext.class) ) {
            for ( final VariantContext vc : dataSource ) {
                final StructuralVariantType svType = vc.getStructuralVariantType();
                if ( svType == null ) continue;
                final String eventName = vc.getID();
                final int contigID = dictionary.getSequenceIndex(vc.getContig());
                if ( contigID < 0 ) {
                    throw new UserException("VCF contig " + vc.getContig() + " does not appear in dictionary.");
                }
                final int start = vc.getStart();
                switch ( svType ) {
                    case DEL:
                    case INV:
                    case CNV:
                        final int end = vc.getEnd();
                        breakpoints.put(new SVInterval(contigID,start-padding, end+padding), eventName);
                        break;
                    case INS:
                    case DUP:
                    case BND:
                        breakpoints.put(new SVInterval(contigID,start-padding, start+padding), eventName);
                        break;
                }
            }
        }
        return breakpoints;
    }
}
