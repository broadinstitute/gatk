package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A utility class that writes out variants to a VCF file.
 * A work that should be improved as GATK SV pipeline matures.
 */
public class SVVCFWriter {

    public static void writeVCF(final List<VariantContext> localVariants, final String vcfFileName,
                                final SAMSequenceDictionary referenceSequenceDictionary, final Logger logger) {

        final List<VariantContext> sortedVariantsList = sortVariantsByCoordinate(localVariants, referenceSequenceDictionary);

        logNumOfVarByTypes(sortedVariantsList, logger);

        writeVariants(vcfFileName, sortedVariantsList, referenceSequenceDictionary);
    }

    private static void logNumOfVarByTypes(final List<VariantContext> sortedVariantsList, final Logger logger) {

        logger.info("Discovered " + sortedVariantsList.size() + " variants.");

        sortedVariantsList.stream()
                .collect(Collectors.groupingBy(vc -> (String)vc.getAttribute(GATKSVVCFConstants.SVTYPE), Collectors.counting()))
                .forEach((key, value) -> logger.info(key + ": " + value));
    }

    // TODO: right now there's an edge case that the "same" inversion events would be called three times on a test sample
    //       such that they have the same start, end and inversion evidence type but differ only in their inserted sequence,
    //       sorting these variants must take into account of such complications.
    // the solution below is hackish
    @VisibleForTesting
    static List<VariantContext> sortVariantsByCoordinate(final List<VariantContext> variants,
                                                         final SAMSequenceDictionary referenceSequenceDictionary) {
        return variants.stream().sorted((VariantContext v1, VariantContext v2) -> {
            final int x = IntervalUtils.compareLocatables(v1, v2, referenceSequenceDictionary);
            if (x == 0) {
                final String s1 = v1.getAttributeAsString(GATKSVVCFConstants.INSERTED_SEQUENCE, "");
                final String s2 = v2.getAttributeAsString(GATKSVVCFConstants.INSERTED_SEQUENCE, "");
                return s1.compareTo(s2);
            } else {
                return x;
            }
        }).collect(SVUtils.arrayListCollector(variants.size()));
    }

    private static void writeVariants(final String fileName, final List<VariantContext> variantsArrayList,
                                      final SAMSequenceDictionary referenceSequenceDictionary) {
        try (final OutputStream outputStream
                     = new BufferedOutputStream(BucketUtils.createFile(fileName))) {

            final VariantContextWriter vcfWriter = getVariantContextWriter(outputStream, referenceSequenceDictionary);

            vcfWriter.writeHeader(getVcfHeader(referenceSequenceDictionary));
            variantsArrayList.forEach(vcfWriter::add);
            vcfWriter.close();

        } catch (final IOException e) {
            throw new GATKException("Could not create output file", e);
        }
    }

    @VisibleForTesting
    static VCFHeader getVcfHeader(final SAMSequenceDictionary referenceSequenceDictionary) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>(GATKSVVCFHeaderLines.getSymbAltAlleleLines());
        headerLines.addAll(GATKSVVCFHeaderLines.getInfoLines());
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        final VCFHeader header = new VCFHeader(new VCFHeader( headerLines ));
        header.setSequenceDictionary(referenceSequenceDictionary);
        return header;
    }

    private static VariantContextWriter getVariantContextWriter(final OutputStream outputStream,
                                                                final SAMSequenceDictionary referenceSequenceDictionary) {
        VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder()
                .clearOptions()
                .setOutputStream(outputStream);

        if (null != referenceSequenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceSequenceDictionary);
        }
        // todo: remove this when things are solid?
        vcWriterBuilder = vcWriterBuilder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        for (final Options opt : new Options[]{}) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }
}
