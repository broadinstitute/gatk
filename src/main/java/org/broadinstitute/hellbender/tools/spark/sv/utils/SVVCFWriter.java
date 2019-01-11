package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
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
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A utility class that writes out variants to a VCF file.
 * A work that should be improved as GATK SV pipeline matures.
 */
public class SVVCFWriter {

    /**
     * {@code referenceSequenceDictionary} is required because 2bit Broadcast references currently order their
     * sequence dictionaries in a scrambled order, see https://github.com/broadinstitute/gatk/issues/2037.
     */
    public static void writeVCF(final List<VariantContext> localVariants, final String vcfFileName,
                                final SAMSequenceDictionary referenceSequenceDictionary,
                                final Set<VCFHeaderLine> defaultToolVCFHeaderLines,
                                final Logger logger) {

        final List<VariantContext> sortedVariantsList = sortVariantsByCoordinate(localVariants, referenceSequenceDictionary);

        if (logger != null)
            logNumOfVarByTypes(sortedVariantsList, logger);

        writeVariants(vcfFileName, sortedVariantsList, referenceSequenceDictionary, defaultToolVCFHeaderLines);
    }

    private static void logNumOfVarByTypes(final List<VariantContext> variants, final Logger logger) {

        logger.info("Discovered a total of " + variants.size() + " records: ");

        final Map<String, Long> variantsCountByType = variants.stream()
                .collect(Collectors.groupingBy(vc -> (String) vc.getAttribute(GATKSVVCFConstants.SVTYPE), Collectors.counting()));

        variantsCountByType.forEach((key, value) -> logger.info("  " + key + ": " + value));

        logger.info("  And none of : " + Sets.difference(SvType.getKnownTypes(), variantsCountByType.keySet()).toString());
    }


    // TODO: 5/31/18 it has been the case since we output VCF files that some records have exactly the same POS and END, yet with slight differences in annotations (e.g. inserted sequence, homology, etc.) pointing to difference variants;
    //       this sort is to make sure such records are sorted. Ultimately we should decide on what to do (squash them into single records when possible?) with such records.
    @VisibleForTesting
    public static List<VariantContext> sortVariantsByCoordinate(final List<VariantContext> variants,
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
                                      final SAMSequenceDictionary referenceSequenceDictionary,
                                      final Set<VCFHeaderLine> defaultToolVCFHeaderLines) {
        try (final OutputStream outputStream
                     = new BufferedOutputStream(BucketUtils.createFile(fileName))) {

            final VariantContextWriter vcfWriter = getVariantContextWriter(outputStream, referenceSequenceDictionary);

            final VCFHeader vcfHeader = getVcfHeader(referenceSequenceDictionary);
            defaultToolVCFHeaderLines.forEach(vcfHeader::addMetaDataLine);
            vcfWriter.writeHeader(vcfHeader);
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
        headerLines.addAll(GATKSVVCFHeaderLines.getFormatLines());
        headerLines.addAll(GATKSVVCFHeaderLines.getFilterLines());
        final VCFHeader header = new VCFHeader(new VCFHeader( headerLines ));
        header.setSequenceDictionary(referenceSequenceDictionary);
        return header;
    }

    // TODO: 8/10/18 see 5083
    private static VariantContextWriter getVariantContextWriter(final OutputStream outputStream,
                                                                final SAMSequenceDictionary referenceSequenceDictionary) {
        VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder()
                .clearOptions()
                .setOutputStream(outputStream);

        if (null != referenceSequenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceSequenceDictionary);
        }
        for (final Options opt : new Options[]{}) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }
}
