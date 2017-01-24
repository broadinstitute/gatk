package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A utility class that writes out variants to a VCF file.
 * A work that should be improved as GATK SV pipeline matures.
 */
class SVVCFWriter {

    /**
     * FASTA and Broadcast references are both required because 2bit Broadcast references currently order their
     * sequence dictionaries in a scrambled order, see https://github.com/broadinstitute/gatk/issues/2037.
     */
    static void writeVCF(final PipelineOptions pipelineOptions, final String outputPath, String vcfFileName,
                         final String fastaReference, final JavaRDD<VariantContext> variantContexts) {

        final SAMSequenceDictionary referenceSequenceDictionary = new ReferenceMultiSource(pipelineOptions, fastaReference, ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);

        final List<VariantContext> sortedVariantsList = sortVariantsByCoordinate(variantContexts.collect(), referenceSequenceDictionary);

        writeVariants(pipelineOptions, outputPath, vcfFileName, sortedVariantsList, referenceSequenceDictionary);
    }

    // TODO: right now there's an edge case that the "same" inversion events would be called three times on a test sample such that they have the same start, end and inversion evidence type
    //     but differ only in their inserted sequence, sorting these variants must take into account of such complications. the solution below is hackish
    private static List<VariantContext> sortVariantsByCoordinate(final List<VariantContext> variants, SAMSequenceDictionary referenceSequenceDictionary) {
        return variants.stream().sorted((VariantContext v1, VariantContext v2) -> {
            final int x = IntervalUtils.compareLocatables(v1, v2, referenceSequenceDictionary);
            if (x==0) {
                final String s1 = v1.getAttributeAsString(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, "");
                final String s2 = v2.getAttributeAsString(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, "");
                return s1.compareTo(s2);
            } else {
                return x;
            }
        }).collect(Collectors.toList());
    }

    private static void writeVariants(final PipelineOptions pipelineOptions, final String outputPath, final String fileName,
                                      final List<VariantContext> variantsArrayList, final SAMSequenceDictionary referenceSequenceDictionary) {
        try (final OutputStream outputStream = new BufferedOutputStream(
                BucketUtils.createFile(outputPath + "/" + fileName, pipelineOptions))) {

            final VariantContextWriter vcfWriter = getVariantContextWriter(outputStream, referenceSequenceDictionary);

            vcfWriter.writeHeader(getVcfHeader(referenceSequenceDictionary));
            variantsArrayList.forEach(vcfWriter::add);
            vcfWriter.close();

        } catch (final IOException e) {
            throw new GATKException("Could not create output file", e);
        }
    }

    private static VCFHeader getVcfHeader(final SAMSequenceDictionary referenceSequenceDictionary) {
        final VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(referenceSequenceDictionary);
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        GATKSVVCFHeaderLines.vcfHeaderLines.values().forEach(header::addMetaDataLine);
        return header;
    }

    private static VariantContextWriter getVariantContextWriter(final OutputStream outputStream, final SAMSequenceDictionary referenceSequenceDictionary) {
        VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder()
                .clearOptions()
                .setOutputStream(outputStream);

        if (null != referenceSequenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceSequenceDictionary);
        }
        // todo: remove this when things are solid?
        vcWriterBuilder = vcWriterBuilder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        for (Options opt : new Options[]{}) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }
}
