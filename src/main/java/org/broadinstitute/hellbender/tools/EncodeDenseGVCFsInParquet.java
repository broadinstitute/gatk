package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.generic.GenericRecordBuilder;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.avro.AvroParquetReader;
import org.apache.parquet.avro.AvroParquetWriter;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.hadoop.ParquetWriter;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.io.*;
import java.util.List;

@CommandLineProgramProperties(
        summary = "FOO",
        oneLineSummary = "BAR",
        programGroup = VariantProgramGroup.class
)
public final class EncodeDenseGVCFsInParquet extends CommandLineProgram{

    @Argument(shortName = "F", fullName = "file", doc = "File or dir")
    public File arg;

    @Argument(shortName = "out", fullName = "O", doc = "file out")
    public File out;

    @Override
    protected Object doWork() {
        final File[] files;
        if (arg.isDirectory()) {
            files = arg.listFiles(fn -> fn.getName().endsWith(".vcf"));
        } else {
            files = new File[]{arg};
        }
        System.out.println("Processing " + files.length + " files");
        final Schema.Parser parser = new Schema.Parser();
        try (final InputStream stream = new FileInputStream(new File("src/main/resources/Alleles.avsc"))) {
            final Schema schema = parser.parse(stream);
            final Path path = new Path(out.getCanonicalPath());

            try (ParquetWriter<GenericRecord> writer = AvroParquetWriter.<GenericRecord>builder(path).withSchema(schema).build()) {
                final GenericRecordBuilder builder = new GenericRecordBuilder(schema);
                for (int i = 0; i < files.length; i++) {
                    System.out.println("File " + i + " of " + files.length + " " + files[i]);
                    encodeGVCF(files[i], writer, builder);
                }
            }

            final AvroParquetReader.Builder<GenericRecord> builder = AvroParquetReader.<GenericRecord>builder(path);
            //I'm not sure why the suppress is needed
            try (@SuppressWarnings("unchecked") final ParquetReader<GenericRecord> reader = builder.build()) {
                GenericRecord genericRecord = reader.read();
                System.out.println(genericRecord);
            }
        } catch (IOException e) {
            throw new GATKException("ERROR", e);
        }
        return null;
    }

    private static void encodeGVCF(final File file, final ParquetWriter<GenericRecord> writer, final GenericRecordBuilder builder) throws IOException {
        try(final VCFFileReader vcf = new VCFFileReader(file)){
        final List<String> namesInOrder = vcf.getFileHeader().getSampleNamesInOrder();
        if (namesInOrder.size() != 1){
            System.err.println("expected 1 sample but got " + namesInOrder.size());
            return;
        }
        final String sampleName= namesInOrder.get(0);
        for (final VariantContext vc : vcf) {
            final String contig = vc.getContig();
            final Allele ref = vc.getReference();
            final int pos = vc.getStart();
            final Genotype genotype = vc.getGenotype(0);
            final int[] ad = genotype.getAD();
            final int refCount;
            final int altCount;
            if (ad != null){
               refCount = ad[0];
               altCount = (int) MathUtils.sum(ad) - refCount;
	    } else {
	       refCount = 0;
               altCount = 0;
	    }
            builder.set("sampleName", sampleName);
            builder.set("contig", contig);
            builder.set("REF", ref.getBaseString());
            builder.set("pos", pos);
            builder.set("NREF", refCount); //UGH, I want to write 1 byte. Not sure how to in avro.
            builder.set("NALT", altCount);
            writer.write(builder.build());
        }
    }
    }
}
