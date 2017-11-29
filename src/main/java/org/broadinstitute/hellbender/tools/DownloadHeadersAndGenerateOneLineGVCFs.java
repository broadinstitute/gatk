package org.broadinstitute.hellbender.tools;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.stream.Collector;
import java.util.stream.Collectors;

@BetaFeature
@CommandLineProgramProperties(summary = "Fix the sample names in a vcf that ran through the GenomicsDBImport when " +
        "the batch ordering bug was present, this will restore the correct sample names provided it is given the exact sample name mapping / vcf ordering" +
        "and batch sized that was used in the initial import. " +
        "See https://github.com/broadinstitute/gatk/issues/3682 for more information",
        oneLineSummary = "fix sample names in a shuffled callset",
        programGroup = VariantProgramGroup.class,
        omitFromCommandLine = true)
public final class DownloadHeadersAndGenerateOneLineGVCFs extends CommandLineProgram {
    public static final String SAMPLE_NAME_KEY = "SN";
    public static final String ORIGINAL_FILE_NAME_KEY = "OF";

    private static final String CHR_1 = "chr1";

    @Argument(fullName = GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME,
            doc = "the same sampleNameMap file which was used to import the callset using GenomicsDBImport",
            optional = false)
    public String sampleNameMapPath;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "where to write a reheadered version of the input VCF with the sample names in the correct order",
            optional = false)
    public String output;


    @Argument(fullName = "threads", doc = "how many threads to use")
    public int threads = 4;

    private static final Allele AREF = Allele.create("A", true);
    private static final Allele C = Allele.create("C");

    @Override
    protected Object doWork() {

        final LinkedHashMap<String, Path> sampleNameMap = GenomicsDBImport.loadSampleNameMapFile(IOUtils.getPath(sampleNameMapPath));
        final ThreadPoolExecutor executorService = initializeExecutorService(threads);

        try {
            Files.createDirectory(IOUtils.getPath(output));
        } catch (IOException e) {
            throw new UserException("Couldn't create output directory", e);
        }

        final List<Future<Path>> futures = sampleNameMap.entrySet().stream()
                .map(e -> executorService.submit(() -> getHeaderAndWriteNewGVCF(e.getValue(), e.getKey(), IOUtils.getPath(output))))
                .collect(Collectors.toList());
        
        executorService.shutdown();
        try {
            while (!executorService.awaitTermination(1, TimeUnit.SECONDS)){
                System.out.println("task count: "+ executorService.getTaskCount() + " running: " + executorService.getActiveCount() + " completed: " + executorService.getCompletedTaskCount());
            }
        } catch (InterruptedException e){
            throw  new UserException("interrupted", e);
        }

        futures.forEach(f -> {
            try {
                f.get();
            } catch (InterruptedException | ExecutionException e) {
                throw new UserException("problem with futures", e);
            }
        });



//        final List<Path> outputPaths = sampleNameMap.entrySet()
//                .stream()
//                .unordered()
//                .parallel()
//                .map(e -> getHeaderAndWriteNewGVCF(e.getValue(), e.getKey(), IOUtils.getPath(output)))
//                .peek(path -> {
//                    final int count = counter.incrementAndGet();
//                    if (count % 10 == 0) {
//                        System.out.println("Done " + count);
//                    }
//                })
//                .collect(Collectors.toList());
        return 0;
    }

    private static Path getHeaderAndWriteNewGVCF(Path gvcf, String renamedSampleName, Path outputDir){
        try ( final AbstractFeatureReader<VariantContext, ?> reader = AbstractFeatureReader.getFeatureReader(gvcf.toAbsolutePath().toUri().toString(), null, new VCFCodec(), false, Function.identity(), Function
                .identity()))
        {

            final VCFFormatHeaderLine sampleNameFormatField = new VCFFormatHeaderLine(SAMPLE_NAME_KEY, 1, VCFHeaderLineType.String,
                                                                            "the name of the sample this genotype came from");
            final VCFFormatHeaderLine originalFileNameFormatField = new VCFFormatHeaderLine(ORIGINAL_FILE_NAME_KEY, 1, VCFHeaderLineType.String,
                                                                                      "the full path to the file this came from");
            final VCFHeader header = (VCFHeader) reader.getHeader();
            final List<String> samples = header.getSampleNamesInOrder();
            if (samples.size() != 1){
                throw new UserException("Expected exactly 1 sample per file");
            }

            final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(header.getMetaDataInInputOrder());
            headerLines.add(sampleNameFormatField);
            headerLines.add(originalFileNameFormatField);
            final VCFHeader newHeader = new VCFHeader(headerLines, samples);

            final Path outputVcf = outputDir.resolve(renamedSampleName + "g.vcf.gz");

            try(VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                    outputVcf.toFile(), header.getSequenceDictionary(), false,
                    Options.INDEX_ON_THE_FLY)){
                writer.writeHeader(newHeader);

                final List<Allele> alleles = Arrays.asList(AREF, C);
                final VariantContext variant = new VariantContextBuilder("invented", CHR_1, 100, 100, alleles)
                        .genotypes(new GenotypeBuilder(samples.get(0), alleles)
                                           .attribute(SAMPLE_NAME_KEY, Utils.nonNull(renamedSampleName))
                                           .attribute(ORIGINAL_FILE_NAME_KEY, Utils.nonNull(gvcf.toString()))
                                           .make())
                        .make();
                writer.add(variant);
            }

            return outputVcf;
        }
        catch ( IOException e ) {
            throw new UserException("Error opening header for gvcf " + gvcf);
        }
    }



    private static ThreadPoolExecutor initializeExecutorService(final int threads) {
        final ThreadFactory threadFactory = new ThreadFactoryBuilder()
                .setNameFormat("readerInitializer-thread-%d")
                .setDaemon(true)
                .build();
        return (ThreadPoolExecutor)Executors.newFixedThreadPool(threads, threadFactory);
    }

}

