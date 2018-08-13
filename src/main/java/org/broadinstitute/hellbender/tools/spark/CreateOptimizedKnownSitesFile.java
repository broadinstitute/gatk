package org.broadinstitute.hellbender.tools.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.*;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Converts a VCF of known sites into an optimized format (serialized Java object format)",
        oneLineSummary = "Converts a VCF of known sites into an optimized format (serialized Kryo object format)",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class CreateOptimizedKnownSitesFile extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger(CreateOptimizedKnownSitesFile.class);

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input VCF file.")
    private String inputFile;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "The output file, in serialized Kryo object format",
            optional = true)
    public String outputFile;

    @Override
    protected Object doWork() {
        IntervalsSkipList<GATKVariant> variants = retrieveVariants(Collections.singletonList(inputFile));
        logger.info("Retrieved variants.");
        try {
            logger.info("Writing file to " + outputFile);
            serialize(variants, Files.newOutputStream(IOUtils.getPath(outputFile)));
            logger.info("Finished writing file");
        } catch (Throwable e) {
            throw new UserException("Problem writing file", e);
        }
        return null;
    }

    private static IntervalsSkipList<GATKVariant> retrieveVariants(List<String> paths) {
        return new IntervalsSkipList<>(paths
                .stream()
                .map(CreateOptimizedKnownSitesFile::loadFromFeatureDataSource)
                .flatMap(Collection::stream)
                .collect(Collectors.toList()));
    }

    private static List<GATKVariant> loadFromFeatureDataSource(String path) {
        int cloudPrefetchBuffer = 40; // only used for GCS
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(path, null, 0, null, cloudPrefetchBuffer, cloudPrefetchBuffer) ) {
            return wrapQueryResults(dataSource.iterator());
        }
    }

    private static List<GATKVariant> wrapQueryResults(final Iterator<VariantContext> queryResults ) {
        final List<GATKVariant> wrappedResults = new ArrayList<>();
        long count = 0;
        while ( queryResults.hasNext() ) {
            if (count++ % 100000 == 0) {
                logger.info("Number of variants read: " + count);
            }
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }

    private static void serialize(Object object, OutputStream out) {
        Kryo kryo = new Kryo();
        try (Output output = new Output(out)) {
            kryo.writeClassAndObject(output, object);
        }
    }
}
