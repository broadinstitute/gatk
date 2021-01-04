package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Calculates copy number posteriors for a given set of structural variants. Supports multiple samples.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF produced by SVCopyNumberPosteriors
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Multiple VCFs containing uniform sample ploidy
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVShardVcfByPloidy
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Shards SV VCF by ploidy",
        oneLineSummary = "Shards SV VCF by ploidy",
        programGroup = CoverageAnalysisProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class SVShardVcfByPloidy extends VariantWalker {
    public static final String OUTPUT_DIR_LONG_NAME = "output-dir";
    public static final String OUTPUT_NAME_LONG_NAME = "output-name";

    @Argument(
            doc = "Output directory",
            fullName = OUTPUT_DIR_LONG_NAME
    )
    private GATKPath outputDir;

    @Argument(
            doc = "Output basename",
            fullName = OUTPUT_NAME_LONG_NAME
    )
    private String outputName;

    public static final String PLOIDY_SHARD_VCF_DESIGNATION_PREFIX = "ploidy_shard_";
    public static final String SAMPLE_SHARD_VCF_DESIGNATION_PREFIX = "sample_set_";

    private PloidySharder sharder;

    @Override
    public void onTraversalStart() {
        final VCFHeader header = getHeaderForVariants();
        if (header == null) {
            throw new UserException("VCF header not found");
        }
        sharder = new PloidySharder(header);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        sharder.add(variant);
    }

    @Override
    public Object onTraversalSuccess() {
        sharder.close();
        return null;
    }

    private final class PloidySharder {
        private final VCFHeader header;
        private int shardCount;
        private String currentContig;
        private ShardWriter currentWriter;
        private final Map<ShardIndex,ShardWriter> shardsMap;

        public PloidySharder(final VCFHeader header) {
            Utils.nonNull(header);
            this.header = header;
            this.shardCount = 0;
            this.shardsMap = new HashMap<>();
        }

        public void add(final VariantContext baseVariant) {
            Utils.nonNull(baseVariant);
            if (baseVariant.getContig() != currentContig) {
                currentContig = baseVariant.getContig();
                currentWriter = checkForAndPossiblyCreateNewShard(baseVariant);
            }
            currentWriter.add(baseVariant);
        }

        private ShardWriter checkForAndPossiblyCreateNewShard(final VariantContext baseVariant) {
            final Collection<Map.Entry<Integer, List<Genotype>>> genotypesByPloidy = baseVariant.getGenotypes().stream()
                    .collect(Collectors.groupingBy(g -> Integer.valueOf((String)g.getExtendedAttribute(GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY)), Collectors.toList()))
                    .entrySet();
            final ShardIndex index = new ShardIndex();
            for (final Map.Entry<Integer, List<Genotype>> entry : genotypesByPloidy) {
                final Integer ploidy = entry.getKey();
                final Set<String> samples = entry.getValue().stream().map(Genotype::getSampleName).collect(Collectors.toSet());
                index.add(samples, ploidy);
            }
            ShardWriter writer = shardsMap.get(index);
            if (writer == null) {
                final Path basePath = Paths.get(outputDir.toString(), outputName + "."
                       + baseVariant.getContig() + "." + (shardCount++));
                writer = new ShardWriter(basePath, header, index.getSampleSetsSet());
                shardsMap.put(index, writer);
            }
            return writer;
        }

        public void close() {
            shardsMap.values().stream().forEach(ShardWriter::close);
        }
    }

    private final class ShardWriter {

        private final Map<Set<String>,VariantContextWriter> writers;

        public ShardWriter(final Path basePath, final VCFHeader baseHeader, final Set<Set<String>> sampleSets) {
            Utils.nonNull(basePath);
            Utils.nonNull(baseHeader);
            Utils.nonNull(sampleSets);
            writers = new HashMap<>(SVUtils.hashMapCapacity(sampleSets.size()));
            int shard = 0;
            for (final Set<String> samples : sampleSets) {
                final VCFHeader header = new VCFHeader(baseHeader.getMetaDataInInputOrder(), samples);
                final Path vcfPath = Paths.get(basePath.toAbsolutePath() + "." + SAMPLE_SHARD_VCF_DESIGNATION_PREFIX + (shard++) + FileExtensions.COMPRESSED_VCF);
                final VariantContextWriter writer = createVCFWriter(vcfPath);
                writer.writeHeader(header);
                writers.put(samples, writer);
            }
        }

        public void add(final VariantContext baseVariant) {
            Utils.nonNull(baseVariant);
            for (final Map.Entry<Set<String>,VariantContextWriter> entry : writers.entrySet()) {
                final VariantContextBuilder builder = new VariantContextBuilder(baseVariant);
                final GenotypesContext fitleredGenotypes = builder.getGenotypes().subsetToSamples(entry.getKey());
                builder.genotypes(fitleredGenotypes);
                entry.getValue().add(builder.make());
            }
        }

        public void close() {
            writers.values().stream().forEach(VariantContextWriter::close);
        }
    }

    private final class ShardIndex {
        private final Map<Set<String>,Integer> samplePloidyMap;

        public ShardIndex() {
            samplePloidyMap = new HashMap<>();
        }

        public void add(final Set<String> samples, final Integer ploidy) {
            Utils.nonNull(ploidy);
            Utils.nonNull(samples);
            samplePloidyMap.put(samples, ploidy);
        }

        public Set<Set<String>> getSampleSetsSet() {
            return samplePloidyMap.keySet();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof ShardIndex)) return false;
            ShardIndex that = (ShardIndex) o;
            return Objects.equals(samplePloidyMap, that.samplePloidyMap);
        }

        @Override
        public int hashCode() {
            return Objects.hash(samplePloidyMap);
        }
    }
}
