package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.file.Path;
import java.nio.file.Paths;

public class ShardingVCFWriter implements VariantContextWriter {

    public static final String SHARD_INDEX_PREFIX = ".shard_";
    public static final String SHARD_INDEX_SUFFIX = FileExtensions.COMPRESSED_VCF;

    private VariantContextWriter writer;
    private VCFHeader header;
    int shardIndex;
    int shardSize;
    private final int maxVariantsPerShard;
    private final Path basePath;
    private final SAMSequenceDictionary dictionary;
    private final boolean createMD5;
    private final Options[] options;

    public ShardingVCFWriter(final Path basePath,
                             final int maxVariantsPerShard,
                             final SAMSequenceDictionary dictionary,
                             final boolean createMD5,
                             final Options... options) {
        Utils.validateArg(maxVariantsPerShard > 0, "maxVariantsPerShard must be positive");
        this.basePath = basePath;
        this.maxVariantsPerShard = maxVariantsPerShard;
        this.dictionary = dictionary;
        this.createMD5 = createMD5;
        this.options = options;

        // Initialize first shard
        this.shardIndex = 0;
        this.shardSize = 0;
        this.writer = createNewWriter();
    }

    private void createNextShard() {
        writer.close();
        shardIndex++;
        shardSize = 0;
        writer = createNewWriter();
        if (header == null) {
            throw new GATKException("Attempted to create new shard before header has been set");
        }
        writer.writeHeader(header);
    }

    private VariantContextWriter createNewWriter() {
        final Path outPath = Paths.get(basePath + SHARD_INDEX_PREFIX + shardIndex + SHARD_INDEX_SUFFIX);
        return GATKVariantContextUtils.createVCFWriter(
                outPath,
                dictionary,
                createMD5,
                options);
    }

    @Override
    public void writeHeader(final VCFHeader header) {
        this.header = header;
        writer.writeHeader(header);
    }

    @Override
    public void close() {
        writer.close();
    }

    @Override
    public boolean checkError() {
        return writer.checkError();
    }

    @Override
    public void add(final VariantContext vc) {
        if (shardSize + 1 > maxVariantsPerShard) {
            createNextShard();
        }
        writer.add(vc);
        shardSize++;
    }

    @Override
    public void setHeader(final VCFHeader header) {
        this.header = header;
        writer.setHeader(header);
    }
}
