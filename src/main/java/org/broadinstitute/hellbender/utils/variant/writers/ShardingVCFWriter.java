package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Variant writer tha splits output to multiple VCFs given the maximum records per file. Before using {@link #add},
 * the header should be set using either {@link #setHeader} or {@link #writeHeader}, which may only be called once
 * and will determine whether headers are written to all shards.
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
public class ShardingVCFWriter implements VariantContextWriter {

    public static final String SHARD_INDEX_PREFIX = ".shard_";
    public static final String SHARD_INDEX_SUFFIX = FileExtensions.COMPRESSED_VCF;

    private VariantContextWriter writer;
    private VCFHeader header;
    private final int maxVariantsPerShard;
    private final Path basePath;
    private final SAMSequenceDictionary dictionary;
    private final boolean createMD5;
    private final Options[] options;

    /** Current shard  */
    private int shardIndex;
    /** Number of records written to current shard  */
    private int shardSize;
     /** Whether to write header, or null if header is undefined */
    private Boolean enableWriteHeader;

    /**
     * Create a new sharding VCF writer
     *
     * @param basePath              base path of the output VCFs. The shard designation and file extension will be added.
     * @param maxVariantsPerShard   max number of records per file (last shard may have less)
     * @param dictionary            sequence dictionary for this writer
     * @param createMD5             enable MD5 file creation
     * @param options               vcf writer options
     */
    public ShardingVCFWriter(final Path basePath,
                             final int maxVariantsPerShard,
                             final SAMSequenceDictionary dictionary,
                             final boolean createMD5,
                             final Options... options) {
        Utils.nonNull(basePath);
        Utils.validateArg(maxVariantsPerShard > 0, "maxVariantsPerShard must be positive");
        this.basePath = IOUtils.removeExtension(basePath, FileExtensions.VCF_LIST);
        this.maxVariantsPerShard = maxVariantsPerShard;
        this.dictionary = dictionary;
        this.createMD5 = createMD5;
        this.options = options;

        // Initialize first shard
        this.shardIndex = 0;
        this.shardSize = 0;
        this.writer = createNewWriter();
    }

    /**
     * Initializes a new sharded file.
     */
    protected void createNextShard() {
        writer.close();
        shardIndex++;
        shardSize = 0;
        writer = createNewWriter();
        Utils.nonNull(header, "Attempted to create new shard before header has been set");
        initializeShardHeader();
    }

    /**
     * Initializes shard header depending on which header function (set or write) was used
     */
    protected void initializeShardHeader() {
        if (enableWriteHeader != null) {
            if (enableWriteHeader.booleanValue()) {
                writer.writeHeader(header);
            } else {
                writer.setHeader(header);
            }
        }
    }

    /**
     * Creates a writer for a new shard
     *
     * @return the new writer
     */
    protected VariantContextWriter createNewWriter() {
        final Path outPath = Paths.get(getShardFilename(basePath, shardIndex));
        return GATKVariantContextUtils.createVCFWriter(
                outPath,
                dictionary,
                createMD5,
                options);
    }

    /**
     * Gets filepath for the given shard and base path
     *
     * @param basePath path without extension
     * @param shardIndex
     * @return path as String
     */
    public static String getShardFilename(final Path basePath, final int shardIndex) {
        return String.format("%s%s%05d%s", basePath, SHARD_INDEX_PREFIX, shardIndex, SHARD_INDEX_SUFFIX);
    }

    /**
     * Defines header for the writer. The header will not be written to any shards. May only be called once and only
     * if {@link #writeHeader} has not been called.
     *
     * @param header header to use
     */
    @Override
    public void setHeader(final VCFHeader header) {
        Utils.validate(this.header == null, "Cannot redefine header");
        this.header = header;
        enableWriteHeader = Boolean.FALSE;
        writer.setHeader(header);
    }

    /**
     * Defines header for the writer that will be written to all shards. May only be called once and only
     * if {@link #setHeader} has not been called.
     *
     * @param header header to use
     */
    @Override
    public void writeHeader(final VCFHeader header) {
        Utils.validate(this.header == null, "Cannot redefine header");
        this.header = header;
        enableWriteHeader = Boolean.TRUE;
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

    /**
     * Adds variant to writer. Note that a header must be assigned first.
     *
     * @param vc variant to write
     */
    @Override
    public void add(final VariantContext vc) {
        if (shardSize + 1 > maxVariantsPerShard) {
            createNextShard();
        }
        writer.add(vc);
        shardSize++;
    }
}
