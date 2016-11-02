package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class ContigAligner implements Closeable {

    static String referencePath;

    final BwaIndex index;
    final BwaMem bwaMem;

    private static final Logger log = LogManager.getLogger(ContigAligner.class);

    static {
        BWANativeLibrary.load();
    }

    public ContigAligner(final String referenceFilePath) throws IOException {
        referencePath = referenceFilePath;

        try {
            index = new BwaIndex(LocalizedReference.INSTANCE);
            log.info("Created BWA index");
        } catch (final IOException e) {
            throw new GATKException("Unable to load reference", e);
        }
        bwaMem = new BwaMem(index);
        log.info("Created BWA MEM");
    }

    /**
     * Takes a collection of assembled contigs and aligns them to the reference with jBWA. Non-canonical
     * (secondary) alignments are filtered out, preserving the primary and supplementary alignments.
     * Within the output list, alignments are sorted first by contig (based on the order in which
     * the contigs were passed in, and then by their start position on the contig).
     *
     * @param assemblyId An identifier for the assembly or set of contigs
     * @param contigsCollection The set of all canonical (primary or supplementary) alignments for the contigs.
     * @return
     */
    @VisibleForTesting
    public List<AlignmentRegion> alignContigs(final String assemblyId, final ContigsCollection contigsCollection) {
        final List<AlignmentRegion> alignedContigs = new ArrayList<>(contigsCollection.getContents().size());
        try {
            for(final LocalAssemblyContig contigInfo : contigsCollection.getContents()) {
                final String contigId = contigInfo.contigID;
                final byte[] sequence = contigInfo.seq.getBytes();
                final AlnRgn[] alnRgns = bwaAlignSequence(bwaMem, contigId, sequence);

                // filter out secondary alignments, convert to AlignmentRegion objects and sort by alignment start pos
                Arrays.stream(alnRgns)
                        .filter(a -> a.getSecondary() < 0)
                        .map(a -> new AlignmentRegion(assemblyId, contigId, a))
                        .sorted(Comparator.comparing(a -> a.startInAssembledContig))
                        .forEach(alignedContigs::add);
            }
        } catch (final IOException e) {
            throw new GATKException("could not execute BWA");
        }

        return alignedContigs;
    }

    /**
     * Wrap a contig sequence in a ShortRead object and pass it to BWA to align
     */
    private AlnRgn[] bwaAlignSequence(final BwaMem bwaMem, final String contigId, final byte[] sequence) throws IOException {
        final ShortRead contigShortRead = new ShortRead(contigId, sequence, qualSequence(sequence.length));
        log.debug("Calling bwaMem.align");
        return bwaMem.align(contigShortRead);
    }

    /**
     * Generate a bogus base quality sequence to pass in for the aligned contig (since the jBWA API requires that reads have qualities)
     */
    private byte[] qualSequence(final int length) {
        final byte[] quals = new byte[length];
        Arrays.fill(quals, (byte)'A');
        return quals;
    }

    @Override
    public void close() {
        log.info("closing BWA mem and index");
        bwaMem.dispose();
        index.close();
    }

    private static class LocalizedReference {
        static File INSTANCE;

        static {
            try {
                INSTANCE = BucketUtils.isHadoopUrl(ContigAligner.referencePath) ? BwaSparkEngine.localizeReferenceAndBwaIndexFiles(ContigAligner.referencePath) : new File(ContigAligner.referencePath);
            } catch (final IOException e) {
                throw new GATKException("unable to localize reference", e);
            }
        }
    }

}
