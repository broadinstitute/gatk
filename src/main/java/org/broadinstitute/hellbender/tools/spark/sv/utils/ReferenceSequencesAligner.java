package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public final class ReferenceSequencesAligner implements AutoCloseable {

    private final File fasta;
    private final File image;
    private final BwaMemAligner aligner;
    private final BwaMemIndex index;
    private final List<String> refNames;

    public static final class DescribedRefContig {
        private final String name;
        private final String description;
        private final byte[] bases;

        public DescribedRefContig(final String name, final String description, final byte[] bases) {
            this.name = name;
            this.description = description;
            this.bases = bases;
        }

        public String getName() {
            return name;
        }

        public String getDescription() {
            return description;
        }

        public byte[] getBases() {
            return bases;
        }
    }

    public ReferenceSequencesAligner(final String name, final byte[] bases) {
        try {
            Utils.nonNull(bases);
            refNames = Collections.singletonList( Utils.nonNull(name) );
            fasta = File.createTempFile("ssvh-temp", ".fasta");
            fasta.deleteOnExit();
            image = File.createTempFile(fasta.getParent(), fasta.toString().replace(".fasta", ".img"));
            image.deleteOnExit();
            FastaReferenceWriter.writeSingleSequenceReference(fasta.toPath(), false, false, name, null, bases);
            BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), image.toString());
            index = new BwaMemIndex(image.toString());
            aligner = new BwaMemAligner(index);
        } catch (final IOException ex) {
            throw new GATKException("could not create index files", ex);
        }
    }

    /**
     * Note: when {@code pathToSaveFasta} is {@code null},
     * the fasta is saved in a temporary directory that will be deleted when the alignment step is finished.
     */
    public ReferenceSequencesAligner(final List<String> refNames, final List<DescribedRefContig> refContigs,
                                     final boolean makeIndex, final boolean makeDict,
                                     final int basesPerLine, final String pathToSaveFasta) {

        this.refNames = Utils.nonNull(refNames);

        try {
            if (pathToSaveFasta == null) {
                fasta = File.createTempFile("ssvh-temp", ".fasta");
                fasta.deleteOnExit();
            } else {
                fasta = new File(pathToSaveFasta);
            }

            image = new File(fasta.toString().replace(".fasta", ".img"));
            image.deleteOnExit();
        } catch (final IOException ioex) {
            throw new GATKException("could not create index files", ioex);
        }

        try (final FastaReferenceWriter fastaReferenceWriter =
                     makeIndex ? new FastaReferenceWriter(fasta.toPath(), basesPerLine, makeIndex, makeDict)
                               : new FastaReferenceWriter(fasta.toPath(), makeIndex, makeDict)
        ) {
            for (final DescribedRefContig contig : refContigs) {
                fastaReferenceWriter.appendSequence(contig.name, contig.description, contig.bases);
            }
            BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), image.toString());
            index = new BwaMemIndex(image.toString());
            aligner = new BwaMemAligner(index);
        } catch (final IOException ioex) {
            throw new GATKException("could not create index files", ioex);
        }
    }

    public final List<List<AlignmentInterval>> align(final List<byte[]> seqs) {
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs);
        final List<List<AlignmentInterval>> result = new ArrayList<>(alignments.size());
        for (int i = 0; i < alignments.size(); i++) {
            final int queryLength = seqs.get(i).length;
            final List<AlignmentInterval> intervals = alignments.get(i).stream()
                    .filter(bwa -> bwa.getRefId() >= 0)
                    .filter(bwa -> SAMFlag.SECONDARY_ALIGNMENT.isUnset(bwa.getSamFlag()))
                    .map(bma -> new AlignmentInterval(bma, refNames, queryLength))
                    .collect(Collectors.toList()); // ignore secondary alignments.
            result.add(intervals);
        }
        return result;
    }

    public final List<List<SAMRecord>> align(final List<SVFastqUtils.FastqRead> reads,
                                             final SAMFileHeader header,
                                             final boolean pairedReadsInterleaved) {

        final int readCount = reads.size();

        final List<String> readNames = new ArrayList<>(readCount);
        final List<byte[]> bases = new ArrayList<>(readCount);
        reads.forEach(fastqRead -> {
            readNames.add(fastqRead.getName());
            bases.add(fastqRead.getBases());
        });

        if (pairedReadsInterleaved) aligner.alignPairs();
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(bases);

        final List<List<SAMRecord>> result = new ArrayList<>(readCount);
        for (int i = 0; i < readNames.size(); ++i) {
            final String readName = readNames.get(i);
            final List<BwaMemAlignment> bwaMemAlignments = alignments.get(i);
            final SAMReadGroupRecord readGroup = getReadGroup(readName);
            final List<SAMRecord> samRecords =
                    BwaMemAlignmentUtils.toSAMStreamForRead(readName, bases.get(i), bwaMemAlignments,
                            header, refNames, readGroup).collect(Collectors.toList());
            result.add(samRecords);
        }
        return result;
    }

    // TODO: 5/8/18 how?
    private static SAMReadGroupRecord getReadGroup(final String readName) {
        return null;
    }

    public SAMSequenceDictionary getDict() {
        return SAMSequenceDictionaryExtractor
                .extractDictionary(IOUtils.getPath(fasta.toString().replace(".fasta", ".dict")));
    }

    @Override
    public void close() throws IOException {
        aligner.close();
        index.close();
    }
}