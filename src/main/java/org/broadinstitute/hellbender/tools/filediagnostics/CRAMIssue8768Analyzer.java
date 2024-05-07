package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.build.CRAMReferenceRegion;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.ref.ReferenceContext;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.*;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.Tuple;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.*;

/**
 * A diagnostic class that analyzes a CRAM input file to look for conditions that trigger issue 8768.
 */
public class CRAMIssue8768Analyzer extends HTSAnalyzer {
    private final ReferenceSource referenceSource;
    private final OutputStream outputStream;
    private final CompressorCache compressorCache = new CompressorCache();
    private final boolean verbose;
    private final boolean echoToConsole;
    private final double mismatchRateThreshold;
    private final GATKPath tsvOutputPath;
    private int retCode = 0;

    private SAMFileHeader samHeader = null;

    // everything we need to record for a bad contig
    private record BadContig(
            String contigName,
            List<ContainerStats> badContainers
    ) { }

    // everything we need to record for a (good or bad) container
    private record ContainerStats(
            int containerOrdinal, // container ordinal # within the contig
            boolean isBad,        // true if this container is bad
            AlignmentContext alignmentContext,  // reference ID, alignment start, alignment span
            long misMatchCount,   // count of mismatched bases
            double misMatchRate   // rate of mismatched bases (mismatches/total bases)
    ) { }

    public CRAMIssue8768Analyzer(
            final GATKPath inputPath,
            final GATKPath textOutputPath,
            final GATKPath tsvOutputPath,
            final GATKPath referencePath,
            final double mismatchRateThreshold,
            final boolean verbose,
            final boolean echoToConsole) {
        super(inputPath, textOutputPath);
        this.verbose = verbose;
        this.echoToConsole = echoToConsole;
        this.tsvOutputPath = tsvOutputPath;
        this.mismatchRateThreshold = mismatchRateThreshold;

        referenceSource = new ReferenceSource(referencePath.toPath());
        try {
            outputStream = Files.newOutputStream(this.outputPath.toPath());
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    @Override
    public void close() throws IOException {
        if (outputStream != null) {
            outputStream.close();
        }
    }

    public int getRetCode() {
        return retCode;
    }
    protected void emitln(final String s) {
        try {
            outputStream.write(s.getBytes());
            outputStream.write('\n');
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Some background:
     *
     * Each CRAM container has an associated htsjdk AlignmentContext that models if/how the contents of the
     * container relates to the reference. An AlignmentContext includes a ReferenceContext, and depending on the
     * type of the ReferenceContext, possibly an alignment start and alignment span. There are 3 possible types of
     * ReferenceContexts:
     *
     * 1. SINGLE_REF: The container contains only reads that are mapped to a single reference contig, in which case
     * the referenceID for that ReferenceContext is the contig ID for the associated contig. This is the most
     * common case, and is the only type of ReferenceContext for which the alignment start and span are meaningful.
     * Call isMappedSingleRef() to determine if the ReferenceContext is SINGLE_REF.
     *
     * Note that it is an error (the code throws) to you attempt to query a ReferenceContext for it's contig ID if
     * the ReferenceContext is not SINGLE_REF.
     *
     * 2. MULTI_REF: The container contains reads that are mapped to more than one reference contig. This is an
     * optimization used primarily when there aren't enough reads mapped to a single reference contig to justify
     * putting the reads into a separate container. Reads in these containers are not reference compressed, and
     * AlignmentContexts for MULTI_REF containers have no meaningful start/span values. Call isMappedMultiRef() to
     * determine if the ReferenceContext is MULTI_REF.
     *
     * 3. UNMAPPED_UNPLACED: The container contains only unmapped unplaced reads. start/span are irrelevant. Call
     * isUnmappedUnplaced() to determine if the ReferenceContext is UNMAPPED_UNPLACED.
     */
    public void doAnalysis() {
        final Map<String, BadContig> badContigs = new LinkedHashMap<>(); // contig name, BadContig
        final List<ContainerStats> goodContainers = new ArrayList<>(); // good containers, for all contigs
        int containerOrdinalForContig = 0;
        final int NUMBER_OF_GOOD_CONTAINERS_PER_CONTIG_TO_REPORT = 4;
        int nGoodContainersReportedForContig = 0;

        // these are HTSJDK CRAM container alignment contexts, not the GATK kind you're thinking of
        AlignmentContext previousAlignmentContext = null;

        try (final SeekablePathStream seekableStream = new SeekablePathStream(this.inputPath.toPath())) {
            List<ContainerStats> badContainersForContig = new ArrayList<>();
            final CramHeader cramHeader = analyzeCRAMHeader(seekableStream);
            samHeader = Container.readSAMFileHeaderContainer(
                    cramHeader.getCRAMVersion(),
                    seekableStream,
                    inputPath.getRawInputString());

            // iterate through the containers looking for ones with alignment spans that trigger the issue
            for (boolean isEOF = false; !isEOF;) {
                final Container container = new Container(
                        cramHeader.getCRAMVersion(),
                        seekableStream,
                        seekableStream.position());
                containerOrdinalForContig++;

                // reject CRAMs with properties that clearly indicate they were not written by GATK/Picard
                if (isForeignCRAM(container)) {
                    return;
                }

                if (previousAlignmentContext == null) {
                    // first container in the whole file can't be bad
                    recordContainerStats(goodContainers, false, container, containerOrdinalForContig);
                    nGoodContainersReportedForContig++;
                } else if (!previousAlignmentContext.getReferenceContext().equals(
                        container.getAlignmentContext().getReferenceContext())) {
                    // this is the first container for a new reference context; emit any bad containers accumulated
                    // for the previous reference context/contig, and reset state for the next one
                    if (badContainersForContig.size() > 0) {
                        recordContigStats(badContigs, badContainersForContig, previousAlignmentContext);
                        badContainersForContig = new ArrayList<>();
                    }
                    containerOrdinalForContig = 1;
                    // the first container for a reference context is never bad, so always add it to the good list
                    recordContainerStats(goodContainers, false, container, containerOrdinalForContig);
                    nGoodContainersReportedForContig = 1;
                } else if (previousAlignmentContext.getReferenceContext().isMappedSingleRef() &&
                        (previousAlignmentContext.getAlignmentStart() == 1)) {
                    // we're on the same reference context as the previous container, and the previous container
                    // was mapped to position 1, so if this container is mapped, it's a candidate for bad, whether
                    // it starts at position 1 (the rare case where there is more than one bad container) or not
                    // (the common case where this is the one bad container for this contig)
                    recordContainerStats(badContainersForContig, true, container, containerOrdinalForContig);
                } else {
                    // we're on the same reference context as the previous container, and the previous container
                    // was NOT mapped to position 1, so this container is not bad - add it to the list of good
                    // containers
                    if (verbose || nGoodContainersReportedForContig < NUMBER_OF_GOOD_CONTAINERS_PER_CONTIG_TO_REPORT) {
                        recordContainerStats(goodContainers, false,container, containerOrdinalForContig);
                        nGoodContainersReportedForContig++;
                    }
                }

                previousAlignmentContext = new AlignmentContext(
                        container.getAlignmentContext().getReferenceContext(),
                        container.getAlignmentContext().getAlignmentStart(),
                        container.getAlignmentContext().getAlignmentSpan());
                isEOF = container.isEOF();
            }
        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }

        retCode = printTextResults(badContigs, goodContainers);
        if (tsvOutputPath != null) {
            printTSVResults(badContigs, goodContainers, tsvOutputPath);
        }
    }

    /**
     * Display metadata for a CRAM file header.
     */
    private CramHeader analyzeCRAMHeader(InputStream is) {
        final CramHeader cramHeader = CramIO.readCramHeader(is);
        emitln("CRAM File Name: " + inputPath.toPath().getFileName());
        emitln("CRAM Version: " + cramHeader.getCRAMVersion().toString());
        emitln("CRAM ID Contents: " + String.format("%s", Base64.getEncoder().encodeToString(cramHeader.getId())));
        return cramHeader;
    }

    // reject any inputs that have containers that are reference-less; have multiple slices per container;
    // or have slices with an embedded reference, since these indicate that the file was not written by GATK/Picard.
    // it is in theory possible that the file could have been written by some other client of htsjdk (i.e., the
    // htsjdk tests can write such file), but analyzing such files is beyond the scope of this tool
    private boolean isForeignCRAM(final Container container) {
        final List<Slice> slices = container.getSlices();
        if (slices.size() > 1 ) {
            emitln("Multi-slice container detected. This file was not written by GATK or Picard.");
            return true;
        } else if (container.getAlignmentContext().getReferenceContext().isMappedSingleRef() &&
                !container.getCompressionHeader().isReferenceRequired()) {
            emitln("Reference-less container detected. This file was not written by GATK or Picard.");
            return true;
        }
        for (final Slice slice : slices) {
            if (slice.getEmbeddedReferenceContentID() != Slice.EMBEDDED_REFERENCE_ABSENT_CONTENT_ID) {
                emitln(String.format("Embedded reference block (ID %d) detected. This file was not written by GATK or Picard.",
                        slice.getEmbeddedReferenceContentID()));
                return true;
            }
        }
        return false;
    }

    private void recordContainerStats(
            final List<ContainerStats> targetList,
            final boolean isBad,
            final Container container,
            final int containerOrdinal) {
        // don't even try to compute stats for unmapped-unplaced containers or multi-ref containers
        if (container.getAlignmentContext().getReferenceContext().isMappedSingleRef()) {
            final Tuple<Long, Double> containerStats = analyzeContainerBaseMismatches(container);
            targetList.add(new ContainerStats(
                    containerOrdinal,
                    isBad,
                    container.getAlignmentContext(),
                    containerStats.a,   // mismatches
                    containerStats.b)); // mismatchPercent
        }
    }

    private void recordContigStats(
            final Map<String, BadContig> badContigs,
            final List<ContainerStats> badContainers,
            final AlignmentContext previousAlignmentContext) {
        if (null != badContigs.putIfAbsent(
                previousAlignmentContext.getReferenceContext().toString(),
                new BadContig(
                        previousAlignmentContext.getReferenceContext().toString(),
                        badContainers))) {
            throw new IllegalStateException(
                    String.format(
                            "Attempt to add a bad contig (%s) more than once",
                            previousAlignmentContext.getReferenceContext().toString()));
        }
    }

    /**
     * Analyze base mismatches CRAM file container.
     * return true if container is EOF container
     */
    private Tuple<Long, Double> analyzeContainerBaseMismatches(final Container container) {
        final SAMSequenceDictionary sequenceDictionary = samHeader.getSequenceDictionary();
        final List<SAMRecord> actualSAMRecords = container.getSAMRecords(
                ValidationStringency.LENIENT,
                new CRAMReferenceRegion(referenceSource, samHeader.getSequenceDictionary()),
                compressorCache,
                samHeader);

        final CRAMReferenceRegion localReferenceRegion = new CRAMReferenceRegion(referenceSource, sequenceDictionary);
        // SequenceUtil.countMismatches wants the full contig's reference bases
        localReferenceRegion.fetchReferenceBases(container.getAlignmentContext().getReferenceContext().getReferenceContextID());
        final byte referenceBases[] = localReferenceRegion.getCurrentReferenceBases();

        long totalBases = 0;
        long misMatches = 0;
        for (SAMRecord samRec : actualSAMRecords) {
            totalBases += (long) samRec.getReadLength();
            misMatches += (long) SequenceUtil.countMismatches(samRec, referenceBases);
        }
        return new Tuple<>(totalBases, misMatches/(double) totalBases);
    }

    private int printTextResults(final Map<String, BadContig> badContigs, final List<ContainerStats> goodContainers) {
        int retCode;
        if (badContigs.isEmpty()) {
            final String NO_CORRUPT_CONTAINERS = "\n**********************NO CORRUPT CONTAINERS DETECTED**********************";
            emitln(NO_CORRUPT_CONTAINERS);
            // always emit the results summary to console
            System.out.println(NO_CORRUPT_CONTAINERS);
            retCode = 0;
        } else {
            final String POSSIBLE_CORRUPT_CONTAINERS = "\n**********************!!!!!Possible CORRUPT CONTAINERS DETECTED!!!!!**********************:\n";
            emitln(POSSIBLE_CORRUPT_CONTAINERS);
            // always emit the results summary to console
            System.out.println(POSSIBLE_CORRUPT_CONTAINERS);
            retCode = 1;
        }

        // before we print out the containers, print out the stats for both the good and the bad containers
        final int totalGoodContainers = goodContainers.size();
        final double sumOfGoodMismatchRates = goodContainers.stream().mapToDouble(c -> c.misMatchRate).sum();
        final double averageGoodMismatchRate = sumOfGoodMismatchRates / totalGoodContainers;
        final String avgGoodMismatchStr = String.format("Average mismatch rate for presumed good containers: %f", averageGoodMismatchRate);
        emitln(avgGoodMismatchStr);
        if (echoToConsole) {
            System.out.println(avgGoodMismatchStr);
        }

        if (!badContigs.isEmpty()) {
            final int totalBadContainers = badContigs.values().stream().mapToInt(bc -> bc.badContainers().size()).sum();
            final double sumOfBadMismatchRates = badContigs.values().stream().mapToDouble(
                    bc -> bc.badContainers().stream().mapToDouble(c -> c.misMatchRate).sum()).sum();
            final double averageBadMismatchRate = sumOfBadMismatchRates / totalBadContainers;
            final String avgBadMismatchStr = String.format("Average mismatch rate for suspected bad containers: %f", averageBadMismatchRate);
            emitln(avgBadMismatchStr);
            if (echoToConsole) {
                System.out.println(avgBadMismatchStr);
            }

            if (averageBadMismatchRate > mismatchRateThreshold) {
                final String exceedThresholdStr = String.format(
                        "The average base mismatch rate of %f for suspected bad containers exceeds the threshold rate of %1.2f, and indicates this file may be corrupt.",
                        averageBadMismatchRate,
                        mismatchRateThreshold);
                emitln(exceedThresholdStr);
                if (echoToConsole) {
                    System.out.println(exceedThresholdStr);
                }
            }

            // now emit the list of corrupt containers for each bad contig
            emitln("\nSuspected CORRUPT Containers:");
            for (final Map.Entry<String, BadContig> entry : badContigs.entrySet()) {
                for (final ContainerStats badContainer : entry.getValue().badContainers()) {
                    final String badStatStr = String.format("  Ordinal: %d (%s) Mismatch Rate/Count: %f/%d",
                            badContainer.containerOrdinal,
                            badContainer.alignmentContext.toString(),
                            badContainer.misMatchRate,
                            badContainer.misMatchCount);
                    emitln(badStatStr);
                    if (echoToConsole) {
                        System.out.println(badStatStr);
                    }
                }
            }
        }

        emitln("\nPresumed GOOD Containers:");
        int lastContig = ReferenceContext.UNINITIALIZED_REFERENCE_ID;
        for (final ContainerStats goodContainer : goodContainers) {
            if (lastContig != ReferenceContext.UNINITIALIZED_REFERENCE_ID &&
                    lastContig != goodContainer.alignmentContext.getReferenceContext().getReferenceContextID()) {
                emitln("");
                if (echoToConsole) {
                    System.out.println("");
                }
            }
            lastContig = goodContainer.alignmentContext.getReferenceContext().getReferenceContextID();
            final String goodDetailStr = String.format("  Ordinal: %d (%s) Mismatch Rate/Count: %f/%d",
                    goodContainer.containerOrdinal,
                    goodContainer.alignmentContext.toString(),
                    goodContainer.misMatchRate,
                    goodContainer.misMatchCount);
            emitln(goodDetailStr);
            if (echoToConsole) {
                System.out.println(goodDetailStr);
            }
        }
        return retCode;
    }

    // write the results out to a machine-readable tsv file
    private void printTSVResults(
            final Map<String, BadContig> badContigs,
            final List<ContainerStats> goodContainers,
            final GATKPath tsvOutputPath) {
        // File name, contig name, container ordinal, good or bad, mismatch rate
        final TableColumnCollection columnNames = new TableColumnCollection(
                "file_name",    // file name
                "contig_name",          // contig name
                "container_ordinal",    // container ordinal
                "container_is_bad",     // good or bad, 1 or 0
                "mismatch_rate",        // mismatch rate (double)
                "alignment_start",      // alignment start (int)
                "alignment_span"        // alignment span (int)
        );

        try (final TableWriter<ContainerStats> tsvWriter = new TableWriter<>(tsvOutputPath.toPath(), columnNames) {
            @Override
            protected void composeLine(final ContainerStats containerStats, final DataLine dataLine) {
                dataLine.set("file_name", inputPath.toPath().getFileName().toString())
                        .set("contig_name", samHeader.getSequenceDictionary().getSequence(containerStats.alignmentContext.getReferenceContext().getReferenceContextID()).getSequenceName())
                        .set("container_ordinal", containerStats.containerOrdinal)
                        .set("container_is_bad", containerStats.isBad ? 1 : 0)
                        .set("mismatch_rate", containerStats.misMatchRate)
                        .set("alignment_start", containerStats.alignmentContext.getAlignmentStart())
                        .set("alignment_span", containerStats.alignmentContext.getAlignmentSpan());
                }
            })
        {
            tsvWriter.writeHeaderIfApplies();
            if (badContigs.isEmpty()) {
                tsvWriter.writeComment("No bad containers detected");
            } else {
                tsvWriter.writeComment("Bad containers:");
                for (final Map.Entry<String, BadContig> entry : badContigs.entrySet()) {
                    tsvWriter.writeAllRecords(entry.getValue().badContainers());
                }
            }
            if (goodContainers.isEmpty()) {
                tsvWriter.writeComment("No good mapped containers detected");
            } else {
                tsvWriter.writeComment("Good containers:");
                tsvWriter.writeAllRecords(goodContainers);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
