package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.util.*;

public final class SVFileUtils {

    private static final String REFERENCE_GAP_INTERVAL_FILE_COMMENT_LINE_PROMPT = "#";

    public static void writeSAMFile(final String outputName, final Iterator<SAMRecord> alignments, final SAMFileHeader header,
                                    final boolean preOrdered) {
        Utils.nonNull(alignments, "provided alignments to write out is null");
        Utils.nonNull(header, "provided header for outputting sam file is null");
        Utils.nonNull(outputName, "provided output name is null");

        try (final SAMFileWriter writer = getSAMFileWriter(outputName, header, preOrdered)) {
            alignments.forEachRemaining(writer::addAlignment);
        } catch ( final UncheckedIOException ie) {
            throw new GATKException("Can't write SAM file to the specified location: " + outputName, ie);
        }
    }

    /**
     * Supported format: BAM and SAM. CRAM is unsupported. Other requested type will default to SAM.
     */
    private static SAMFileWriter getSAMFileWriter(final String outputName, final SAMFileHeader header, final boolean preOrdered) {
        final int idx = outputName.lastIndexOf(".");
        if ( idx < 0) {
            throw new IllegalArgumentException("Provided path doesn't have a proper extension: " + outputName);
        }
        final SAMFileWriterFactory factory = new SAMFileWriterFactory()
                .setCreateIndex(preOrdered && header.getSortOrder() == SAMFileHeader.SortOrder.coordinate);

        final String fileExtension = outputName.substring(idx + 1, outputName.length());
        if (fileExtension.endsWith(SamReader.Type.BAM_TYPE.fileExtension())) {
            return factory.makeBAMWriter(header, preOrdered, BucketUtils.createFile(outputName));
        } else if (fileExtension.endsWith(SamReader.Type.SAM_TYPE.fileExtension())) {
            return factory.makeSAMWriter(header, preOrdered, BucketUtils.createFile(outputName));
        } else if (fileExtension.endsWith(SamReader.Type.CRAM_TYPE.fileExtension())) {
            throw new UnsupportedOperationException("We currently don't support CRAM output");
        } else {
            return factory.makeSAMWriter(header, preOrdered, BucketUtils.createFile(outputName));
        }
    }

    /**
     * Read a file of kmers.
     * Each line must be exactly
     * {@link org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection#KMER_SIZE}
     * characters long, and must match [ACGT]*.
     */
    public static Set<SVKmer> readKmersFile(final String kmersFilePath, final int kSize) {
        Utils.nonNull(kmersFilePath, "provided path for file containing kmers is null");
        Utils.validateArg(kSize > 0, "provided k-size is non positive: " + kSize);

        final Set<SVKmer> kmers;

        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(kmersFilePath))) ) {
            final long fileLength = BucketUtils.fileSize(kmersFilePath);
            kmers = new HopscotchSet<>((int)(fileLength/(kSize+1)));
            String line;
            while ( (line = rdr.readLine()) != null ) {
                if ( line.length() != kSize ) {
                    throw new GATKException("SVKmer kill set contains a line of length " + line.length() +
                            " but we were expecting K = " + kSize);
                }

                final SVKmerizer kmerizer = new SVKmerizer(line, kSize, 1, new SVKmerLong(kSize));
                if ( !kmerizer.hasNext() ) {
                    throw new GATKException("Unable to kmerize the kmer kill set string '" + line + "'.");
                }

                kmers.add(kmerizer.next());
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read kmers from " + kmersFilePath, ioe);
        }

        return kmers;
    }

    /** Write kmers to file. */
    public static <KType extends SVKmer> void writeKmersFile(final String kmersFilePath, final int kSize,
                                                             final Collection<KType> kmers) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(kmersFilePath))) ) {
            for ( final KType kmer : kmers ) {
                writer.write(kmer.toString(kSize));
                writer.write('\n');
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to write kmers to " + kmersFilePath, ioe);
        }
    }

    /** Read intervals from file. */
    public static List<SVInterval> readIntervalsFile(final String intervalsFilePath,
                                                     final Map<String, Integer> contigNameMap ) {
        Utils.nonNull(intervalsFilePath, "provided intervals file path is null");
        Utils.nonNull(contigNameMap, "provided map for contig index lookup is null");

        final List<SVInterval> intervals;
        try ( final BufferedReader rdr =
                      new BufferedReader(new InputStreamReader(BucketUtils.openFile(intervalsFilePath))) ) {
            final long sizeGuess = BucketUtils.fileSize(intervalsFilePath)/25; // 25 is a guess on file line length
            intervals = new ArrayList<>((int)sizeGuess);
            String line;
            int lineNo = 0;
            while ( (line = rdr.readLine()) != null ) {
                ++lineNo;
                if (line.startsWith(REFERENCE_GAP_INTERVAL_FILE_COMMENT_LINE_PROMPT)) {
                    continue;
                }
                final List<String> tokens = Utils.split(line, '\t');
                if ( tokens.size() != 3 ) {
                    throw new GATKException("Interval file " + intervalsFilePath + " line " +
                            lineNo + " did not contain 3 columns: " + line);
                }
                try {
                    final Integer contigId = contigNameMap.get(tokens.get(0));
                    if ( contigId == null ) throw new GATKException("contig name " + tokens.get(0) + " not in dictionary");
                    final int start = Integer.valueOf(tokens.get(1));
                    final int end = Integer.valueOf(tokens.get(2));
                    intervals.add(new SVInterval(contigId, start, end));
                }
                catch ( final Exception e ) {
                    throw new GATKException("Unable to parse interval file " + intervalsFilePath + " line " + lineNo + ": " + line, e);
                }
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to read intervals from " + intervalsFilePath, ioe);
        }
        return intervals;
    }

    /** Write intervals to a file. */
    public static void writeIntervalsFile( final String intervalsFilePath,
                                           final Collection<SVInterval> intervals, final List<String> contigNames ) {
        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(intervalsFilePath)))) {
            for (final SVInterval interval : intervals) {
                final String seqName = contigNames.get(interval.getContig());
                writer.write(seqName + "\t" + interval.getStart() + "\t" + interval.getEnd() + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write intervals file " + intervalsFilePath, ioe);
        }
    }
}
