package org.broadinstitute.hellbender.tools.utilities;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.encoding.EncodingFactory;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.structure.*;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Base64;
import java.util.HashMap;
import java.util.Map;

//TODO: Add javadoc

/**
 * Analyzer for CRAM files. Displays metadata for each CRAM container,
 * slice, and block.
 *
 * Note: the analyzer does not require a reference for the since it only
 * enumerates metadata and doesn't attempt to dereference the reads.
 */
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Print information about a CRAM file",
        oneLineSummary = "Print information about a CRAM index",
        programGroup = OtherProgramGroup.class
)
public class AnalyzeCRAM extends CommandLineProgram {

    @Argument(
            shortName=StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName=StandardArgumentDefinitions.INPUT_LONG_NAME,
            doc="File to be analyzed",
            optional=false)
    private GATKPath inputPath; // File due to htsjdk signature

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "File to which to write file analysis information (if not specified, prints to standard output)",
            optional = false)
    private GATKPath outputPath;

    final Map<DataSeries, Long> externalDataSeriesDataSizes = new HashMap<>();
    final Map<Integer, Long> externalTagDataSizes = new HashMap<>();
    long coreBlocksDataSize = 0L;
//    long containerHeaderSize = 0L;
//    long compressionHeaderSize = 0L;
//    long sliceHeaderSize = 0L;
    long recordCount = 0;

    /**
     * Run the analyzer for the file.
     */
    protected Object doWork() {
        int containerCount = 0;
        try (final SeekableStream seekableStream = new SeekablePathStream(inputPath.toPath());
             final PrintStream outputStream = new PrintStream(new BufferedOutputStream(
                     outputPath == null ?
                     System.out :
                     outputPath.getOutputStream()))) {

            outputStream.println("\nCRAM File: " + inputPath.getRawInputString());
            final CramHeader cramHeader = analyzeCRAMHeader(seekableStream, outputStream);
            boolean isEOF = false;
            while (!isEOF) {
                long containerOffset = seekableStream.position();
                final Container container = new Container(
                        cramHeader.getCRAMVersion(),
                        seekableStream,
                        containerOffset);
                 isEOF = analyzeContainer(seekableStream, outputStream, container, ++containerCount) ;
            }
            outputStream.println("\nTotal Record Count: " + recordCount);
            emitDataDistribution(outputStream);
        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        return 1;
    }

    /**
     * Display metadata for a CRAM file header.
     *
     */
    public CramHeader analyzeCRAMHeader(final InputStream is, final PrintStream os) {
        final CramHeader cramHeader = CramIO.readCramHeader(is);
        os.println("CRAM Version: " + cramHeader.getCRAMVersion().toString());
        os.println("CRAM ID Contents: " + String.format("%s", Base64.getEncoder().encodeToString(cramHeader.getId())));

        final SAMFileHeader samHeader = Container.readSAMFileHeaderContainer(cramHeader.getCRAMVersion(), is, inputPath.getRawInputString());
        os.println("\n" + samHeader.toString());
        final SAMSequenceDictionary dict = samHeader.getSequenceDictionary();
        os.println(dict.toString());
        dict.getSequences().forEach(e -> os.println(e.toString()));
        return cramHeader;
    }

    /**
     * Display metadata for a CRAM file container.
     * return true if container is EOF container
     */
    public boolean analyzeContainer(final InputStream is, final PrintStream os, Container container, int containerCount) {
        final ContainerHeader containerHeader = container.getContainerHeader();
        os.println(String.format(
                "\n***Container #:%d %s byteOffset=%d",
                containerCount, containerHeader.toString(), container.getContainerByteOffset()));
        if (container.isEOF()) {
            return true;
        }
        // contentID to DataSeries (excludes tags ?)
        final Map<Integer, DataSeries> dataSeriesByContentID = analyzeCompressionHeader(os, container.getCompressionHeader());
        int sliceCount = 0;
        for (final Slice slice : container.getSlices()) {
            analyzeSlice(os, slice, ++sliceCount, dataSeriesByContentID);
        }
        return false;
    }

    public Map<Integer, DataSeries> analyzeCompressionHeader(final PrintStream os, final CompressionHeader compressionHeader) {
        //preservation map, data series encoding map, and tag encoding map
        analyzePreservationMap(os, compressionHeader);
        final Map<Integer, DataSeries> dataSeriesByContentID = analyzeDataSeriesEncodingMap(os, compressionHeader);
        analyzeTagEncodingMap(os, compressionHeader);
        return dataSeriesByContentID;
    }

    public void analyzePreservationMap(final PrintStream os, final CompressionHeader compressionHeader) {
        os.println(String.format(
                "Requires reference (%b); Preserved read names (%b); APDelta (%b)",
                    compressionHeader.isReferenceRequired(),
                    compressionHeader.isPreserveReadNames(),
                    compressionHeader.isAPDelta()));
    }

    public Map<Integer, DataSeries> analyzeDataSeriesEncodingMap(final PrintStream os, final CompressionHeader compressionHeader) {
        // Since each container can in theory use a different set of ids for the data series in that
        // container, we need to reconstruct the mapping for each container, and return it so that as
        // we process Slice blocks and track data, we can calculate the data series to which we should
        // accrue that block's data
        final Map<Integer, DataSeries> dataSeriesByContentID = new HashMap<>();

        os.println("\nData Series Encodings:\n");
        final CompressionHeaderEncodingMap encodingMap = compressionHeader.getEncodingMap();
        for (final DataSeries dataSeries: DataSeries.values()) {
            final EncodingDescriptor encodingDescriptor = encodingMap.getEncodingDescriptorForDataSeries(dataSeries);
            if (encodingDescriptor == null) {
                os.println(String.format("%-50s not present",
                        String.format("DataSeries (%s/%s)",
                                dataSeries.getCanonicalName(),
                                dataSeries.name())));
            } else {
                // the encoding map has an encoding descriptor for this data series; determine if its external
                // and if so, update the contentID to data series map so we can track how much data is used for
                // the blocks for this data series
                final Integer externalContentID = contentIDFromExternalDescriptor(dataSeries.getType(), encodingDescriptor);
                if (externalContentID != null) {
                    dataSeriesByContentID.put(externalContentID,  dataSeries);
                }
                os.println(String.format("%-50s %s",
                        String.format("DataSeries (%s/%s)",
                                dataSeries.getCanonicalName(),
                                dataSeries.name()),
                        encodingDescriptorAsString(
                                dataSeries.getType(),
                                encodingDescriptor)));
            }
        }
        return dataSeriesByContentID;
    }

    public void analyzeTagEncodingMap(final PrintStream os, final CompressionHeader compressionHeader) {
        os.println("\nTag Encodings:");
        for (final Map.Entry<Integer, EncodingDescriptor> entry : compressionHeader.getTagEncodingMap().entrySet()) {
            final Integer contentID = entry.getKey(); // is this content ID ?
            final EncodingDescriptor ep = entry.getValue();
            os.println(String.format("%-50s %s",
                    String.format("Content ID/Tag (%s/%s)",
                            contentID,
                            decomposeTagNameAndType(contentID)),
                     encodingDescriptorAsString(DataSeriesType.BYTE_ARRAY, ep)));
       }
    }

    /**
     * Display metadata for a CRAM container slice.
     *
     */
    public void analyzeSlice(
            final PrintStream os,
            final Slice slice,
            final int sliceCount,
            final Map<Integer, DataSeries> dataSeriesByContentID) {
        os.println(String.format("\n******Slice #: %d %s",
                sliceCount,
                slice.toString()));
        os.println(String.format("%-50s %s",
                "Header block ",
                slice.getSliceHeaderBlock()));
        os.println(String.format("%-50s %s",
                "Core block ",
                slice.getSliceBlocks().getCoreBlock()));

        if (slice.getEmbeddedReferenceContentID() != Slice.EMBEDDED_REFERENCE_ABSENT_CONTENT_ID) {
            os.println(String.format("Embedded reference block ID %d", slice.getEmbeddedReferenceContentID()));
        }

        slice.getSliceBlocks().getExternalContentIDs().forEach((id) -> {
            final String blockID = dataSeriesByContentID.get(id) == null ?
                    Integer.toString(id) : // not a fixed data series (i.e., tag block)
                    dataSeriesByContentID.get(id).getCanonicalName();
            os.println(String.format("%-50s %s",
                    String.format("External Block (%s):", blockID),
                    slice.getSliceBlocks().getExternalBlock(id).toString()));
        });

        updateDataDistribution(slice, dataSeriesByContentID);
        recordCount += slice.getNumberOfRecords();
    }

    final void updateDataDistribution(final Slice slice, final Map<Integer, DataSeries> dataSeriesByContentID) {
        final SliceBlocks sliceBlocks = slice.getSliceBlocks();

        coreBlocksDataSize += sliceBlocks.getCoreBlock().getCompressedContentSize();
        final Map<Integer, EncodingDescriptor> tagContentIDs = slice.getCompressionHeader().getTagEncodingMap();
        for (final Integer contentID : sliceBlocks.getExternalContentIDs()) {
            //if (tagContentIDs.containsKey(contentID)) {
            final DataSeries ds = dataSeriesByContentID.get(contentID);
            if (ds == null) {
                // accrue to tag data
                externalTagDataSizes.merge(
                        contentID,
                        new Long(slice.getSliceBlocks().getExternalBlock(contentID).getCompressedContentSize()),
                        (oldValue, increment) -> oldValue + increment);
            } else {
                // accrue to fixed DataSeries ID
                externalDataSeriesDataSizes.merge(
                        ds,
                        new Long(sliceBlocks.getExternalBlock(contentID).getCompressedContentSize()),
                        (oldValue, increment) -> oldValue + increment);
            }
        }
    }

    public void emitDataDistribution(final PrintStream os) {
        os.println("\nCore Block(s) Total: " + String.format("%,d\n", coreBlocksDataSize));
        os.println("External Data Series Totals (external block resolution - all core data encodings accrue to the core block):\n");
        for (final DataSeries ds : externalDataSeriesDataSizes.keySet()) {
            os.println(String.format("%s: %,d", ds, externalDataSeriesDataSizes.get(ds)));
        }

        os.println("\nTag Series Distribution:\n");
        for (final Map.Entry<Integer, Long> externalEntry : externalTagDataSizes.entrySet()) {
            final Integer contentID = externalEntry.getKey();
            final String tagName = String.format("%d (%s)", contentID, decomposeTagNameAndType(externalEntry.getKey()));
            os.println(String.format("%s: %,d", tagName, externalEntry.getValue()));
        }
    }

    private String encodingDescriptorAsString(final DataSeriesType dsType, final EncodingDescriptor descriptor) {
        final String encodingIDString = EncodingID.values()[descriptor.getEncodingID().getId()].toString();
        final CRAMEncoding<?> encoding = EncodingFactory.createCRAMEncoding(
                dsType, descriptor.getEncodingID(), descriptor.getEncodingParameters());
        return String.format("%s (%s)", encodingIDString, encoding.toString());
    }

    private Integer contentIDFromExternalDescriptor(final DataSeriesType dsType, final EncodingDescriptor descriptor) {
        final CRAMEncoding<?> encoding = EncodingFactory.createCRAMEncoding(
                dsType, descriptor.getEncodingID(), descriptor.getEncodingParameters());
        switch (descriptor.getEncodingID()) {
            case EXTERNAL:
                return ITF8.readUnsignedITF8(encoding.toSerializedEncodingParams());

            case BYTE_ARRAY_STOP:
                final ByteBuffer buf = ByteBuffer.wrap(encoding.toSerializedEncodingParams());
                buf.order(ByteOrder.LITTLE_ENDIAN);
                final byte stopByte = buf.get(); // discard it
                return ITF8.readUnsignedITF8(buf);

            // Everything else is either not external, or is hybrid (BYTE_ARRAY_LEN), so we can't
            // track the data for these at the block level since the data foes either partially or
            // fully into the core block
            case BYTE_ARRAY_LEN:
            case NULL:
            case GOLOMB:
            case HUFFMAN:
            case BETA:
            case SUBEXPONENTIAL:
            case GOLOMB_RICE:
            case GAMMA:
            default:
                return null;
        }
    }

    private String decomposeTagNameAndType(final int contentID) {
        return ReadTag.intToNameType4Bytes(contentID);
    }

}

