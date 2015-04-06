package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclReader;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.util.*;
import java.util.zip.GZIPOutputStream;

public abstract class FileFaker {

    int size;
    List<Integer> tiles;

    protected abstract void fakeFile(ByteBuffer buffer);

    protected abstract boolean addLeadingZeros();

    protected abstract int bufferSize();

    public void fakeFile(final File base, final int tile, final int lane, final String extension) throws IOException {
        fakeFile(base, Collections.singletonList(tile), lane, extension);
    }

    public void fakeFile(final File base, final List<Integer> expectedTiles, final int lane, final String extension)
            throws IOException {
        if (base.exists() || base.mkdirs()) {
            this.tiles = expectedTiles;
            final File fakeFile;
            if (expectedTiles.size() == 1) {
                String longTileName = String.valueOf(tiles.get(0));
                if (addLeadingZeros()) {
                    while (longTileName.length() < 4) {
                        longTileName = "0" + longTileName;
                    }
                }
                fakeFile = new File(base, String.format("s_%d_%s%s", lane, longTileName, extension));
            } else {
                fakeFile = new File(base, String.format("s_%s%s", lane, extension));
            }

            fakeFile(fakeFile, bufferSize());
        }

    }

    public void fakeFile(final File cycleFile, Integer size) throws IOException {
        if (size == null) {
            size = 1;
        }
        this.size = size;

        final OutputStream outputStream;
        if (BclReader.isGzipped(cycleFile)) outputStream = new GZIPOutputStream(new FileOutputStream(cycleFile));
        else if (BclReader.isBlockGzipped(cycleFile)) outputStream = new BlockCompressedOutputStream(cycleFile);
        else outputStream = new FileOutputStream(cycleFile);

        final WritableByteChannel channel = Channels.newChannel(outputStream);
        final ByteBuffer buffer = ByteBuffer.allocate(this.bufferSize());
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        fakeFile(buffer);

        buffer.flip();

        channel.write(buffer);

        channel.close();
        outputStream.close();
    }
}
