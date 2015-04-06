package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * Created by jcarey on 3/14/14.
 */
public class BciFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        for (final Integer tile : tiles) {
            buffer.putInt(tile);
            buffer.putInt(1);
        }
    }

    @Override
    protected boolean addLeadingZeros() {
        return false;
    }

    @Override
    protected int bufferSize() {
        return 8 * tiles.size();
    }

    public void fakeBciFile(final File bci, final List<Integer> expectedTiles) throws IOException {
        tiles = expectedTiles;
        final FileOutputStream fileOutputStream = new FileOutputStream(bci);
        final FileChannel channel = fileOutputStream.getChannel();
        final ByteBuffer buffer = ByteBuffer.allocate(8 * expectedTiles.size());
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        fakeFile(buffer);
        buffer.flip();

        channel.write(buffer);
        channel.force(true);

        CloserUtil.close(channel);
        CloserUtil.close(fileOutputStream);
    }
}
