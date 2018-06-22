package org.broadinstitute.hellbender.testutils;

import com.google.common.annotations.VisibleForTesting;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.NonWritableChannelException;
import java.nio.channels.SeekableByteChannel;
import java.util.function.Function;

/**
 * A SeekableByteChannel wrapper for testing.
 */
public class XorWrapper implements SeekableByteChannel {

  private final byte KEY;
  private final SeekableByteChannel in;

  public static Function<SeekableByteChannel,SeekableByteChannel> forKey(byte key) {
    return (SeekableByteChannel c) -> new XorWrapper(c, key);
  }

  public XorWrapper(SeekableByteChannel in, byte key) {
    this.in = in;
    this.KEY = key;
  }

  @VisibleForTesting
  public static void xor(byte[] buffer, byte key, int start, int count) {
    for (int i = start; i < start + count; i++) {
      buffer[i] = (byte)(buffer[i] ^ key);
    }
  }

  @Override
  public int read(ByteBuffer buf) throws IOException {
    if (!buf.hasArray()) {
      throw new RuntimeException("Sorry, support for readonly arrays hasn't been written yet.");
    }
    int pos = buf.position();
    int count = in.read(buf);
    xor(buf.array(), KEY, pos, count);
    return count;
  }

  @Override
  public int write(ByteBuffer src) throws IOException {
    throw new NonWritableChannelException();
  }

  @Override
  public long position() throws IOException {
    return in.position();
  }

  @Override
  public SeekableByteChannel position(long newPosition) throws IOException {
    in.position(newPosition);
    return this;
  }

  @Override
  public long size() throws IOException {
    return in.size();
  }

  @Override
  public SeekableByteChannel truncate(long size) throws IOException {
    throw new NonWritableChannelException();
  }

  @Override
  public boolean isOpen() {
    return in.isOpen();
  }

  @Override
  public void close() throws IOException {
    in.close();
  }
}
