package org.broadinstitute.hellbender.utils.nio;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import org.testng.Assert;
import org.testng.annotations.Test;

public class LineIteratorTest {

  @Test
  public void testLineIteratorInForEach() throws IOException {
    try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
      final Path path = jimfs.getPath("test.txt");
      try (BufferedWriter bufferedWriter = Files.newBufferedWriter(path)) {
        bufferedWriter.write("Hello world\n");
        bufferedWriter.write("What's new?");
      }
      ArrayList<String> got = new ArrayList<>();
      for (String s: new LineIterator(path)) {
        got.add(s);
      }
      Assert.assertEquals(got.size(), 2);
      Assert.assertEquals(got.get(0), "Hello world");
      Assert.assertEquals(got.get(1), "What's new?");
    }
  }

  @Test
  public void testLineIteratorAsResource() throws IOException {
    try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
      final Path path = jimfs.getPath("test.txt");
      try (BufferedWriter bufferedWriter = Files.newBufferedWriter(path)) {
        bufferedWriter.write("Hello world\n");
        bufferedWriter.write("What's new?");
      }
      ArrayList<String> got = new ArrayList<>();
      try (LineIterator lines = new LineIterator(path)) {
        lines.forEach(s -> got.add(s));
      }
      Assert.assertEquals(got.size(), 2);
      Assert.assertEquals(got.get(0), "Hello world");
      Assert.assertEquals(got.get(1), "What's new?");
    }
  }

}