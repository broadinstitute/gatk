package org.broadinstitute.hellbender.utils.nio;

import org.broadinstitute.hellbender.testutils.XorWrapper;
import org.testng.annotations.Test;
import org.testng.Assert;

public class XorWrapperTest {

  @Test
  public void testXor() throws Exception {
    byte key = 42;
    byte[] a = new byte[256];
    byte[] b = new byte[256];
    for (int i=0; i < a.length; i++) {
      a[i] = (byte)i;
      b[i] = (byte)i;
    }
    XorWrapper.xor(b, key, 0, 256);
    for (int i=0; i < a.length; i++) {
      Assert.assertNotEquals( a[i], b[i] );
    }
    XorWrapper.xor(b, key, 0, 256);
    for (int i=0; i < a.length; i++) {
      Assert.assertEquals( a[i], b[i] );
    }
  }

}