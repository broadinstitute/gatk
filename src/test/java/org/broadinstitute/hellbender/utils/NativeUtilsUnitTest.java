package org.broadinstitute.hellbender.utils;

import org.apache.commons.lang3.SystemUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class NativeUtilsUnitTest extends GATKBaseTest {

    @Test
    public void testArchitectureDetectorsAreMutuallyConsistentWithRuntime() {
        final String arch = SystemUtils.OS_ARCH;
        final boolean x86 = NativeUtils.runningOn64BitX86Architecture();
        final boolean aarch64 = NativeUtils.runningOnAarch64Architecture();

        // On the architectures GATK actually targets, exactly one of the two detectors must be true,
        // and it must agree with the JVM-reported os.arch.
        if ("x86_64".equals(arch) || "amd64".equals(arch)) {
            Assert.assertTrue(x86, "x86 detector should be true on os.arch=" + arch);
            Assert.assertFalse(aarch64, "aarch64 detector should be false on os.arch=" + arch);
        } else if ("aarch64".equals(arch) || "arm64".equals(arch)) {
            Assert.assertTrue(aarch64, "aarch64 detector should be true on os.arch=" + arch);
            Assert.assertFalse(x86, "x86 detector should be false on os.arch=" + arch);
        }

        // The two families are always mutually exclusive regardless of platform.
        Assert.assertFalse(x86 && aarch64, "x86 and aarch64 detectors must never both be true");
    }

    @Test
    public void testRunningOnArchitectureMatchesAarch64Helper() {
        final boolean viaHelper = NativeUtils.runningOnAarch64Architecture();
        final boolean viaRaw = NativeUtils.runningOnArchitecture("aarch64")
                || NativeUtils.runningOnArchitecture("arm64");
        Assert.assertEquals(viaHelper, viaRaw);
    }

    @Test
    public void testRunningOnArchitectureRejectsNonMatching() {
        Assert.assertFalse(NativeUtils.runningOnArchitecture("not-a-real-architecture"));
    }
}
