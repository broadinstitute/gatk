package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link Target}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetUnitTest {

    @Test()
    public void testConstructorWithOnlyName() {
        final Target subject = new Target("my-name");
        Assert.assertEquals(subject.getName(), "my-name");
        Assert.assertNull(subject.getInterval());
    }

    @Test()
    public void testGetName() {
        final Target subject = new Target("my-name");
        Assert.assertEquals(subject.getName(), "my-name");
    }

    @Test()
    public void testConstructorWithInterval() {
        final Target subject1 = new Target("my-name", null);
        Assert.assertEquals(subject1.getName(), "my-name");
        Assert.assertNull(subject1.getInterval());
        final Target subject2 = new Target("my-name", new SimpleInterval("1", 1, 2));
        Assert.assertEquals(subject2.getName(), "my-name");
        Assert.assertEquals(subject2.getInterval(), new SimpleInterval("1", 1, 2));
    }

    @Test()
    public void testGetContigWithInterval() {
        final Target subject = new Target("my-name", new SimpleInterval("1", 1, 2));
        Assert.assertEquals(subject.getContig(), "1");
    }

    @Test()
    public void testGetStartWithInterval() {
        final Target subject = new Target("my-name", new SimpleInterval("1", 1, 2));
        Assert.assertEquals(subject.getStart(), 1);
    }

    @Test()
    public void testGetEndWithInterval() {
        final Target subject = new Target("my-name", new SimpleInterval("1", 1, 2));
        Assert.assertEquals(subject.getEnd(), 2);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testGetContigWithoutInterval() {
        final Target subject = new Target("my-name", null);
        subject.getContig();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testGetStartWithoutInterval() {
        final Target subject = new Target("my-name", null);
        subject.getStart();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testGetEndWithoutInterval() {
        final Target subject = new Target("my-name", null);
        Assert.assertEquals(subject.getEnd(), 2);
    }

    @Test
    public void testLength() {
        final Target subject1 = new Target("target", new SimpleInterval("chr",1,1));
        final Target subject2 = new Target("target", new SimpleInterval("chr",10,100));
        Assert.assertEquals(subject1.length(), 1);
        Assert.assertEquals(subject2.length(), 91);
    }

    @Test()
    public void testEquals() {
        final Target subject1 = new Target("my-name");
        final Target subject1bis = new Target("my-name", new SimpleInterval("1", 1, 2));
        final Target subject1tris = new Target("my-name", new SimpleInterval("2", 1, 2));
        final Target subject2 = new Target("other-name");

        Assert.assertTrue(subject1.equals(subject1bis));
        Assert.assertTrue(subject1.equals(subject1tris));
        Assert.assertTrue(subject1bis.equals(subject1));
        Assert.assertTrue(subject1bis.equals(subject1tris));
        Assert.assertTrue(subject1tris.equals(subject1));
        Assert.assertTrue(subject1tris.equals(subject1bis));
        Assert.assertTrue(subject1.equals(subject1));
        Assert.assertTrue(subject1bis.equals(subject1bis));
        Assert.assertTrue(subject1tris.equals(subject1tris));
        Assert.assertTrue(subject2.equals(subject2));

        Assert.assertFalse(subject1.equals(null));
        Assert.assertFalse(subject1.equals(subject2));
        Assert.assertFalse(subject2.equals(subject1));
        Assert.assertFalse(subject1bis.equals(subject2));
        Assert.assertFalse(subject2.equals(subject1bis));
        Assert.assertFalse(subject1tris.equals(subject2));
        Assert.assertFalse(subject2.equals(subject1tris));
    }
}
