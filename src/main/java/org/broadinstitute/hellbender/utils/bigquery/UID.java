package org.broadinstitute.hellbender.utils.bigquery;
import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.io.Serializable;
import java.security.SecureRandom;

public final class UID implements Serializable {

    private static int hostUnique;
    private static boolean hostUniqueSet = false;

    private static final Object lock = new Object();
    private static long lastTime = System.currentTimeMillis();
    private static short lastCount = Short.MIN_VALUE;

    /** indicate compatibility with JDK 1.1.x version of class */
    private static final long serialVersionUID = 1086053664494604050L;

    /**
     * number that uniquely identifies the VM that this UID
     * was generated in with respect to its host and at the given time
     * @serial
     */
    private final int unique;

    /**
     * a time (as returned by {@link System#currentTimeMillis()})
     */
    private final long time;

    /**
     * 16-bit number to distinguish UID instances created
     * in the same VM with the same time value
     */
    private final short count;

    /**
     * Generates a UID that is unique over time with
     * respect to the host that it was generated on.
     */
    public UID() {

        synchronized (lock) {
            if (!hostUniqueSet) {
                hostUnique = (new SecureRandom()).nextInt();
                hostUniqueSet = true;
            }
            unique = hostUnique;
            if (lastCount == Short.MAX_VALUE) {
                boolean interrupted = Thread.interrupted();
                boolean done = false;
                while (!done) {
                    long now = System.currentTimeMillis();
                    if (now == lastTime) {
                        // wait for time to change
                        try {
                            Thread.sleep(1);
                        } catch (InterruptedException e) {
                            interrupted = true;
                        }
                    } else {
                        // If system time has gone backwards increase
                        // original by 1ms to maintain uniqueness
                        lastTime = (now < lastTime) ? lastTime+1 : now;
                        lastCount = Short.MIN_VALUE;
                        done = true;
                    }
                }
                if (interrupted) {
                    Thread.currentThread().interrupt();
                }
            }
            time = lastTime;
            count = lastCount++;
        }
    }

    /**
     * Creates a "well-known" UID.
     *
     * There are two 16 possible such well-known ids.
     *
     * UID created via this constructor will not
     * clash with any UIDs generated via the no-arg
     * constructor.
     *
     * @param   num number for well-known UID
     */
    public UID(short num) {
        unique = 0;
        time = 0;
        count = num;
    }

    /**
     * Constructs a UID given data read from a stream.
     */
    private UID(int unique, long time, short count) {
        this.unique = unique;
        this.time = time;
        this.count = count;
    }

    /**
     * Returns the hash code value for this UID.
     *
     * @return  the hash code value for this UID
     */
    public int hashCode() {
        return (int) time + (int) count;
    }

    /**
     * Compares the specified object with this UID for
     * equality.
     *
     * This method returns true if and only if the
     * specified object is a UID instance with the same
     * unique, time, and count
     * values as this one.
     *
     * @param   obj the object to compare this UID to
     *
     * @return  true if the given object is equivalent to
     * this one, and false otherwise
     */
    public boolean equals(Object obj) {
        if (obj instanceof UID) {
            UID uid = (UID) obj;
            return (unique == uid.unique &&
                    count == uid.count &&
                    time == uid.time);
        } else {
            return false;
        }
    }

    /**
     * Returns a string representation of this UID.
     *
     * @return  a string representation of this UID
     */
    public String toString() {
        return Integer.toString(unique,16) +
                Long.toString(time,16) +
                Integer.toString(Math.abs(count) ,16);
    }

    /**
     * Marshals a binary representation of this UIDto
     * a DataOutput instance.
     *
     * @param   out the DataOutput instance to write
     * this UID to
     *
     * @throws  IOException if an I/O error occurs while performing
     * this operation
     */
    public void write(DataOutput out) throws IOException {
        out.writeInt(unique);
        out.writeLong(time);
        out.writeShort(count);
    }

    /**
     * Constructs and returns a new UID instance by
     * unmarshalling a binary representation from an DataInput instance.
     *
     * @param   in the DataInputinstance to read UID from
     *
     * @return  unmarshalled UID instance
     *
     * @throws  IOException if an I/O error occurs while performing
     * this operation
     */
    public static UID read(DataInput in) throws IOException {
        int unique = in.readInt();
        long time = in.readLong();
        short count = in.readShort();
        return new UID(unique, time, count);
    }
}
