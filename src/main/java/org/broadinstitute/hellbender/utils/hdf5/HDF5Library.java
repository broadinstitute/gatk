package org.broadinstitute.hellbender.utils.hdf5;

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.exceptions.HDF5LibraryException;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Represents the HDF5 IO library.
 *
 * <p>
 *     This is a singleton class whose only instance acts a token for the existence of a working and initialized
 *     HDF support library.
 * </p>
 *
 * <p>
 *     Other components that require of a initialized HDF support library to function must acquire a non-null reference
 *     to this class singleton instance to assert that they can go ahead in performing HDF5 dependent functionality.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5Library {

    /**
     * Reference to the singleton instance of this class.
     */
    private static HDF5Library global;

    /**
     * Instantiates the singleton instance and with that loads and initializes the HDF5 access native library
     * provided by {@link H5}.
     *
     * @throws GATKException if the library cannot be instantiated.
     */
    private HDF5Library() {
        final int code;
        try {
            code = H5.H5open();
        } catch (HDF5LibraryException e) {
            throw new GATKException("could not instantiate the HDF5 Library, due to an exception.",e);
        } catch (UnsatisfiedLinkError e) {
            throw new GATKException(
                    String.format("could not instantiate the HDF5 Library, due to a library linking exception: %s",e.getMessage()), e);
        }
        if (code < 0) {
            throw new GATKException("could not instantiate the HDF5 library. H5open returned a negative value: " + code);
        }
    }

    /**
     * Returns the VM wide HDF5 library singleton instance.
     *
     * <p>
     *     This method we throw an exception when HDF5 is not supported in the current host.
     * </p>
     *
     * <p>
     *     You can check whether this is the case without the prospect of an exception using {@link #isSupported()} instead.
     * </p>
     *
     * @throws GATKException if HDF5 is not supported.
     */
    public static HDF5Library getLibrary() {
        if (global == null) {
            global = new HDF5Library();
        }
        return global;
    }

    /**
     * Checks whether HDF5 is supported in the current system.
     * <p>
     *     This method won't result in an exception if HDF5 is not currently supported, just return {@code false}.
     * </p>
     * <p>
     *     This method will load the corresponding HDF5 library for further use if available.
     * </p>
     *
     * @return {@code true} iff supported.
     */
    public static boolean isSupported() {

        if (global != null) {
            return true;
        } else {
            try {
                global = new HDF5Library();
                return true;
            } catch (final GATKException ex) {
                return false;
            }
        }
    }
}
