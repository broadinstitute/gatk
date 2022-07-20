package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * A class to hold the mappings of sample names to VCF / VCF index paths. Used by GenomicsDBImport.
 *
 * This class can be constructed from a textual file containing lines in the format:
 *
 * Sample\tVCF
 * or:
 * Sample\tVCF\tIndex
 *
 * The sample names may have internal whitespace, but not leading/trailing whitespace.
 * The VCF and Index URIs may have leading/trailing whitespace, which is ignored.
 *
 * The third Index column is optional. It is permitted to specify the index for some samples
 * and not others. If an index is not specified for a sample, its location is inferred from
 * the VCF URI.
 *
 * It is also possible to construct an empty SampleNameMap using the no-arg constructor, and
 * add sample mappings one at a time using addSample().
 */
public final class SampleNameMap {
    // Sorted mapping between sample names and corresponding GVCF file name
    //
    // IMPORTANT: This must be sorted or it will result in sample name swaps in the output database.
    // This happens because the callset json is generated independently from the import process
    // each imported batch is then sorted, so if we have an unsorted list we'll end up with different
    // global vs batch sorting.
    // We preemptively sort here so we will have consistent sorting.
    private SortedMap<String, URI> sampleNameToVcfPath;

    // Mapping between sample names and corresponding VCF index path
    //
    // This Map contains only indices specified explicitly via the sample name map file.
    // If an explicit index is not specified for a given sample, it will not have an
    // entry in this Map, and the index path will be automatically inferred based on
    // the location of the VCF.
    //
    // The ordering of the entries in this Map does not actually matter, since it's not
    // directly exposed, and is used only for individual lookups via getVCFIndexForSample()
    private SortedMap<String, URI> sampleNameToVcfIndexPath;

    /**
     * Create an empty SampleNameMap. Samples can be added later using addSample()
     */
    public SampleNameMap() {
        sampleNameToVcfPath = new TreeMap<>();
        sampleNameToVcfIndexPath = new TreeMap<>();
    }

    /**
     * Create a SampleNameMap from a textual file containing the sample mappings. The
     * lines in this file must be in the format:
     *
     * Sample\tVCF
     * or:
     * Sample\tVCF\tIndex
     *
     * The sample names may have internal whitespace, but not leading/trailing whitespace.
     * The VCF and Index URIs may have leading/trailing whitespace, which is ignored.
     *
     * The third Index column is optional. It is permitted to specify the index for some samples
     * and not others. If an index is not specified for a sample, its location is inferred from
     * the VCF URI.
     *
     * @param sampleMapFilePath Path to the file containing the sample name mappings to load
     */
    public SampleNameMap(final Path sampleMapFilePath) {
        this(sampleMapFilePath, false);
    }

    /**
     * Create a SampleNameMap from a textual file containing the sample mappings. The
     * lines in this file must be in the format:
     *
     * SampleName1\tVCF
     * or:
     * SampleName1\tVCF\tIndex
     *
     * The sample names may have internal whitespace, but not leading/trailing whitespace.
     * The VCF and Index URIs may have leading/trailing whitespace, which is ignored.
     *
     * The third Index column is optional. It is permitted to specify the index for some samples
     * and not others. If an index is not specified for a sample, its location is inferred from
     * the VCF URI.
     *
     * @param sampleMapFilePath Path to the file containing the sample name mappings to load
     * @param checkVcfIsCompressedAndIndexed If true, check each VCF to make sure it's compressed and indexed
     */
    public SampleNameMap(final Path sampleMapFilePath, final boolean checkVcfIsCompressedAndIndexed) {
        sampleNameToVcfPath = new TreeMap<>();
        sampleNameToVcfIndexPath = new TreeMap<>();

        loadSampleNameMapFile(sampleMapFilePath, checkVcfIsCompressedAndIndexed);
    }

    private void loadSampleNameMapFile(final Path sampleToFileMapPath, final boolean checkVcfIsCompressedAndIndexed) {
        try {
            final List<String> lines = Files.readAllLines(sampleToFileMapPath);
            if (lines.isEmpty()) {
                throw new UserException.BadInput( "At least 1 sample is required but none were found in the sample mapping file");
            }

            for (final String line : lines) {
                final String[] split = line.split("\\t",-1);
                if (split.length != 2 && split.length != 3) {
                    throw new UserException.BadInput("Sample name map file must have 2 or 3 fields per line in the format:\nSample\tFile\nor:\nSample\tFile\tIndex\nbut found line: \""
                            + line +"\" with "+split.length+" fields");
                }
                if ( ! sampleNameIsLegal(split[0]) || split[1].trim().isEmpty()) {
                    throw new UserException.BadInput("Sample name map file must have lines in the format:\nSample\tFile\nor:\nSample\tFile\tIndex\n but found line: '" + line + "'\nValid sample names must be non-empty strings that cannot begin or end with whitespace and valid file names must be non-empty and not all whitespace");
                }
                final String sample = split[0];
                final String vcfPath = split[1].trim();

                String vcfIndexPath = null;
                if ( split.length == 3 ) {
                    vcfIndexPath = split[2].trim();

                    if ( vcfIndexPath.isEmpty() ) {
                        throw new UserException.BadInput("Found a line in the sample name map file with an empty or all-whitespace value for the index:\n" + "\"" + line + "\"");
                    }
                }

                try {
                    final URI existingVCFPath = sampleNameToVcfPath.put(sample, new URI(vcfPath));
                    if (existingVCFPath != null){
                        throw new UserException.BadInput("Found two mappings for the same sample: " + sample + "\n" + vcfPath + "\n" + existingVCFPath);
                    }

                    if ( vcfIndexPath != null ) {
                        final URI existingVCFIndexPath = sampleNameToVcfIndexPath.put(sample, new URI(vcfIndexPath));
                        if (existingVCFIndexPath != null) {
                            throw new UserException.BadInput("Found two indices for the same sample: " + sample + "\n" + vcfIndexPath + "\n" + existingVCFIndexPath);
                        }
                    }

                    if (checkVcfIsCompressedAndIndexed) {
                        GATKGenomicsDBUtils.assertVariantFileIsCompressedAndIndexed(IOUtils.getPath(vcfPath), vcfIndexPath == null ? null : IOUtils.getPath(vcfIndexPath));
                    }
                }
                catch(final URISyntaxException e) {
                    throw new UserException("Malformed URI: " + e.toString());
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(sampleToFileMapPath, "exception while reading sample->filename mapping file",  e);
        }
    }

    /**
     * Tests whether the sample name is legal. Sample names must be non-empty, and
     * may have internal whitespace but not leading/trailing whitespace.
     *
     * @param sampleName sample name to test
     * @return true if sampleName is legal, otherwise false
     */
    private boolean sampleNameIsLegal(final String sampleName) {
        return sampleName != null &&
                ! sampleName.trim().isEmpty() &&
                sampleName.trim().equals(sampleName);
    }

    /**
     * Add a new sample mapping
     *
     * @param sampleName name of the sample
     * @param vcfPath path to the VCF for the sample
     */
    public void addSample(final String sampleName, final URI vcfPath) {
        addSample(sampleName, vcfPath, null);
    }

    /**
     * Add a new sample mapping
     *
     * @param sampleName name of the sample
     * @param vcfPath path to the VCF for the sample (not null)
     * @param vcfIndexPath path to the index for the sample (may be null)
     */
    public void addSample(final String sampleName, final URI vcfPath, final URI vcfIndexPath) {
        if ( ! sampleNameIsLegal(sampleName) ) {
            throw new UserException.BadInput("Sample name " + sampleName + " is not legal. Sample names must be non-empty and not contain leading or trailing whitespace");
        }
        if ( vcfPath == null ) {
            throw new UserException.BadInput("VCF path for sample " + sampleName + " was null");
        }

        final URI previousPath = sampleNameToVcfPath.put(sampleName, vcfPath);
        if (previousPath != null) {
            throw new UserException.BadInput("Duplicate sample: " + sampleName + ". Sample was found in both "
                    + vcfPath + " and " + previousPath + ".");
        }

        if (vcfIndexPath != null) {
            final URI previousIndexPath = sampleNameToVcfIndexPath.put(sampleName, vcfIndexPath);
            if (previousIndexPath != null) {
                throw new UserException.BadInput("For sample " + sampleName + ", attempted to specify multiple indices: " + vcfIndexPath + " and " + previousIndexPath);
            }
        }
    }

    /**
     * @return The full mapping of sample names -> VCF paths, with the sample names in sorted order
     */
    public SortedMap<String, URI> getSampleNameToVcfPath() {
        return sampleNameToVcfPath;
    }

    /**
     * @param sample sample name
     * @return the VCF associated with that sample name, as a URI
     */
    public URI getVCFForSample(final String sample) {
        return sampleNameToVcfPath.get(sample);
    }

    /**
     * @param sample sample name
     * @return the VCF associated with that sample name, as a Path
     */
    public Path getVCFForSampleAsPath(final String sample) {
        final URI vcfURI = sampleNameToVcfPath.get(sample);
        return vcfURI == null ? null : IOUtils.getPath(vcfURI.toString());
    }

    /**
     * @param sample sample name
     * @return the VCF index associated with that sample name, as a URI, or null if no index
     */
    public URI getVCFIndexForSample(final String sample) {
        return sampleNameToVcfIndexPath.get(sample);
    }

    /**
     * @param sample sample name
     * @return the VCF index associated with that sample name, as a Path, or null if no index
     */
    public Path getVCFIndexForSampleAsPath(final String sample) {
        final URI vcfIndexURI = sampleNameToVcfIndexPath.get(sample);
        return vcfIndexURI == null ? null : IOUtils.getPath(vcfIndexURI.toString());
    }

    /**
     * @return number of samples in this Map
     */
    public int getNumSamples() {
        return sampleNameToVcfPath.size();
    }

    /**
     * @return a List of the sample names in this Map in sorted order
     */
    public List<String> getSampleNamesInSortedOrder() {
        return new ArrayList<>(sampleNameToVcfPath.keySet());
    }

    /**
     * @return true if an index was specified for at least one sample, otherwise false
     */
    public boolean indicesSpecified() {
        return ! sampleNameToVcfIndexPath.isEmpty();
    }
}
