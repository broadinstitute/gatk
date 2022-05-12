package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.util.*;

public class PostFilterOnRamp extends OnRampBase {

    static private final int HAPLOTYPE_START_COLUMN = 2;
    static private final boolean IGNORE_DANGLING_READS = true;

    public PostFilterOnRamp(final String filename) throws IOException {
        super(filename);
    }

    public AlleleLikelihoods<GATKRead, Haplotype> getAlleleLikelihoods(final Locatable loc, final String name, final SampleList sampleList,
                                                                       final Map<String,List<GATKRead>> readsBySample, final Map<String,List<GATKRead>> readsBySample2) throws IOException {

        // get haplotypes
        final List<Haplotype> haplotypes = getHaplotypes(loc, name + ".haplotypes");

        // collect reads by name (incl supp)
        final Map<String, GATKRead> readsByName = collectReadNames(readsBySample);
        final Map<String, GATKRead> readsByName2 = collectReadNames(readsBySample2);

        // walk the samples
        final List<List<GATKRead>>        evidenceBySampleList = new LinkedList<>();
        final double[][][]                values = new double[sampleList.numberOfSamples()][][];
        for ( int sampleIndex = 0 ; sampleIndex < sampleList.numberOfSamples() ; sampleIndex++ ) {

            // establish context
            final String                  sampleName = sampleList.getSample(sampleIndex);
            final String                  baseName = name + ".samples." + sampleName;

            // enrich reads
            final List<GATKRead>          enrichedReads = enrichReads(loc,baseName + ".reads", readsByName, readsByName2);
            evidenceBySampleList.add(enrichedReads);

            // read matrix
            final double[][]              matrix = getSampleMatrix(loc, baseName + ".matrix", haplotypes, enrichedReads);
            values[sampleIndex] = matrix;
        }

        // create returned object
        final IndexedAlleleList<Haplotype> allelList = new IndexedAlleleList<>(haplotypes);
        @SuppressWarnings({"unchecked", "rawtypes"})
        final AlleleLikelihoods<GATKRead, Haplotype>      alleleLikelihoods = AlleleLikelihoods.<GATKRead, Haplotype>createAlleleLikelihoods(
                                                        allelList, sampleList, evidenceBySampleList, null, values);
        for ( int sampleIndex = 0 ; sampleIndex < sampleList.numberOfSamples() ; sampleIndex++ ) {
            alleleLikelihoods.sampleMatrix(sampleIndex);
        }

        return alleleLikelihoods;
    }

    private Map<String, GATKRead> collectReadNames(final Map<String,List<GATKRead>> readsBySample) {

        final Map<String, GATKRead> readsByName = new LinkedHashMap<>();
        for ( List<GATKRead> readList : readsBySample.values() ) {
            for (GATKRead read : readList) {
                readsByName.put(getReadSuppName(read), read);
            }
        }

        return readsByName;
    }

    private List<GATKRead> enrichReads(final Locatable loc, final String name, final Map<String, GATKRead> readsByName, final Map<String, GATKRead> readsByName2) throws IOException {

        // open csv file
        final CSVReader           reader = new CSVReader(getEntry(loc, name));
        final int                 nameColumn = reader.getColumn("name");
        final int                 suppColumn = reader.getColumn("supp");

        // read rows, enrich
        final  List<GATKRead>      enrichedReads = new LinkedList<>();
        while ( reader.next() ) {

            // locate read
            final String          readName = getReadSuppName(reader.get(nameColumn), reader.getBoolean(suppColumn));
            GATKRead        read = readsByName.get(readName);
            if ( read == null ) {
                read = readsByName2.get(readName);
                if ( read == null ) {
                    if (IGNORE_DANGLING_READS)
                        continue;
                    else
                        throw new IOException("no such read: " + readName);
                }
            }

            // add to list
            enrichedReads.add(read);
        }

        // close up
        reader.close();

        return enrichedReads;
    }

    private double[][] getSampleMatrix(final Locatable loc, final String name, final List<Haplotype> haplotypes, List<GATKRead> reads) throws IOException {

        // open csv file
        final CSVReader           reader = new CSVReader(getEntry(loc, name));
        if ( reader.getColumnCount() != (haplotypes.size() + HAPLOTYPE_START_COLUMN) ) {
            throw new IOException(String.format("column count %d does not match number of haplotypes %d",
                    reader.getColumnCount(), haplotypes.size()));
        }

        // allocate matrix
        // for some (unknown) reason the matrix has one more row at the end with value dup to the one before it???
        final double[][]          matrix = new double[haplotypes.size()][];
        for ( int i = 0 ; i < haplotypes.size() ; i++ ) {
            matrix[i] = new double[reads.size() + 1];
        }
        int                 row = 0;
        while ( reader.next() ) {

            // establish read
            final String      readName = reader.get(0);
            if ( !readName.equals(reads.get(row).getName()) ) {
                if (IGNORE_DANGLING_READS) {
                    logger.warn("ignoring dangling read " + readName + " on " + loc);
                    continue;
                    // TODO: investigate further
                }
                else
                    throw new IOException(String.format("Matrix read %s does not match expected %s",
                            readName, reads.get(row).getName()));
            }

            // read data
            for ( int i = 0 ; i < haplotypes.size() ; i++ ) {
                matrix[i][row] = reader.getDouble(i + HAPLOTYPE_START_COLUMN);
            }

            // advance
            row++;
        }
        if ( row != reads.size() ) {
            throw new IOException(String.format("%d rows found, but there are %d reads", row, reads.size()));
        }

        // add the extra row
        for ( int i = 0 ; i < haplotypes.size() ; i++ ) {
            matrix[i][row] = matrix[i][row - 1];
        }

        // return
        return matrix;
    }
}

