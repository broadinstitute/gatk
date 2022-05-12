package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.json.JSONArray;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class OnRampBase extends RampBase {

    final private ZipFile inputZip;
    private Map<String,Integer> regionIndex = new LinkedHashMap<>();

    public OnRampBase(final String filename) throws IOException {
        super(filename, Type.OnRamp);

        // open zip for reading
        inputZip = new ZipFile(this.file);

        // read info object
        final InputStream is = getEntry(null,"info.json");
        info = new JSONObject(new JSONTokener(is));

        // build region index
        final JSONArray regions = info.getJSONArray("regions");
        final int       regionCount = regions.length();
        for ( int i = 0 ; i < regionCount ; i++ )
            regionIndex.put(regionKey(regions.getJSONObject(i)), i);
    }

    public void close() throws IOException {

        // close file
        inputZip.close();

        super.close();
    }

    protected InputStream getEntry(final Locatable loc, final String nameParam) throws IOException {

        final String      prefix = loc != null ? getLocFilenameSuffix(loc) + "/" : "";
        final String      name = prefix + nameParam;

        Enumeration<? extends ZipEntry> entries = inputZip.entries();
        while ( entries.hasMoreElements() ) {
            ZipEntry entry = entries.nextElement();
            if ( entry.getName().equals(name) ) {
                return inputZip.getInputStream(entry);
            }
        }

        // if here, not found
        throw new IOException("no such: " + name);
    }

    public boolean hasRegion(final Locatable loc) {

        final JSONArray       regions = info.getJSONArray("regions");
        final String          locKey = loc.toString();

        return regionIndex.containsKey(locKey);
    }

    public JSONObject getRegion(final Locatable loc, final String name) {

        // get top level object
        final JSONArray       regions = info.getJSONArray("regions");
        final String          locKey = loc.toString();
        if ( !regionIndex.containsKey(locKey) )
            return null;
        JSONObject      obj = regions.getJSONObject(regionIndex.get(locKey));
        if ( obj == null || name == null ) {
            return obj;
        }

        // access nested object
        for ( String key : name.split("\\.") ) {
            if ( !obj.has(key) )
                return null;
            obj = obj.getJSONObject(key);
            if ( obj == null ) {
                return obj;
            }
        }

        return obj;
    }

    public List<Haplotype> getHaplotypes(final Locatable loc, final String name) throws IOException {

        // open csv file
        final CSVReader           reader = new CSVReader(getEntry(loc, name));
        final int                 contigColumn = reader.getColumn("contig");
        final int                 startColumn = reader.getColumn("start");
        final int                 endColumn = reader.getColumn("end");
        final int                 refColumn = reader.getColumn("ref");
        final int                 cigarColumn = reader.getColumn("cigar");
        final int                 basesColumn = reader.getColumn("bases");
        final int                 scoreColumn = reader.getColumn("score");
        final int                 alignmentStartHapwrtRefColumn = reader.getColumn("alignmentStartHapwrtRef");

        // read haplotypes
        final List<Haplotype>     haplotypes = new LinkedList<>();
        while ( reader.next() ) {

            final Haplotype       haplotype = new Haplotype(reader.get(basesColumn).getBytes(), reader.getBoolean(refColumn));
            haplotype.setGenomeLocation(reader.getLocatable(contigColumn, startColumn, endColumn));
            haplotype.setScore(reader.getDouble(scoreColumn));
            haplotype.setCigar(reader.getCigar(cigarColumn));
            haplotype.setAlignmentStartHapwrtRef(reader.getInteger(alignmentStartHapwrtRefColumn));

            /**
             * what to do about contigs? can they be regenerted?
             */
            // haplotype.contigs = h.contigs;

            haplotypes.add(haplotype);

        }

        // close up
        reader.close();

        return haplotypes;
    }

    private String regionKey(final JSONObject json) {
        return String.format("%s:%d-%d", json.get("contig"), json.getInt("start"), json.getInt("end"));
    }

    protected class CSVReader {

        final private BufferedReader    reader;
        final Map<String,Integer>       columns = new LinkedHashMap<>();
        String[]                        row;

        CSVReader(final InputStream is) throws  IOException {

            // open reader
            reader = new BufferedReader(new InputStreamReader(is));

            // read column headings
            final String      toks[] = reader.readLine().split(",");
            for ( int n = 0 ; n < toks.length ; n++ ) {
                columns.put(toks[n], n);
            }
        }

        void close() throws IOException {
            reader.close();
        }

        boolean next() throws IOException {

            final String  line = reader.readLine();
            if ( line == null ) {
                return false;
            }

            row = line.split(",");
            return true;
        }

        int getColumnCount() {
            return columns.size();
        }

        int getColumn(final String name) throws IOException {

            Integer         column = columns.get(name);
            if ( column == null ) {
                throw new IOException("no such column: " + name);
            }

            return column;
        }

        String get(final int column) throws IOException {
            return row[column];
        }

        int getInteger(final int column) throws IOException {
            return Integer.parseInt(get(column));
        }

        double getDouble(final int column) throws IOException {
            return Double.parseDouble(get(column));
        }

        boolean getBoolean(final int column) throws IOException {
            return getInteger(column) != 0;
        }

        Cigar getCigar(final int column) throws IOException {
            return TextCigarCodec.decode(get(column));
        }

        Locatable getLocatable(final int contigColumn, final int startColumn, final int endColumn) throws IOException {
            return new SimpleInterval(get(contigColumn), getInteger(startColumn), getInteger(endColumn));
        }
    }

}
