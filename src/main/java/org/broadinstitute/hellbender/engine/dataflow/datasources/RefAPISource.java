package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.ListBasesResponse;
import com.google.api.services.genomics.model.Reference;
import com.google.api.services.genomics.model.ReferenceSet;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.GenomicsFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.dev.BunnyLog;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * RefAPISource makes calls to the Google Genomics API to get ReferenceBases. This is designed so it also works
 * at workers.
 * It also has helper methods for parsing Ref API URLs (intended for the Hellbender command line).
 */
public class RefAPISource implements Serializable {
    private static final long serialVersionUID = 1L;
    private static final int pageSize = 1000000; // The number results per request (the default is 200k).
    public static final String URL_PREFIX = "gg://reference/";
    private final static Logger logger = LogManager.getLogger(RefAPISource.class);

    // With our current design, there are three steps required (1) Create a map from reference name to Id,
    // (2) instantiate RefAPISource, and (3) call getReferenceBases, which then calls the Google Genomics API.
    // A one-off analysis showed that (1) takes ~seconds (2) ~1ms and (3) ~100ms for a small read.
    // So, creating the table is expensive and instantiating RefAPISource is cheap relative to calling getReferenceBases.
    //
    // genomicsService is created on each worker.
    private transient Genomics genomicsService;

    // We use the singleton pattern because we need only once instance.
    private static RefAPISource singletonSource = null;

    /**
     * Sets the singleton source for the reference API. This is intended for testing use only. This will not work
     * using the Dataflow or Spark runners as the static values will be cleared when the class is reconstituted on
     * each worker.
     * @param source the (mock) Ref API source to use.
     */
    public static void setInstance(RefAPISource source) {
        singletonSource = source;
    }

    public static RefAPISource getInstance() {
        if (singletonSource == null) {
            singletonSource = new RefAPISource();
        }
        return singletonSource;
    }
    private RefAPISource() { /* To prevent instantiation */ }

    public static boolean isApiSourceUrl(String url) {
        return url.startsWith(URL_PREFIX);
    }

    public static String getReferenceSetID(String url) {
        if (!isApiSourceUrl(url)) throw new IllegalArgumentException("Not a reference API source URL: "+url);
        return url.substring(URL_PREFIX.length());
    }

    /**
     * buildReferenceNameToIdTable produces a table from reference name (sadly not standard names like GRCh37),
     * to IDs need for Google Genomics API calls. The table is produced via a query to get the mapping. The query
     * currently takes ~10 seconds, so this table should be cached and no produced per query to get the reference bases.
     * This is clumsy and should be refactored with better Genomics API use (issue 643).
     * @param pipelineOptions options needed for credentials to produce a Genomics class for the API call.
     * @param referenceSetID the ID of the reference set, e.g., EOSt9JOVhp3jkwE.
     * @return returns a mapping from reference name to String ID.
     */
    public Map<String, String> buildReferenceNameToIdTable(final PipelineOptions pipelineOptions, final String referenceSetID) {
        Map<String, Reference> referenceMap = buildReferenceNameToReferenceTable(pipelineOptions, referenceSetID);
        return buildReferenceNameToIdTableFromMap(referenceMap);
    }

    /**
     * Given a reference set ID, get a SAMSequenceDictionary that lists its contigs and their length.
     * Note that this does NOT use the referenceContig order from the API. Instead, we'll sort them
     * and try to match the read sequence dictionary order, if given.
     *
     * @param pipelineOptions - are used to get the credentials necessary to call the Genomics API
     * @param referenceSetID - the ID of the reference set to use
     * @param optReadSequenceDictionaryToMatch - (optional) the sequence dictionary of the reads, we'll match its order if possible.
     * @return a SAMSequenceDictionary that lists the referenceset's contigs and their length.
     */
    public SAMSequenceDictionary getReferenceSequenceDictionary(final PipelineOptions pipelineOptions, String referenceSetID, final SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        Map<String, Reference> referenceMap = buildReferenceNameToReferenceTable(pipelineOptions, referenceSetID);
        return getReferenceSequenceDictionaryFromMap(referenceMap, optReadSequenceDictionaryToMatch);
    }

    /**
     * Given a reference set ID, get a SAMSequenceDictionary that lists its contigs and their length.
     * Note that this does NOT use the referenceContig order from the API. Instead, we'll sort them
     * and try to match the read sequence dictionary order, if given.
     *
     * Also produces a table from reference name (sadly not standard names like GRCh37),
     * to IDs need for Google Genomics API calls. The table is produced via a query to get the mapping. The query
     * currently takes ~10 seconds, so this table should be cached and no produced per query to get the reference bases.
     * This is clumsy and should be refactored with better Genomics API use (issue 643).
     * @param pipelineOptions options needed for credentials to produce a Genomics class for the API call.
     * @param referenceSetID the ID of the reference set, e.g., EOSt9JOVhp3jkwE.
     * @param optReadSequenceDictionaryToMatch - (optional) the sequence dictionary of the reads, we'll match its order if possible.
     * @return returns KV(a SAMSequenceDictionary that lists the referenceset's contigs and their length, a mapping from reference name to String ID)
     */
    public KV<SAMSequenceDictionary, Map<String, String>> getReferenceSequenceDictionaryAndBuildReferenceNameToIdTable(final PipelineOptions pipelineOptions, String referenceSetID, SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        Map<String, Reference> referenceMap = buildReferenceNameToReferenceTable(pipelineOptions, referenceSetID);
        final SAMSequenceDictionary referenceSequenceDictionary = getReferenceSequenceDictionaryFromMap(referenceMap, optReadSequenceDictionaryToMatch);
        final Map<String, String> referenceNameToIdTable = buildReferenceNameToIdTableFromMap(referenceMap);
        return KV.of(referenceSequenceDictionary, referenceNameToIdTable);
    }


    /**
     * Query the Google Genomics API for reference bases spanning the specified interval from the specified
     * reference name.
     *
     * Footnote: queries larger than pageSize will be truncated.
     *
     * @param pipelineOptions -- are used to get the credentials necessary to call the Genomics API
     * @param apiData - contains the hashmap that maps from reference name to Id needed for the API call.
     * @param interval - the range of bases to retrieve.
     * @return the reference bases specified by interval and apiData (using the Google Genomics API).
     */
    public ReferenceBases getReferenceBases(final PipelineOptions pipelineOptions, final RefAPIMetadata apiData, final SimpleInterval interval) {
        Utils.nonNull(pipelineOptions);
        Utils.nonNull(apiData);
        Utils.nonNull(interval);

        if (genomicsService == null) {
            genomicsService = createGenomicsService(pipelineOptions);
        }
        if ( !apiData.getReferenceNameToIdTable().containsKey(interval.getContig()) ) {
            throw new UserException("Contig " + interval.getContig() + " not in our set of reference names for this reference source");
        }

        try {
            final Genomics.References.Bases.List listRequest = genomicsService.references().bases().list(apiData.getReferenceNameToIdTable().get(interval.getContig())).setPageSize(pageSize);
            // We're subtracting 1 with the start but not the end because GA4GH is zero-based (inclusive,exclusive)
            // for its intervals.
            listRequest.setStart((long) interval.getGA4GHStart());
            listRequest.setEnd((long) interval.getGA4GHEnd());

            final ListBasesResponse result = listRequest.execute();
            if ( result.getSequence() == null ) {
                throw new UserException("No reference bases returned in query for interval " + interval + ". Is the interval valid for this reference?");
            }
            byte[] bases = result.getSequence().getBytes();
            if (bases.length != interval.size()) {
                throw new GATKException("Did not get all bases for query, is the query longer than pageSize?");
            }

            return new ReferenceBases(bases, interval);
        }
        catch ( IOException e ) {
            throw new UserException("Query to genomics service failed for reference interval " + interval, e);
        }
    }

    // --------------------------------------------------------------------------------------------
    // non-public methods

    private void fillGenomicsService(final PipelineOptions pipelineOptions) {
        if (null==genomicsService) genomicsService = createGenomicsService(pipelineOptions);
    }

    private Genomics createGenomicsService(final PipelineOptions pipelineOptions) {
        try {
            final GenomicsFactory.OfflineAuth auth = GenomicsOptions.Methods.getGenomicsAuth(pipelineOptions.as(GCSOptions.class));
            return auth.getGenomics(auth.getDefaultFactory());
        }
        catch ( GeneralSecurityException e ) {
            throw new UserException("Authentication failed for Google genomics service", e);
        }
        catch ( IOException e ) {
            throw new UserException("Unable to access Google genomics service", e);
        }
    }

    private Map<String, String> buildReferenceNameToIdTableFromMap(final Map<String, Reference> referenceMap) {
        Map<String, String> ret = new HashMap<>();
        for (Map.Entry<String,Reference> e : referenceMap.entrySet()) {
            ret.put(e.getKey(), e.getValue().getId());
        }
        return ret;
    }

    private SAMSequenceDictionary getReferenceSequenceDictionaryFromMap(final Map<String, Reference> referenceMap, final SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        SAMSequenceDictionary refDictionary = new SAMSequenceDictionary();
        ArrayList<SAMSequenceRecord> refContigs = new ArrayList<>();

        for (Map.Entry<String, Reference> e : referenceMap.entrySet()) {
            if (e.getKey()!=null && e.getValue().getLength()!=null) {
                refContigs.add(new SAMSequenceRecord(e.getKey(), e.getValue().getLength().intValue()));
            }
        }

        HashMap<String,Integer> indexBuilder = null;
        if (null!=optReadSequenceDictionaryToMatch) {
            indexBuilder = new HashMap<>();
            for (int i=0; i<optReadSequenceDictionaryToMatch.size(); i++) {
                final SAMSequenceRecord sequence = optReadSequenceDictionaryToMatch.getSequence(i);
                indexBuilder.put(sequence.getSequenceName(), i);
            }
        }
        final HashMap<String,Integer> optReadSequenceDictionaryToMatchIndex = indexBuilder;

        // GATK requires human contigs in karyotypic order: 1, 2, ..., 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs.
        // So we sort them.
        Collections.sort(refContigs, new Comparator<SAMSequenceRecord>() {
            @Override
            public int compare(SAMSequenceRecord o1, SAMSequenceRecord o2) {

                // if those are ordered in the readDictionary, then match that order
                if (null!=optReadSequenceDictionaryToMatchIndex) {
                    if (optReadSequenceDictionaryToMatchIndex.containsKey(o1.getSequenceName()) && optReadSequenceDictionaryToMatchIndex.containsKey(o2.getSequenceName())) {
                        return optReadSequenceDictionaryToMatchIndex.get(o1.getSequenceName()).compareTo(optReadSequenceDictionaryToMatchIndex.get(o2.getSequenceName()));
                    }
                }

                // otherwise, order them in karyotypic order.
                int r1 = getRank(o1.getSequenceName());
                int r2 = getRank(o2.getSequenceName());
                if (r1<r2) return -1;
                if (r2<r1) return 1;
                return o1.getSequenceName().compareTo(o2.getSequenceName());
            }

            private int getRank(String name) {
                if (name.equalsIgnoreCase("x")) return 23;
                if (name.equalsIgnoreCase("y")) return 24;
                if (name.equalsIgnoreCase("m")) return 25;
                if (name.equalsIgnoreCase("mt")) return 25;
                StringBuilder b = new StringBuilder();
                for (char c : name.toCharArray()) {
                    if (c>='0' && c<='9') {
                        b.append(c);
                    }
                }
                String numsOnly = b.toString();
                if (numsOnly.length()==0) return 0;
                return Integer.parseInt(numsOnly);
            }
        });

        for (SAMSequenceRecord s : refContigs) {
            refDictionary.addSequence(s);
        }
        return refDictionary;
    }

    // fetch the data from the Google Genomics API
    private Map<String, Reference> buildReferenceNameToReferenceTable(final PipelineOptions pipelineOptions, final String referenceSetID) {
        Utils.nonNull(pipelineOptions);
        Utils.nonNull(referenceSetID);
        fillGenomicsService(pipelineOptions);
        final ReferenceSet referenceSet;
        try {
            referenceSet = genomicsService.referencesets().get(referenceSetID).execute();
        }
        catch ( IOException e ) {
            throw new UserException("Could not load reference set for reference set ID " + referenceSetID, e);
        }
        final Map<String, Reference> ret = new HashMap<>();

        try {
            for ( final String referenceID : referenceSet.getReferenceIds() ) {
                final Reference reference = genomicsService.references().get(referenceID).execute();
                ret.put(reference.getName(), reference);
            }
        }

        catch ( IOException e ) {
            throw new UserException("Error while looking up reference contigs for reference set ID " + referenceSet.getId(), e);
        }

        return ret;
    }

}
