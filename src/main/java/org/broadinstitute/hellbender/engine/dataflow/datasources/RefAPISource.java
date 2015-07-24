package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.ListBasesResponse;
import com.google.api.services.genomics.model.Reference;
import com.google.api.services.genomics.model.ReferenceSet;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.GenomicsFactory;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;
import java.security.GeneralSecurityException;
import java.util.HashMap;
import java.util.Map;


/**
 * RefAPISource makes calls to the Google Genomics API to get ReferenceBases. This is designed so it also works
 * at workers.
 */
public class RefAPISource implements Serializable {
    private static final long serialVersionUID = 1L;
    private static final int pageSize = 1000000; // The number results per request (the default is 200k).

    // With our current design, there are three steps required (1) Create a map from reference name to Id,
    // (2) instantiate RefAPISource, and (3) call getReferenceBases, which then calls the Google Genomics API.
    // A one-off analysis showed that (1) takes ~seconds (2) ~1ms and (3) ~100ms for a small read.
    // So, creating the table is expensive and instantiating RefAPISource is cheap relative to calling getReferenceBases.
    // We've benchmarked this and found that class instantiation is much cheaper.
    //
    // genomicsService is created on each worker.
    private transient Genomics genomicsService;

    // We use the singleton pattern here because Dataflow probably isn't smart enough to prevent multiple instantiations
    // of this object (when it's not needed).
    private static RefAPISource singletonSource = null;

    /**
     * Sets the singleton source for the reference API. This is intended for testing use only. This will not work
     * using the Dataflow or Spark runners as the static values will be cleared when the class is reconstituted on
     * each worker.
     * @param source the (mock) Ref API source to use.
     */
    public static void setRefAPISource(RefAPISource source) {
        singletonSource = source;
    }

    public static RefAPISource getRefAPISource() {
        if (singletonSource == null) {
            singletonSource = new RefAPISource();
        }
        return singletonSource;
    }
    private RefAPISource() { /* To prevent instantiation */ }

    /**
     * buildReferenceNameToIdTable produces a table from reference name (sadly not standard names like GRCh37),
     * to Ids need for Google Genomics API calls. The table is produced via a query to get the mapping. The query
     * currently takes ~10 seconds, so this table should be cached and no produced per query to get the reference bases.
     * This is clumsy and should be refactored with better Genomics API use (issue 643).
     * @param pipelineOptions options needed for credentials to produce a Genomics class for the API call.
     * @param referenceName the "name of the reference, e.g., EOSt9JOVhp3jkwE.
     * @return returns a mapping from reference name to String id.
     */
    public static Map<String, String> buildReferenceNameToIdTable(final PipelineOptions pipelineOptions, final String referenceName) {
        Genomics genomicsService = createGenomicsService(pipelineOptions);
        final ReferenceSet referenceSet;
        try {
            referenceSet = genomicsService.referencesets().get(referenceName).execute();
        }
        catch ( IOException e ) {
            throw new UserException("Could not load reference set for reference name " + referenceName, e);
        }
        final Map<String, String> referenceNameToIdTable = new HashMap<>();

        try {
            for ( final String referenceID : referenceSet.getReferenceIds() ) {
                final Reference reference = genomicsService.references().get(referenceID).execute();
                referenceNameToIdTable.put(reference.getName(), reference.getId());
            }
        }

        catch ( IOException e ) {
            throw new UserException("Error while looking up references for reference set " + referenceSet.getId(), e);
        }

        return referenceNameToIdTable;
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

    private static Genomics createGenomicsService( final PipelineOptions pipelineOptions ) {
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
}
