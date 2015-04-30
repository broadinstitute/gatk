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
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.security.GeneralSecurityException;
import java.util.HashMap;
import java.util.Map;

public class ReferenceSource {

    public static final int REFERENCE_SHARD_SIZE = 100000;

    private final String referenceName;
    private final Genomics genomicsService;
    private final Map<String, String> referenceNameToIdTable;

    public ReferenceSource( final String referenceName, final PipelineOptions pipelineOptions ) {
        this.referenceName = referenceName;
        this.genomicsService = createGenomicsService(pipelineOptions);
        this.referenceNameToIdTable = buildReferenceNameToIdTable(getReferenceSet());
    }

    public static int getShardIDForInterval( final Locatable l ) {
        return l.getStart() / REFERENCE_SHARD_SIZE;
    }

    public ReferenceBases getReferenceBases( final SimpleInterval interval ) {
        if ( ! referenceNameToIdTable.containsKey(interval.getContig()) ) {
            throw new IllegalArgumentException("Contig " + interval.getContig() + " not in our set of reference names for this reference source");
        }

        try {
            final Genomics.References.Bases.List listRequest = genomicsService.references().bases().list(referenceNameToIdTable.get(interval.getContig()));
            listRequest.setStart((long)interval.getStart() - 1);
            listRequest.setEnd((long)interval.getEnd() - 1);

            final ListBasesResponse result = listRequest.execute();
            return new ReferenceBases(result.getSequence().getBytes(), interval);
        }
        catch ( IOException e ) {
            throw new UserException("Query to genomics service failed for reference interval " + interval, e);
        }
    }

    private Genomics createGenomicsService( final PipelineOptions pipelineOptions ) {
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

    private ReferenceSet getReferenceSet() {
        try {
            return genomicsService.referencesets().get(referenceName).execute();
        }
        catch ( IOException e ) {
            throw new UserException("Could not load reference set for reference name " + referenceName, e);
        }
    }

    private Map<String, String> buildReferenceNameToIdTable( final ReferenceSet referenceSet ) {
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

}
