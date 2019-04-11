package org.broadinstitute.hellbender.utils.Bigquery;

// NOTE:
// File adapted from: https://github.com/google/google-api-java-client-samples/bigquery-appengine-sample/src/main/java/com/google/api/client/sample/bigquery/appengine/dashboard

// Copyright 2011 Google Inc. All Rights Reserved.

import com.google.api.client.util.Preconditions;
import com.google.api.services.bigquery.Bigquery;
import com.google.api.services.bigquery.model.*;
import com.google.common.base.Joiner;

import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

/**
 * Utility methods for beginning jobs, waiting for jobs, and instantiating Bigqueries.
 *
 * @author lparkinson@google.com (Laura Parkinson)
 */
public class BigqueryUtils {

    private static final Logger log       = Logger.getLogger(BigqueryUtils.class.getName());
    private static final String projectId = System.getProperty("com.google.api.client.sample.bigquery.appengine.dashboard.projectId");

    private final String   userId;
    private final Bigquery bigquery;
    private       Job      job;

    public BigqueryUtils(final String userId) throws IOException {
        this(userId, null);
    }

    private BigqueryUtils(final String userId, final String jobId) throws IOException {
        this.userId = userId;

        bigquery = ServiceUtils.loadBigqueryClient(userId);

        if (jobId != null) {
            job = tryToDo(() -> bigquery.jobs().get(projectId, jobId).execute());

            if (job == null) {
                throw new SampleDashboardException(HttpServletResponse.SC_INTERNAL_SERVER_ERROR,
                        "Wasn't able to get a job for jobId " + jobId);
            }
        }
    }

    public void beginQuery() throws SampleDashboardException {
        final Job queryJob = makeJob(buildExampleQuery());

        job = tryToDo(() -> bigquery.jobs().insert(projectId, queryJob).execute());

        Preconditions.checkNotNull(job);
        enqueueWaitingTask();
    }

    public boolean jobSucceeded() {
        return (job != null && job.getStatus().getErrorResult() == null);
    }

    public String getJobErrorMessage() {
        if (job != null && job.getStatus().getErrorResult() != null) {
            return job.getStatus().getErrorResult().getMessage();
        }
        return "";
    }

    public boolean jobIsDone() {
        final String status = getJobStatus();
        return (status != null && ("DONE").equalsIgnoreCase(status));
    }

    private String getJobStatus() {
        return (job != null) ? job.getStatus().getState() : null;
    }

    public List<TableFieldSchema> getSchemaFieldNames() throws SampleDashboardException {
        if (job != null) {
            final TableReference tableReference = job.getConfiguration().getQuery().getDestinationTable();

            final Table table = tryToDo(
                    () -> bigquery
                            .tables()
                            .get(tableReference.getProjectId(),tableReference.getDatasetId(),tableReference.getTableId())
                            .execute()
            );

            Preconditions.checkNotNull(table);
            Preconditions.checkNotNull(table.getSchema());
            Preconditions.checkNotNull(table.getSchema().getFields());
            return table.getSchema().getFields();
        }
        return null;
    }

    public List<TableRow> getTableData() throws SampleDashboardException {
        if (job != null) {
            final TableReference tableReference = job.getConfiguration().getQuery().getDestinationTable();

            final TableDataList tableDataList = tryToDo(
                    () -> bigquery
                            .tabledata()
                            .list(tableReference.getProjectId(),tableReference.getDatasetId(),tableReference.getTableId())
                            .execute()
            );

            Preconditions.checkNotNull(tableDataList);
            Preconditions.checkNotNull(tableDataList.getRows());
            return tableDataList.getRows();
        }
        return null;
    }

    /**
     * Constructs a task with necessary parameters and options and puts it in App Engine's default
     * task queue.
     */
    private void enqueueWaitingTask() {
//        final TaskOptions options = TaskOptions.Builder.withDefaults();
//        options.param("jobId", job.getJobReference().getJobId());
//        options.param("userId", userId);
//        options.url("/task");
//        options.countdownMillis(1000);
//        options.retryOptions(RetryOptions.Builder.withTaskRetryLimit(0));
//
//        final Queue queue = QueueFactory.getDefaultQueue();
//        queue.add(options);
    }

    private static String buildExampleQuery() {
        final String[] west = {"WA", "OR", "CA", "AK", "HI", "ID", "MT", "WY", "NV", "UT", "CO", "AZ", "NM"};
        final String[] south = {"OK", "TX", "AR", "LA", "TN", "MS", "AL", "KY", "GA", "FL", "SC", "NC", "VA",
                "WV", "MD", "DC", "DE"};
        final String[] midwest   = {"ND", "SD", "NE", "KS", "MN", "IA", "MO", "WI", "IL", "IN", "MI", "OH"};
        final String[] northeast = {"NY", "PA", "NJ", "CT", "RI", "MA", "VT", "NH", "ME"};

        final Joiner joiner = Joiner.on("', '");

        return "SELECT IF (state IN ('" + joiner.join(west) + "'), 'West', \n\t"
                + "IF (state IN ('" + joiner.join(south) + "'), 'South', \n\t" + "IF (state IN ('"
                + joiner.join(midwest) + "'), 'Midwest', \n\t" + "IF (state IN ('" + joiner.join(northeast)
                + "'), 'Northeast', 'None')))) "
                + "as region, \naverage_mother_age, \naverage_father_age, \nstate, \nyear \n"
                + "FROM (SELECT year, \n\t\tstate, \n\t\tSUM(mother_age)/COUNT(mother_age) as "
                + "average_mother_age, \n\t\tSUM(father_age)/COUNT(father_age) as average_father_age \n\t"
                + "FROM publicdata:samples.natality \n\tWHERE father_age < 99 \n\tGROUP BY year, state) \n"
                + "ORDER BY year, region;";
    }

    /**
     * Instantiates an example job and sets required fields.
     */
    private Job makeJob(final String query) {
        final JobConfigurationQuery jobconfigurationquery = new JobConfigurationQuery();

        jobconfigurationquery.setQuery(query);
        jobconfigurationquery.setCreateDisposition("CREATE_IF_NEEDED");

        final JobConfiguration jobconfiguration = new JobConfiguration();
        jobconfiguration.setQuery(jobconfigurationquery);

        final JobReference jobreference = new JobReference();
        jobreference.setProjectId(projectId);

        final Job newJob = new Job();
        newJob.setConfiguration(jobconfiguration);
        newJob.setJobReference(jobreference);

        return newJob;
    }

    /**
     * Attempts to run the given callback with a number of retries. If the callback responds with
     * SC_UNAUTHORIZED, the tokens are refreshed.
     *
     * @throws SampleDashboardException
     */
    private <T> T tryToDo(final Callable<T> callback) throws SampleDashboardException {
        final int                retries    = 3;
        int                      currentTry = 0;
        SampleDashboardException sdex       = null;
        while (currentTry < retries) {
            currentTry++;
            try {
                return callback.call();
            } catch ( final Exception ex) {
                sdex = new SampleDashboardException(ex);
                log.warning("Caught exception (" + sdex.getStatusCode() + "): " + ex);
            }
        }
        throw Preconditions.checkNotNull(sdex);
    }
}