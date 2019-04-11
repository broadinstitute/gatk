package org.broadinstitute.hellbender.utils.Bigquery;

// NOTE:
// File adapted from: https://github.com/google/google-api-java-client-samples/bigquery-appengine-sample/src/main/java/com/google/api/client/sample/bigquery/appengine/dashboard


// Copyright 2011 Google Inc. All Rights Reserved.

import org.apache.http.client.HttpResponseException;
import javax.servlet.http.HttpServletResponse;

/**
 * Exception to wrap an arbitrary exception as a HttpResponseException.
 *
 * @author lparkinson@google.com (Laura Parkinson)
 */
class SampleDashboardException extends HttpResponseException {

    private static final long serialVersionUID = 1L;

    SampleDashboardException(final int statusCode, final String s) {
        super(statusCode, s);
    }

    SampleDashboardException(final Exception ex) {
        super(getStatusFromException(ex), getMessageFromException(ex));
    }

    private static String getMessageFromException(final Exception ex) {
        if (ex instanceof com.google.api.client.http.HttpResponseException) {
            final com.google.api.client.http.HttpResponseException hrex =
                    (com.google.api.client.http.HttpResponseException) ex;
            return "The server encountered an exception: " + hrex.getStatusMessage();
        }
        return "The server encountered an exception: " + ex.getMessage();
    }

    private static int getStatusFromException(final Exception ex) {
        if (ex instanceof com.google.api.client.http.HttpResponseException) {
            final com.google.api.client.http.HttpResponseException hrex =
                    (com.google.api.client.http.HttpResponseException) ex;
            return hrex.getStatusCode();
        }
        return HttpServletResponse.SC_INTERNAL_SERVER_ERROR;
    }
}
