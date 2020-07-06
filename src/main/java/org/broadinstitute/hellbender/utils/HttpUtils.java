package org.broadinstitute.hellbender.utils;

import org.apache.http.HttpResponse;
import org.apache.http.client.ServiceUnavailableRetryStrategy;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.impl.conn.PoolingHttpClientConnectionManager;
import org.apache.http.protocol.HttpContext;

/**
 * Utilities for using an HTTP server
 */
public class HttpUtils {
    private static CloseableHttpClient client = HttpUtils.makeClient();

    private HttpUtils() {}

    /**
     * Returns a singleton thread safe HTTP client that can be shared
     * @return HTTP client
     */
    public static CloseableHttpClient getClient() {
        return HttpUtils.client;
    }

    private static CloseableHttpClient makeClient() {
        return HttpUtils.client = HttpClientBuilder.create()
            .setConnectionManager(new PoolingHttpClientConnectionManager())
            .setServiceUnavailableRetryStrategy(new ServiceUnavailableRetryStrategy() {
                private int interval = 1000;

                @Override
                // retry at most 4 times if a 5xx status code is received, or no status line is present
                public boolean retryRequest(final HttpResponse resp, final int executionCount, final HttpContext context) {
                    if (executionCount > 4) {
                        return false;
                    }
                    if (resp.getStatusLine() == null) {
                        return true;
                    }
                    final int statusCode = resp.getStatusLine().getStatusCode();
                    return 500 <= statusCode && statusCode < 600;
                }

                @Override
                public long getRetryInterval() {
                    final int retryInterval = interval;
                    interval *= 2;
                    return retryInterval;
                }
            })
            .build();
    }
}
