package org.broadinstitute.hellbender.utils;

import org.apache.http.HttpResponse;
import org.apache.http.client.ServiceUnavailableRetryStrategy;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.impl.conn.PoolingHttpClientConnectionManager;
import org.apache.http.protocol.HttpContext;

public class HttpUtils {
    private static CloseableHttpClient client;

    private HttpUtils() {}

    public static synchronized CloseableHttpClient getClient() {
        if (HttpUtils.client == null) {
            HttpUtils.client = HttpClientBuilder.create()
                .setConnectionManager(new PoolingHttpClientConnectionManager())
                .setServiceUnavailableRetryStrategy(new ServiceUnavailableRetryStrategy() {
                    private int interval = 1;

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
        return HttpUtils.client;
    }
}
