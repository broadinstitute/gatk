package org.broadinstitute.hellbender.utils.bigquery;

// NOTE:
// File adapted from: https://github.com/google/google-api-java-client-samples/bigquery-appengine-sample/src/main/java/com/google/api/client/sample/bigquery/appengine/dashboard

/*
 * Copyright (c) 2012 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */

import com.google.api.client.auth.oauth2.Credential;
import com.google.api.client.extensions.appengine.datastore.AppEngineDataStoreFactory;
import com.google.api.client.extensions.appengine.http.UrlFetchTransport;
import com.google.api.client.googleapis.auth.oauth2.GoogleAuthorizationCodeFlow;
import com.google.api.client.googleapis.auth.oauth2.GoogleClientSecrets;
import com.google.api.client.http.GenericUrl;
import com.google.api.client.http.HttpTransport;
import com.google.api.client.json.JsonFactory;
import com.google.api.client.json.jackson2.JacksonFactory;
import com.google.api.client.util.Preconditions;
import com.google.api.client.util.store.DataStoreFactory;
import com.google.api.services.bigquery.Bigquery;
import com.google.api.services.bigquery.BigqueryScopes;

import javax.servlet.http.HttpServletRequest;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;

/**
 * Utility class for Google service related tasks, for example JDO persistence, OAuth flow helpers,
 * and others.
 *
 * @author Matthias Linder (mlinder)
 */
class ServiceUtils {

    /** Global instance of the HTTP transport. */
    private static final HttpTransport HTTP_TRANSPORT = new UrlFetchTransport();

    /** Global instance of the JSON factory. */
    private static final JsonFactory JSON_FACTORY = JacksonFactory.getDefaultInstance();

    /**
     * Global instance of the {@link DataStoreFactory}. The best practice is to make it a single
     * globally shared instance across your application.
     */
    private static final AppEngineDataStoreFactory DATA_STORE_FACTORY =
            AppEngineDataStoreFactory.getDefaultInstance();

    private static GoogleClientSecrets clientSecrets = null;

    private static GoogleClientSecrets getClientCredential() throws IOException {
        if (clientSecrets == null) {
            clientSecrets = GoogleClientSecrets.load(JSON_FACTORY,
                    new InputStreamReader(ServiceUtils.class.getResourceAsStream("/client_secrets.json")));
            Preconditions.checkArgument(!clientSecrets.getDetails().getClientId().startsWith("Enter ")
                            && !clientSecrets.getDetails().getClientSecret().startsWith("Enter "),
                    "Enter Client ID and Secret from https://code.google.com/apis/console/?api=bigquery "
                            + "into bigquery-appengine-sample/src/main/resources/client_secrets.json");
        }
        return clientSecrets;
    }

    static String getRedirectUri(final HttpServletRequest req) {
        final GenericUrl url = new GenericUrl(req.getRequestURL().toString());
        url.setRawPath("/oauth2callback");
        return url.build();
    }

    static void deleteCredentials(final String userId) throws IOException {
        final GoogleAuthorizationCodeFlow flow       = newFlow();
        final Credential                  credential = flow.loadCredential(userId);
        if (credential != null) {
            flow.getCredentialDataStore().delete(userId);
        }
    }

    private static GoogleAuthorizationCodeFlow newFlow() throws IOException {
        return new GoogleAuthorizationCodeFlow.Builder(HTTP_TRANSPORT, JSON_FACTORY,
                getClientCredential(), Collections.singleton(BigqueryScopes.BIGQUERY)).setDataStoreFactory(
                DATA_STORE_FACTORY).setAccessType("offline").build();
    }

    static Bigquery loadBigqueryClient(final String userId) throws IOException {
        final Credential credential = newFlow().loadCredential(userId);
        return new Bigquery.Builder(HTTP_TRANSPORT, JSON_FACTORY, credential).build();
    }

    private ServiceUtils() {
    }
}
