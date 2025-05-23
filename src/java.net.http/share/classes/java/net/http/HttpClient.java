/*
 * Copyright (c) 2015, 2025, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.  Oracle designates this
 * particular file as subject to the "Classpath" exception as provided
 * by Oracle in the LICENSE file that accompanied this code.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 */

package java.net.http;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.InetAddress;
import java.net.http.HttpResponse.BodyHandlers;
import java.net.http.HttpResponse.BodySubscriber;
import java.net.http.HttpResponse.BodySubscribers;
import java.nio.channels.Selector;
import java.net.Authenticator;
import java.net.CookieHandler;
import java.net.InetSocketAddress;
import java.net.Proxy;
import java.net.ProxySelector;
import java.time.Duration;
import java.util.Objects;
import java.util.Optional;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLParameters;
import java.net.http.HttpResponse.BodyHandler;
import java.net.http.HttpResponse.PushPromiseHandler;
import java.util.concurrent.Flow.Subscription;

import jdk.internal.net.http.HttpClientBuilderImpl;

/**
 * An HTTP Client.
 *
 * <p> An {@code HttpClient} can be used to send {@linkplain HttpRequest
 * requests} and retrieve their {@linkplain HttpResponse responses}. An {@code
 * HttpClient} is created through a {@link HttpClient.Builder builder}.
 * The {@link #newBuilder() newBuilder} method returns a builder that creates
 * instances of the default {@code HttpClient} implementation.
 * The builder can be used to configure per-client state, like: the preferred
 * protocol version ( HTTP/1.1 or HTTP/2 ), whether to follow redirects, a
 * proxy, an authenticator, etc. Once built, an {@code HttpClient} is immutable,
 * and can be used to send multiple requests.
 *
 * <p> An {@code HttpClient} provides configuration information, and resource
 * sharing, for all requests sent through it. An {@code HttpClient} instance
 * typically manages its own pools of connections, which it may then reuse
 * as and when necessary. Connection pools are  typically not shared between
 * {@code HttpClient} instances. Creating a new client for each operation,
 * though possible, will usually prevent reusing such connections.
 *
 * <p> A {@link BodyHandler BodyHandler} must be supplied for each {@link
 * HttpRequest} sent. The {@code BodyHandler} determines how to handle the
 * response body, if any. Once an {@link HttpResponse} is received, the
 * headers, response code, and body (typically) are available. Whether the
 * response body bytes have been read or not depends on the type, {@code T}, of
 * the response body.
 *
 * <p> Requests can be sent either synchronously or asynchronously:
 * <ul>
 *     <li>{@link HttpClient#send(HttpRequest, BodyHandler)} blocks
 *     until the request has been sent and the response has been received.</li>
 *
 *     <li>{@link HttpClient#sendAsync(HttpRequest, BodyHandler)} sends the
 *     request and receives the response asynchronously. The {@code sendAsync}
 *     method returns immediately with a {@link CompletableFuture
 *     CompletableFuture}&lt;{@link HttpResponse}&gt;. The {@code
 *     CompletableFuture} completes when the response becomes available. The
 *     returned {@code CompletableFuture} can be combined in different ways to
 *     declare dependencies among several asynchronous tasks.</li>
 * </ul>
 *
 * <p><b>Synchronous Example</b>
 * {@snippet :
 *   HttpClient client = HttpClient.newBuilder()
 *        .version(Version.HTTP_1_1)
 *        .followRedirects(Redirect.NORMAL)
 *        .connectTimeout(Duration.ofSeconds(20))
 *        .proxy(ProxySelector.of(new InetSocketAddress("proxy.example.com", 80)))
 *        .authenticator(Authenticator.getDefault())
 *        .build();
 *   HttpResponse<String> response = client.send(request, BodyHandlers.ofString());
 *   System.out.println(response.statusCode());
 *   System.out.println(response.body());  }
 *
 * <p><b>Asynchronous Example</b>
 * {@snippet :
 *    HttpRequest request = HttpRequest.newBuilder()
 *        .uri(URI.create("https://foo.com/"))
 *        .timeout(Duration.ofMinutes(2))
 *        .header("Content-Type", "application/json")
 *        .POST(BodyPublishers.ofFile(Paths.get("file.json")))
 *        .build();
 *   client.sendAsync(request, BodyHandlers.ofString())
 *        .thenApply(HttpResponse::body)
 *        .thenAccept(System.out::println);  }
 *
 * @apiNote
 * Resources allocated by the {@code HttpClient} may be
 * reclaimed early by {@linkplain #close() closing} the client.
 *
 * @implNote
 *  <p id="streaming">
 *  The {@link BodyHandlers} and {@link BodySubscribers}
 *  classes provide some {@linkplain BodySubscribers##streaming-body streaming
 *  or publishing {@code BodyHandler} and {@code BodySubscriber}
 *  implementations} which allow to stream body data back to the caller.
 *  In order for the resources associated with these streams to be
 *  reclaimed, and for the HTTP request to be considered completed,
 *  a caller must eventually {@linkplain  HttpResponse#body()
 *  obtain the streaming response body} and close, cancel, or
 *  read the returned streams to exhaustion. Likewise, a custom
 *  {@link BodySubscriber} implementation should either {@linkplain
 *  Subscription#request(long) request} all data until {@link
 *  BodySubscriber#onComplete() onComplete} or {@link
 *  BodySubscriber#onError(Throwable) onError} is signalled, or eventually
 *  {@linkplain Subscription#cancel() cancel} its subscription.
 *
 * <p id="closing">
 * The JDK built-in implementation of the {@code HttpClient} overrides
 * {@link #close()}, {@link #shutdown()}, {@link #shutdownNow()},
 * {@link #awaitTermination(Duration)}, and {@link #isTerminated()} to
 * provide a best effort implementation. Failing to close, cancel, or
 * read {@link ##streaming streaming or publishing bodies} to exhaustion
 * may stop delivery of data while leaving the request open, and
 * {@linkplain #awaitTermination(Duration) stall an
 * orderly shutdown}. The {@link #shutdownNow()} method, if called, will
 * attempt to cancel any such non-completed requests, but may cause
 * abrupt termination of any on going operation.
 *
 * <p id="gc">
 * If not {@linkplain ##closing explicitly closed}, the JDK
 * built-in implementation of the {@code HttpClient} releases
 * its resources when an {@code HttpClient} instance is no longer
 * strongly reachable, and all operations started on that instance have
 * eventually completed. This relies both on the garbage collector
 * to notice that the instance is no longer reachable, and on all
 * requests started on the client to eventually complete. Failure
 * to properly close {@linkplain ##streaming streaming or publishing bodies}
 * may prevent the associated requests from running to completion, and
 * prevent the resources allocated by the associated client from
 * being reclaimed by the garbage collector.
 *
 * @since 11
 */
public abstract class HttpClient implements AutoCloseable {

    /**
     * Creates an HttpClient.
     */
    protected HttpClient() {}

    /**
     * Returns a new {@code HttpClient} with default settings.
     *
     * <p> Equivalent to {@code newBuilder().build()}.
     *
     * <p> The default settings include: the "GET" request method, a preference
     * of {@linkplain HttpClient.Version#HTTP_2 HTTP/2}, a redirection policy of
     * {@linkplain Redirect#NEVER NEVER}, the {@linkplain
     * ProxySelector#getDefault() default proxy selector}, and the {@linkplain
     * SSLContext#getDefault() default SSL context}.
     *
     * @implNote The system-wide default values are retrieved at the time the
     * {@code HttpClient} instance is constructed. Changing the system-wide
     * values after an {@code HttpClient} instance has been built, for
     * instance, by calling {@link ProxySelector#setDefault(ProxySelector)}
     * or {@link SSLContext#setDefault(SSLContext)}, has no effect on already
     * built instances.
     *
     * @return a new HttpClient
     * @throws UncheckedIOException if necessary underlying IO resources required to
     * {@linkplain Builder#build() build a new HttpClient} cannot be allocated.
     */
    public static HttpClient newHttpClient() {
        return newBuilder().build();
    }

    /**
     * Creates a new {@code HttpClient} builder.
     *
     * <p> Builders returned by this method create instances
     * of the default {@code HttpClient} implementation.
     *
     * @return an {@code HttpClient.Builder}
     */
    public static Builder newBuilder() {
        return new HttpClientBuilderImpl();
    }

    /**
     * A builder of {@linkplain HttpClient HTTP Clients}.
     *
     * <p> Builders are created by invoking {@link HttpClient#newBuilder()
     * newBuilder}. Each of the setter methods modifies the state of the builder
     * and returns the same instance. Builders are not thread-safe and should not be
     * used concurrently from multiple threads without external synchronization.
     *
     * @since 11
     */
    public interface Builder {

        /**
         * A proxy selector that always return {@link Proxy#NO_PROXY} implying
         * a direct connection.
         *
         * <p> This is a convenience object that can be passed to
         * {@link #proxy(ProxySelector)} in order to build an instance of
         * {@link HttpClient} that uses no proxy.
         */
        public static final ProxySelector NO_PROXY = ProxySelector.of(null);


        /**
         * Sets a cookie handler.
         *
         * @param cookieHandler the cookie handler
         * @return this builder
         */
        public Builder cookieHandler(CookieHandler cookieHandler);

        /**
         * Sets the connect timeout duration for this client.
         *
         * <p> In the case where a new connection needs to be established, if
         * the connection cannot be established within the given {@code
         * duration}, then {@link HttpClient#send(HttpRequest,BodyHandler)
         * HttpClient::send} throws an {@link HttpConnectTimeoutException}, or
         * {@link HttpClient#sendAsync(HttpRequest,BodyHandler)
         * HttpClient::sendAsync} completes exceptionally with an
         * {@code HttpConnectTimeoutException}. If a new connection does not
         * need to be established, for example if a connection can be reused
         * from a previous request, then this timeout duration has no effect.
         *
         * @param duration the duration to allow the underlying connection to be
         *                 established
         * @return this builder
         * @throws IllegalArgumentException if the duration is non-positive
         */
        public Builder connectTimeout(Duration duration);

        /**
         * Sets an {@code SSLContext}.
         *
         * <p> If this method is not invoked prior to {@linkplain #build()
         * building}, then newly built clients will use the {@linkplain
         * SSLContext#getDefault() default context}, which is normally adequate
         * for client applications that do not need to specify protocols, or
         * require client authentication.
         *
         * @param sslContext the SSLContext
         * @return this builder
         */
        public Builder sslContext(SSLContext sslContext);

        /**
         * Sets an {@code SSLParameters}.
         *
         * <p> If this method is not invoked prior to {@linkplain #build()
         * building}, then newly built clients will use a default,
         * implementation specific, set of parameters.
         *
         * <p> Some parameters which are used internally by the HTTP Client
         * implementation (such as the application protocol list) should not be
         * set by callers, as they may be ignored. The contents of the given
         * object are copied.
         *
         * @param sslParameters the SSLParameters
         * @return this builder
         */
        public Builder sslParameters(SSLParameters sslParameters);

        /**
         * Sets the executor to be used for asynchronous and dependent tasks.
         *
         * <p> If this method is not invoked prior to {@linkplain #build()
         * building}, a default executor is created for each newly built {@code
         * HttpClient}.
         *
         * @implNote The default executor uses a thread pool, with a custom
         * thread factory.
         *
         * @param executor the Executor
         * @return this builder
         */
        public Builder executor(Executor executor);

        /**
         * Specifies whether requests will automatically follow redirects issued
         * by the server.
         *
         * <p> If this method is not invoked prior to {@linkplain #build()
         * building}, then newly built clients will use a default redirection
         * policy of {@link Redirect#NEVER NEVER}.
         *
         * @param policy the redirection policy
         * @return this builder
         */
        public Builder followRedirects(Redirect policy);

        /**
         * Requests a specific HTTP protocol version where possible.
         *
         * <p> If this method is not invoked prior to {@linkplain #build()
         * building}, then newly built clients will prefer {@linkplain
         * Version#HTTP_2 HTTP/2}.
         *
         * <p> If set to {@linkplain Version#HTTP_2 HTTP/2}, then each request
         * will attempt to upgrade to HTTP/2. If the upgrade succeeds, then the
         * response to this request will use HTTP/2 and all subsequent requests
         * and responses to the same
         * <a href="https://tools.ietf.org/html/rfc6454#section-4">origin server</a>
         * will use HTTP/2. If the upgrade fails, then the response will be
         * handled using HTTP/1.1
         *
         * @implNote Constraints may also affect the selection of protocol version.
         * For example, if HTTP/2 is requested through a proxy, and if the implementation
         * does not support this mode, then HTTP/1.1 may be used
         *
         * @param version the requested HTTP protocol version
         * @return this builder
         */
        public Builder version(HttpClient.Version version);

        /**
         * Sets the default priority for any HTTP/2 requests sent from this
         * client. The value provided must be between {@code 1} and {@code 256}
         * (inclusive).
         *
         * @param priority the priority weighting
         * @return this builder
         * @throws IllegalArgumentException if the given priority is out of range
         */
        public Builder priority(int priority);

        /**
         * Sets a {@link java.net.ProxySelector}.
         *
         * @apiNote {@link ProxySelector#of(InetSocketAddress) ProxySelector::of}
         * provides a {@code ProxySelector} which uses a single proxy for all
         * requests. The system-wide proxy selector can be retrieved by
         * {@link ProxySelector#getDefault()}.
         *
         * @implNote
         * If this method is not invoked prior to {@linkplain #build() building},
         * then newly built clients will use the {@linkplain
         * ProxySelector#getDefault() default proxy selector}, which is usually
         * adequate for client applications. The default proxy selector supports
         * a set of system properties related to
         * <a href="{@docRoot}/java.base/java/net/doc-files/net-properties.html#Proxies">
         * proxy settings</a>. This default behavior can be disabled by
         * supplying an explicit proxy selector, such as {@link #NO_PROXY} or
         * one returned by {@link ProxySelector#of(InetSocketAddress)
         * ProxySelector::of}, before {@linkplain #build() building}.
         *
         * @param proxySelector the ProxySelector
         * @return this builder
         */
        public Builder proxy(ProxySelector proxySelector);

        /**
         * Sets an authenticator to use for HTTP authentication.
         *
         * @implNote
         * In the JDK built-in implementation of the {@code HttpClient},
         * if a {@link HttpRequest} has an {@code Authorization} or {@code
         * Proxy-Authorization} header set then its value is used and
         * the {@link Authenticator} is not invoked for the corresponding
         * authentication. In this case, any authentication errors are returned
         * to the user and requests are not automatically retried.
         * Additionally, the JDK built-in implementation currently only supports HTTP
         * {@code Basic} authentication.
         *
         * @param authenticator the Authenticator
         * @return this builder
         */
        public Builder authenticator(Authenticator authenticator);

        /**
         * Binds the socket to this local address when creating
         * connections for sending requests.
         *
         * <p> If no local address is set or {@code null} is passed
         * to this method then sockets created by the
         * HTTP client will be bound to an automatically
         * assigned socket address.
         *
         * <p> Common usages of the {@code HttpClient} do not require
         * this method to be called. Setting a local address, through this
         * method, is only for advanced usages where users of the {@code HttpClient}
         * require specific control on which network interface gets used
         * for the HTTP communication. Callers of this method are expected to
         * be aware of the networking configurations of the system where the
         * {@code HttpClient} will be used and care should be taken to ensure the
         * correct {@code localAddr} is passed. Failure to do so can result in
         * requests sent through the {@code HttpClient} to fail.
         *
         * @implSpec The default implementation of this method throws
         * {@code UnsupportedOperationException}. {@code Builder}s obtained
         * through {@link HttpClient#newBuilder()} provide an implementation
         * of this method that allows setting the local address.
         *
         * @param localAddr The local address of the socket. Can be null.
         * @return this builder
         * @throws UnsupportedOperationException if this builder doesn't support
         *         configuring a local address or if the passed {@code localAddr}
         *         is not supported by this {@code HttpClient} implementation.
         * @since 19
         */
        default Builder localAddress(InetAddress localAddr) {
            throw new UnsupportedOperationException();
        }

        /**
         * Returns a new {@link HttpClient} built from the current state of this
         * builder.
         *
         * @return a new {@code HttpClient}
         *
         * @throws UncheckedIOException may be thrown if underlying IO resources required
         * by the implementation cannot be allocated. For instance,
         * if the implementation requires a {@link Selector}, and opening
         * one fails due to {@linkplain Selector#open() lack of necessary resources}.
         */
        public HttpClient build();
    }


    /**
     * Returns an {@code Optional} containing this client's {@link
     * CookieHandler}. If no {@code CookieHandler} was set in this client's
     * builder, then the {@code Optional} is empty.
     *
     * @return an {@code Optional} containing this client's {@code CookieHandler}
     */
    public abstract Optional<CookieHandler> cookieHandler();

    /**
     * Returns an {@code Optional} containing the <i>connect timeout duration</i>
     * for this client. If the {@linkplain Builder#connectTimeout(Duration)
     * connect timeout duration} was not set in the client's builder, then the
     * {@code Optional} is empty.
     *
     * @return an {@code Optional} containing this client's connect timeout
     *         duration
     */
     public abstract Optional<Duration> connectTimeout();

    /**
     * Returns the follow redirects policy for this client. The default value
     * for client's built by builders that do not specify a redirect policy is
     * {@link HttpClient.Redirect#NEVER NEVER}.
     *
     * @return this client's follow redirects setting
     */
    public abstract Redirect followRedirects();

    /**
     * Returns an {@code Optional} containing the {@code ProxySelector}
     * supplied to this client. If no proxy selector was set in this client's
     * builder, then the {@code Optional} is empty.
     *
     * <p> Even though this method may return an empty optional, the {@code
     * HttpClient} may still have a non-exposed {@linkplain
     * Builder#proxy(ProxySelector) default proxy selector} that is
     * used for sending HTTP requests.
     *
     * @return an {@code Optional} containing the proxy selector supplied
     *        to this client.
     */
    public abstract Optional<ProxySelector> proxy();

    /**
     * Returns this client's {@code SSLContext}.
     *
     * <p> If no {@code SSLContext} was set in this client's builder, then the
     * {@linkplain SSLContext#getDefault() default context} is returned.
     *
     * @return this client's SSLContext
     */
    public abstract SSLContext sslContext();

    /**
     * Returns a copy of this client's {@link SSLParameters}.
     *
     * <p> If no {@code SSLParameters} were set in the client's builder, then an
     * implementation specific default set of parameters, that the client will
     * use, is returned.
     *
     * @return this client's {@code SSLParameters}
     */
    public abstract SSLParameters sslParameters();

    /**
     * Returns an {@code Optional} containing the {@link Authenticator} set on
     * this client. If no {@code Authenticator} was set in the client's builder,
     * then the {@code Optional} is empty.
     *
     * @return an {@code Optional} containing this client's {@code Authenticator}
     */
    public abstract Optional<Authenticator> authenticator();

    /**
     * Returns the preferred HTTP protocol version for this client. The default
     * value is {@link HttpClient.Version#HTTP_2}
     *
     * @implNote Constraints may also affect the selection of protocol version.
     * For example, if HTTP/2 is requested through a proxy, and if the
     * implementation does not support this mode, then HTTP/1.1 may be used
     *
     * @return the HTTP protocol version requested
     */
    public abstract HttpClient.Version version();

    /**
     * Returns an {@code Optional} containing this client's {@link
     * Executor}. If no {@code Executor} was set in the client's builder,
     * then the {@code Optional} is empty.
     *
     * <p> Even though this method may return an empty optional, the {@code
     * HttpClient} may still have an non-exposed {@linkplain
     * HttpClient.Builder#executor(Executor) default executor} that is used for
     * executing asynchronous and dependent tasks.
     *
     * @return an {@code Optional} containing this client's {@code Executor}
     */
    public abstract Optional<Executor> executor();

    /**
     * The HTTP protocol version.
     *
     * @since 11
     */
    public enum Version {

        /**
         * HTTP version 1.1
         */
        HTTP_1_1,

        /**
         * HTTP version 2
         */
        HTTP_2
    }

    /**
     * Defines the automatic redirection policy.
     *
     * <p> The automatic redirection policy is checked whenever a {@code 3XX}
     * response code is received. If redirection does not happen automatically,
     * then the response, containing the  {@code 3XX} response code, is returned,
     * where it can be handled manually.
     *
     * <p> {@code Redirect} policy is set through the {@linkplain
     * HttpClient.Builder#followRedirects(Redirect) Builder.followRedirects}
     * method.
     *
     * @implNote When automatic redirection occurs, the request method of the
     * redirected request may be modified depending on the specific {@code 30X}
     * status code, as specified in <a href="https://tools.ietf.org/html/rfc7231">
     * RFC 7231</a>. In addition, the {@code 301} and {@code 302} status codes
     * cause a {@code POST} request to be converted to a {@code GET} in the
     * redirected request.
     *
     * @since 11
     */
    public enum Redirect {

        /**
         * Never redirect.
         */
        NEVER,

        /**
         * Always redirect.
         */
        ALWAYS,

        /**
         * Always redirect, except from HTTPS URLs to HTTP URLs.
         */
        NORMAL
    }

    /**
     * Sends the given request using this client, blocking if necessary to get
     * the response. The returned {@link HttpResponse}{@code <T>} contains the
     * response status, headers, and body ( as handled by given response body
     * handler ).
     *
     * <p> If the operation is interrupted, the default {@code HttpClient}
     * implementation attempts to cancel the HTTP exchange and
     * {@link InterruptedException} is thrown.
     * No guarantee is made as to exactly <em>when</em> the cancellation request
     * may be taken into account. In particular, the request might still get sent
     * to the server, as its processing might already have started asynchronously
     * in another thread, and the underlying resources may only be released
     * asynchronously.
     * <ul>
     *     <li>With HTTP/1.1, an attempt to cancel may cause the underlying
     *         connection to be closed abruptly.
     *     <li>With HTTP/2, an attempt to cancel may cause the stream to be reset,
     *         or in certain circumstances, may also cause the connection to be
     *         closed abruptly, if, for instance, the thread is currently trying
     *         to write to the underlying socket.
     * </ul>
     *
     * @param <T> the response body type
     * @param request the request
     * @param responseBodyHandler the response body handler
     * @return the response
     * @throws IOException if an I/O error occurs when sending or receiving, or
     *         the client has {@linkplain ##closing shut down}
     * @throws InterruptedException if the operation is interrupted
     * @throws IllegalArgumentException if the {@code request} argument is not
     *         a request that could have been validly built as specified by {@link
     *         HttpRequest.Builder HttpRequest.Builder}.
     */
    public abstract <T> HttpResponse<T>
    send(HttpRequest request, HttpResponse.BodyHandler<T> responseBodyHandler)
        throws IOException, InterruptedException;

    /**
     * Sends the given request asynchronously using this client with the given
     * response body handler.
     *
     * <p> Equivalent to: {@code sendAsync(request, responseBodyHandler, null)}.
     *
     * @param <T> the response body type
     * @param request the request
     * @param responseBodyHandler the response body handler
     * @return a {@code CompletableFuture<HttpResponse<T>>}
     * @throws IllegalArgumentException if the {@code request} argument is not
     *         a request that could have been validly built as specified by {@link
     *         HttpRequest.Builder HttpRequest.Builder}.
     */
    public abstract <T> CompletableFuture<HttpResponse<T>>
    sendAsync(HttpRequest request,
              BodyHandler<T> responseBodyHandler);

    /**
     * Sends the given request asynchronously using this client with the given
     * response body handler and push promise handler.
     *
     * <p> The returned completable future, if completed successfully, completes
     * with an {@link HttpResponse}{@code <T>} that contains the response status,
     * headers, and body ( as handled by given response body handler ).
     *
     * <p> {@linkplain PushPromiseHandler Push promises} received, if any, are
     * handled by the given {@code pushPromiseHandler}. A {@code null} valued
     * {@code pushPromiseHandler} rejects any push promises.
     *
     * <p> The returned completable future completes exceptionally with:
     * <ul>
     * <li>{@link IOException} - if an I/O error occurs when sending or receiving,
     *      or the client has {@linkplain ##closing shut down}.</li>
     * </ul>
     *
     * <p id="cancel"> The default {@code HttpClient} implementation returns
     * {@code CompletableFuture} objects that are <em>cancelable</em>.
     * {@code CompletableFuture} objects {@linkplain CompletableFuture#newIncompleteFuture()
     * derived} from cancelable futures are themselves <em>cancelable</em>.
     * Invoking {@linkplain CompletableFuture#cancel(boolean) cancel(true)}
     * on a cancelable future that is not completed, attempts to cancel the HTTP exchange
     * in an effort to release underlying resources as soon as possible.
     * No guarantee is made as to exactly <em>when</em> the cancellation request
     * may be taken into account. In particular, the request might still get sent
     * to the server, as its processing might already have started asynchronously
     * in another thread, and the underlying resources may only be released
     * asynchronously.
     * <ul>
     *     <li>With HTTP/1.1, an attempt to cancel may cause the underlying connection
     *         to be closed abruptly.
     *     <li>With HTTP/2, an attempt to cancel may cause the stream to be reset.
     * </ul>
     *
     * @param <T> the response body type
     * @param request the request
     * @param responseBodyHandler the response body handler
     * @param pushPromiseHandler push promise handler, may be null
     * @return a {@code CompletableFuture<HttpResponse<T>>}
     * @throws IllegalArgumentException if the {@code request} argument is not
     *         a request that could have been validly built as specified by {@link
     *         HttpRequest.Builder HttpRequest.Builder}.
     */
    public abstract <T> CompletableFuture<HttpResponse<T>>
    sendAsync(HttpRequest request,
              BodyHandler<T> responseBodyHandler,
              PushPromiseHandler<T> pushPromiseHandler);

    /**
     * Creates a new {@code WebSocket} builder (optional operation).
     *
     * <p> <b>Example</b>
     * {@snippet :
     *   HttpClient client = HttpClient.newHttpClient();
     *   CompletableFuture<WebSocket> ws = client.newWebSocketBuilder()
     *      .buildAsync(URI.create("ws://websocket.example.com"), listener);  }
     *
     * <p> Finer control over the WebSocket Opening Handshake can be achieved
     * by using a custom {@code HttpClient}.
     *
     * <p> <b>Example</b>
     * {@snippet :
     *   InetSocketAddress addr = new InetSocketAddress("proxy.example.com", 80);
     *   HttpClient client = HttpClient.newBuilder()
     *           .proxy(ProxySelector.of(addr))
     *           .build();
     *
     *   CompletableFuture<WebSocket> ws = client.newWebSocketBuilder()
     *           .buildAsync(URI.create("ws://websocket.example.com"), listener);  }
     *
     * @implSpec The default implementation of this method throws
     * {@code UnsupportedOperationException}. Clients obtained through
     * {@link HttpClient#newHttpClient()} or {@link HttpClient#newBuilder()}
     * return a {@code WebSocket} builder.
     *
     * @implNote Both builder and {@code WebSocket}s created with it operate in
     * a non-blocking fashion. That is, their methods do not block before
     * returning a {@code CompletableFuture}. Asynchronous tasks are executed in
     * this {@code HttpClient}'s executor.
     *
     * <p> When a {@code CompletionStage} returned from
     * {@link WebSocket.Listener#onClose Listener.onClose} completes,
     * the {@code WebSocket} will send a Close message that has the same code
     * the received message has and an empty reason.
     *
     * @return a {@code WebSocket.Builder}
     * @throws UnsupportedOperationException
     *         if this {@code HttpClient} does not provide WebSocket support
     */
    public WebSocket.Builder newWebSocketBuilder() {
        throw new UnsupportedOperationException();
    }

    /**
     * Initiates an orderly shutdown in which  requests previously
     * submitted with {@code send} or {@code sendAsync}
     * are run to completion, but no new request will be accepted.
     * Running a request to completion may involve running several
     * operations in the background, including {@linkplain ##closing
     * waiting for responses to be delivered}, which will all have to
     * run to completion until the request is considered completed.
     *
     * Invocation has no additional effect if already shut down.
     *
     * <p>This method does not wait for previously submitted request
     * to complete execution.  Use {@link #awaitTermination(Duration)
     * awaitTermination} or {@link #close() close} to do that.
     *
     * @implSpec
     * The default implementation of this method does nothing. Subclasses should
     * override this method to implement the appropriate behavior.
     *
     * @see ##closing Implementation Note on closing the HttpClient
     *
     * @since 21
     */
    public void shutdown() { }

    /**
     * Blocks until all operations have completed execution after a shutdown
     * request, or the {@code duration} elapses, or the current thread is
     * {@linkplain Thread#interrupt() interrupted}, whichever happens first.
     * Operations are any tasks required to run a request previously
     * submitted with {@code send} or {@code sendAsync} to completion.
     *
     * <p> This method does not wait if the duration to wait is less than or
     * equal to zero. In this case, the method just tests if the thread has
     * terminated.
     *
     * @implSpec
     * The default implementation of this method checks for null arguments, but
     * otherwise does nothing and returns true.
     * Subclasses should override this method to implement the proper behavior.
     *
     * @param duration the maximum time to wait
     * @return {@code true} if this client terminated and
     *         {@code false} if the timeout elapsed before termination
     * @throws InterruptedException if interrupted while waiting
     *
     * @see ##closing Implementation Note on closing the HttpClient
     *
     * @since 21
     */
    public boolean awaitTermination(Duration duration) throws InterruptedException {
        Objects.requireNonNull(duration);
        return true;
    }

    /**
     * Returns {@code true} if all operations have completed following
     * a shutdown.
     * Operations are any tasks required to run a request previously
     * submitted with {@code send} or {@code sendAsync} to completion.
     * <p> Note that {@code isTerminated} is never {@code true} unless
     * either {@code shutdown} or {@code shutdownNow} was called first.
     *
     * @implSpec
     * The default implementation of this method does nothing and returns false.
     * Subclasses should override this method to implement the proper behavior.
     *
     * @return {@code true} if all tasks have completed following a shutdown
     *
     * @see ##closing Implementation Note on closing the HttpClient
     *
     * @since 21
     */
    public boolean isTerminated() {
        return false;
    }

    /**
     * This method attempts to initiate an immediate shutdown.
     * An implementation of this method may attempt to
     * interrupt operations that are actively running.
     * Operations are any tasks required to run a request previously
     * submitted with {@code send} or {@code sendAsync} to completion.
     * The behavior of actively running operations when interrupted
     * is undefined. In particular, there is no guarantee that
     * interrupted operations will terminate, or that code waiting
     * on these operations will ever be notified.
     *
     * @implSpec
     * The default implementation of this method simply calls {@link #shutdown()}.
     * Subclasses should override this method to implement the appropriate
     * behavior.
     *
     * @see ##closing Implementation Note on closing the HttpClient
     *
     * @since 21
     */
    public void shutdownNow() {
        shutdown();
    }

    /**
     * Initiates an orderly shutdown in which  requests previously
     * submitted to {@code send} or {@code sendAsync}
     * are run to completion, but no new request will be accepted.
     * Running a request to completion may involve running several
     * operations in the background, including {@linkplain ##closing
     * waiting for responses to be delivered}.
     * This method waits until all operations have completed execution
     * and the client has terminated.
     *
     * <p> If interrupted while waiting, this method may attempt to stop all
     * operations by calling {@link #shutdownNow()}. It then continues to wait
     * until all actively executing operations have completed.
     * The interrupt status will be re-asserted before this method returns.
     *
     * <p> If already terminated, invoking this method has no effect.
     *
     * @implSpec
     * The default implementation invokes {@code shutdown()} and waits for tasks
     * to complete execution with {@code awaitTermination}.
     *
     * @see ##closing Implementation Note on closing the HttpClient
     *
     * @since 21
     */
    @Override
    public void close() {
        boolean terminated = isTerminated();
        if (!terminated) {
            shutdown();
            boolean interrupted = false;
            while (!terminated) {
                try {
                    terminated = awaitTermination(Duration.ofDays(1L));
                } catch (InterruptedException e) {
                    if (!interrupted) {
                        interrupted = true;
                        shutdownNow();
                        if (isTerminated()) break;
                    }
                }
            }
            if (interrupted) {
                Thread.currentThread().interrupt();
            }
        }
    }

}
