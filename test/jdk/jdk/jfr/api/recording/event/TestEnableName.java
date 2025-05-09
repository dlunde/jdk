/*
 * Copyright (c) 2018, 2025, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.
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

package jdk.jfr.api.recording.event;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Iterator;

import jdk.jfr.Recording;
import jdk.jfr.consumer.RecordedEvent;
import jdk.test.lib.Asserts;
import jdk.test.lib.jfr.EventNames;
import jdk.test.lib.jfr.Events;

/**
 * @test
 * @summary Simple enable Event class.
 * @requires vm.flagless
 * @requires vm.hasJFR
 * @library /test/lib
 * @run main/othervm jdk.jfr.api.recording.event.TestEnableName
 */
public class TestEnableName {

    public static void main(String[] args) throws Throwable {
        Recording r = new Recording();
        final String eventType = EventNames.FileWrite;
        r.enable(eventType).withoutStackTrace();
        r.start();
        Files.write(Paths.get(".", "dummy.txt"), "DummyFile".getBytes());
        r.stop();

        Iterator<RecordedEvent> iterator = Events.fromRecording(r).iterator();
        Asserts.assertTrue(iterator.hasNext(), "No events found");
        System.out.printf("Event:%n%s%n", iterator.next());
        r.close();
    }

}
