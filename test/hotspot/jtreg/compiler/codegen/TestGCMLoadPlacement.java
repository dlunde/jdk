/*
 * Copyright (c) 2023, Oracle and/or its affiliates. All rights reserved.
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

package compiler.codegen;

/**
 * @test
 * @bug 8333393
 * @summary Test that a load is not scheduled too late.
 * @run main/othervm -XX:CompileCommand=compileonly,*Test*::test
 *                   -XX:-TieredCompilation -Xbatch
 *                   -XX:PerMethodTrapLimit=0
 *                   compiler.codegen.TestGCMLoadPlacement
 * @run main compiler.codegen.TestGCMLoadPlacement
 */

public class TestGCMLoadPlacement {
    volatile byte volFld;
    int iFld;

    int test() {
        for (int i = 0; i < 50; ++i)
            for (int j = 0; j < 50; ++j) {
                iFld = 0;
                for (int k = 0; k < 1; ++k) {

                }
        }
        int res = iFld; // This load needs to schedule before the loop below ...
        for (int i = 0; i < 50; ++i) {
            volFld = 0;
            iFld -= 42;
        }
        // ... and was previously scheduled here.
        return res;
    }

    public static void main(String[] args) {
        TestGCMLoadPlacement t = new TestGCMLoadPlacement();
        for (int i = 0; i < 10; i++) {
            int res = t.test();
            if (res != 0) {
                throw new RuntimeException("Unexpected result: " + res);
            }
        }
    }
}
