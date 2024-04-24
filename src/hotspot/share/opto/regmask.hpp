/*
 * Copyright (c) 1997, 2024, Oracle and/or its affiliates. All rights reserved.
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
 *
 */

#ifndef SHARE_OPTO_REGMASK_HPP
#define SHARE_OPTO_REGMASK_HPP

#include "code/vmreg.hpp"
#include "opto/optoreg.hpp"
#include "utilities/count_leading_zeros.hpp"
#include "utilities/count_trailing_zeros.hpp"
#include "utilities/globalDefinitions.hpp"
#include "memory/resourceArea.hpp"

class LRG;
class RegMaskStatic;

//-------------Non-zero bit search methods used by RegMask---------------------
// Find lowest 1, undefined if empty/0
static unsigned int find_lowest_bit(uintptr_t mask) {
  return count_trailing_zeros(mask);
}
// Find highest 1, undefined if empty/0
static unsigned int find_highest_bit(uintptr_t mask) {
  return count_leading_zeros(mask) ^ (BitsPerWord - 1U);
}

//------------------------------RegMask----------------------------------------
// The ADL file describes how to print the machine-specific registers, as well
// as any notion of register classes.  We provide a register mask, which is
// just a collection of Register numbers.

// The ADLC defines 2 macros, RM_SIZE and FORALL_BODY.
// RM_SIZE is the size of a register mask in 32-bit words.
// FORALL_BODY replicates a BODY macro once per word in the register mask.
// The usage is somewhat clumsy and limited to the regmask.[h,c]pp files.
// However, it means the ADLC can redefine the unroll macro and all loops
// over register masks will be unrolled by the correct amount.

class RegMask {

  friend class RegMaskIterator;

  // The RM_SIZE is aligned to 64-bit - assert that this holds
  LP64_ONLY(STATIC_ASSERT(is_aligned(RM_SIZE, 2)));

  protected:

  static const unsigned int _WordBitMask = BitsPerWord - 1U;
  static const unsigned int _LogWordBits = LogBitsPerWord;
  static const unsigned int _RM_SIZE     = LP64_ONLY(RM_SIZE >> 1) NOT_LP64(RM_SIZE);
  static const unsigned int _RM_MAX      = _RM_SIZE - 1U;

  union {
    // Array of Register Mask bits.  This array is large enough to cover
    // all the machine registers and all parameters that need to be passed
    // on the stack (stack registers) up to some interesting limit.  Methods
    // that need more parameters will NOT be compiled.  On Intel, the limit
    // is something like 90+ parameters.
    int*       _RM_I;
    uintptr_t* _RM_UP;
  };

  // Size in words
  unsigned int _rm_size;

  // The maximum word index
  unsigned int _rm_max() const { return _rm_size - 1U; }

  // We can offset register masks to present different views of the register
  // space. We define the mask not to contain registers before _offset (i.e.,
  // they are zero).
  // NOTE We should only offset by full words
  unsigned int _offset_up;
  unsigned int _offset_bits() const { return _offset_up * BitsPerWord; };

  bool _all_stack = false;

  // The low and high water marks represents the lowest and highest word
  // that might contain set register mask bits, respectively. We guarantee
  // that there are no bits in words outside this range, but any word at
  // and between the two marks can still be 0.
  unsigned int _lwm;
  unsigned int _hwm;

  // Protected constructor
  RegMask(int _rm_size, int _offset_up = 0)
    : _rm_size(_rm_size), _offset_up(_offset_up),
      _lwm(_rm_max()), _hwm(0) {}

  // Common copying functionality
  static void _copy(const RegMask& src, RegMask& dst) {
    assert(src._offset_up == dst._offset_up, "");
    assert(src._rm_size == dst._rm_size, "");
    dst._hwm = src._hwm;
    dst._lwm = src._lwm;
    for (unsigned i = 0; i < dst._rm_size; i++) {
      dst._RM_UP[i] = src._RM_UP[i];
    }
    dst.set_AllStack(src.is_AllStack());
    assert(dst.valid_watermarks(), "post-condition");
  }

  public:
  enum { CHUNK_SIZE = _RM_SIZE * BitsPerWord }; // TODO Move to private and expose _rm_size instead

  unsigned int rm_size() const { return _rm_size; }
  unsigned int rm_size_bits() const { return _rm_size * BitsPerWord; }

  // SlotsPerLong is 2, since slots are 32 bits and longs are 64 bits.
  // Also, consider the maximum alignment size for a normally allocated
  // value.  Since we allocate register pairs but not register quads (at
  // present), this alignment is SlotsPerLong (== 2).  A normally
  // aligned allocated register is either a single register, or a pair
  // of adjacent registers, the lower-numbered being even.
  // See also is_aligned_Pairs() below, and the padding added before
  // Matcher::_new_SP to keep allocated pairs aligned properly.
  // If we ever go to quad-word allocations, SlotsPerQuad will become
  // the controlling alignment constraint.  Note that this alignment
  // requirement is internal to the allocator, and independent of any
  // particular platform.
  enum { SlotsPerLong = 2,
         SlotsPerVecA = 4,
         SlotsPerVecS = 1,
         SlotsPerVecD = 2,
         SlotsPerVecX = 4,
         SlotsPerVecY = 8,
         SlotsPerVecZ = 16,
         SlotsPerRegVectMask = X86_ONLY(2) NOT_X86(1)
         };

  // Pure virtual copy assignment
  virtual RegMask& operator= (const RegMask& rm) = 0;

  // Check for register being in mask
  bool Member(OptoReg::Name reg) const {
    int reg_offset = reg - _offset_bits();
    if(reg_offset < 0) { return false; }
    unsigned r = (unsigned)reg_offset;
    if (r >= _rm_size * BitsPerWord) { return is_AllStack(); }
    return _RM_UP[r >> _LogWordBits] & (uintptr_t(1) << (r & _WordBitMask));
  }

  // The last bit in the register mask indicates that the mask should repeat
  // indefinitely with ONE bits.  Returns TRUE if mask is infinite or
  // unbounded in size.  Returns FALSE if mask is finite size.
  bool is_AllStack() const { return _all_stack; }
  void set_AllStack(bool value = true) { _all_stack = value; }

  // Test for being a not-empty mask.
  bool is_NotEmpty() const {
    assert(valid_watermarks(), "sanity");
    uintptr_t tmp = 0;
    for (unsigned i = _lwm; i <= _hwm; i++) {
      tmp |= _RM_UP[i];
    }
    assert(tmp || (!tmp && !is_AllStack()), ""); // !tmp => !is_AllStack()
    return tmp;
  }

  // Find lowest-numbered register from mask, or BAD if mask is empty.
  OptoReg::Name find_first_elem() const {
    assert(valid_watermarks(), "sanity");
    for (unsigned i = _lwm; i <= _hwm; i++) {
      uintptr_t bits = _RM_UP[i];
      if (bits) {
        return OptoReg::Name(_offset_bits() + (i << _LogWordBits) + find_lowest_bit(bits));
      }
    }
    return OptoReg::Name(OptoReg::Bad);
  }

  // Get highest-numbered register from mask, or BAD if mask is empty.
  OptoReg::Name find_last_elem() const {
    assert(valid_watermarks(), "sanity");
    if(is_AllStack()) { return OptoReg::Name(INT_MAX); }
    // Careful not to overflow if _lwm == 0
    unsigned i = _hwm + 1;
    while (i > _lwm) {
      uintptr_t bits = _RM_UP[--i];
      if (bits) {
        return OptoReg::Name(_offset_bits() + (i << _LogWordBits) + find_highest_bit(bits));
      }
    }
    return OptoReg::Name(OptoReg::Bad);
  }

  // Clear out partial bits; leave only aligned adjacent bit pairs.
  void clear_to_pairs();

#ifdef ASSERT
  // Verify watermarks are sane, i.e., within bounds and that no
  // register words below or above the watermarks have bits set.
  bool valid_watermarks() const {
    assert(_hwm < _rm_size, "_hwm out of range: %d", _hwm);
    assert(_lwm < _rm_size, "_lwm out of range: %d", _lwm);
    for (unsigned i = 0; i < _lwm; i++) {
      assert(_RM_UP[i] == 0, "_lwm too high: %d regs at: %d", _lwm, i);
    }
    for (unsigned i = _hwm + 1; i < _rm_size; i++) {
      assert(_RM_UP[i] == 0, "_hwm too low: %d regs at: %d", _hwm, i);
    }
    return true;
  }
#endif // !ASSERT

  // Test that the mask contains only aligned adjacent bit pairs
  bool is_aligned_pairs() const;

  // mask is a pair of misaligned registers
  bool is_misaligned_pair() const;
  // Test for single register
  bool is_bound1() const;
  // Test for a single adjacent pair
  bool is_bound_pair() const;
  // Test for a single adjacent set of ideal register's size.
  bool is_bound(uint ireg) const;

  // Check that whether given reg number with size is valid
  // for current regmask, where reg is the highest number.
  bool is_valid_reg(OptoReg::Name reg, const int size) const;

  // Find the lowest-numbered register set in the mask.  Return the
  // HIGHEST register number in the set, or BAD if no sets.
  // Assert that the mask contains only bit sets.
  OptoReg::Name find_first_set(LRG &lrg, const int size) const;

  // Clear out partial bits; leave only aligned adjacent bit sets of size.
  void clear_to_sets(const unsigned int size);
  // Smear out partial bits to aligned adjacent bit sets.
  void smear_to_sets(const unsigned int size);
  // Test that the mask contains only aligned adjacent bit sets
  bool is_aligned_sets(const unsigned int size) const;

  // Test for a single adjacent set
  bool is_bound_set(const unsigned int size) const;

  static bool is_vector(uint ireg);
  static int num_registers(uint ireg);
  static int num_registers(uint ireg, LRG &lrg);

  // Fast overlap test.  Non-zero if any registers in common.
  bool overlap(const RegMask &rm) const {
    assert(_offset_up == rm._offset_up, "");
    assert(valid_watermarks() && rm.valid_watermarks(), "sanity");
    unsigned hwm = MIN2(_hwm, rm._hwm);
    unsigned lwm = MAX2(_lwm, rm._lwm);
    uintptr_t result = 0;
    for (unsigned i = lwm; i <= hwm; i++) {
      result |= _RM_UP[i] & rm._RM_UP[i];
    }
    bool overlap_all_stack = is_AllStack() && rm.is_AllStack();
    return result || overlap_all_stack;
  }

  // Special test for register pressure based splitting
  // UP means register only, Register plus stack, or stack only is DOWN
  bool is_UP() const;

  // Clear a register mask
  void Clear() {
    _lwm = _rm_max();
    _hwm = 0;
    memset(_RM_UP, 0, sizeof(uintptr_t) * _rm_size);
    set_AllStack(false);
    assert(valid_watermarks(), "sanity");
  }

  // Fill a register mask with 1's
  void Set_All() {
    assert(_offset_up == 0, "");
    _lwm = 0;
    _hwm = _rm_max();
    memset(_RM_UP, 0xFF, sizeof(uintptr_t) * _rm_size);
    set_AllStack(true);
    assert(valid_watermarks(), "sanity");
  }

  virtual void Set_From(OptoReg::Name reg) {
    assert(reg != OptoReg::Bad, "sanity");
    assert(reg != OptoReg::Special, "sanity");
    int reg_offset = reg - _offset_bits();
    assert(reg_offset >= 0, "");
    unsigned r = (unsigned)reg_offset;
    assert(r < _rm_size * BitsPerWord, "sanity");
    assert(valid_watermarks(), "pre-condition");
    unsigned index = r >> _LogWordBits;
    _RM_UP[index] |= (uintptr_t(-1) << (r & _WordBitMask));
    if (index < _rm_max()) {
      memset(_RM_UP + index + 1, 0xFF, sizeof(uintptr_t) * (_rm_max() - index));
    }
    if (index < _lwm) _lwm = index;
    _hwm = _rm_max();
    set_AllStack();
    assert(valid_watermarks(), "post-condition");
  }

  // Insert register into mask
  void Insert(OptoReg::Name reg) {
    assert(reg != OptoReg::Bad, "sanity");
    assert(reg != OptoReg::Special, "sanity");
    int reg_offset = reg - _offset_bits();
    assert(reg_offset >= 0, "");
    unsigned r = (unsigned)reg_offset;
    assert(r < _rm_size * BitsPerWord, "sanity");
    assert(valid_watermarks(), "pre-condition");
    unsigned index = r >> _LogWordBits;
    if (index > _hwm) _hwm = index;
    if (index < _lwm) _lwm = index;
    _RM_UP[index] |= (uintptr_t(1) << (r & _WordBitMask));
    assert(valid_watermarks(), "post-condition");
  }

  // Remove register from mask
  void Remove(OptoReg::Name reg) {
    int reg_offset = reg - _offset_bits();
    assert(reg_offset >= 0, "");
    unsigned r = (unsigned)reg_offset;
    assert(r < _rm_size * BitsPerWord, "");
    _RM_UP[r >> _LogWordBits] &= ~(uintptr_t(1) << (r & _WordBitMask));
  }

  // OR 'rm' into 'this'
  virtual void OR(const RegMask &rm) {
    assert(_offset_up == rm._offset_up, "");
    assert(_rm_size >= rm._rm_size, "");
    assert(valid_watermarks() && rm.valid_watermarks(), "sanity");
    // OR widens the live range
    if (_lwm > rm._lwm) _lwm = rm._lwm;
    if (_hwm < rm._hwm) _hwm = rm._hwm;
    for (unsigned i = _lwm; i <= _hwm && i < rm._rm_size; i++) {
      _RM_UP[i] |= rm._RM_UP[i];
    }
    if (rm.is_AllStack() && rm._rm_size < _rm_size ) {
      memset(_RM_UP + rm._rm_size, 0xFF, sizeof(uintptr_t) * (_rm_size - rm._rm_size));
      _hwm = _rm_max();
    }
    set_AllStack(is_AllStack() || rm.is_AllStack());
    assert(valid_watermarks(), "sanity");
  }

  // AND 'rm' into 'this'
  virtual void AND(const RegMask &rm) {
    assert(_offset_up == rm._offset_up, "");
    assert(_rm_size >= rm._rm_size, "");
    assert(valid_watermarks() && rm.valid_watermarks(), "sanity");
    // Do not evaluate words outside the current watermark range, as they are
    // already zero and an &= would not change that
    for (unsigned i = _lwm; i <= _hwm && i < rm._rm_size; i++) {
      _RM_UP[i] &= rm._RM_UP[i];
    }
    // Narrow the watermarks if &rm spans a narrower range.
    // Update after to ensure non-overlapping words are zeroed out.
    if (_lwm < rm._lwm) _lwm = rm._lwm;
    if (_hwm > rm._hwm) _hwm = rm._hwm;
    if (!rm.is_AllStack() && rm._rm_size < _rm_size ) {
      memset(_RM_UP + rm._rm_size, 0, sizeof(uintptr_t) * (_rm_size - rm._rm_size));
      _hwm = rm._rm_max();
    }
    set_AllStack(is_AllStack() && rm.is_AllStack());
  }

  // Subtract 'rm' from 'this'
  virtual void SUBTRACT(const RegMask &rm) {
    assert(_offset_up == rm._offset_up, "");
    assert(_rm_size >= rm._rm_size, "");
    assert(valid_watermarks() && rm.valid_watermarks(), "sanity");
    unsigned hwm = MIN2(_hwm, rm._hwm);
    unsigned lwm = MAX2(_lwm, rm._lwm);
    for (unsigned i = lwm; i <= hwm; i++) {
      _RM_UP[i] &= ~rm._RM_UP[i];
    }
    if (rm.is_AllStack() && rm._rm_size < _rm_size ) {
      memset(_RM_UP + rm._rm_size, 0, sizeof(uintptr_t) * (_rm_size - rm._rm_size));
      _hwm = rm._rm_max();
    }
    set_AllStack(is_AllStack() && !rm.is_AllStack());
  }

  void rollover() {
    assert(is_AllStack_only(),"");
    _offset_up += _rm_size;
    Set_All();
  }

  // Compute size of register mask: number of bits
  uint Size() const;

#ifndef PRODUCT
  void print() const { dump(); }
  void dump(outputStream *st = tty) const; // Print a mask

  bool is_AllStack_only() const {
    assert(valid_watermarks(), "sanity");
    uintptr_t tmp = 0;
    for (unsigned i = _lwm; i <= _hwm; i++) {
      tmp |= _RM_UP[i];
    }
    return !tmp && is_AllStack();
  }
#endif

  static const RegMaskStatic Empty;   // Common empty mask
  static const RegMaskStatic All;     // Common all mask

};

class RegMaskIterator {
 private:
  uintptr_t _current_bits;
  unsigned int _next_index;
  OptoReg::Name _reg;
  const RegMask& _rm;
 public:
  RegMaskIterator(const RegMask& rm) : _current_bits(0), _next_index(rm._lwm), _reg(OptoReg::Bad), _rm(rm) {
    // Calculate the first element
    next();
  }

  bool has_next() {
    return _reg != OptoReg::Bad;
  }

  // Get the current element and calculate the next
  OptoReg::Name next() {
    OptoReg::Name r = _reg;

    // This bit shift scheme, borrowed from IndexSetIterator,
    // shifts the _current_bits down by the number of trailing
    // zeros - which leaves the "current" bit on position zero,
    // then subtracts by 1 to clear it. This quirk avoids the
    // undefined behavior that could arise if trying to shift
    // away the bit with a single >> (next_bit + 1) shift when
    // next_bit is 31/63. It also keeps number of shifts and
    // arithmetic ops to a minimum.

    // We have previously found bits at _next_index - 1, and
    // still have some left at the same index.
    if (_current_bits != 0) {
      unsigned int next_bit = find_lowest_bit(_current_bits);
      assert(_reg != OptoReg::Bad, "can't be in a bad state");
      assert(next_bit > 0, "must be");
      assert(((_current_bits >> next_bit) & 0x1) == 1, "lowest bit must be set after shift");
      _current_bits = (_current_bits >> next_bit) - 1;
      _reg = OptoReg::add(_reg, next_bit);
      return r;
    }

    // Find the next word with bits
    while (_next_index <= _rm._hwm) {
      _current_bits = _rm._RM_UP[_next_index++];
      if (_current_bits != 0) {
        // Found a word. Calculate the first register element and
        // prepare _current_bits by shifting it down and clearing
        // the lowest bit
        unsigned int next_bit = find_lowest_bit(_current_bits);
        assert(((_current_bits >> next_bit) & 0x1) == 1, "lowest bit must be set after shift");
        _current_bits = (_current_bits >> next_bit) - 1;
        _reg = OptoReg::Name(_rm._offset_bits() + ((_next_index - 1) << RegMask::_LogWordBits) + next_bit);
        return r;
      }
    }

    // No more bits
    _reg = OptoReg::Name(OptoReg::Bad);
    return r;
  }
};

class RegMaskStatic : public RegMask {

  // Array of Register Mask bits.  This array is large enough to cover
  // all the machine registers and all parameters that need to be passed
  // on the stack (stack registers) up to some interesting limit.  Methods
  // that need more parameters will NOT be compiled.  On Intel, the limit
  // is something like 90+ parameters.
  int _RM_STORAGE[RM_SIZE];

  public:

  // A constructor only used by the ADLC output.  All mask fields are filled
  // in directly.  Calls to this look something like RM(1,2,3,4);
  RegMaskStatic(
#   define BODY(I) int a##I,
    FORALL_BODY
#   undef BODY
    int dummy = 0) : RegMask(_RM_SIZE) {
    _RM_I = _RM_STORAGE;
#if defined(VM_LITTLE_ENDIAN) || !defined(_LP64)
#   define BODY(I) _RM_I[I] = a##I;
#else
    // We need to swap ints.
#   define BODY(I) _RM_I[I ^ 1] = a##I;
#endif
    FORALL_BODY
#   undef BODY
    _lwm = 0;
    _hwm = _RM_MAX;
    while (_hwm > 0      && _RM_UP[_hwm] == 0) _hwm--;
    while ((_lwm < _hwm) && _RM_UP[_lwm] == 0) _lwm++;
    assert(valid_watermarks(), "post-condition");
  }

  // Construct an empty mask
  RegMaskStatic() : RegMask(_RM_SIZE), _RM_STORAGE() {
    _RM_I = _RM_STORAGE;
    assert(valid_watermarks(), "post-condition");
  }

  // Construct a mask with a single bit
  RegMaskStatic(OptoReg::Name reg) : RegMaskStatic() {
    Insert(reg);
  }

  RegMaskStatic(const RegMaskStatic& rm) : RegMask(_RM_SIZE) {
    _RM_I = _RM_STORAGE;
    _copy(rm,*this);
  }

  RegMaskStatic(const RegMask& rm) : RegMask(_RM_SIZE) {
    _RM_I = _RM_STORAGE;
    _copy(rm,*this);
  }

  RegMaskStatic& operator= (const RegMaskStatic& rm) {
    _copy(rm,*this);
    return *this;
  }

  RegMaskStatic& operator= (const RegMask& rm) {
    _copy(rm,*this);
    return *this;
  }

};

class RegMaskGrowable : public RegMask {

  Arena* _arena;

  void _grow(unsigned int min_size) {
    if(min_size > _rm_size) {
      unsigned int old_size = _rm_size;
      _rm_size = min_size;
      _RM_UP = REALLOC_ARENA_ARRAY(_arena, uintptr_t, _RM_UP, old_size, _rm_size);
      int fill = 0;
      if(is_AllStack()) {
        fill = 0xFF;
        _hwm = _rm_max();
      }
      memset(_RM_UP + old_size, fill, sizeof(uintptr_t) * (_rm_size - old_size));
    }
  }

  public:

  RegMaskGrowable();
  RegMaskGrowable(Arena* arena) : RegMask(_RM_SIZE), _arena(arena) {
    _RM_UP = NEW_ARENA_ARRAY(_arena, uintptr_t, _RM_SIZE);
    memset(_RM_UP, 0, sizeof(uintptr_t) * _RM_SIZE);
  }

  RegMaskGrowable(const RegMask& rm);
  RegMaskGrowable(const RegMask& rm, Arena* arena)
    : RegMask(rm.rm_size()), _arena(arena) {
      _RM_UP = NEW_ARENA_ARRAY(_arena, uintptr_t, rm.rm_size());
      _copy(rm,*this);
    }

  RegMaskGrowable(const RegMaskGrowable& rm);

  void Insert(OptoReg::Name reg) {
    unsigned index = reg >> _LogWordBits;
    unsigned int min_size = index + 1;
    _grow(min_size);
    RegMask::Insert(reg);
  }

  RegMaskGrowable& operator= (const RegMaskGrowable& rm) {
    _grow(rm.rm_size());
    _copy(rm,*this);
    return *this;
  }

  RegMaskGrowable& operator= (const RegMask& rm) {
    _grow(rm.rm_size());
    _copy(rm,*this);
    return *this;
  }

  void OR(const RegMask &rm) {
    _grow(rm.rm_size());
    RegMask::OR(rm);
  }

  void AND(const RegMask &rm) {
    _grow(rm.rm_size());
    RegMask::AND(rm);
  }

  void SUBTRACT(const RegMask &rm) {
    _grow(rm.rm_size());
    RegMask::SUBTRACT(rm);
  }

  void Set_From(OptoReg::Name reg) {
    unsigned index = reg >> _LogWordBits;
    unsigned int min_size = index + 1;
    _grow(min_size);
    RegMask::Set_From(reg);
  }

};

// Do not use this constant directly in client code!
#undef RM_SIZE

#endif // SHARE_OPTO_REGMASK_HPP
