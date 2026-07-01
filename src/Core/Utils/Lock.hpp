/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_LOCK_HPP
#define BOGUS_LOCK_HPP

#include <mutex>

namespace bogus {

//! Simple RAII mutex wrapper, used to guard concurrent writes to shared
//! block-matrix state from OpenMP-parallelized loops.
/*!
        Copying a Lock does not copy its locked state: it default-constructs
        a fresh, unlocked mutex. This is required as Lock is a plain member
        of copyable classes such as SparseBlockMatrixBase.
        */
class Lock {
 public:
  Lock() = default;
  Lock(const Lock&) {}
  Lock& operator=(const Lock&) { return *this; }

  void set() const { m_mutex.lock(); }
  void unset() const { m_mutex.unlock(); }

  //! RAII lock guard. With DoLock = false, locking is skipped entirely
  //! (resolved at compile-time), so call sites that don't need thread-safety
  //! pay no runtime cost.
  template <bool DoLock = true>
  struct Guard {
    explicit Guard(const Lock& lock) : m_lock(lock) {
      if constexpr (DoLock) m_lock.set();
    }

    ~Guard() {
      if constexpr (DoLock) m_lock.unset();
    }

    Guard(const Guard&) = delete;
    Guard& operator=(const Guard&) = delete;

   private:
    const Lock& m_lock;
  };

 private:
  mutable std::mutex m_mutex;
};

}  // namespace bogus

#endif
