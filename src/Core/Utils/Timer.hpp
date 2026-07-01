/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_TIMER_HPP
#define BOGUS_TIMER_HPP

#include <chrono>

namespace bogus {

//! Simple timer class. Starts when constructed.
class Timer {
 public:
  Timer() { reset(); }

  //! Returns the elapsed time, in seconds, since the last call to reset()
  double elapsed() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - m_start).count(); }

  //! Restarts the timer
  void reset() { m_start = std::chrono::steady_clock::now(); }

 private:
  std::chrono::steady_clock::time_point m_start;
};

}  // namespace bogus

#endif
