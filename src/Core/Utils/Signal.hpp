/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_SIGNAL_HPP
#define BOGUS_SIGNAL_HPP

#include <functional>
#include <vector>

namespace bogus {

//! Signal class, to which an arbitrary number of listeners can be connected
/*!
        Each time the Signal::trigger() method is called with arguments ( Arg 1, ..., Arg n ),
        the listener functions are called with those same arguments.
        The number and types of arguments are determined by the template parameters of the Signal class.
        */
template <typename... Args>
class Signal {
 public:
  //! Connects the signal to a free function, lambda, or any Args-compatible callable
  void connect(std::function<void(Args...)> callback) { m_callees.emplace_back(std::move(callback)); }

  //! Connects the signal to a member function
  /*! Its signature should be T::member_func( Arg 1, ..., Arg n ) ; */
  template <typename T>
  void connect(T& object, void (T::*member_func)(Args...)) {
    connect([&object, member_func](Args... args) { (object.*member_func)(args...); });
  }

  //! Connects the signal to another Signal
  /*! It should have the same template parameters */
  void connect(const Signal& other) {
    connect([&other](Args... args) { other.trigger(args...); });
  }

  //! Disconnects all listeners
  void disconnectAll() { m_callees.clear(); }

  //! Triggers the signal
  void trigger(Args... args) const {
    for (const auto& callee : m_callees) callee(args...);
  }

 private:
  std::vector<std::function<void(Args...)>> m_callees;
};

}  // namespace bogus

#endif
