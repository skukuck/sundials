/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Base classes for C++ implementations wrappers (views) of SUNDIALS objects.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CLASSVIEW_HPP
#define _SUNDIALS_CLASSVIEW_HPP

#include <memory>
#include <type_traits>
#include <utility>

#include <sundials/sundials_convertibleto.hpp>

namespace sundials {
namespace experimental {


template<class T, class Creator, class Deleter, class... Args>
std::shared_ptr<T> make_our_shared(Args&&... args)
{
  return std::shared_ptr<T>(Creator{}(std::forward<Args>(args)...), Deleter{});
}

template<class T, class Deleter>
std::shared_ptr<T> make_our_shared(T* ptr)
{
  return std::shared_ptr<T>(ptr, Deleter{});
}

template<class T, class Deleter>
class ClassView : public sundials::ConvertibleTo<T>
{
public:
  ClassView() : object_(nullptr, Deleter{})
  {
    fprintf(stderr, ">>>> Creating ClassView:%p holding object.get()=%p\n",
            this, object_.get());
  }

  ClassView(T& object) : object_(std::forward<T>(object), Deleter{})
  {
    fprintf(stderr, ">>>> Creating ClassView:%p holding object.get()=%p\n",
            this, object_.get());
  }

  ClassView(T&& object) : object_(std::forward<T>(object), Deleter{})
  {
    fprintf(stderr, ">>>> Creating ClassView:%p holding object.get()=%p\n",
            this, object_.get());
  }

  ClassView(const ClassView&) = delete;

  ClassView(ClassView&& other)
  {
    fprintf(stderr, ">>>> Move creating ClassView:%p holding object.get()=%p\n",
            this, other.object_.get());
    this->object_ = std::move(other.object_);
  };

  ClassView& operator=(const ClassView&) = delete;

  ClassView& operator=(ClassView&& rhs) = default;

  ~ClassView()
  {
    fprintf(stderr, ">>>> deleting ClassView:%p holding object.get()=%p with use count = %lu\n",
            this, object_.get(), object_.use_count());
    object_.reset();
  };

  // Override ConvertibleTo functions
  T get() override { return object_.get(); }

  T get() const override { return object_.get(); }

  operator T() override { return object_.get(); }

  operator T() const override { return object_.get(); }

protected:
  std::shared_ptr<std::remove_pointer_t<T>> object_;
};

template<class Deleter>
class ClassView<void*, Deleter> : public sundials::ConvertibleTo<void*>
{
public:
  ClassView() : object_(nullptr) {}

  ClassView(void* object) : object_(object) {}

  ClassView(const ClassView&) = delete;

  ClassView(ClassView&& other) noexcept
    : object_(std::exchange(other.object_, nullptr))
  {}

  ClassView& operator=(const ClassView&) = delete;

  ClassView& operator=(ClassView&& rhs) noexcept
  {
    if (object_) { Deleter{}(this->get()); }
    this->object_ = std::exchange(rhs.object_, nullptr);
    return *this;
  };

  ~ClassView()
  {
    if (object_) { Deleter{}(this->get()); }
  };

  void* get() override { return object_; }

  void* get() const override { return object_; }

  operator void*() override { return object_; }

  operator void*() const override { return object_; }

protected:
  void* object_;
};

} // namespace experimental
} // namespace sundials

#endif
