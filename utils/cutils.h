/***************************************************************************
 * Copyright (C) 2023 Francesco Florian
 * This file is part of SkeletAlg.
 *
 * SkeletAlg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SkeletAlg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SkeletAlg.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The copyright holders give you permission to combine SkeletAlg
 * with code included in the standard release of Netgen (from Joachim
 * Schöberl), METIS (from George Karypis at the University of
 * Minnesota), OpenCASCADE (from Open CASCADE S.A.S) and ParaView
 * (from Kitware, Inc.) under their respective licenses. You may copy
 * and distribute such a system following the terms of the GNU GPL for
 * Gmsh and the licenses of the other code concerned, provided that
 * you include the source code of that other code when and as the GNU
 * GPL requires distribution of source code.
 *
 * Note that people who make modified versions of SkeletAlg are not
 * obligated to grant this special exception for their modified
 * versions; it is their choice whether to do so. The GNU General
 * Public License gives permission to release a modified version
 * without this exception; this exception also makes it possible to
 * release a modified version which carries forward this exception.
 *
 * Additional permission under GNU GPL version 3 section 7
 * If you modify this Program, or any covered work, by linking or combining it
 * with H2Lib (https://github.com/H2Lib/H2Lib), (or a modified version of that
 * library), containing parts covered by the terms all right reserved,the
 * licensors of this Program grant you additional permission to convey the
 * resulting work.
 ***************************************************************************/
#ifndef CUTILS_H
#define CUTILS_H

#include <type_traits>

#define DECLAREcsTRUCT(cstruct) \
  struct _ ## cstruct; \
  typedef struct _ ## cstruct cstruct; \
  extern "C" void del_ ## cstruct ( cstruct * );

DECLAREcsTRUCT(amatrix)
DECLAREcsTRUCT(avector)
DECLAREcsTRUCT(cluster)
DECLAREcsTRUCT(clustergeometry)
DECLAREcsTRUCT(dblock)
DECLAREcsTRUCT(dcluster)
DECLAREcsTRUCT(dclusterbasis)
DECLAREcsTRUCT(dclusteroperator)
DECLAREcsTRUCT(helmholtz3d)
DECLAREcsTRUCT(leveldir)
DECLAREcsTRUCT(reporter)
DECLAREcsTRUCT(sparsematrix)
DECLAREcsTRUCT(stopwatch)
DECLAREcsTRUCT(surface3d)
#undef DECLAREcsTRUCT

struct _diradmdata;
typedef struct _diradmdata diradmdata;

void del_diradmdata(diradmdata* object);
diradmdata* new_diradmdata();

#include "h2lib/basic.h"
#include "h2lib/dh2matrix.h"

namespace MainNamespace {
  namespace CFunctions {
    // the following lines is an "aya" snippet, to create more overloads
    // inline void freeCPointer(~type *~object) { del_type(~object); }
    // DECLAREcsTRUCT(~type)
    inline void freeCPointer(amatrix *object) { del_amatrix(object); }
    inline void freeCPointer(avector *object) { del_avector(object); }
    inline void freeCPointer(cluster *object) { del_cluster(object); }
    inline void freeCPointer(clustergeometry *object) { del_clustergeometry(object); }
    inline void freeCPointer(dblock *object) { del_dblock(object); }
    inline void freeCPointer(dcluster *object) { del_dcluster(object); }
    inline void freeCPointer(dclusteroperator *object) { del_dclusteroperator(object); }
    inline void freeCPointer(diradmdata *object) { del_diradmdata(object); }
    inline void freeCPointer(helmholtz3d *object) { del_helmholtz3d(object); }
    inline void freeCPointer(leveldir *object) { del_leveldir(object); }
    inline void freeCPointer(reporter *object) { del_reporter(object); }
    inline void freeCPointer(sparsematrix *object) { del_sparsematrix(object); }
    inline void freeCPointer(stopwatch *object) { del_stopwatch(object); }
    inline void freeCPointer(surface3d *object) { del_surface3d(object); }

    inline void freeCPointer(double *object) { freemem(static_cast<void *>(object)); }
    inline void freeCPointer(dh2matrix *object) {
      del_dclusterbasis(const_cast<dclusterbasis*>(object->cb));
      del_dclusterbasis(const_cast<dclusterbasis*>(object->rb));
      del_dh2matrix(object);
    }
  }
  template<typename T>
  inline void maybeFree(T *object) {
    if(object != nullptr)
      CFunctions::freeCPointer(const_cast<std::remove_const_t<T> *>(object));
  }

  template<class BaseData>
  class DumbPointer {
  public:
    inline DumbPointer() {}
    inline DumbPointer(BaseData * const pointer) : data_(pointer) {}
    DumbPointer(const DumbPointer &other) = delete;
    inline DumbPointer& operator=(DumbPointer &other) = delete;
    inline DumbPointer& operator=(DumbPointer &&other) noexcept {
      if ( this != &other ) {
        maybeFree(data_);
        data_ = other.data_;
        other.data_ = nullptr;
      }
      return *this;
    }
    inline DumbPointer& operator=(BaseData * const pointer) noexcept {
      maybeFree(data_);
      data_ = pointer;
      return *this;
    }
    inline DumbPointer(DumbPointer &&other) noexcept {
      data_ = other.data_;
      other.data_ = nullptr;
    }
    inline ~DumbPointer() {
      maybeFree(data_);
    }

    inline BaseData *data() {
      return data_;
    };
    inline bool operator==(std::nullptr_t nullPointer) const { return data_ == nullPointer; }
    inline BaseData * operator*() { return data_; }
    inline BaseData const * operator*() const { return data_; }
    inline BaseData * operator->() { return data_; }
    inline BaseData const * operator->() const { return data_; }
  private:
    BaseData *data_ = nullptr;
  };
}

#endif /* CUTILS_H */
