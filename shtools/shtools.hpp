/**
 * Copyright (c) 2019 Melown Technologies SE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * *  Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * *  Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef shtools_shtools_hpp_included_
#define shtools_shtools_hpp_included_

#include <vector>

#include <boost/optional.hpp>

#include "math/geometry_core.hpp"

namespace shtools {

/** DH Output grid.
 */
class DHGrid {
public:
    DHGrid(int width, int height)
        : size_(width, height), data_(width * height) {}

    double operator()(int x, int y) const {
        return data_[size_.height * x + y];
    }

    int width() const { return size_.width; }
    int height() const { return size_.height; }

    double* raw() { return data_.data(); }

private:
    math::Size2 size_;
    std::vector<double> data_;
};

enum class Norm : int {
    geodesy = 1
    , schmidt = 2
    , unnormalized = 3
    , orthonormal = 4
};

enum class Sampling : int {
    equallySampled = 1
    , equallySpaced = 2
};

struct GridParams {
    /** Normalization.
     */
    Norm norm = Norm::geodesy;

    /** Condon-Shortley phase factor application.
     */
    bool csPhase = false;

    /** Maximum degree used in MakeGridDH.
     */
    boost::optional<int> lmaxCalc = boost::none;

    /** Deallocate local data if true.
     */
    bool deallocate = false;
};

/** Generates DH grid from sphedical harmonics.
 *
 *  \param sh list of parametrs in the form to build both matrices.
 */
DHGrid makeGridDH(const std::vector<double> &sh
                  , Sampling sampling = Sampling::equallySampled
                  , const GridParams &gridParams = {});

std::vector<double> expand(const std::vector<double> &sh
                           , const math::Points2 &lonlat
                           , const GridParams &gridParams = {});

} // namespace shtools

#endif // shtools_shtools_hpp_included_
