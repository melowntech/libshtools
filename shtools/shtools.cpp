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

#include "dbglog/dbglog.hpp"

#include "utility/streams.hpp"

#include "shtools.hpp"

extern "C" {

void pymakegriddh_(int *exitstatus
                   , double *griddh, int *n
                   , double *cilm, const int *lmax
                   , const int *norm, const int *sampling
                   , const int *csphase, const int *lmax_calc
                   , const int *cilm_d0, const int *cilm_d1, const int *cilm_d2
                   , const int *griddh_d0, const int *griddh_d1);

double pymakegridpoint_(double *cilm, const int *lmax
                        , const double *lat, const double *lon
                        , const int *norm, const int *csphase
                        , const int *dealloc
                        , const int *cilm_d0, const int *cilm_d1
                        , const int *cilm_d2);

} // extern "C"

namespace shtools {

namespace {

std::vector<double> prepareCilm(int &lmax, const std::vector<double> &sh)
{
    auto tmp(std::sqrt(sh.size()));
    int size(tmp);
    if (double(size) != tmp) {
        LOGTHROW(err2, std::logic_error)
            << "Number of spherical harmonics must be power of two, not "
            << sh.size() << ".";
    }

    if (size == 1) {
        LOGTHROW(err2, std::logic_error)
            << "Too little data in harmonics.";
    }

    // left matrix
    const auto matrixSize(size * size);
    std::vector<double> cilm(2 * matrixSize, 0);

    auto ish(sh.begin());

    auto base(cilm.data());
    for (int j(0); j < size; ++j) {
        for (int i(0); i <= j; ++i) {
            base[2 * (i * size + j)] = *ish++;
        }
    }

    ++base;
    for (int j(1); j < size; ++j) {
        for (int i(1); i <= j; ++i) {
            base[2 * (i * size + j)] = *ish++;
        }
    }

    lmax = size - 1;
    return cilm;
}

} // namespace

/** Generates DH longlat grid from spherical harmonics.
 */
DHGrid makeGridDH(const std::vector<double> &sh, Sampling sampling
                  , const GridParams &gridParams)
{
    int lmax;
    auto cilm(prepareCilm(lmax, sh));

    const auto size(lmax + 1);
    int n(2 * lmax + 2);

    const auto sampling_(static_cast<int>(sampling));
    const auto norm_(static_cast<int>(gridParams.norm));
    const auto csPhase_(gridParams.csPhase ? -1 : 1);
    const auto lmaxCalc_(gridParams.lmaxCalc ? *gridParams.lmaxCalc : lmax);

    DHGrid grid(sampling_ * n, n);

    const int cilm_d0(2);
    const int cilm_d1(size);
    const int cilm_d2(size);

    const int griddh_d0(grid.height());
    const int griddh_d1(grid.width());

    int exitstatus(-1);
    pymakegriddh_(&exitstatus
                  , grid.raw(), &n
                  , cilm.data(), &lmax
                  , &norm_, &sampling_, &csPhase_, &lmaxCalc_
                  , &cilm_d0, &cilm_d1, &cilm_d2, &griddh_d0, &griddh_d1);

    switch (exitstatus) {
    case 0: break;
    case 1:
        LOGTHROW(err2, std::runtime_error)
            << "makeGridDH: Improper dimensions of input array.";
    case 2:
        LOGTHROW(err2, std::runtime_error)
            << "makeGridDH: Improper bounds for input variable.";
    case 3:
        LOGTHROW(err2, std::runtime_error)
            << "makeGridDH: Error allocating memory.";
    case 4:
        LOGTHROW(err2, std::runtime_error)
            << "makeGridDH: File IO error.";
    default:
        LOGTHROW(err2, std::runtime_error)
            << "makeGridDH: Invalid exit status <" << exitstatus << ">.";
    }

    return grid;
}

std::vector<double> expand(const std::vector<double> &sh
                           , const math::Points2 &lonlat
                           , const GridParams &gridParams)
{
    int lmax;
    auto cilm(prepareCilm(lmax, sh));

    const auto norm_(static_cast<int>(gridParams.norm));
    const auto csPhase_(gridParams.csPhase ? -1 : 1);
    const int &dealloc_(gridParams.deallocate);

    const auto size(lmax + 1);
    const int cilm_d0(2);
    const int cilm_d1(size);
    const int cilm_d2(size);

    std::vector<double> values;
    values.reserve(lonlat.size());
    for (const auto &p : lonlat) {
        values.push_back
            (pymakegridpoint_
             (cilm.data(), &lmax, &p(0), &p(1), &norm_, &csPhase_, &dealloc_
              , &cilm_d0, &cilm_d1, &cilm_d2));
    }

    return values;
}

} // namespace shtools
