/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, 2024 Peter Caspers

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file laplaceinterpolation.hpp
    \brief Laplace interpolation of missing values
*/

#include <ql/experimental/math/laplaceinterpolation.hpp>
#include <ql/math/matrixutilities/bicgstab.hpp>
#include <ql/math/matrixutilities/sparsematrix.hpp>
#include <ql/methods/finitedifferences/meshers/fdm1dmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/meshers/predefined1dmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/secondderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/triplebandlinearop.hpp>

namespace QuantLib {

    LaplaceInterpolation::LaplaceInterpolation(std::function<Real(const std::vector<Size>&)> y,
                                               std::vector<std::vector<Real>> x,
                                               const Real relTol)
    : y_(std::move(y)), x_(std::move(x)), relTol_(relTol) {

        // set up the mesher

        std::vector<Size> dim(x.size());
        std::transform(x.begin(), x.end(), dim.begin(),
                       [](const std::vector<Real>& x) { return x.size(); });
        for (Size i = 0; i < dim.size(); ++i) {
            QL_REQUIRE(
                dim[i] >= 2,
                "LaplaceInterpolation: at least two points are required in each dimension, have "
                    << dim[i] << " for #" << i);
        }
        layout_ = boost::make_shared<FdmLinearOpLayout>(dim);
        std::vector<boost::shared_ptr<Fdm1dMesher>> meshers(x.size());
        std::transform(x.begin(), x.end(), meshers.begin(), [](const std::vector<Real>& x) {
            return boost::make_shared<Predefined1dMesher>(x);
        });
        auto mesher = boost::make_shared<FdmMesherComposite>(layout_, meshers);

        // set up the Laplace operator and convert it to matrix

        struct LaplaceOp : public FdmLinearOpComposite {
            LaplaceOp(const boost::shared_ptr<FdmMesher>& mesher)
            : map_(SecondDerivativeOp(0, mesher)) {
                for (Size direction = 0; direction < mesher->layout()->dim().size(); ++direction) {
                    map_ = map_.add(SecondDerivativeOp(direction, mesher));
                }
            }
            TripleBandLinearOp map_;

            Size size() const override { QL_FAIL("no impl"); }
            void setTime(Time t1, Time t2) override { QL_FAIL("no impl"); }
            Array apply(const array_type& r) const override { QL_FAIL("no impl"); }
            Array apply_mixed(const Array& r) const override { QL_FAIL("no impl"); }
            Array apply_direction(Size direction, const Array& r) const override {
                QL_FAIL("no impl");
            }
            Array solve_splitting(Size direction, const Array& r, Real s) const override {
                QL_FAIL("no impl");
            }
            Array preconditioner(const Array& r, Real s) const override { QL_FAIL("no impl"); }
        };

        SparseMatrix op = LaplaceOp(mesher).toMatrix();

        // set up the linear system to solve

        Size N = layout_->size();

        SparseMatrix g(N, N, 5 * N);
        Array rhs(N, 0.0), guess(N, 0.0);
        Real guessTmp = 0.0;

        struct f_A {
            const SparseMatrix& g;
            explicit f_A(const SparseMatrix& g) : g(g) {}
            Array operator()(const Array& x) const { return prod(g, x); }
        };

        auto rowit = op.begin1();
        Size count = 0;
        std::vector<Real> corner_h(dim.size());
        std::vector<Size> corner_neighbour_index(dim.size());
        for (auto const& pos : *layout_) {
            auto coord = pos.coordinates();
            Real val = y(coord);
            QL_REQUIRE(rowit != op.end1() && rowit.index1() == count,
                       "LaplaceInterpolation: op matrix row iterator ("
                           << (rowit != op.end1() ? std::to_string(rowit.index1()) : "na")
                           << ") does not match expected row count (" << count << ")");
            if (val == Null<Real>()) {
                bool isCorner = true;
                for (Size d = 0; d < dim.size() && isCorner; ++d) {
                    if (coord[d] == 0) {
                        corner_h[d] = meshers[d]->dplus(0);
                        corner_neighbour_index[d] = 1;
                    } else if (coord[d] == layout_->dim()[d] - 1) {
                        corner_h[d] = meshers[d]->dminus(dim[d] - 1);
                        corner_neighbour_index[d] = dim[d] - 2;
                    } else {
                        isCorner = false;
                    }
                }
                if (isCorner) {
                    // shandling of the "corners", all second derivs are zero in the op
                    // this directly generalizes Numerical Recipes, 3rd ed, eq 3.8.6
                    Real sum_corner_h =
                        std::accumulate(corner_h.begin(), corner_h.end(), 0.0, std::plus<Real>());
                    for (Size j = 0; j < dim.size(); ++j) {
                        std::vector<Size> coord_j(coord);
                        coord_j[j] = corner_neighbour_index[j];
                        Real weight = 0.0;
                        for (Size i = 0; i < dim.size(); ++i) {
                            if (i != j)
                                weight += corner_h[i];
                        }
                        weight = 1.0 - weight / sum_corner_h;
                        g(count, layout_->index(coord_j)) = -weight;
                    }
                } else {
                    // point with at least one dimension with non-trivial second derivative
                    for (auto colit = rowit.begin(); colit != rowit.end(); ++colit)
                        g(count, colit.index2()) = *colit;
                }
                rhs[count] = 0.0;
                guess[count] = guessTmp;
            } else {
                g(count, count) = 1;
                rhs[count] = val;
                guess[count] = guessTmp = val;
            }
            ++count;
        }

        interpolatedValues_ = BiCGstab(f_A(g), 10 * N, relTol_).solve(rhs, guess).x;
    }


    Real LaplaceInterpolation::operator()(const std::vector<Size>& coordinates) const {
        return interpolatedValues_[layout_->index(coordinates)];
    }

    template <class M>
    void laplaceInterpolation(M& A,
                              Real relTol,
                              const std::vector<Real>& x,
                              const std::vector<Real>& y) {

        std::vector<std::vector<Real>> tmp;
        tmp.push_back(y);
        tmp.push_back(x);

        if (y.empty()) {
            tmp[0].resize(A.rows());
            std::iota(tmp[0].begin(), tmp[0].end(), 0.0);
        }

        if (x.empty()) {
            tmp[1].resize(A.columns());
            std::iota(tmp[1].begin(), tmp[1].end(), 0.0);
        }

        LaplaceInterpolation interpolation(
            [&A](const std::vector<Size>& coordinates) {
                return A(coordinates[0], coordinates[1]);
            },
            tmp, relTol);

        for (Size i = 0; i < A.rows(); ++i) {
            for (Size j = 0; j < A.columns(); ++j) {
                if (A(i, j) == Null<Real>())
                    A(i, j) = interpolation({i, j});
            }
        }
    }

    // template instantiations for matrix classes we want to support

    template void laplaceInterpolation(Matrix& A,
                                       Real relTol,
                                       const std::vector<Real>& x,
                                       const std::vector<Real>& y);

} // namespace QuantLib
