/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2007 StatPro Italia srl

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
#include <ql/exercise.hpp>
#include <ql/pricingengines/blackcalculator.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <utility>

namespace QuantLib {

    AnalyticEuropeanEngine::AnalyticEuropeanEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process)
    : process_(std::move(process)) {
        registerWith(process_);
    }

    AnalyticEuropeanEngine::AnalyticEuropeanEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process,
                                                   Handle<YieldTermStructure> discountCurve,
                                                   ext::optional<unsigned int> spotDays,
                                                   ext::optional<Calendar> spotCalendar)
        : process_(std::move(process)), discountCurve_(std::move(discountCurve)), spotDays_(spotDays), spotCalendar_(spotCalendar) {
        registerWith(process_);
        registerWith(discountCurve_);
    }

    void AnalyticEuropeanEngine::calculate() const {

        // if the discount curve is not specified, we default to the
        // risk free rate curve embedded within the GBM process
        ext::shared_ptr<YieldTermStructure> discountPtr = 
            discountCurve_.empty() ? 
            process_->riskFreeRate().currentLink() :
            discountCurve_.currentLink();

        QL_REQUIRE(arguments_.exercise->type() == Exercise::European,
                   "not an European option");

        ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::dynamic_pointer_cast<StrikedTypePayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-striked payoff given");

        Real variance =
            process_->blackVolatility()->blackVariance(
                                              arguments_.exercise->lastDate(),
                                              payoff->strike());
        const unsigned int spotDays = spotDays_.get_value_or(0);
        const Calendar spotCalendar = spotCalendar_.get_value_or(NullCalendar());
        Date expirySpotDate = spotDays > 0 ? spotCalendar.advance(arguments_.exercise->lastDate(), spotDays * Days)
                                           : arguments_.exercise->lastDate();

        DiscountFactor dividendDiscount = process_->dividendYield()->discount(expirySpotDate);
        DiscountFactor df = discountPtr->discount(arguments_.exercise->lastDate());
        DiscountFactor riskFreeDiscountForFwdEstimation = process_->riskFreeRate()->discount(expirySpotDate);

        Date refDate = process_->dividendYield()->referenceDate();

        const Date spotDate = spotDays > 0 ? spotCalendar.advance(refDate, spotDays * Days) : refDate;
        DiscountFactor dividendDiscountSpotDate = spotDate <= process_->dividendYield()->referenceDate()
                                                      ? 1.0
                                                      : process_->dividendYield()->discount(spotDate);
        DiscountFactor riskFreeDiscountSpotDate =
            spotDate <= process_->riskFreeRate()->referenceDate() ? 1.0 : process_->riskFreeRate()->discount(spotDate);

        Real spot = process_->stateVariable()->value();
        // In case spot is not today, compute underlying price as today
        Real s0 = spot * riskFreeDiscountSpotDate / dividendDiscountSpotDate;
        QL_REQUIRE(spot > 0.0, "negative or null underlying given");
        Real forwardPrice = s0 * dividendDiscount / riskFreeDiscountForFwdEstimation;

        BlackCalculator black(payoff, forwardPrice, std::sqrt(variance), df);

        results_.value = black.value();
        results_.delta = black.delta(spot);
        results_.deltaForward = black.deltaForward();
        results_.elasticity = black.elasticity(spot);
        results_.gamma = black.gamma(spot);

        DayCounter rfdc  = discountPtr->dayCounter();
        DayCounter divdc = process_->dividendYield()->dayCounter();
        DayCounter voldc = process_->blackVolatility()->dayCounter();
        Time t = rfdc.yearFraction(process_->riskFreeRate()->referenceDate(),
                                   arguments_.exercise->lastDate());
        results_.rho = black.rho(t);

        t = divdc.yearFraction(process_->dividendYield()->referenceDate(),
                               arguments_.exercise->lastDate());
        results_.dividendRho = black.dividendRho(t);

        t = voldc.yearFraction(process_->blackVolatility()->referenceDate(),
                               arguments_.exercise->lastDate());
        results_.vega = black.vega(t);
        try {
            results_.theta = black.theta(spot, t);
            results_.thetaPerDay =
                black.thetaPerDay(spot, t);
        } catch (Error&) {
            results_.theta = Null<Real>();
            results_.thetaPerDay = Null<Real>();
        }

        results_.strikeSensitivity  = black.strikeSensitivity();
        results_.itmCashProbability = black.itmCashProbability();

        Real tte = process_->blackVolatility()->timeFromReference(arguments_.exercise->lastDate());
        results_.additionalResults["spot"] = spot;
        results_.additionalResults["dividendDiscount"] = dividendDiscount;
        results_.additionalResults["riskFreeDiscount"] = riskFreeDiscountForFwdEstimation;
        results_.additionalResults["forward"] = forwardPrice;
        results_.additionalResults["strike"] = payoff->strike();
        results_.additionalResults["volatility"] = Real(std::sqrt(variance / tte));
        results_.additionalResults["timeToExpiry"] = tte;
        results_.additionalResults["discountFactor"] = df;
    }

}

