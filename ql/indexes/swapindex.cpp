/*
 Copyright (C) 2006, 2009 Ferdinando Ametrano
 Copyright (C) 2006, 2007, 2009 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.


 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE. See the license for more details. */

#include <ql/indexes/iborindex.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/makeois.hpp>
#include <ql/instruments/makevanillaswap.hpp>
#include <ql/time/schedule.hpp>
#include <sstream>
#include <utility>

namespace QuantLib {

    SwapIndex::SwapIndex(const std::string& familyName,
                         const Period& tenor,
                         Natural settlementDays,
                         const Currency& currency,
                         const Calendar& fixingCalendar,
                         const Period& fixedLegTenor,
                         BusinessDayConvention fixedLegConvention,
                         const DayCounter& fixedLegDayCounter,
                         ext::shared_ptr<IborIndex> iborIndex)
    : InterestRateIndex(
          familyName, tenor, settlementDays, currency, fixingCalendar, fixedLegDayCounter),
      tenor_(tenor), iborIndex_(std::move(iborIndex)), fixedLegTenor_(fixedLegTenor),
      fixedLegConvention_(fixedLegConvention), exogenousDiscount_(false) {
        registerWith(iborIndex_);
    }

    SwapIndex::SwapIndex(const std::string& familyName,
                         const Period& tenor,
                         Natural settlementDays,
                         const Currency& currency,
                         const Calendar& fixingCalendar,
                         const Period& fixedLegTenor,
                         BusinessDayConvention fixedLegConvention,
                         const DayCounter& fixedLegDayCounter,
                         ext::shared_ptr<IborIndex> iborIndex,
                         Handle<YieldTermStructure> discount)
    : InterestRateIndex(
          familyName, tenor, settlementDays, currency, fixingCalendar, fixedLegDayCounter),
      tenor_(tenor), iborIndex_(std::move(iborIndex)), fixedLegTenor_(fixedLegTenor),
      fixedLegConvention_(fixedLegConvention), exogenousDiscount_(true),
      discount_(std::move(discount)) {
        registerWith(iborIndex_);
        registerWith(discount_);
    }

    Handle<YieldTermStructure> SwapIndex::forwardingTermStructure() const {
        return iborIndex_->forwardingTermStructure();
    }

    Handle<YieldTermStructure> SwapIndex::discountingTermStructure() const {
        return discount_;  // empty if not exogenous
    }

    Rate SwapIndex::forecastFixing(const Date& fixingDate) const {
        return underlyingSwap(fixingDate)->fairRate();
    }

    ext::shared_ptr<VanillaSwap>
    SwapIndex::underlyingSwap(const Date& fixingDate) const {

        QL_REQUIRE(fixingDate!=Date(), "null fixing date");

        // caching mechanism
        if (lastFixingDate_!=fixingDate) {
            Rate fixedRate = 0.0;
            if (exogenousDiscount_)
                lastSwap_ = MakeVanillaSwap(tenor_, iborIndex_, fixedRate)
                    .withEffectiveDate(valueDate(fixingDate))
                    .withFixedLegCalendar(fixingCalendar())
                    .withFixedLegDayCount(dayCounter_)
                    .withFixedLegTenor(fixedLegTenor_)
                    .withFixedLegConvention(fixedLegConvention_)
                    .withFixedLegTerminationDateConvention(fixedLegConvention_)
                    .withDiscountingTermStructure(discount_);
            else
                lastSwap_ = MakeVanillaSwap(tenor_, iborIndex_, fixedRate)
                    .withEffectiveDate(valueDate(fixingDate))
                    .withFixedLegCalendar(fixingCalendar())
                    .withFixedLegDayCount(dayCounter_)
                    .withFixedLegTenor(fixedLegTenor_)
                    .withFixedLegConvention(fixedLegConvention_)
                    .withFixedLegTerminationDateConvention(fixedLegConvention_);
            lastFixingDate_ = fixingDate;
        }
        return lastSwap_;
    }

    Date SwapIndex::maturityDate(const Date& valueDate) const {
        Date fixDate = fixingDate(valueDate);
        return underlyingSwap(fixDate)->maturityDate();
    }

    ext::shared_ptr<SwapIndex>
    SwapIndex::clone(const Handle<YieldTermStructure>& forwarding) const {

        if (exogenousDiscount_)
            return ext::make_shared<SwapIndex>(familyName(),
                          tenor(),
                          fixingDays(),
                          currency(),
                          fixingCalendar(),
                          fixedLegTenor(),
                          fixedLegConvention(),
                          dayCounter(),
                          iborIndex_->clone(forwarding),
                          discount_);
        else
            return ext::make_shared<SwapIndex>(familyName(),
                          tenor(),
                          fixingDays(),
                          currency(),
                          fixingCalendar(),
                          fixedLegTenor(),
                          fixedLegConvention(),
                          dayCounter(),
                          iborIndex_->clone(forwarding));
    }

    ext::shared_ptr<SwapIndex>
    SwapIndex::clone(const Handle<YieldTermStructure>& forwarding,
                     const Handle<YieldTermStructure>& discounting) const {
        return ext::make_shared<SwapIndex>(familyName(),
                       tenor(),
                       fixingDays(),
                       currency(),
                       fixingCalendar(),
                       fixedLegTenor(),
                       fixedLegConvention(),
                       dayCounter(),
                       iborIndex_->clone(forwarding),
                       discounting);
    }

    ext::shared_ptr<SwapIndex>
    SwapIndex::clone(const Period& tenor) const {

        if (exogenousDiscount_)
            return ext::make_shared<SwapIndex>(familyName(),
                          tenor,
                          fixingDays(),
                          currency(),
                          fixingCalendar(),
                          fixedLegTenor(),
                          fixedLegConvention(),
                          dayCounter(),
                          iborIndex(),
                          discountingTermStructure());
        else
            return ext::make_shared<SwapIndex>(familyName(),
                          tenor,
                          fixingDays(),
                          currency(),
                          fixingCalendar(),
                          fixedLegTenor(),
                          fixedLegConvention(),
                          dayCounter(),
                          iborIndex());

    }

    OvernightIndexedSwapIndex::OvernightIndexedSwapIndex(
        const std::string& familyName,
        const Period& tenor,
        Natural settlementDays,
        const Currency& currency,
        const ext::shared_ptr<OvernightIndex>& overnightIndex,
        bool telescopicValueDates,
        RateAveraging::Type averagingMethod,
        const Period& fixedLegTenor,
        Handle<YieldTermStructure> discount)
    : SwapIndex(familyName,
                tenor,
                settlementDays,
                currency,
                overnightIndex->fixingCalendar(),
                fixedLegTenor,
                ModifiedFollowing,
                overnightIndex->dayCounter(),
                overnightIndex,
                discount),
      overnightIndex_(overnightIndex), telescopicValueDates_(telescopicValueDates),
      averagingMethod_(averagingMethod) {}


    ext::shared_ptr<OvernightIndexedSwap>
    OvernightIndexedSwapIndex::underlyingSwap(const Date& fixingDate) const {

        QL_REQUIRE(fixingDate!=Date(), "null fixing date");

        // caching mechanism
        if (lastFixingDate_!=fixingDate) {
            Rate fixedRate = 0.0;
            if (exogenousDiscount_) {
                lastSwap_ = MakeOIS(tenor_, overnightIndex_, fixedRate)
                                .withEffectiveDate(valueDate(fixingDate))
                                .withFixedLegDayCount(dayCounter_)
                                .withTelescopicValueDates(telescopicValueDates_)
                                .withAveragingMethod(averagingMethod_)
                                .withPaymentFrequency(fixedLegTenor_.frequency())
                                .withDiscountingTermStructure(discount_);
            } else {
                lastSwap_ = MakeOIS(tenor_, overnightIndex_, fixedRate)
                                .withEffectiveDate(valueDate(fixingDate))
                                .withFixedLegDayCount(dayCounter_)
                                .withTelescopicValueDates(telescopicValueDates_)
                                .withAveragingMethod(averagingMethod_)
                                .withPaymentFrequency(fixedLegTenor_.frequency())
                                .withDiscountingTermStructure(discount_);
            }
            lastFixingDate_ = fixingDate;
        }
        return lastSwap_;
    }

    Rate OvernightIndexedSwapIndex::forecastFixing(const Date& fixingDate) const {
        return underlyingSwap(fixingDate)->fairRate();
    }

    ext::shared_ptr<SwapIndex>
    OvernightIndexedSwapIndex::clone(const Handle<YieldTermStructure>& forwarding) const {

        if (exogenousDiscount_)
            return ext::shared_ptr<SwapIndex>(new OvernightIndexedSwapIndex(
                familyName(), tenor(), fixingDays(), currency(),
                ext::dynamic_pointer_cast<OvernightIndex>(overnightIndex()->clone(forwarding)),
                telescopicValueDates(), averagingMethod(), fixedLegTenor(), discount_));
        else
            return ext::shared_ptr<SwapIndex>(new OvernightIndexedSwapIndex(
                familyName(), tenor(), fixingDays(), currency(),
                ext::dynamic_pointer_cast<OvernightIndex>(overnightIndex()->clone(forwarding)),
                telescopicValueDates(), averagingMethod(), fixedLegTenor()));
    }

    ext::shared_ptr<SwapIndex>
    OvernightIndexedSwapIndex::clone(const Handle<YieldTermStructure>& forwarding,
                                     const Handle<YieldTermStructure>& discounting) const {
        return ext::shared_ptr<SwapIndex>(new OvernightIndexedSwapIndex(
            familyName(), tenor(), fixingDays(), currency(),
            ext::dynamic_pointer_cast<OvernightIndex>(overnightIndex()->clone(forwarding)),
            telescopicValueDates(), averagingMethod(), fixedLegTenor(), discounting));
    }

    ext::shared_ptr<SwapIndex> OvernightIndexedSwapIndex::clone(const Period& tenor) const {

        if (exogenousDiscount_)
            return ext::shared_ptr<SwapIndex>(new OvernightIndexedSwapIndex(
                familyName(), tenor, fixingDays(), currency(),
                ext::dynamic_pointer_cast<OvernightIndex>(overnightIndex()), telescopicValueDates(),
                averagingMethod(), fixedLegTenor(), discount_));
        else
            return ext::shared_ptr<SwapIndex>(new OvernightIndexedSwapIndex(
                familyName(), tenor, fixingDays(), currency(),
                ext::dynamic_pointer_cast<OvernightIndex>(overnightIndex()), telescopicValueDates(),
                averagingMethod(), fixedLegTenor()));
    }
}
