
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2018 Matthias Lungwitz

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

#ifndef quantlib_lazy_object_i
#define quantlib_lazy_object_i

%include common.i
%include types.i
%include observer.i

%{
using QuantLib::LazyObject;
%}

%shared_ptr(LazyObject)
class LazyObject : public Observable {
  private:
    LazyObject();
  public:
    void recalculate();
    void freeze();
    void unfreeze();
    %extend {
        static void forwardFirstNotificationOnly() {
            LazyObject::Defaults::instance().forwardFirstNotificationOnly();
        }
        static void alwaysForwardNotifications() {
            LazyObject::Defaults::instance().alwaysForwardNotifications();
        }
        static bool forwardsAllNotifications() {
            return LazyObject::Defaults::instance().forwardsAllNotifications();
        }
    }
};


#endif
