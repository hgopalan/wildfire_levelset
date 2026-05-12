// sim_datetime.H
// Shared helper for computing human-readable simulation date/time strings
// from a calendar start date + elapsed simulation time.
//
// Used by main.cpp (per-step log output) and fire_report.H (HTML header).
//
// Notes:
//   - Gregorian leap-year rule is applied (divisible by 4, except centuries,
//     except those divisible by 400).
//   - Result is a "YYYY-MM-DD HH:MM" string.
//   - Returns empty string when year == 0 (no start date configured).

#ifndef SIM_DATETIME_H
#define SIM_DATETIME_H

#include <cmath>
#include <cstdio>
#include <string>

namespace WildfireDatetime {

/// Returns true when `year` is a Gregorian leap year.
inline bool is_leap_year(int year)
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

/// Returns the number of days in `month` of `year`.
inline int days_in_month(int month, int year)
{
    static const int dom_common[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
    if (month == 2 && is_leap_year(year)) return 29;
    if (month >= 1 && month <= 12)        return dom_common[month];
    return 30; // fallback (should not reach)
}

/// Build a "YYYY-MM-DD HH:MM" string.
///
/// \param year        Calendar year at simulation t=0 (0 = no date configured)
/// \param month       Calendar month at t=0 (1–12)
/// \param day         Calendar day at t=0 (1–31)
/// \param start_hour  Decimal clock hour at t=0 (e.g. 8.5 = 08:30)
/// \param elapsed_s   Elapsed simulation time [s]
/// \return Formatted datetime string, or "" when year == 0.
inline std::string elapsed_to_datetime(int year, int month, int day,
                                        double start_hour, double elapsed_s)
{
    if (year <= 0) return "";

    double total_h = start_hour + elapsed_s / 3600.0;
    int extra_days = static_cast<int>(std::floor(total_h / 24.0));
    double frac_h  = total_h - extra_days * 24.0;

    int hh = static_cast<int>(frac_h);
    int mm = static_cast<int>((frac_h - hh) * 60.0 + 0.5);
    if (mm >= 60) { mm -= 60; ++hh; }
    if (hh >= 24) { hh -= 24; ++extra_days; }

    // Advance date by extra_days with proper month lengths and leap years
    int d = day + extra_days;
    int m = month;
    int y = year;

    while (d > days_in_month(m, y)) {
        d -= days_in_month(m, y);
        ++m;
        if (m > 12) { m = 1; ++y; }
    }

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d",
                  y, m, d, hh, mm);
    return buf;
}

} // namespace WildfireDatetime

#endif // SIM_DATETIME_H
