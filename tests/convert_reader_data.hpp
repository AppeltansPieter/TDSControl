#pragma once

#include <tdscontrol/tds.hpp>
#include <tdsreader/system.hpp>

/**
 * @brief Utility function to convert a tdsreader tds to a tdscontrol tds.
 *
 * @param data The tdsreader data to convert
 * @return tds::tds The equivalent time-delay system.
 */
inline tds::tds create_tds_from_data(tdsreader::tds data)
{
    return tds::tds{std::move(data.A), std::move(data.hA)};
}