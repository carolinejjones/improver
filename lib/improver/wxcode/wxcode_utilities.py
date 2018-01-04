# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# (C) British Crown Copyright 2017 Met Office.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""This module defines the utilities required for wxcode plugin """

from collections import OrderedDict
import numpy as np
import math
import cartopy.crs as ccrs

from improver.utilities.temporal import iris_time_to_datetime
from improver.utilities.spatial import lat_lon_determine


_WX_DICT_IN = {0: 'Clear_Night',
               1: 'Sunny_Day',
               2: 'Partly_Cloudy_Night',
               3: 'Partly_Cloudy_Day',
               4: 'Dust',
               5: 'Mist',
               6: 'Fog',
               7: 'Cloudy',
               8: 'Overcast',
               9: 'Light_Shower_Night',
               10: 'Light_Shower_Day',
               11: 'Drizzle',
               12: 'Light_Rain',
               13: 'Heavy_Shower_Night',
               14: 'Heavy_Shower_Day',
               15: 'Heavy_Rain',
               16: 'Sleet_Shower_Night',
               17: 'Sleet_Shower_Day',
               18: 'Sleet',
               19: 'Hail_Shower_Night',
               20: 'Hail_Shower_Day',
               21: 'Hail',
               22: 'Light_Snow_Shower_Night',
               23: 'Light_Snow_Shower_Day',
               24: 'Light_Snow',
               25: 'Heavy_Snow_Shower_Night',
               26: 'Heavy_Snow_Shower_Day',
               27: 'Heavy_Snow',
               28: 'Thunder_Shower_Night',
               29: 'Thunder_Shower_Day',
               30: 'Thunder'}

WX_DICT = OrderedDict(sorted(_WX_DICT_IN.items(), key=lambda t: t[0]))

DAYNIGHT_CODES = [1, 3, 10, 14, 17, 20, 23, 26, 29]


def add_wxcode_metadata(cube):
    """ Add weather code metadata to a cube
    Args:
        cube (Iris.cube.Cube):
            Cube which needs weather code metadata added.
    Returns:
        cube (Iris.cube.Cube):
            Cube with weather code metadata added.
    """
    cube.long_name = "weather_code"
    cube.standard_name = None
    cube.var_name = None
    cube.units = "1"
    wx_keys = np.array(WX_DICT.keys())
    cube.attributes.update({'weather_code': wx_keys})
    wxstring = " ".join(WX_DICT.values())
    cube.attributes.update({'weather_code_meaning': wxstring})
    return cube


def expand_nested_lists(query, key):
    """
    Produce flat lists from list and nested lists.

    Args:
        query (dict):
            A single query from the decision tree.
        key (string):
            A string denoting the field to be taken from the dict.

    Returns:
        items (list):
            A 1D list containing all the values for a given key.
    """
    items = []
    for item in query[key]:
        if isinstance(item, list):
            items.extend(item)
        else:
            items.extend([item])
    return items


def calc_day_number(day, month, year):
    """
    Calculate the day of the year.

    Args:
        day (int):
            Day of the month.
        month (int):
            Integer value for month, January = 1, December = 12 etc.
        year (int):
            Year YYYY

    Returns:
        day_number:
             Day of the year, 1 to 365 or 366 in a leap year.
    """
    month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if year % 4 == 0:
        if year % 100 != 0:
            month_days[1] = 29
        elif year % 400 != 0:
            month_days[1] = 28
        else:
            month_days[1] = 29
    day_number = 0
    for i in range(0, month-1):
        day_number += month_days[i]
    day_number += day
    return day_number


def solar_pos(latitude, longitude, day_number, utc_hours):
    """
    Calculate solar elevation for a latitude and longitude
    for a day of the year and an hour (in UTC)

    Args:
        latitude (float):
            Latitude of the point.
        longitude (float):
            Longitude of the point.
        day_number (int):
            Day of the year.
        utc_hour (float):
            Hour of the day in UTC

    Returns:
        solar_elevation (float):
            Solar elevation in radians

    """
    deg2rad = math.pi/180.0

    thetao = 2*math.pi*day_number/365.0

    eqt = (0.000075 + 0.001868 * math.cos(thetao) -
           0.032077 * math.sin(thetao) - 0.014615 * math.cos(2*thetao) -
           0.040849 * math.sin(2*thetao))
    # Declination (degrees):
    decl = -23.5 * math.cos((0.9856*day_number+9.3)*deg2rad)
    # Longitudinal Correction from the Grenwich Meridian
    lon_correction = 24.0*longitude/360.0
    # Solar time (hours):
    solar_time = utc_hours + lon_correction + eqt*12/math.pi
    # Hour angle (degrees):
    degree_hour = (solar_time - 12.0) * 15.0

    # Calculate solar position:

    solar_elevation = (math.asin(math.sin(decl * deg2rad) *
                                 math.sin(latitude * deg2rad) +
                                 math.cos(decl*deg2rad) *
                                 math.cos(latitude*deg2rad) *
                                 math.cos(degree_hour*deg2rad)))

    return solar_elevation


def night_wxcode(wxcode):
    """
    Update weather code from Day value to Night value

    Args:
        wxcode (int):
            Weather symbol code.

    Returns:
        night_wxcode(int):
            Night version of Weather symbol code.
    """
    night_wxcode = wxcode
    if wxcode in DAYNIGHT_CODES:
        night_wxcode = wxcode - 1
    return night_wxcode


def night_day_wxcode(wxcube):
    """
    Modify weather symbols within a cube to account for
    whether it is night or day.

    Please note this method is very slow

    Args:
        wxcube (iris.cube.Cube):
            Weather code cube.

    Returns:
        updated_wxcube (iris.cube.Cube):
            Updated Weather code cube.
    """
    updated_wxcube = wxcube
    dtval = iris_time_to_datetime(wxcube.coord('time'))[0]
    day_number = calc_day_number(dtval.day, dtval.month, dtval.year)
    utc_hours = (dtval.hour * 60.0 + dtval.minute) / 60.0
    trg_latlon = ccrs.PlateCarree()
    trg_crs = lat_lon_determine(wxcube)
    y_points = wxcube.coord('projection_y_coordinate').points
    x_points = wxcube.coord('projection_x_coordinate').points
    for i, yval in enumerate(y_points):
        for j, xval in enumerate(x_points):
            lon, lat = trg_latlon.transform_point(xval, yval, trg_crs)
            solar_elevation = solar_pos(lat, lon, day_number, utc_hours)
            sin_solar_elevation = math.sin(solar_elevation)
            if sin_solar_elevation < 0.0:
                updated_wxcube.data[i, j] = (
                    night_wxcode(updated_wxcube.data[i, j]))
    return updated_wxcube
