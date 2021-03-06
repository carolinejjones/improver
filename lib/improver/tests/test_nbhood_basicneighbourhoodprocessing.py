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
"""Unit tests for the nbhood.BasicNeighbourhoodProcessing plugin."""


import unittest

from cf_units import Unit
from iris.coords import AuxCoord, DimCoord
from iris.coord_systems import OSGB
from iris.cube import Cube
from iris.tests import IrisTest
import numpy as np


from improver.grids.osgb import OSGBGRID
from improver.nbhood import BasicNeighbourhoodProcessing as NBHood


SINGLE_POINT_RANGE_3_CENTROID = np.array([
    [0.992, 0.968, 0.96, 0.968, 0.992],
    [0.968, 0.944, 0.936, 0.944, 0.968],
    [0.96, 0.936, 0.928, 0.936, 0.96],
    [0.968, 0.944, 0.936, 0.944, 0.968],
    [0.992, 0.968, 0.96, 0.968, 0.992]
])

SINGLE_POINT_RANGE_X_6_Y_3_CENTROID = np.array([
    [1.0, 1.0, 0.99799197, 0.99598394, 0.99799197, 1.0, 1.0],
    [1.0, 0.98995984, 0.98393574, 0.98192771, 0.98393574, 0.98995984, 1.0],
    [0.98995984, 0.97991968, 0.97389558, 0.97188755, 0.97389558, 0.97991968,
     0.98995984],
    [0.98393574, 0.97389558, 0.96787149, 0.96586345, 0.96787149, 0.97389558,
     0.98393574],
    [0.98192771, 0.97188755, 0.96586345, 0.96385542, 0.96586345, 0.97188755,
     0.98192771],
    [0.98393574, 0.97389558, 0.96787149, 0.96586345, 0.96787149, 0.97389558,
     0.98393574],
    [0.98995984, 0.97991968, 0.97389558, 0.97188755, 0.97389558, 0.97991968,
     0.98995984],
    [1.0, 0.98995984, 0.98393574, 0.98192771, 0.98393574, 0.98995984, 1.0],
    [1.0, 1.0, 0.99799197, 0.99598394, 0.99799197, 1.0, 1.0]
])

SINGLE_POINT_RANGE_2_CENTROID_FLAT = np.array([
    [1.0, 1.0, 0.92307692, 1.0, 1.0],
    [1.0, 0.92307692, 0.92307692, 0.92307692, 1.0],
    [0.92307692, 0.92307692, 0.92307692, 0.92307692, 0.92307692],
    [1.0, 0.92307692, 0.92307692, 0.92307692, 1.0],
    [1.0, 1.0, 0.92307692, 1.0, 1.0]
])

SINGLE_POINT_RANGE_5_CENTROID = np.array([
    [1.0, 1.0, 0.99486125, 0.99177801, 0.99075026, 0.99177801,
     0.99486125, 1.0, 1.0],
    [1.0, 0.99280576, 0.98766701, 0.98458376, 0.98355601,
     0.98458376, 0.98766701, 0.99280576, 1.0],
    [0.99486125, 0.98766701, 0.98252826, 0.97944502, 0.97841727,
     0.97944502, 0.98252826, 0.98766701, 0.99486125],
    [0.99177801, 0.98458376, 0.97944502, 0.97636177, 0.97533402,
     0.97636177, 0.97944502, 0.98458376, 0.99177801],
    [0.99075026, 0.98355601, 0.97841727, 0.97533402, 0.97430627,
     0.97533402, 0.97841727, 0.98355601, 0.99075026],
    [0.99177801, 0.98458376, 0.97944502, 0.97636177, 0.97533402,
     0.97636177, 0.97944502, 0.98458376, 0.99177801],
    [0.99486125, 0.98766701, 0.98252826, 0.97944502, 0.97841727,
     0.97944502, 0.98252826, 0.98766701, 0.99486125],
    [1.0, 0.99280576, 0.98766701, 0.98458376, 0.98355601,
     0.98458376, 0.98766701, 0.99280576, 1.0],
    [1.0, 1.0, 0.99486125, 0.99177801, 0.99075026, 0.99177801,
     0.99486125, 1.0, 1.0]
])


def set_up_cube(zero_point_indices=((0, 7, 7),), num_time_points=1,
                num_grid_points=16):
    """Set up a normal OSGB UK National Grid cube."""
    data = np.ones((num_time_points, num_grid_points, num_grid_points))
    for time_index, lat_index, lon_index in zero_point_indices:
        data[time_index][lat_index][lon_index] = 0
    cube = Cube(data, standard_name="precipitation_amount",
                units="kg m^-2 s^-1")
    coord_system = OSGB()
    scaled_y_coord = OSGBGRID.coord('projection_y_coordinate')
    cube.add_dim_coord(
        DimCoord(
            scaled_y_coord.points[:num_grid_points],
            'projection_y_coordinate',
            units='m', coord_system=coord_system
        ),
        1
    )
    scaled_x_coord = OSGBGRID.coord('projection_x_coordinate')
    cube.add_dim_coord(
        DimCoord(
            scaled_x_coord.points[:num_grid_points],
            'projection_x_coordinate',
            units='m', coord_system=coord_system
        ),
        2
    )
    tunit = Unit("hours since 1970-01-01 00:00:00", "gregorian")
    time_points = [402192.5 + _ for _ in range(num_time_points)]
    cube.add_aux_coord(AuxCoord(time_points,
                                "time", units=tunit), 0)
    return cube


def set_up_cube_lat_long(zero_point_indices=((0, 7, 7),), num_time_points=1,
                         num_grid_points=16):
    """Set up a lat-long coord cube."""
    data = np.ones((num_time_points, num_grid_points, num_grid_points))
    for time_index, lat_index, lon_index in zero_point_indices:
        data[time_index][lat_index][lon_index] = 0
    cube = Cube(data, standard_name="precipitation_amount",
                units="kg m^-2 s^-1")
    cube.add_dim_coord(
        DimCoord(np.linspace(0.0, float(num_grid_points - 1),
                             num_grid_points),
                 'latitude',
                 units='degrees'),
        1
    )
    cube.add_dim_coord(
        DimCoord(np.linspace(0.0, float(num_grid_points - 1),
                             num_grid_points),
                 'longitude',
                 units='degrees'),
        2
    )
    tunit = Unit("hours since 1970-01-01 00:00:00", "gregorian")
    time_points = [402192.5 + _ for _ in range(num_time_points)]
    cube.add_aux_coord(AuxCoord(time_points,
                                "time", units=tunit), 0)
    return cube


class Test_operation_radius_to_grid_cells(IrisTest):

    """Test conversion of kernel radius in kilometres to grid cells."""

    RADIUS_IN_KM = 6.1

    def test_basic_radius_to_grid_cells(self):
        """Test the lat-long radius-to-grid-cell conversion."""
        cube = set_up_cube()
        plugin = NBHood(self.RADIUS_IN_KM)
        result = plugin.get_grid_x_y_kernel_ranges(cube)
        self.assertEqual(result, (3, 3))

    def test_basic_radius_to_grid_cells_km_grid(self):
        """Test the radius-to-grid-cell conversion, grid in km."""
        cube = set_up_cube()
        cube.coord("projection_x_coordinate").convert_units("kilometres")
        cube.coord("projection_y_coordinate").convert_units("kilometres")
        plugin = NBHood(self.RADIUS_IN_KM)
        result = plugin.get_grid_x_y_kernel_ranges(cube)
        self.assertEqual(result, (3, 3))


class Test_operation_neighbourhooding(IrisTest):

    """Test neighbourhood processing plugin on the OS National Grid."""

    RADIUS_IN_KM = 6.3  # Gives 3 grid cells worth.

    def test_basic(self):
        """Test that the plugin returns an iris.cube.Cube."""
        cube = set_up_cube()
        plugin = NBHood(self.RADIUS_IN_KM)
        result = plugin.process(cube)
        self.assertIsInstance(result, Cube)

    def test_single_point(self):
        """Test behaviour for a single non-zero grid cell."""
        cube = set_up_cube()
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
            expected[0][5 + index][5:10] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_flat(self):
        """Test behaviour for a single non-zero grid cell, flat weighting.

        Note that this gives one more grid cell range than weighted! As the
        affected area is one grid cell more in each direction, an equivalent
        range of 2 was chosen for this test.

        """
        cube = set_up_cube()
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_2_CENTROID_FLAT):
            expected[0][5 + index][5:10] = slice_
        radius_in_km = 4.2  # Equivalent to a range of 2.
        result = NBHood(radius_in_km, unweighted_mode=True).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_multi_point_multitimes(self):
        """Test behaviour for points over multiple times."""
        cube = set_up_cube(
            zero_point_indices=[(0, 10, 10), (1, 7, 7)],
            num_time_points=2
        )
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
            expected[0][8 + index][8:13] = slice_
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
            expected[1][5 + index][5:10] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_nan(self):
        """Test behaviour for a single NaN grid cell."""
        cube = set_up_cube()
        cube.data[0][6][7] = np.NAN
        msg = "NaN detected in input cube data"
        with self.assertRaisesRegexp(ValueError, msg):
            NBHood(self.RADIUS_IN_KM).process(cube)

    def test_single_point_lat_long(self):
        """Test behaviour for a single grid cell on lat long grid."""
        cube = set_up_cube_lat_long()
        msg = "Invalid grid: projection_x/y coords required"
        expected = np.zeros_like(cube.data)
        with self.assertRaisesRegexp(ValueError, msg):
            NBHood(self.RADIUS_IN_KM).process(cube)

    def test_single_point_masked_to_null(self):
        """Test behaviour with a masked non-zero point.

        The behaviour here is not right, as the mask is ignored.
        This comes directly from the scipy.ndimage.correlate base
        behaviour.

        """
        cube = set_up_cube()
        expected = np.ones_like(cube.data)
        mask = np.zeros_like(cube.data)
        mask[0][7][7] = 1
        cube.data = np.ma.masked_array(cube.data, mask=mask)
        for time_index in range(len(expected)):
            for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
                expected[time_index][5 + index][5:10] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_masked_other_point(self):
        """Test behaviour with a non-zero point next to a masked point.

        The behaviour here is not right, as the mask is ignored.

        """
        cube = set_up_cube()
        expected = np.ones_like(cube.data)
        mask = np.zeros_like(cube.data)
        mask[0][6][7] = 1
        cube.data = np.ma.masked_array(cube.data, mask=mask)
        for time_index in range(len(expected)):
            for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
                expected[time_index][5 + index][5:10] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_range_negative(self):
        """Test behaviour with a non-zero point with negative range."""
        cube = set_up_cube()
        msg = "negative dimensions are not allowed"
        with self.assertRaisesRegexp(ValueError, msg):
            NBHood(-1.0 * self.RADIUS_IN_KM).process(cube)

    def test_single_point_range_0(self):
        """Test behaviour with a non-zero point with zero range."""
        cube = set_up_cube()
        msg = "radius of 0.005 km gives zero cell extent"
        with self.assertRaisesRegexp(ValueError, msg):
            expected = np.zeros_like(cube.data)
            NBHood(0.005).process(cube)

    def test_single_point_range_1(self):
        """Test behaviour with a non-zero point with unit range."""
        cube = set_up_cube()
        expected = np.ones_like(cube.data)
        expected[0][7][7] = 0.0
        radius_in_km = 2.1  # Equivalent to a range of 1 grid cell.
        result = NBHood(radius_in_km).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_range_5(self):
        """Test behaviour with a non-zero point with a large range."""
        cube = set_up_cube()
        expected = np.ones_like(cube.data)
        for time_index in range(len(expected)):
            for index, slice_ in enumerate(SINGLE_POINT_RANGE_5_CENTROID):
                expected[time_index][3 + index][3:12] = slice_
        radius_in_km = 10.5  # Equivalent to a range of 5 grid cells.
        result = NBHood(radius_in_km).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_range_5_small_domain(self):
        """Test behaviour - non-zero point, small domain, large range.

        This exhibits the undesirable edge reflection behaviour.

        """
        cube = set_up_cube(zero_point_indices=((0, 1, 1),), num_grid_points=4)
        expected = np.array([
            [[0.97636177, 0.97533402, 0.97636177, 0.97944502],
             [0.97533402, 0.97430627, 0.97533402, 0.97841727],
             [0.97636177, 0.97533402, 0.97636177, 0.97944502],
             [0.97944502, 0.97841727, 0.97944502, 0.98252826]]
        ])
        radius_in_km = 10.5  # Equivalent to a range of 5 grid cells.
        result = NBHood(radius_in_km).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_range_lots(self):
        """Test behaviour with a non-zero point with unhandleable range."""
        cube = set_up_cube()
        msg = "radius of 500000.0 km exceeds maximum grid cell extent"
        with self.assertRaisesRegexp(ValueError, msg):
            expected = np.zeros_like(cube.data)
            NBHood(500000.0).process(cube)

    def test_point_pair(self):
        """Test behaviour for two nearby non-zero grid cells."""
        cube = set_up_cube(
            zero_point_indices=[(0, 7, 6), (0, 7, 8)])
        expected_snippet = np.array([
            [0.992, 0.968, 0.952, 0.936, 0.952, 0.968, 0.992],
            [0.968, 0.944, 0.904, 0.888, 0.904, 0.944, 0.968],
            [0.96, 0.936, 0.888, 0.872, 0.888, 0.936, 0.96],
            [0.968, 0.944, 0.904, 0.888, 0.904, 0.944, 0.968],
            [0.992, 0.968, 0.952, 0.936, 0.952, 0.968, 0.992]
        ])
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(expected_snippet):
            expected[0][5 + index][4:11] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_almost_edge(self):
        """Test behaviour for a non-zero grid cell quite near the edge."""
        cube = set_up_cube(
            zero_point_indices=[(0, 7, 2)])  # Just within range of the edge.
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
            expected[0][5 + index][0:5] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_adjacent_edge(self):
        """Test behaviour for a single non-zero grid cell near the edge."""
        cube = set_up_cube(
            zero_point_indices=[(0, 7, 1)])  # Range 3 goes over the edge.
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
            expected[0][5 + index][0:4] = slice_[1:]
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_on_edge(self):
        """Test behaviour for a non-zero grid cell on the edge.

        Note that this behaviour is 'wrong' and is a result of
        scipy.ndimage.correlate 'nearest' mode. We need to fix
        this in the future.

        """
        cube = set_up_cube(
            zero_point_indices=[(0, 7, 0)])  # On the (y) edge.
        expected = np.ones_like(cube.data)
        expected_centroid = np.array([
            [0.92, 0.96, 0.992],
            [0.848, 0.912, 0.968],
            [0.824, 0.896, 0.96],
            [0.848, 0.912, 0.968],
            [0.92, 0.96, 0.992],
        ])
        for index, slice_ in enumerate(expected_centroid):
            expected[0][5 + index][0:3] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_almost_corner(self):
        """Test behaviour for a non-zero grid cell quite near a corner."""
        cube = set_up_cube(
            zero_point_indices=[(0, 2, 2)])  # Just within corner range.
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
            expected[0][index][0:5] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_adjacent_corner(self):
        """Test behaviour for a non-zero grid cell near the corner."""
        cube = set_up_cube(
            zero_point_indices=[(0, 1, 1)])  # Kernel goes over the corner.
        expected = np.ones_like(cube.data)
        for index, slice_ in enumerate(SINGLE_POINT_RANGE_3_CENTROID):
            if index == 0:
                continue
            expected[0][index - 1][0:4] = slice_[1:]
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_single_point_on_corner(self):
        """Test behaviour for a single non-zero grid cell on the corner.

        Note that this behaviour is 'wrong' and is a result of
        scipy.ndimage.correlate 'nearest' mode. We need to fix
        this in the future.

        """
        cube = set_up_cube(
            zero_point_indices=[(0, 0, 0)])  # Point is right on the corner.
        expected = np.ones_like(cube.data)
        expected_centroid = np.array([
            [0.592, 0.768, 0.92],
            [0.768, 0.872, 0.96],
            [0.92, 0.96, 0.992],
        ])
        for index, slice_ in enumerate(expected_centroid):
            expected[0][index][0:3] = slice_
        result = NBHood(self.RADIUS_IN_KM).process(cube)
        self.assertArrayAlmostEqual(result.data, expected)

    def test_fail_multiple_realisations(self):
        """Test failing when the array has a realisation dimension."""
        data = np.ones((14, 1, 16, 16))
        data[0][0][7][7] = 0.0

        cube = Cube(data, standard_name="precipitation_amount",
                    units="kg m^-2 s^-1")
        num_grid_points = 16
        coord_system = OSGB()
        scaled_y_coord = OSGBGRID.coord('projection_y_coordinate')
        cube.add_dim_coord(
            DimCoord(
                scaled_y_coord.points[:num_grid_points],
                'projection_y_coordinate',
                units='m', coord_system=coord_system
            ),
            2
        )
        scaled_x_coord = OSGBGRID.coord('projection_x_coordinate')
        cube.add_dim_coord(
            DimCoord(
                scaled_x_coord.points[:num_grid_points],
                'projection_x_coordinate',
                units='m', coord_system=coord_system
            ),
            3
        )
        time_origin = "hours since 1970-01-01 00:00:00"
        calendar = "gregorian"
        tunit = Unit(time_origin, calendar)
        cube.add_aux_coord(AuxCoord([402192.5],
                                    "time", units=tunit), 1)
        cube.add_aux_coord(AuxCoord(np.array(range(14)),
                                    standard_name="realization"), 0)
        msg = "Does not operate across realizations"
        with self.assertRaisesRegexp(ValueError, msg):
            plugin = NBHood(self.RADIUS_IN_KM).process(cube)


if __name__ == '__main__':
    unittest.main()
